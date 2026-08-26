"""ROI-aware ECM analysis pipeline for annotated QuPath coverage exports.

Input: 15 ROI-annotated QuPath coverage workbooks placed in
``inputs/ecm_coverage``. The pipeline does not modify these input workbooks.

Analysis hierarchy:
    tile -> area-weighted ROI coverage -> section/region summary

The violin/box plots use ROI-level values. Statistical inference uses a
mixed-effects model with Sample_ID as the grouping variable, so ROIs from the
same section are not treated as independent biological samples.
"""

import glob
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.formula.api as smf
from scipy.stats import norm, pearsonr
from statsmodels.stats.multitest import multipletests


warnings.filterwarnings('ignore', category=FutureWarning)

# ==========================================
# 1. Paths and plot settings
# ==========================================
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUT_DIR = os.path.join(PROJECT_DIR, 'inputs', 'ecm_coverage')
OUT_DIR = os.path.join(PROJECT_DIR, 'outputs', 'ecm')
os.makedirs(OUT_DIR, exist_ok=True)

print(f'Loading ROI-annotated coverage workbooks from: {INPUT_DIR}')
print('Tile-level data will be summarised to area-weighted ROI values before analysis.')

sns.set_theme(style='ticks')
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 300

GROUP_ORDER = ['Normal Skin', 'Wound Centre', 'Distant Edge', 'Proximal Edge']
COMPARISON_ORDER = ['Wound Centre', 'Distant Edge', 'Proximal Edge']
REGION_MAP = {
    'normal': 'Normal Skin',
    'normal skin': 'Normal Skin',
    'center': 'Wound Centre',
    'wound center': 'Wound Centre',
    'wound centre': 'Wound Centre',
    'wound skin': 'Wound Centre',
    'edge_distant': 'Distant Edge',
    'distant edge': 'Distant Edge',
    'edge_proximal': 'Proximal Edge',
    'proximal edge': 'Proximal Edge',
}

METRICS_INFO = {
    'Collagen': {
        'file_keyword': 'collagen coverage',
        'coverage_suffix': 'collagen %',
        'component_area_label': 'collagen area',
        'other_area_label': 'non-collagen area',
    },
    'Mucin': {
        'file_keyword': 'mucin coverage',
        'coverage_suffix': 'mucin %',
        'component_area_label': 'mucin area',
        'other_area_label': 'background area',
    },
    'Necrosis': {
        'file_keyword': 'necrosis coverage',
        'coverage_suffix': 'necrotic tissue %',
        'component_area_label': 'necrotic tissue area',
        'other_area_label': 'background area',
    },
}


# ==========================================
# 2. Data loading and ROI-level aggregation
# ==========================================
def find_column_by_suffix(columns, required_suffix):
    """Find a QuPath export column by its text following the colon.

    Area columns are matched by a stable text prefix so that Excel/Windows
    encoding variations in the micrometre-squared unit do not affect loading.
    """
    for column in columns:
        suffix = str(column).split(':')[-1].strip().lower()
        required = required_suffix.lower()
        if suffix == required or suffix.startswith(f'{required} '):
            return column
    return None


def standardise_region(value):
    return REGION_MAP.get(str(value).strip().lower(), str(value).strip())


def metric_from_filename(filename):
    lowered = filename.lower()
    for metric, info in METRICS_INFO.items():
        if info['file_keyword'] in lowered:
            return metric
    return None


coverage_files = sorted(glob.glob(os.path.join(INPUT_DIR, '* coverage *.xlsx')))
coverage_files = [
    file_path for file_path in coverage_files
    if not os.path.basename(file_path).startswith('~$')
    and 'summary' not in os.path.basename(file_path).lower()
]

if len(coverage_files) != 15:
    print(f'Warning: expected 15 coverage workbooks, found {len(coverage_files)}.')

tile_data_by_metric = {metric: [] for metric in METRICS_INFO}

for file_path in coverage_files:
    filename = os.path.basename(file_path)
    metric = metric_from_filename(filename)
    if metric is None:
        continue

    info = METRICS_INFO[metric]
    data = pd.read_excel(file_path, engine='openpyxl')
    required_base_columns = {'Sample_ID', 'region', 'ROI_ID'}
    missing_base_columns = required_base_columns - set(data.columns)
    if missing_base_columns:
        raise ValueError(
            f'{filename} is not ROI-annotated. Missing columns: {sorted(missing_base_columns)}. '
            'Run annotate_roi_ids.py first.'
        )

    component_area_column = find_column_by_suffix(data.columns, info['component_area_label'])
    other_area_column = find_column_by_suffix(data.columns, info['other_area_label'])
    coverage_column = find_column_by_suffix(data.columns, info['coverage_suffix'])

    if not all([component_area_column, other_area_column, coverage_column]):
        raise ValueError(f'{filename}: could not find the expected area/coverage columns for {metric}.')

    selected = data[[
        'Sample_ID', 'region', 'ROI_ID', coverage_column,
        component_area_column, other_area_column,
    ]].copy()
    selected.columns = [
        'Sample_ID', 'region', 'ROI_ID', 'Tile_Coverage',
        'Component_Area', 'Other_Area',
    ]
    selected['Source_File'] = filename
    selected['Metric'] = metric
    selected['region'] = selected['region'].map(standardise_region)

    for column in ['Tile_Coverage', 'Component_Area', 'Other_Area']:
        selected[column] = pd.to_numeric(selected[column], errors='coerce')

    selected['Tissue_Area'] = selected['Component_Area'] + selected['Other_Area']
    selected = selected.dropna(subset=['Sample_ID', 'region', 'ROI_ID', 'Tissue_Area'])
    selected = selected[selected['Tissue_Area'] > 0].copy()
    tile_data_by_metric[metric].append(selected)

if not any(tile_data_by_metric.values()):
    raise ValueError('No valid tile-level data were loaded.')

roi_summaries = []
section_region_summaries = []
for metric, frames in tile_data_by_metric.items():
    if not frames:
        continue

    metric_tiles = pd.concat(frames, ignore_index=True)
    roi_summary = (
        metric_tiles.groupby(['Metric', 'Sample_ID', 'region', 'ROI_ID'], as_index=False)
        .agg(
            Tile_Count=('Tissue_Area', 'size'),
            Component_Area_Total=('Component_Area', 'sum'),
            Tissue_Area_Total=('Tissue_Area', 'sum'),
        )
    )
    roi_summary['Coverage'] = (
        100 * roi_summary['Component_Area_Total'] / roi_summary['Tissue_Area_Total']
    )
    roi_summaries.append(roi_summary)

    # Normal sections contain three ROI_ID values. They are equally weighted here.
    # Wound sections have one ROI per region, so their section/region mean is unchanged.
    section_summary = (
        roi_summary.groupby(['Metric', 'Sample_ID', 'region'], as_index=False)
        .agg(
            ROI_Count=('ROI_ID', 'nunique'),
            Coverage=('Coverage', 'mean'),
            Coverage_SD_between_ROIs=('Coverage', 'std'),
        )
    )
    section_region_summaries.append(section_summary)

df_roi = pd.concat(roi_summaries, ignore_index=True)
df_section = pd.concat(section_region_summaries, ignore_index=True)
df_roi.to_csv(os.path.join(OUT_DIR, '0_ROI_AreaWeighted_Summary.csv'), index=False)
df_section.to_csv(os.path.join(OUT_DIR, '0_SectionRegion_Summary.csv'), index=False)

print('Loaded ROI-level records by metric:')
for metric in METRICS_INFO:
    metric_roi = df_roi[df_roi['Metric'] == metric]
    metric_section = df_section[df_section['Metric'] == metric]
    if not metric_roi.empty:
        print(
            f'  {metric}: {len(metric_roi)} ROI values from '
            f'{metric_section.Sample_ID.nunique()} sections.'
        )


# ==========================================
# 3. Descriptive statistics and outlier report
# ==========================================
print('\nGenerating descriptive statistics and outlier reports...')
descriptive_statistics = (
    df_roi.groupby(['Metric', 'Sample_ID', 'region', 'ROI_ID'])['Coverage']
    .describe()
    .reset_index()
)
descriptive_statistics.to_csv(os.path.join(OUT_DIR, '1_Descriptive_Statistics.csv'), index=False)


def detect_outliers(group):
    if len(group) < 4:
        return group.iloc[0:0].copy()
    q1 = group['Coverage'].quantile(0.25)
    q3 = group['Coverage'].quantile(0.75)
    iqr = q3 - q1
    return group[
        (group['Coverage'] < q1 - 1.5 * iqr)
        | (group['Coverage'] > q3 + 1.5 * iqr)
    ]


outliers = (
    df_roi.groupby(['Metric', 'region'], group_keys=False)
    .apply(detect_outliers)
    .reset_index(drop=True)
)
outliers.to_csv(os.path.join(OUT_DIR, '2_Potential_Outliers.csv'), index=False)


# ==========================================
# 4. Mixed-effects model and Holm contrasts
# ==========================================
def fit_mixed_effects_model(metric_roi_data):
    """Fit ROI coverage ~ region with Sample_ID as the repeated-measures group."""
    model_data = metric_roi_data.copy()
    model_data['region'] = pd.Categorical(model_data['region'], categories=GROUP_ORDER, ordered=True)
    model_data = model_data.dropna(subset=['Coverage', 'region', 'Sample_ID'])

    if model_data['region'].nunique() < 2:
        raise ValueError('At least two regions are required for a mixed-effects model.')

    formula = "Coverage ~ C(region, Treatment(reference='Normal Skin'))"
    model = smf.mixedlm(formula, data=model_data, groups=model_data['Sample_ID'])
    try:
        result = model.fit(reml=False, method='lbfgs', maxiter=500, disp=False)
    except Exception:
        # Powell is slower but can resolve some boundary/convergence cases in small datasets.
        result = model.fit(reml=False, method='powell', maxiter=1000, disp=False)

    return result, model_data


def normal_vs_region_contrasts(result, comparison_regions):
    """Calculate model-based Normal Skin contrasts and apply Holm correction."""
    fe_names = result.fe_params.index.tolist()
    fixed_covariance = result.cov_params().loc[fe_names, fe_names]
    records = []

    for region in comparison_regions:
        matching_terms = [name for name in fe_names if name.endswith(f'[T.{region}]')]
        if not matching_terms:
            records.append({
                'Comparison': f'Normal Skin vs {region}',
                'Estimate': np.nan,
                'SE': np.nan,
                'Z': np.nan,
                'P_Unadjusted': np.nan,
            })
            continue

        term = matching_terms[0]
        estimate = float(result.fe_params[term])
        standard_error = float(np.sqrt(fixed_covariance.loc[term, term]))
        z_value = estimate / standard_error if standard_error > 0 else np.nan
        p_value = 2 * norm.sf(abs(z_value)) if not np.isnan(z_value) else np.nan
        records.append({
            'Comparison': f'Normal Skin vs {region}',
            'Estimate': estimate,
            'SE': standard_error,
            'Z': z_value,
            'P_Unadjusted': p_value,
        })

    contrasts = pd.DataFrame(records)
    valid = contrasts['P_Unadjusted'].notna()
    contrasts['P_Holm'] = np.nan
    if valid.any():
        contrasts.loc[valid, 'P_Holm'] = multipletests(
            contrasts.loc[valid, 'P_Unadjusted'], method='holm'
        )[1]
    return contrasts


def significance_label(p_value):
    if pd.isna(p_value) or p_value > 0.05:
        return 'ns'
    if p_value <= 0.0001:
        return '****'
    if p_value <= 0.001:
        return '***'
    if p_value <= 0.01:
        return '**'
    return '*'


def add_significance_line(axis, left_index, right_index, y_position, line_height, label):
    axis.plot(
        [left_index, left_index, right_index, right_index],
        [y_position, y_position + line_height, y_position + line_height, y_position],
        lw=1.2,
        c='black',
    )
    axis.text(
        (left_index + right_index) * 0.5,
        y_position + line_height,
        label,
        ha='center',
        va='bottom',
        color='black',
        fontsize=10,
        fontweight='bold',
    )


print('\nFitting mixed-effects models and generating violin/box plots...')
model_results = []
all_contrasts = []

for metric in METRICS_INFO:
    metric_roi = df_roi[df_roi['Metric'] == metric].copy()
    if metric_roi.empty:
        continue

    existing_groups = [group for group in GROUP_ORDER if group in metric_roi['region'].unique()]
    metric_roi['region'] = pd.Categorical(metric_roi['region'], categories=existing_groups, ordered=True)
    metric_roi = metric_roi.dropna(subset=['Coverage'])

    model_result, model_data = fit_mixed_effects_model(metric_roi)
    contrasts = normal_vs_region_contrasts(model_result, COMPARISON_ORDER)
    contrasts.insert(0, 'Metric', metric)
    all_contrasts.append(contrasts)

    model_results.append({
        'Metric': metric,
        'Model': 'Mixed-effects: ROI coverage ~ region + (1 | Sample_ID)',
        'ROI_Observations': len(model_data),
        'Independent_Sections': model_data['Sample_ID'].nunique(),
        'Converged': bool(model_result.converged),
        'Log_Likelihood': model_result.llf,
        'Random_Intercept_Variance': float(np.asarray(model_result.cov_re)[0, 0]),
    })

    palette_colors = sns.color_palette('pastel', len(existing_groups))
    box_palette = sns.color_palette('muted', len(existing_groups))
    plt.figure(figsize=(7, 6))
    axis = sns.violinplot(
        x='region', y='Coverage', data=metric_roi, order=existing_groups,
        inner=None, palette=palette_colors, cut=0, linewidth=1.2, alpha=0.5,
    )
    sns.boxplot(
        x='region', y='Coverage', data=metric_roi, order=existing_groups,
        width=0.15, palette=box_palette, showcaps=True, showfliers=False,
        whiskerprops={'color': 'black', 'linewidth': 1.2},
        boxprops={'edgecolor': 'black', 'linewidth': 1.2, 'zorder': 3},
    )

    for index, group in enumerate(existing_groups):
        group_values = metric_roi.loc[metric_roi['region'] == group, 'Coverage'].dropna()
        if not group_values.empty:
            median_value = group_values.median()
            axis.hlines(median_value, index - 0.075, index + 0.075,
                        color='black', linewidth=2.2, zorder=4)

    y_min, y_max = metric_roi['Coverage'].min(), metric_roi['Coverage'].max()
    y_range = y_max - y_min
    y_step = y_range * 0.15 if y_range > 0 else 0.5
    line_height = y_step * 0.15
    p_by_region = {
        comparison.replace('Normal Skin vs ', ''): p_value
        for comparison, p_value in zip(contrasts['Comparison'], contrasts['P_Holm'])
    }
    for index, group in enumerate(existing_groups[1:], start=1):
        add_significance_line(
            axis, 0, index, y_max + index * y_step, line_height,
            significance_label(p_by_region.get(group, np.nan)),
        )

    axis.set_ylim(y_min - y_range * 0.1 if y_range > 0 else y_min - 0.5,
                  y_max + (len(existing_groups) + 1) * y_step)
    # Figure title intentionally omitted; use the thesis figure caption instead.
    plt.ylabel(f'{metric} Coverage (%)', fontsize=12, fontweight='bold')
    plt.xlabel('', fontsize=12, fontweight='bold')
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, f'3_ViolinBoxPlot_{metric}.png'), dpi=300)
    plt.close()

pd.DataFrame(model_results).to_csv(os.path.join(OUT_DIR, '5_Main_MixedEffects_Results.csv'), index=False)
if all_contrasts:
    contrast_results = pd.concat(all_contrasts, ignore_index=True)
    contrast_results.to_csv(os.path.join(OUT_DIR, '5_Posthoc_Holm_Contrasts.csv'), index=False)


# ==========================================
# 5. Correlation scatter plots using section/region summaries
# ==========================================
print('\nGenerating correlation scatter plots with regression lines...')
correlation_wide = (
    df_section.pivot(index=['Sample_ID', 'region'], columns='Metric', values='Coverage')
    .reset_index()
)
correlation_wide.to_csv(os.path.join(OUT_DIR, '6_Aggregated_Means_for_Correlation.csv'), index=False)

available_metrics = [metric for metric in METRICS_INFO if metric in correlation_wide.columns]
correlation_results = []
for first_index, metric_x in enumerate(available_metrics):
    for metric_y in available_metrics[first_index + 1:]:
        plot_data = correlation_wide[[metric_x, metric_y]].dropna()
        if len(plot_data) < 3:
            continue

        r_value, p_value = pearsonr(plot_data[metric_x], plot_data[metric_y])
        correlation_results.append({
            'Metric_X': metric_x,
            'Metric_Y': metric_y,
            'N_SectionRegion_Values': len(plot_data),
            'Pearson_r': r_value,
            'P_Value': p_value,
        })

        plt.figure(figsize=(6, 6))
        sns.regplot(
            data=plot_data, x=metric_x, y=metric_y,
            color='#2b5c8f', marker='o',
            scatter_kws={'s': 50, 'alpha': 0.8, 'edgecolor': 'black', 'linewidths': 0.8},
            line_kws={'color': '#d95f02', 'linewidth': 2},
        )
        # Figure title intentionally omitted; use the thesis figure caption instead.
        plt.xlabel(f'{metric_x} Coverage (%)', fontsize=11, fontweight='bold')
        plt.ylabel(f'{metric_y} Coverage (%)', fontsize=11, fontweight='bold')

        if metric_x == 'Collagen' and metric_y == 'Mucin':
            plt.ylim(bottom=0)

        label = f"r = {r_value:.3f}\np {'< 0.001' if p_value < 0.001 else f'= {p_value:.3f}'}"
        plt.gca().text(
            0.05, 0.92, label, transform=plt.gca().transAxes,
            fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'),
        )
        sns.despine()
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, f'8_Correlation_{metric_x}_vs_{metric_y}.png'), dpi=300)
        plt.close()

pd.DataFrame(correlation_results).to_csv(os.path.join(OUT_DIR, '8_Correlation_Results.csv'), index=False)
print(f'\nAll analyses and figure generation completed successfully. Results are saved in: {OUT_DIR}')
