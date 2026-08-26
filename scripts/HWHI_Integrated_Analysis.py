"""Build an exploratory Histological Wound Healing Index (HWHI).

This is a portable, reproducible pipeline. It reads the outputs created by
the ECM and FracLac pipelines in this repository:

outputs/
├── ecm/0_SectionRegion_Summary.csv
└── fraclac/0_FracLac_QC_Passed_Data.csv

It never modifies the two input files. It joins spatially matched
section-region observations and creates a new output folder.

Primary index, prespecified with equal component weights:
HWHI = [z(Collagen) + z(Fractal Dimension) - z(Mucin)
        - z(Necrosis) - z(Lacunarity)] / 5

Higher values indicate a more normal-like / structurally mature ECM phenotype.
This is an exploratory research composite, not a clinical diagnostic score.
"""

import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.formula.api as smf
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


warnings.filterwarnings('ignore', category=FutureWarning)

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ECM_FILE = os.path.join(PROJECT_DIR, 'outputs', 'ecm', '0_SectionRegion_Summary.csv')
FRACLAC_FILE = os.path.join(PROJECT_DIR, 'outputs', 'fraclac', '0_FracLac_QC_Passed_Data.csv')
OUT_DIR = os.path.join(PROJECT_DIR, 'outputs', 'hwhi')
os.makedirs(OUT_DIR, exist_ok=True)

REGION_ORDER = ['Normal Skin', 'Wound Centre', 'Distant Edge', 'Proximal Edge']
COMPARISON_ORDER = ['Wound Centre', 'Distant Edge', 'Proximal Edge']
INDEX_DEFINITION = (
    'HWHI = [z(Collagen) + z(Fractal Dimension) - z(Mucin) '
    '- z(Necrosis) - z(Lacunarity)] / 5'
)


def require_columns(table, expected_columns, table_name):
    missing = [column for column in expected_columns if column not in table.columns]
    if missing:
        raise ValueError(f'{table_name} is missing required columns: {missing}')


def normalise_region(value):
    """Make ECM and FracLac labels identical before joining."""
    mapping = {
        'Wound Skin': 'Wound Centre',
        'Wound Center': 'Wound Centre',
        'Center': 'Wound Centre',
    }
    return mapping.get(str(value).strip(), str(value).strip())


def z_score(series):
    """Population-standardise a component across the integrated observations."""
    standard_deviation = series.std(ddof=0)
    if pd.isna(standard_deviation) or standard_deviation == 0:
        raise ValueError(f'Cannot standardise {series.name}: standard deviation is zero or missing.')
    return (series - series.mean()) / standard_deviation


def fit_mixed_model(data):
    """Fit HWHI ~ Region with Sample_ID as the repeated-measures grouping variable."""
    model_data = data.dropna(subset=['HWHI', 'Region', 'Sample_ID']).copy()
    formula = "HWHI ~ C(Region, Treatment(reference='Normal Skin'))"
    model = smf.mixedlm(formula, data=model_data, groups=model_data['Sample_ID'])
    try:
        result = model.fit(reml=False, method='lbfgs', maxiter=500, disp=False)
    except Exception:
        result = model.fit(reml=False, method='powell', maxiter=1000, disp=False)
    return result, model_data


def normal_vs_region_contrasts(result):
    fixed_names = result.fe_params.index.tolist()
    fixed_covariance = result.cov_params().loc[fixed_names, fixed_names]
    rows = []
    for region in COMPARISON_ORDER:
        candidates = [name for name in fixed_names if name.endswith(f'[T.{region}]')]
        if not candidates:
            rows.append({'Comparison': f'Normal Skin vs {region}', 'Estimate': np.nan,
                         'SE': np.nan, 'Z': np.nan, 'P_Unadjusted': np.nan})
            continue
        term = candidates[0]
        estimate = float(result.fe_params[term])
        standard_error = float(np.sqrt(fixed_covariance.loc[term, term]))
        z_value = estimate / standard_error if standard_error > 0 else np.nan
        p_value = 2 * norm.sf(abs(z_value)) if not np.isnan(z_value) else np.nan
        rows.append({'Comparison': f'Normal Skin vs {region}', 'Estimate': estimate,
                     'SE': standard_error, 'Z': z_value, 'P_Unadjusted': p_value})
    contrasts = pd.DataFrame(rows)
    contrasts['P_Holm'] = np.nan
    valid = contrasts['P_Unadjusted'].notna()
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


def add_significance_line(axis, left, right, y_position, height, label):
    axis.plot([left, left, right, right],
              [y_position, y_position + height, y_position + height, y_position],
              lw=1.2, c='black')
    axis.text((left + right) / 2, y_position + height, label, ha='center', va='bottom',
              color='black', fontsize=10, fontweight='bold')


def plot_hwhi(data, contrasts):
    """Create a title-free violin/box/point plot for the thesis figure."""
    plt.figure(figsize=(8, 6))
    axis = sns.violinplot(x='Region', y='HWHI', data=data, order=REGION_ORDER,
                          palette='pastel', inner=None, cut=0, linewidth=1.2)
    sns.boxplot(x='Region', y='HWHI', data=data, order=REGION_ORDER,
                palette='muted', width=0.15, showfliers=False,
                whiskerprops={'color': 'black', 'linewidth': 1.2},
                boxprops={'edgecolor': 'black', 'linewidth': 1.2, 'zorder': 3})
    sns.stripplot(x='Region', y='HWHI', data=data, order=REGION_ORDER,
                  color='black', alpha=0.6, jitter=0.15, size=4, zorder=4)

    y_min, y_max = data['HWHI'].min(), data['HWHI'].max()
    y_range = y_max - y_min
    line_step = y_range * 0.14 if y_range > 0 else 0.1
    line_height = line_step * 0.18
    p_by_region = {
        row['Comparison'].replace('Normal Skin vs ', ''): row['P_Holm']
        for _, row in contrasts.iterrows()
    }
    for index, region in enumerate(COMPARISON_ORDER, start=1):
        add_significance_line(axis, 0, index, y_max + index * line_step, line_height,
                              significance_label(p_by_region.get(region, np.nan)))
    axis.set_ylim(y_min - y_range * 0.10 if y_range > 0 else y_min - 0.1,
                  y_max + (len(COMPARISON_ORDER) + 1) * line_step)
    plt.ylabel('HWHI (standardised composite score)', fontsize=12, fontweight='bold')
    plt.xlabel('')
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, '3_HWHI_ViolinBoxPlot.png'), dpi=300, bbox_inches='tight')
    plt.close()


if not os.path.isfile(ECM_FILE):
    raise FileNotFoundError(
        'ECM summary table not found. Run HistologyData_Pipeline_MixedEffects.py first. '\
        f'Expected: {ECM_FILE}'
    )
if not os.path.isfile(FRACLAC_FILE):
    raise FileNotFoundError(
        'FracLac QC-passed table not found. Run FracLac_Spatial_Analysis_MixedEffects_HWHI_Ready.py first. '\
        f'Expected: {FRACLAC_FILE}'
    )

print('Reading ECM and FracLac summaries...')
ecm_long = pd.read_csv(ECM_FILE)
fraclac_raw = pd.read_csv(FRACLAC_FILE)
require_columns(ecm_long, ['Metric', 'Sample_ID', 'region', 'Coverage'], 'ECM summary table')
require_columns(fraclac_raw, ['Sample_ID', 'Region', 'Fractal_Dimension', 'Lacunarity',
                              'Low_Foreground_Fraction_Flag'], 'FracLac table')

ecm_long['Region'] = ecm_long['region'].map(normalise_region)
ecm_wide = ecm_long.pivot(index=['Sample_ID', 'Region'], columns='Metric', values='Coverage').reset_index()
required_ecm_metrics = ['Collagen', 'Mucin', 'Necrosis']
require_columns(ecm_wide, required_ecm_metrics, 'pivoted ECM summary table')

# Three normal FracLac masks are averaged to one section-region value;
# each wound section-region has one FracLac mask.
fraclac_raw['Region'] = fraclac_raw['Region'].map(normalise_region)
fraclac_section = (
    fraclac_raw.groupby(['Sample_ID', 'Region'], as_index=False)
    .agg(
        Fractal_Dimension=('Fractal_Dimension', 'mean'),
        Lacunarity=('Lacunarity', 'mean'),
        FracLac_Mask_Count=('Mask_ID', 'count'),
        Low_Foreground_Fraction_Flag=('Low_Foreground_Fraction_Flag', 'max'),
        Minimum_Foreground_Fraction=('Foreground_Fraction', 'min'),
    )
)

master = ecm_wide.merge(fraclac_section, on=['Sample_ID', 'Region'], how='outer', indicator=True)
master['Integration_Status'] = master['_merge'].map({
    'both': 'Matched ECM and FracLac',
    'left_only': 'ECM only - investigate',
    'right_only': 'FracLac only - investigate',
})
master = master.drop(columns='_merge')
master['Region'] = pd.Categorical(master['Region'], categories=REGION_ORDER, ordered=True)
master = master.sort_values(['Region', 'Sample_ID']).reset_index(drop=True)
master.to_csv(os.path.join(OUT_DIR, '0_HWHI_Integration_QC.csv'), index=False)

unmatched = master[master['Integration_Status'] != 'Matched ECM and FracLac']
if not unmatched.empty:
    raise ValueError('ECM and FracLac records did not fully match. See 0_HWHI_Integration_QC.csv.')
if len(master) != 22:
    raise ValueError(f'Expected 22 matched section-region observations, found {len(master)}.')

for component in required_ecm_metrics + ['Fractal_Dimension', 'Lacunarity']:
    master[f'Z_{component}'] = z_score(master[component])

master['HWHI'] = (
    master['Z_Collagen'] + master['Z_Fractal_Dimension']
    - master['Z_Mucin'] - master['Z_Necrosis'] - master['Z_Lacunarity']
) / 5
master['Index_Definition'] = INDEX_DEFINITION
master['Interpretation'] = np.where(
    master['HWHI'] >= 0,
    'More normal-like / structurally mature composite phenotype',
    'More wound-associated / structurally disrupted composite phenotype',
)
master.to_csv(os.path.join(OUT_DIR, '0_HWHI_Master_Table.csv'), index=False)

print(f'Matched {len(master)} section-region observations from {master.Sample_ID.nunique()} sections.')
print(f'Low-foreground FracLac flags retained in primary HWHI: '
      f'{master.Low_Foreground_Fraction_Flag.sum()}')

descriptive = master.groupby('Region', observed=True)['HWHI'].agg(
    Count='count', Mean='mean', SD='std', Median='median', Min='min', Max='max'
).reset_index()
descriptive.to_csv(os.path.join(OUT_DIR, '1_HWHI_Descriptive_Statistics.csv'), index=False)

print('Fitting the primary mixed-effects model and generating the figure...')
result, model_data = fit_mixed_model(master)
contrasts = normal_vs_region_contrasts(result)
pd.DataFrame([{
    'Metric': 'HWHI',
    'Index_Definition': INDEX_DEFINITION,
    'Model': 'Mixed-effects: HWHI ~ region + (1 | Sample_ID)',
    'Section_Region_Observations': len(model_data),
    'Independent_Sections': model_data['Sample_ID'].nunique(),
    'Converged': bool(result.converged),
    'Log_Likelihood': result.llf,
    'Random_Intercept_Variance': float(np.asarray(result.cov_re)[0, 0]),
}]).to_csv(os.path.join(OUT_DIR, '5_HWHI_MixedEffects_Results.csv'), index=False)
contrasts.to_csv(os.path.join(OUT_DIR, '5_HWHI_Posthoc_Holm_Contrasts.csv'), index=False)
plot_hwhi(model_data, contrasts)

# Prespecified transparency output: descriptive sensitivity analysis after omitting
# section-regions whose FracLac mask had <1% foreground pixels. The primary HWHI
# retains them, because a collagen-poor wound centre is itself biologically relevant.
sensitivity = master[~master['Low_Foreground_Fraction_Flag']].copy()
sensitivity.groupby('Region', observed=True)['HWHI'].agg(
    Count='count', Mean='mean', SD='std', Median='median', Min='min', Max='max'
).reset_index().to_csv(
    os.path.join(OUT_DIR, '6_HWHI_Sensitivity_Descriptive_Excluding_LowForeground.csv'), index=False
)

with open(os.path.join(OUT_DIR, 'README_HWHI_Method.txt'), 'w', encoding='utf-8') as handle:
    handle.write('Exploratory Histological Wound Healing Index (HWHI)\n\n')
    handle.write(INDEX_DEFINITION + '\n\n')
    handle.write('Positive directions: collagen coverage and fractal dimension.\n')
    handle.write('Negative directions: mucin coverage, necrosis coverage, and lacunarity.\n')
    handle.write('All five components have equal weight after z-standardisation across the 22 matched section-region observations.\n')
    handle.write('ECM and FracLac metrics were spatially matched by section and anatomical region, but were not necessarily derived from identical ROI boundaries.\n')
    handle.write('The primary index retains low-foreground FracLac masks and reports a descriptive sensitivity table excluding them.\n')
    handle.write('HWHI is exploratory and must not be presented as a validated clinical or diagnostic score.\n')

print(f'Analysis complete. Results saved in: {OUT_DIR}')
