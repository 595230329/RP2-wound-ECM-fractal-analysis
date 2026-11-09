"""Analyse FracLac Box Count Summary workbooks at the ROI level.

Place FracLac exports in the two project folders:
    inputs/fraclac/normal_skin_boxcount
    inputs/fraclac/wound_skin_boxcount

It reads the 30 original FracLac exports, creates an auditable master table,
and compares Normal Skin, Wound Centre, Proximal Edge and Distant Edge using
mixed-effects models with Sample_ID as the repeated-measures grouping variable.
"""

import glob
import os
import re
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
INPUT_DIR = os.path.join(PROJECT_DIR, 'inputs', 'fraclac')
NORMAL_DIR = os.path.join(INPUT_DIR, 'normal_skin_boxcount')
WOUND_DIR = os.path.join(INPUT_DIR, 'wound_skin_boxcount')
OUT_DIR = os.path.join(PROJECT_DIR, 'outputs', 'fraclac')
os.makedirs(OUT_DIR, exist_ok=True)

REGION_ORDER = ['Normal Skin', 'Wound Centre', 'Distant Edge', 'Proximal Edge']
COMPARISON_ORDER = ['Wound Centre', 'Distant Edge', 'Proximal Edge']
R2_QC_THRESHOLD = 0.95
LOW_FOREGROUND_FRACTION_THRESHOLD = 0.01  # QC flag only; no masks are excluded by this rule.


def numeric_value(row, zero_based_column_index, source_file):
    """Read a numeric value from a known FracLac summary position."""
    value = pd.to_numeric(row.iloc[zero_based_column_index], errors='coerce')
    if pd.isna(value):
        raise ValueError(
            f'{source_file}: FracLac value at column {zero_based_column_index + 1} is not numeric.'
        )
    return float(value)


def parse_file_metadata(path, is_normal):
    """Parse sample, region and mask ID from the user-defined file naming scheme."""
    filename = os.path.basename(path)
    upper_name = filename.upper()
    sample_match = re.search(r'9M(0501|0502|06|07|08|09)', upper_name)
    if not sample_match:
        raise ValueError(f'Could not identify the sample code in {filename}.')

    sample_code = sample_match.group(1)
    if is_normal:
        sample_id = f'Normal_{sample_code}'
        roi_match = re.search(r'NORMAL\s+SKIN\s*0?(\d+)', upper_name)
        if not roi_match:
            raise ValueError(f'Could not identify the normal ROI number in {filename}.')
        roi_number = int(roi_match.group(1))
        region = 'Normal Skin'
        mask_id = f'{sample_id}_ROI_{roi_number}'
    else:
        # The wound-image file names are explicitly mapped to the Sample_ID values
        # used by the ECM coverage dataset, enabling section-region HWHI integration.
        wound_sample_map = {
            '0501': 'Wound_05_1',
            '0502': 'Wound_05_2',
            '06': 'Wound_06',
            '07': 'Wound_07',
            '08': 'Wound_08',
            '09': 'Wound_09',
        }
        sample_id = wound_sample_map[sample_code]
        if 'DISTANT EDGE' in upper_name:
            region = 'Distant Edge'
        elif 'PROXIMAL EDGE' in upper_name:
            region = 'Proximal Edge'
        elif 'WOUND CENTER' in upper_name or 'WOUND CENTRE' in upper_name:
            region = 'Wound Centre'
        else:
            raise ValueError(f'Could not identify the wound region in {filename}.')
        mask_id = f"{sample_id}_{region.replace(' ', '_')}"

    return sample_id, region, mask_id


def read_fraclac_workbook(path, is_normal):
    """Extract the prespecified FracLac outputs from a Box Count Summary workbook."""
    source_file = os.path.basename(path)
    table = pd.read_excel(path, header=None)
    if table.shape[0] < 2 or table.shape[1] < 120:
        raise ValueError(f'{source_file}: unexpected FracLac summary layout.')

    values = table.iloc[1]
    sample_id, region, mask_id = parse_file_metadata(path, is_normal)

    # Prespecified FracLac fields, by spreadsheet column number (1-based):
    # 6 = mean D; 11 = star D r2; 87 = lacunarity for image;
    # 118 = mean foreground pixels; 120 = total pixels.
    # The indices passed below are the same columns, zero-based.
    mean_foreground_pixels = numeric_value(values, 117, source_file)
    total_pixels = numeric_value(values, 119, source_file)

    return {
        'Source_File': source_file,
        'Sample_ID': sample_id,
        'Region': region,
        'Mask_ID': mask_id,
        'Fractal_Dimension': numeric_value(values, 5, source_file),
        'R2_Fit': numeric_value(values, 10, source_file),
        'Lacunarity': numeric_value(values, 86, source_file),
        'Mean_Foreground_Pixels': mean_foreground_pixels,
        'Total_Pixels': total_pixels,
        'Foreground_Fraction': mean_foreground_pixels / total_pixels,
        'Scan_Positions': numeric_value(values, 1, source_file),
    }


if not os.path.isdir(NORMAL_DIR) or not os.path.isdir(WOUND_DIR):
    raise FileNotFoundError(
        'Required folders were not found. See inputs/README.md for the required '\
        "'normal_skin_boxcount' and 'wound_skin_boxcount' folder structure."
    )

normal_files = sorted(glob.glob(os.path.join(NORMAL_DIR, '*.xlsx')))
wound_files = sorted(glob.glob(os.path.join(WOUND_DIR, '*.xlsx')))
if not normal_files or not wound_files:
    raise FileNotFoundError('No .xlsx FracLac workbooks were found in one or both input folders.')

print(f'Found {len(normal_files)} normal and {len(wound_files)} wound FracLac workbooks.')
records = [read_fraclac_workbook(path, is_normal=True) for path in normal_files]
records.extend(read_fraclac_workbook(path, is_normal=False) for path in wound_files)
df_all = pd.DataFrame(records)
duplicate_masks = df_all.loc[df_all['Mask_ID'].duplicated(keep=False), ['Mask_ID', 'Source_File']]
if not duplicate_masks.empty:
    raise ValueError(
        'More than one FracLac workbook maps to the same mask:\n'
        f'{duplicate_masks.sort_values("Mask_ID").to_string(index=False)}'
    )
if len(normal_files) != 12 or len(wound_files) != 18:
    print(f'Warning: expected 12 normal and 18 wound FracLac workbooks, found '
          f'{len(normal_files)} and {len(wound_files)}.')
df_all['Region'] = pd.Categorical(df_all['Region'], categories=REGION_ORDER, ordered=True)
df_all['Passes_R2_QC'] = df_all['R2_Fit'] >= R2_QC_THRESHOLD
df_all['Low_Foreground_Fraction_Flag'] = (
    df_all['Foreground_Fraction'] < LOW_FOREGROUND_FRACTION_THRESHOLD
)
df_all.to_csv(os.path.join(OUT_DIR, '0_FracLac_Mask_Master_Table.csv'), index=False)

print('\nMask counts by region:')
print(df_all.groupby('Region', observed=True)['Mask_ID'].nunique())
print(f'R2 QC threshold: {R2_QC_THRESHOLD}; masks passing: {df_all.Passes_R2_QC.sum()} / {len(df_all)}')
print(f'Low-foreground QC flags (< {LOW_FOREGROUND_FRACTION_THRESHOLD:.0%}): '
      f'{df_all.Low_Foreground_Fraction_Flag.sum()}')

failed_r2 = df_all.loc[~df_all['Passes_R2_QC'], ['Mask_ID', 'R2_Fit']]
if not failed_r2.empty:
    print('Masks excluded by the R2 QC (not used in models or HWHI):')
    print(failed_r2.to_string(index=False))

df_clean = df_all[df_all['Passes_R2_QC']].copy()
df_clean.to_csv(os.path.join(OUT_DIR, '0_FracLac_QC_Passed_Data.csv'), index=False)


def fit_mixed_model(data, response):
    """Fit response ~ Region with Sample_ID as the repeated-measures group."""
    model_data = data.dropna(subset=[response, 'Region', 'Sample_ID']).copy()
    formula = f"{response} ~ C(Region, Treatment(reference='Normal Skin'))"
    model = smf.mixedlm(formula, data=model_data, groups=model_data['Sample_ID'])
    try:
        result = model.fit(reml=False, method='lbfgs', maxiter=500, disp=False)
    except Exception:
        result = model.fit(reml=False, method='powell', maxiter=1000, disp=False)
    return result, model_data


def normal_vs_region_contrasts(result):
    """Calculate model-based Normal Skin contrasts with Holm correction."""
    fixed_effect_names = result.fe_params.index.tolist()
    fixed_covariance = result.cov_params().loc[fixed_effect_names, fixed_effect_names]
    records = []

    for region in COMPARISON_ORDER:
        terms = [name for name in fixed_effect_names if name.endswith(f'[T.{region}]')]
        if not terms:
            records.append({
                'Comparison': f'Normal Skin vs {region}',
                'Estimate': np.nan, 'SE': np.nan, 'Z': np.nan, 'P_Unadjusted': np.nan,
            })
            continue

        term = terms[0]
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


def add_significance_line(axis, left_index, right_index, y_position, height, label):
    axis.plot(
        [left_index, left_index, right_index, right_index],
        [y_position, y_position + height, y_position + height, y_position],
        lw=1.2, c='black',
    )
    axis.text((left_index + right_index) / 2, y_position + height, label,
              ha='center', va='bottom', color='black', fontsize=10, fontweight='bold')


def plot_metric(data, metric_column, y_label, filename, contrasts):
    """Create a title-free violin/box/ROI-point plot with adjusted significance lines."""
    plt.figure(figsize=(8, 6))
    axis = sns.violinplot(
        x='Region', y=metric_column, data=data, order=REGION_ORDER,
        palette='pastel', inner=None, cut=0, linewidth=1.2,
    )
    sns.boxplot(
        x='Region', y=metric_column, data=data, order=REGION_ORDER,
        palette='muted', width=0.15, showfliers=False,
        whiskerprops={'color': 'black', 'linewidth': 1.2},
        boxprops={'edgecolor': 'black', 'linewidth': 1.2, 'zorder': 3},
    )
    sns.stripplot(
        x='Region', y=metric_column, data=data, order=REGION_ORDER,
        color='black', alpha=0.6, jitter=0.15, size=4, zorder=4,
    )

    y_min, y_max = data[metric_column].min(), data[metric_column].max()
    y_range = y_max - y_min
    line_step = y_range * 0.14 if y_range > 0 else 0.1
    line_height = line_step * 0.18
    p_by_region = {
        row['Comparison'].replace('Normal Skin vs ', ''): row['P_Holm']
        for _, row in contrasts.iterrows()
    }
    for index, region in enumerate(COMPARISON_ORDER, start=1):
        add_significance_line(
            axis, 0, index, y_max + index * line_step, line_height,
            significance_label(p_by_region.get(region, np.nan)),
        )

    axis.set_ylim(
        y_min - y_range * 0.1 if y_range > 0 else y_min - 0.1,
        y_max + (len(COMPARISON_ORDER) + 1) * line_step,
    )
    # Figure title and x-axis title intentionally omitted; use the thesis caption instead.
    plt.ylabel(y_label, fontsize=12, fontweight='bold')
    plt.xlabel('')
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, filename), dpi=300, bbox_inches='tight')
    plt.close()


print('\nFitting mixed-effects models and generating figures...')
model_records = []
contrast_tables = []
for metric_column, label, image_filename in [
    ('Fractal_Dimension', 'Fractal Dimension', '3_Fractal_Dimension_ViolinBoxPlot.png'),
    ('Lacunarity', 'Lacunarity', '3_Lacunarity_ViolinBoxPlot.png'),
]:
    result, model_data = fit_mixed_model(df_clean, metric_column)
    contrasts = normal_vs_region_contrasts(result)
    contrasts.insert(0, 'Metric', metric_column)
    contrast_tables.append(contrasts)

    model_records.append({
        'Metric': metric_column,
        'Model': 'Mixed-effects: metric ~ region + (1 | Sample_ID)',
        'ROI_Masks': len(model_data),
        'Independent_Sections': model_data['Sample_ID'].nunique(),
        'Converged': bool(result.converged),
        'Log_Likelihood': result.llf,
        'Random_Intercept_Variance': float(np.asarray(result.cov_re)[0, 0]),
    })
    plot_metric(model_data, metric_column, label, image_filename, contrasts)

pd.DataFrame(model_records).to_csv(os.path.join(OUT_DIR, '5_FracLac_MixedEffects_Results.csv'), index=False)
pd.concat(contrast_tables, ignore_index=True).to_csv(
    os.path.join(OUT_DIR, '5_FracLac_Posthoc_Holm_Contrasts.csv'), index=False
)

descriptive_statistics = (
    df_clean.groupby(['Region'])[ ['Fractal_Dimension', 'Lacunarity', 'R2_Fit', 'Foreground_Fraction']]
    .agg(['count', 'mean', 'std', 'median', 'min', 'max'])
)
descriptive_statistics.to_csv(os.path.join(OUT_DIR, '1_FracLac_Descriptive_Statistics.csv'))

print(f'\nAnalysis complete. Results saved in: {OUT_DIR}')
