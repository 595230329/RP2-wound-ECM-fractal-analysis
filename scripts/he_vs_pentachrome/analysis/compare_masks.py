"""Compare H&E and pentachrome collagen-selective binary masks.

Expects these folders in DATA_DIR (see config.py):
    Penta                  binary masks exported for the original pentachrome analysis
    H&E                    binary masks from the registered H&E images (aligned export)
    NORMAL SKIN BOXCOUNT   pentachrome FracLac summaries (optional, Part B)
    WOUND SKIN BOXCOUNT    pentachrome FracLac summaries (optional, Part B)
    HE BOXCOUNT            H&E FracLac summaries (optional, Part B)

Mask file names must follow the naming scheme in naming.py (e.g. contain '9M05',
'NORMAL SKIN 1', 'WOUND CENTER', 'PROXIMAL EDGE'), so that each H&E mask can be
paired with its pentachrome mask by Mask_ID.

Part A: pixel-level agreement between paired masks (quality control).
    exact overlap (Dice, Jaccard, sensitivity, precision)
    registration-tolerant overlap (tolerance in micrometres)
    chance-level Dice baseline (same collagen amount, no spatial match)
    Dice before and after residual-shift correction
    (after a shift, both masks are cropped to the area with real data)

Part B: ROI-level agreement of FracLac outputs (FD, lacunarity, foreground fraction).

With few ROIs, results are reported descriptively; mixed models are only fitted
when there are enough paired masks and samples.
"""

import glob
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.formula.api as smf
from PIL import Image
from scipy import ndimage
from scipy.stats import norm, wilcoxon
from statsmodels.stats.multitest import multipletests

from config import (HE_BOX_DIRS, HE_MASK_DIR, OUT_DIR, PENTA_BOX_DIRS,
                    PENTA_MASK_DIR, PIXEL_SIZE_UM)
from naming import index_mask_folder, parse_file_metadata

warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)
Image.MAX_IMAGE_PIXELS = None

# ----------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------
OVERLAY_DIR = os.path.join(OUT_DIR, 'Overlays')
os.makedirs(OVERLAY_DIR, exist_ok=True)

INVERT_PENTA = False   # set True if a mask set has black collagen on white background
INVERT_HE = False

TOLERANCE_RADII_UM = [2, 5, 10]   # tolerance for registration-tolerant overlap (µm)

ESTIMATE_RESIDUAL_SHIFT = True    # estimate residual translation between masks
APPLY_SHIFT_CORRECTION = True     # correct small translations before scoring
MAX_SHIFT_UM = 25                 # larger estimated shifts are treated as unreliable
MIN_FG_FOR_SHIFT = 0.005          # skip shift estimation on very sparse masks
N_CHANCE_SHIFTS = 20              # random shifts used for the chance-level Dice

MIN_PAIRS_FOR_MODEL = 12          # minimum paired masks before fitting mixed models
MIN_SAMPLES_FOR_MODEL = 3         # minimum distinct Sample_IDs for mixed models

REGION_ORDER = ['Normal Skin', 'Wound Centre', 'Distant Edge', 'Proximal Edge']
COMPARISON_ORDER = ['Wound Centre', 'Distant Edge', 'Proximal Edge']
R2_QC_THRESHOLD = 0.95
LOW_FOREGROUND_FRACTION_THRESHOLD = 0.01

MAX_SHIFT_PX = max(1, int(round(MAX_SHIFT_UM / PIXEL_SIZE_UM)))
TOLERANT_F_COLUMNS = [f'Tolerant_F_{r}um' for r in TOLERANCE_RADII_UM]


# ----------------------------------------------------------------------------
# Part A helpers: mask loading, alignment and overlap metrics
# ----------------------------------------------------------------------------
def load_binary_mask(path, invert=False):
    """Load a mask as boolean, with white (high) pixels as foreground."""
    with Image.open(path) as image:
        if image.mode in ('1', 'P', 'L', 'LA', 'RGB', 'RGBA', 'CMYK'):
            array = np.asarray(image.convert('L')).astype(float)
        else:  # 16-/32-bit images
            array = np.asarray(image).astype(float)
    if array.ndim == 3:
        array = array[..., :3].max(axis=2)
    high = array.max()
    mask = array > high / 2 if high > 0 else np.zeros(array.shape, dtype=bool)
    return ~mask if invert else mask


def resize_mask(mask, shape):
    image = Image.fromarray(mask.astype(np.uint8) * 255)
    resized = image.resize((shape[1], shape[0]), resample=Image.NEAREST)
    return np.asarray(resized) > 127


def estimate_translation(reference, moving):
    """Phase correlation; returns (dy, dx) to apply to `moving` to match `reference`."""
    a = reference.astype(float) - reference.mean()
    b = moving.astype(float) - moving.mean()
    cross = np.fft.fft2(a) * np.conj(np.fft.fft2(b))
    cross /= np.abs(cross) + 1e-12
    correlation = np.fft.ifft2(cross).real
    peak = np.unravel_index(np.argmax(correlation), correlation.shape)
    return tuple(int(p if p <= s // 2 else p - s) for p, s in zip(peak, correlation.shape))


def shift_mask(mask, dy, dx):
    """Integer translation with zero fill (no wrap-around)."""
    height, width = mask.shape
    shifted = np.zeros_like(mask)
    src_y = slice(max(0, -dy), min(height, height - dy))
    dst_y = slice(max(0, dy), min(height, height + dy))
    src_x = slice(max(0, -dx), min(width, width - dx))
    dst_x = slice(max(0, dx), min(width, width + dx))
    shifted[dst_y, dst_x] = mask[src_y, src_x]
    return shifted


def crop_to_valid(penta, he, dy, dx):
    """After shifting H&E by (dy, dx), keep only the area where both masks have real data.

    Without this, the zero-filled strip left by the shift is scored as
    'pentachrome collagen missed by H&E', which lowers Dice and sensitivity.
    """
    height, width = penta.shape
    ys = slice(max(0, dy), min(height, height + dy))
    xs = slice(max(0, dx), min(width, width + dx))
    return penta[ys, xs], he[ys, xs]


def disk(radius):
    y, x = np.ogrid[-radius:radius + 1, -radius:radius + 1]
    return x ** 2 + y ** 2 <= radius ** 2


def safe_divide(numerator, denominator):
    return float(numerator) / float(denominator) if denominator > 0 else np.nan


def dice(a, b):
    tp = np.logical_and(a, b).sum()
    return safe_divide(2 * tp, a.sum() + b.sum())


def overlap_metrics(penta, he):
    """Pentachrome mask is the reference; H&E mask is the test."""
    tp = np.logical_and(penta, he).sum()
    fp = np.logical_and(~penta, he).sum()
    fn = np.logical_and(penta, ~he).sum()
    tn = penta.size - tp - fp - fn
    metrics = {
        'TP_Pixels': int(tp), 'FP_Pixels': int(fp), 'FN_Pixels': int(fn), 'TN_Pixels': int(tn),
        'Penta_Foreground_Fraction': penta.mean(),
        'HE_Foreground_Fraction': he.mean(),
        'Foreground_Fraction_Difference_HE_minus_Penta': he.mean() - penta.mean(),
        'Foreground_Ratio_HE_over_Penta': safe_divide(he.sum(), penta.sum()),
        'Dice': safe_divide(2 * tp, 2 * tp + fp + fn),
        'Jaccard': safe_divide(tp, tp + fp + fn),
        'Sensitivity': safe_divide(tp, tp + fn),   # pentachrome collagen captured by H&E
        'Precision': safe_divide(tp, tp + fp),     # H&E collagen that is pentachrome collagen
        'Specificity': safe_divide(tn, tn + fp),
    }
    for radius_um in TOLERANCE_RADII_UM:
        radius = max(1, int(round(radius_um / PIXEL_SIZE_UM)))
        element = disk(radius)
        penta_dilated = ndimage.binary_dilation(penta, structure=element)
        he_dilated = ndimage.binary_dilation(he, structure=element)
        tol_precision = safe_divide(np.logical_and(he, penta_dilated).sum(), he.sum())
        tol_recall = safe_divide(np.logical_and(penta, he_dilated).sum(), penta.sum())
        tol_f = (2 * tol_precision * tol_recall / (tol_precision + tol_recall)
                 if (tol_precision + tol_recall) > 0 else np.nan)
        metrics[f'Tolerant_Precision_{radius_um}um'] = tol_precision
        metrics[f'Tolerant_Recall_{radius_um}um'] = tol_recall
        metrics[f'Tolerant_F_{radius_um}um'] = tol_f
    return metrics


def chance_dice(penta, he, n=N_CHANCE_SHIFTS, seed=0):
    """Mean Dice after large random circular shifts of the H&E mask.

    Keeps the amount of foreground identical but destroys spatial correspondence,
    so it shows how much Dice you get from collagen amount alone.
    """
    if penta.sum() + he.sum() == 0:
        return np.nan
    rng = np.random.default_rng(seed)
    height, width = penta.shape
    values = []
    for _ in range(n):
        dy = int(rng.integers(height // 4, 3 * height // 4 + 1))
        dx = int(rng.integers(width // 4, 3 * width // 4 + 1))
        values.append(dice(penta, np.roll(he, (dy, dx), axis=(0, 1))))
    return float(np.nanmean(values))


def save_overlay(penta, he, path):
    """White = both, magenta = pentachrome only, green = H&E only, black = neither."""
    overlay = np.zeros(penta.shape + (3,), dtype=np.uint8)
    overlay[np.logical_and(penta, he)] = (255, 255, 255)
    overlay[np.logical_and(penta, ~he)] = (255, 0, 255)
    overlay[np.logical_and(~penta, he)] = (0, 200, 0)
    Image.fromarray(overlay).save(path)


# ----------------------------------------------------------------------------
# Part B helpers: FracLac reading and agreement statistics
# ----------------------------------------------------------------------------
def numeric_value(row, index, source_file):
    value = pd.to_numeric(row.iloc[index], errors='coerce')
    if pd.isna(value):
        raise ValueError(f'{source_file}: FracLac value at column {index + 1} is not numeric.')
    return float(value)


def read_fraclac_table(path):
    if path.lower().endswith('.txt'):
        for encoding in ('utf-8-sig', 'cp1252', 'latin-1'):
            try:
                return pd.read_csv(path, sep='\t', header=None, dtype=str,
                                   encoding=encoding, quoting=3)
            except UnicodeDecodeError:
                continue
    return pd.read_excel(path, header=None)


def read_fraclac_file(path):
    """Read every ROI row in a FracLac summary (.txt or .xlsx).

    Columns (0-based): 5 fractal dimension, 10 R², 86 lacunarity,
    117 foreground pixels, 119 total pixels.
    """
    source_file = os.path.basename(path)
    table = read_fraclac_table(path)
    if table.shape[0] < 2 or table.shape[1] < 120:
        raise ValueError('unexpected FracLac summary layout')
    n_rows = table.shape[0] - 1
    records = []
    for row_number in range(1, table.shape[0]):
        values = table.iloc[row_number]
        label = '' if pd.isna(values.iloc[0]) else str(values.iloc[0]).strip()
        if not label:
            continue
        candidates = [path, label] if n_rows == 1 else [label]
        for candidate in candidates:
            try:
                sample_id, region, mask_id = parse_file_metadata(candidate)
                break
            except ValueError:
                continue
        else:
            print(f'  Skipped row {row_number} in {source_file}: cannot identify ROI '
                  f'from "{label[:60]}"')
            continue
        try:
            foreground = numeric_value(values, 117, source_file)
            total = numeric_value(values, 119, source_file)
            fractal_dimension = numeric_value(values, 5, source_file)
            lacunarity = numeric_value(values, 86, source_file)
        except ValueError as error:
            print(f'  Skipped {mask_id}: {error}')
            continue
        r2 = pd.to_numeric(values.iloc[10], errors='coerce')
        if pd.isna(r2):
            print(f'  NOTE {mask_id}: R2 is "{values.iloc[10]}"; kept with R2 = NaN')
        records.append({
            'Mask_ID': mask_id, 'Sample_ID': sample_id, 'Region': region,
            'Source_File': source_file, 'Image_Label': label[:150],
            'Fractal_Dimension': fractal_dimension,
            'R2_Fit': float(r2),
            'Lacunarity': lacunarity,
            'Foreground_Fraction': foreground / total,
        })
    return records


def read_fraclac_set(folders):
    records = []
    for folder in folders:
        if not os.path.isdir(folder):
            continue
        for path in sorted(glob.glob(os.path.join(folder, '**', '*'), recursive=True)):
            base = os.path.basename(path)
            lower = base.lower()
            if base.startswith('~$'):
                continue
            is_summary_txt = lower.endswith('.txt') and 'summary' in lower
            if not (lower.endswith(('.xlsx', '.xls')) or is_summary_txt):
                continue
            try:
                records.extend(read_fraclac_file(path))
            except Exception as error:
                print(f'  Skipped {base}: {error}')
    table = pd.DataFrame(records)
    if not table.empty and table['Mask_ID'].duplicated().any():
        duplicates = table.loc[table['Mask_ID'].duplicated(keep=False), 'Mask_ID'].unique()
        print(f'  WARNING: duplicate FracLac results for {list(duplicates)}; keeping the last.')
        table = table.drop_duplicates('Mask_ID', keep='last')
    return table


def lins_ccc(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    covariance = np.mean((x - x.mean()) * (y - y.mean()))
    return 2 * covariance / (x.var() + y.var() + (x.mean() - y.mean()) ** 2)


def icc_2_1(x, y):
    """Two-way random effects, absolute agreement, single measurement (Shrout & Fleiss)."""
    data = np.column_stack([x, y]).astype(float)
    n, k = data.shape
    grand = data.mean()
    ss_rows = k * ((data.mean(axis=1) - grand) ** 2).sum()
    ss_cols = n * ((data.mean(axis=0) - grand) ** 2).sum()
    ss_error = ((data - grand) ** 2).sum() - ss_rows - ss_cols
    ms_rows = ss_rows / (n - 1)
    ms_cols = ss_cols / (k - 1)
    ms_error = ss_error / ((n - 1) * (k - 1))
    return (ms_rows - ms_error) / (ms_rows + (k - 1) * ms_error + k * (ms_cols - ms_error) / n)


def agreement_summary(penta, he, metric, subset_label):
    penta, he = np.asarray(penta, float), np.asarray(he, float)
    difference = he - penta
    bias, sd = difference.mean(), difference.std(ddof=1)
    try:
        p_wilcoxon = wilcoxon(he, penta).pvalue
    except ValueError:
        p_wilcoxon = np.nan
    return {
        'Metric': metric, 'Subset': subset_label, 'N_Pairs': len(difference),
        'Mean_Pentachrome': penta.mean(), 'Mean_HE': he.mean(),
        'Bias_HE_minus_Penta': bias, 'SD_Difference': sd,
        'LoA_Lower': bias - 1.96 * sd, 'LoA_Upper': bias + 1.96 * sd,
        'Pearson_r': np.corrcoef(penta, he)[0, 1],
        'Lins_CCC': lins_ccc(penta, he), 'ICC_2_1': icc_2_1(penta, he),
        'Wilcoxon_P': p_wilcoxon,
        'Small_N_Warning': len(difference) < 10,
    }


# ----------------------------------------------------------------------------
# Shared modelling and plotting helpers
# ----------------------------------------------------------------------------
def enough_for_model(data):
    return (len(data) >= MIN_PAIRS_FOR_MODEL
            and data['Sample_ID'].nunique() >= MIN_SAMPLES_FOR_MODEL
            and data['Region'].nunique() >= 2)


def fit_mixed_model(data, response):
    model_data = data.dropna(subset=[response, 'Region', 'Sample_ID']).copy()
    model_data['Region'] = model_data['Region'].astype(str)
    formula = f"{response} ~ C(Region, Treatment(reference='Normal Skin'))"
    model = smf.mixedlm(formula, data=model_data, groups=model_data['Sample_ID'])
    try:
        result = model.fit(reml=False, method='lbfgs', maxiter=500, disp=False)
    except Exception:
        result = model.fit(reml=False, method='powell', maxiter=1000, disp=False)
    return result


def normal_vs_region_contrasts(result, metric):
    names = result.fe_params.index.tolist()
    covariance = result.cov_params().loc[names, names]
    rows = [{
        'Metric': metric, 'Comparison': 'Intercept (Normal Skin)',
        'Estimate': float(result.fe_params['Intercept']),
        'SE': float(np.sqrt(covariance.loc['Intercept', 'Intercept'])),
    }]
    for region in COMPARISON_ORDER:
        terms = [n for n in names if n.endswith(f'[T.{region}]')]
        if not terms:
            continue
        rows.append({
            'Metric': metric, 'Comparison': f'{region} vs Normal Skin',
            'Estimate': float(result.fe_params[terms[0]]),
            'SE': float(np.sqrt(covariance.loc[terms[0], terms[0]])),
        })
    table = pd.DataFrame(rows)
    table['Z'] = table['Estimate'] / table['SE']
    table['P_Unadjusted'] = 2 * norm.sf(table['Z'].abs())
    table['P_Holm'] = np.nan
    is_contrast = table['Comparison'] != 'Intercept (Normal Skin)'
    if is_contrast.any():
        table.loc[is_contrast, 'P_Holm'] = multipletests(
            table.loc[is_contrast, 'P_Unadjusted'], method='holm')[1]
    table['Converged'] = bool(result.converged)
    return table


def plot_by_region(data, column, label, filename, chance_column=None):
    present = [r for r in REGION_ORDER if r in set(data['Region'].astype(str))]
    plt.figure(figsize=(7, 5))
    sns.boxplot(x='Region', y=column, data=data, order=present,
                color='white', width=0.4, showfliers=False)
    sns.stripplot(x='Region', y=column, data=data, order=present,
                  color='black', alpha=0.8, jitter=0.12, size=6)
    if chance_column is not None:
        sns.stripplot(x='Region', y=chance_column, data=data, order=present,
                      color='red', marker='_', size=18, linewidth=2, jitter=0.12)
        plt.title('Black = observed; red = chance level', fontsize=10)
    plt.ylabel(label, fontsize=12, fontweight='bold')
    plt.xlabel('')
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, filename), dpi=300, bbox_inches='tight')
    plt.close()


def plot_foreground_pairs(data, filename):
    """Paired pentachrome vs H&E foreground fraction per mask."""
    palette = dict(zip(REGION_ORDER, sns.color_palette('colorblind', len(REGION_ORDER))))
    plt.figure(figsize=(6, 5))
    for _, row in data.iterrows():
        colour = palette[str(row['Region'])]
        plt.plot([0, 1], [row['Penta_Foreground_Fraction'], row['HE_Foreground_Fraction']],
                 color=colour, marker='o', alpha=0.8)
    for region in REGION_ORDER:
        if region in set(data['Region'].astype(str)):
            plt.plot([], [], color=palette[region], marker='o', label=region)
    plt.xticks([0, 1], ['Pentachrome', 'H&E'])
    plt.ylabel('Collagen foreground fraction', fontweight='bold')
    plt.legend(fontsize=8, frameon=False)
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, filename), dpi=300, bbox_inches='tight')
    plt.close()


def plot_bland_altman_and_identity(paired, metric, label, filename):
    penta, he = paired[f'{metric}_Penta'], paired[f'{metric}_HE']
    mean, difference = (penta + he) / 2, he - penta
    bias, sd = difference.mean(), difference.std(ddof=1)
    flagged = paired['Low_FG_Flag_Any']
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    palette = dict(zip(REGION_ORDER, sns.color_palette('colorblind', len(REGION_ORDER))))
    for region in REGION_ORDER:
        for is_flagged, marker in [(False, 'o'), (True, 'X')]:
            sel = (paired['Region'].astype(str) == region) & (flagged == is_flagged)
            if not sel.any():
                continue
            legend_label = region + (' (<1% foreground)' if is_flagged else '')
            axes[0].scatter(mean[sel], difference[sel], color=palette[region],
                            marker=marker, s=45, label=legend_label)
            axes[1].scatter(penta[sel], he[sel], color=palette[region], marker=marker, s=45)
    axes[0].axhline(bias, color='black', lw=1.2)
    if len(difference) >= 3:
        for limit in (bias - 1.96 * sd, bias + 1.96 * sd):
            axes[0].axhline(limit, color='black', lw=1, ls='--')
    axes[0].set_xlabel(f'Mean of H&E and pentachrome {label}', fontweight='bold')
    axes[0].set_ylabel(f'H&E − pentachrome {label}', fontweight='bold')
    axes[0].legend(fontsize=8, frameon=False)
    low = min(penta.min(), he.min())
    high = max(penta.max(), he.max())
    axes[1].plot([low, high], [low, high], color='grey', ls='--', lw=1)
    axes[1].set_xlabel(f'Pentachrome {label}', fontweight='bold')
    axes[1].set_ylabel(f'H&E {label}', fontweight='bold')
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, filename), dpi=300, bbox_inches='tight')
    plt.close()


# ============================================================================
# PART A: pixel-level mask agreement
# ============================================================================
if os.path.isdir(PENTA_MASK_DIR) and os.path.isdir(HE_MASK_DIR):
    print('PART A: pixel-level mask comparison')
    print(f'  Pixel size: {PIXEL_SIZE_UM} µm; tolerance radii: {TOLERANCE_RADII_UM} µm; '
          f'max reliable shift: {MAX_SHIFT_PX} px')
    penta_index = index_mask_folder(PENTA_MASK_DIR)
    he_index = index_mask_folder(HE_MASK_DIR)
    shared_ids = sorted(set(penta_index) & set(he_index))
    print(f'  Pentachrome masks: {len(penta_index)}; H&E masks: {len(he_index)}; '
          f'paired: {len(shared_ids)}')
    for missing in sorted(set(he_index) - set(penta_index)):
        print(f'  H&E mask without pentachrome partner: {missing}')

    pixel_records = []
    for mask_id in shared_ids:
        penta_path, sample_id, region = penta_index[mask_id]
        he_path = he_index[mask_id][0]
        penta = load_binary_mask(penta_path, INVERT_PENTA)
        he = load_binary_mask(he_path, INVERT_HE)

        record = {'Mask_ID': mask_id, 'Sample_ID': sample_id, 'Region': region,
                  'Penta_File': os.path.basename(penta_path),
                  'HE_File': os.path.basename(he_path),
                  'Penta_Shape': str(penta.shape), 'HE_Shape_Original': str(he.shape),
                  'HE_Resized': penta.shape != he.shape,
                  'HE_Scale_Factor': he.shape[1] / penta.shape[1]}

        if penta.shape != he.shape:
            ratio_penta = penta.shape[1] / penta.shape[0]
            ratio_he = he.shape[1] / he.shape[0]
            if abs(ratio_penta - ratio_he) / ratio_penta > 0.02:
                print(f'  WARNING {mask_id}: aspect ratios differ; check ROI export.')
            print(f'  NOTE {mask_id}: H&E mask resized from {he.shape} to {penta.shape}. '
                  f'Export both at the same µm/pixel before running FracLac.')
            he = resize_mask(he, penta.shape)

        for name, mask in [('Penta', penta), ('HE', he)]:
            if mask.mean() > 0.6:
                print(f'  WARNING {mask_id}: {name} foreground is {mask.mean():.0%}; '
                      f'check mask polarity (INVERT_{name.upper()}).')

        # Agreement before any shift correction (full mask)
        record['Dice_Before_Shift'] = dice(penta, he)
        record.update({'Shift_dy_px': np.nan, 'Shift_dx_px': np.nan,
                       'Shift_um': np.nan, 'Shift_Reliable': False, 'Shift_Applied': False})

        if ESTIMATE_RESIDUAL_SHIFT and min(penta.mean(), he.mean()) >= MIN_FG_FOR_SHIFT:
            dy, dx = estimate_translation(penta, he)
            reliable = max(abs(dy), abs(dx)) <= MAX_SHIFT_PX
            record.update({'Shift_dy_px': dy, 'Shift_dx_px': dx,
                           'Shift_um': float(np.hypot(dy, dx) * PIXEL_SIZE_UM),
                           'Shift_Reliable': reliable})
            if not reliable:
                print(f'  WARNING {mask_id}: estimated shift ({dy}, {dx}) px exceeds '
                      f'{MAX_SHIFT_UM} µm; not applied. Check registration.')
            if APPLY_SHIFT_CORRECTION and reliable and (dy or dx):
                he = shift_mask(he, dy, dx)
                # drop the zero-filled strip so it is not scored as missed collagen
                penta, he = crop_to_valid(penta, he, dy, dx)
                record['Shift_Applied'] = True

        record['Scored_Shape'] = str(penta.shape)
        record.update(overlap_metrics(penta, he))
        record['Dice_Chance'] = chance_dice(penta, he)
        record['Dice_Above_Chance'] = record['Dice'] - record['Dice_Chance']
        record['Low_FG_Flag_Any'] = (min(penta.mean(), he.mean())
                                     < LOW_FOREGROUND_FRACTION_THRESHOLD)
        pixel_records.append(record)

        save_overlay(penta, he, os.path.join(OVERLAY_DIR, f'{mask_id}_overlay.png'))
        print(f"  {mask_id}: shift=({record['Shift_dy_px']}, {record['Shift_dx_px']}) px | "
              f"scored {record['Scored_Shape']} | "
              f"Dice {record['Dice_Before_Shift']:.3f} -> {record['Dice']:.3f} "
              f"(chance {record['Dice_Chance']:.3f}) | "
              f"Tolerant F {TOLERANCE_RADII_UM[-1]} µm = {record[TOLERANT_F_COLUMNS[-1]]:.3f} | "
              f"FG penta {record['Penta_Foreground_Fraction']:.3f}, "
              f"H&E {record['HE_Foreground_Fraction']:.3f}")

    if pixel_records:
        pixel_df = pd.DataFrame(pixel_records)
        pixel_df['Region'] = pd.Categorical(pixel_df['Region'], REGION_ORDER, ordered=True)
        pixel_df = pixel_df.sort_values(['Region', 'Sample_ID', 'Mask_ID'])
        pixel_df.to_csv(os.path.join(OUT_DIR, 'A1_Pixel_Agreement_Per_Mask.csv'), index=False)

        summary_columns = (['Dice_Before_Shift', 'Dice', 'Dice_Chance', 'Dice_Above_Chance',
                            'Jaccard', 'Sensitivity', 'Precision',
                            'Penta_Foreground_Fraction', 'HE_Foreground_Fraction',
                            'Foreground_Ratio_HE_over_Penta'] + TOLERANT_F_COLUMNS)
        (pixel_df.groupby('Region', observed=True)[summary_columns]
         .agg(['count', 'mean', 'median', 'min', 'max'])
         .to_csv(os.path.join(OUT_DIR, 'A2_Pixel_Agreement_By_Region.csv')))

        plot_by_region(pixel_df, 'Dice', 'Dice coefficient',
                       'A3_Dice_By_Region.png', chance_column='Dice_Chance')
        plot_by_region(pixel_df, TOLERANT_F_COLUMNS[-1],
                       f'Tolerant F ({TOLERANCE_RADII_UM[-1]} µm)',
                       f'A3_{TOLERANT_F_COLUMNS[-1]}_By_Region.png')
        plot_foreground_pairs(pixel_df, 'A3_Foreground_Fraction_Paired.png')

        if enough_for_model(pixel_df):
            model_tables = []
            for column in ['Dice', TOLERANT_F_COLUMNS[-1]]:
                try:
                    model_tables.append(
                        normal_vs_region_contrasts(fit_mixed_model(pixel_df, column), column))
                except Exception as error:
                    print(f'  Mixed model for {column} failed: {error}')
            if model_tables:
                pd.concat(model_tables).to_csv(
                    os.path.join(OUT_DIR, 'A4_Pixel_Agreement_MixedModel_Contrasts.csv'),
                    index=False)
        else:
            print(f'  Mixed model skipped: {len(pixel_df)} paired masks from '
                  f'{pixel_df["Sample_ID"].nunique()} samples. Report Part A descriptively.')
else:
    print('PART A skipped: mask folders not found.')


# ============================================================================
# PART B: ROI-level agreement of FracLac outputs
# ============================================================================
print('\nFracLac folders:')
for folder in PENTA_BOX_DIRS + HE_BOX_DIRS:
    print(f'  {"FOUND  " if os.path.isdir(folder) else "MISSING"} {folder}')

penta_box_ok = any(os.path.isdir(folder) for folder in PENTA_BOX_DIRS)
he_box_ok = any(os.path.isdir(folder) for folder in HE_BOX_DIRS)

if penta_box_ok and he_box_ok:
    print('\nPART B: FracLac output comparison')
    penta_fl = read_fraclac_set(PENTA_BOX_DIRS)
    he_fl = read_fraclac_set(HE_BOX_DIRS)
    if penta_fl.empty or he_fl.empty:
        paired = pd.DataFrame()
    else:
        paired = penta_fl.merge(he_fl, on=['Mask_ID', 'Sample_ID', 'Region'],
                                suffixes=('_Penta', '_HE'), how='inner')
    print(f'  Pentachrome: {len(penta_fl)}; H&E: {len(he_fl)}; paired: {len(paired)}')

    if len(paired) == 0:
        print('  No paired FracLac results - check file naming.')
        print('  Pentachrome IDs:', sorted(penta_fl['Mask_ID']) if not penta_fl.empty else 'none')
        print('  H&E IDs:', sorted(he_fl['Mask_ID']) if not he_fl.empty else 'none')
    else:
        if not he_fl.empty:
            unpaired_he = sorted(set(he_fl['Mask_ID']) - set(paired['Mask_ID']))
            for missing in unpaired_he:
                print(f'  H&E FracLac result without pentachrome partner: {missing}')

        paired['Low_FG_Flag_Any'] = (
            paired[['Foreground_Fraction_Penta', 'Foreground_Fraction_HE']].min(axis=1)
            < LOW_FOREGROUND_FRACTION_THRESHOLD)
        # skipna=False: a missing R2 on either stain fails the QC instead of being ignored.
        paired['R2_QC_Both'] = (paired[['R2_Fit_Penta', 'R2_Fit_HE']].min(axis=1, skipna=False)
                                >= R2_QC_THRESHOLD)
        metrics = [('Fractal_Dimension', 'fractal dimension'),
                   ('Lacunarity', 'lacunarity'),
                   ('Foreground_Fraction', 'foreground fraction')]
        for metric, _ in metrics:
            paired[f'{metric}_Diff'] = paired[f'{metric}_HE'] - paired[f'{metric}_Penta']
        paired['Region'] = pd.Categorical(paired['Region'], REGION_ORDER, ordered=True)
        paired = paired.sort_values(['Region', 'Sample_ID', 'Mask_ID'])
        paired.to_csv(os.path.join(OUT_DIR, 'B1_FracLac_Paired_Table.csv'), index=False)

        subsets = {
            'All paired masks': paired,
            'R2 QC passed in both': paired[paired['R2_QC_Both']],
            'Excluding <1% foreground': paired[~paired['Low_FG_Flag_Any']],
        }
        agreement_rows = []
        for metric, _ in metrics:
            for subset_label, subset in subsets.items():
                if len(subset) >= 3:
                    agreement_rows.append(agreement_summary(
                        subset[f'{metric}_Penta'], subset[f'{metric}_HE'],
                        metric, subset_label))
            for region in REGION_ORDER:
                subset = paired[paired['Region'].astype(str) == region]
                if len(subset) >= 3:
                    agreement_rows.append(agreement_summary(
                        subset[f'{metric}_Penta'], subset[f'{metric}_HE'], metric,
                        f'Region: {region}'))
        if agreement_rows:
            pd.DataFrame(agreement_rows).to_csv(
                os.path.join(OUT_DIR, 'B2_FracLac_Agreement_Statistics.csv'), index=False)

        bias_tables = []
        for metric, label in metrics:
            plot_bland_altman_and_identity(paired, metric, label,
                                           f'B3_{metric}_BlandAltman_Identity.png')
            if enough_for_model(paired):
                try:
                    result = fit_mixed_model(paired, f'{metric}_Diff')
                    bias_tables.append(normal_vs_region_contrasts(result, f'{metric}_Diff'))
                except Exception as error:
                    print(f'  Bias model for {metric} failed: {error}')
        if bias_tables:
            pd.concat(bias_tables).to_csv(
                os.path.join(OUT_DIR, 'B4_Region_Dependent_Bias_MixedModel.csv'), index=False)
        elif not enough_for_model(paired):
            print('  Region-dependent bias model skipped: too few pairs. '
                  'Use B1 and the Bland-Altman plots descriptively.')
else:
    print('\nPART B skipped: one or more BOXCOUNT folders not found '
          '(run FracLac on the H&E masks first).')

print(f'\nDone. Results saved in: {OUT_DIR}')
