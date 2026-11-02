"""Figure: H&E vs pentachrome fractal analysis.

Top:    representative image + collagen mask for normal skin and wound,
        H&E block (left) and pentachrome block (right).
Bottom: per-ROI scatter plots (pentachrome on x, H&E on y) for
        fractal dimension and lacunarity, with the line of identity.

Run compare_masks.py first; this script reads B1_FracLac_Paired_Table.csv.
"""

import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from PIL import Image

from config import CHECK_DIR, HE_MASK_DIR, OUT_DIR, PENTA_MASK_DIR, PIXEL_SIZE_UM
from naming import parse_file_metadata

Image.MAX_IMAGE_PIXELS = None

# ===================== Settings =====================
B1_CSV = os.path.join(OUT_DIR, 'B1_FracLac_Paired_Table.csv')

# Which ROIs to show (values from the Mask_ID column of B1_FracLac_Paired_Table.csv)
NORMAL_ID = 'Normal_06_ROI_1'
WOUND_ID = 'Wound_05_2_Wound_Centre'

SCALE_BAR_UM = 100
OUT_NAME = 'Figure_HE_vs_Pentachrome_Fractal'

COLOUR_NORMAL = '#0072B2'
COLOUR_WOUND = '#D55E00'
# ====================================================

plt.rcParams.update({'font.family': 'Arial', 'font.size': 8,
                     'axes.linewidth': 0.8, 'savefig.facecolor': 'white'})


def find_file(folder, mask_id, endings):
    hits = []
    for path in glob.glob(os.path.join(folder, '**', '*'), recursive=True):
        name = os.path.basename(path)
        if not name.lower().endswith(endings):
            continue
        try:
            if parse_file_metadata(name, 'NORMAL' in name.upper())[2] == mask_id:
                hits.append(path)
        except (ValueError, KeyError):
            continue
    if len(hits) != 1:
        raise SystemExit(f'Expected one file for {mask_id} ending {endings} in\n  {folder}\n'
                         f'found {len(hits)}: {hits}')
    print(f'  {mask_id}: {os.path.basename(hits[0])}')
    return hits[0]


def load_rgb(path):
    with Image.open(path) as im:
        return np.asarray(im.convert('RGB'))


def load_mask(path):
    with Image.open(path) as im:
        return np.asarray(im.convert('L')) > 127


def add_scale_bar(ax, width_px):
    length = SCALE_BAR_UM / PIXEL_SIZE_UM
    x1 = width_px * 0.95
    y = width_px * 0.94
    ax.plot([x1 - length, x1], [y, y], color='black', lw=2.5, solid_capstyle='butt')


def label_mask(ax, row):
    if row is None:
        return
    text = (f"D = {row['fd']:.3f}\n"
            f"Λ = {row['lac']:.2f}\n"
            f"FG = {row['fg']:.0%}")
    ax.text(0.03, 0.03, text, transform=ax.transAxes, color='white', fontsize=6.5,
            va='bottom', ha='left', linespacing=1.2,
            bbox=dict(facecolor='black', alpha=0.6, edgecolor='none', pad=2))


def metric_values(df, mask_id, stain):
    r = df[df['Mask_ID'] == mask_id]
    if r.empty:
        print(f'  WARNING: {mask_id} not in B1 table; no values printed on mask.')
        return None
    r = r.iloc[0]
    return {'fd': r[f'Fractal_Dimension_{stain}'],
            'lac': r[f'Lacunarity_{stain}'],
            'fg': r[f'Foreground_Fraction_{stain}']}


def scatter_panel(ax, df, metric, axis_label, show_legend):
    x, y = df[f'{metric}_Penta'], df[f'{metric}_HE']
    normal = df['Region'] == 'Normal Skin'
    low = df['Low_FG_Flag_Any'].astype(str).str.lower() == 'true'

    lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
    pad = 0.06 * (hi - lo)
    lims = (lo - pad, hi + pad)
    ax.plot(lims, lims, ls='--', color='grey', lw=0.8, zorder=0)

    ax.scatter(x[normal], y[normal], s=28, color=COLOUR_NORMAL, edgecolor='black',
               lw=0.4, label='Normal skin', zorder=3)
    ax.scatter(x[~normal & ~low], y[~normal & ~low], s=28, facecolor='white',
               edgecolor=COLOUR_WOUND, lw=1.1, label='Wound', zorder=3)
    ax.scatter(x[~normal & low], y[~normal & low], s=30, marker='x',
               color=COLOUR_WOUND, lw=1.1, label='Wound (<1% foreground)', zorder=3)

    shown = df['Mask_ID'].isin([NORMAL_ID, WOUND_ID])
    ax.scatter(x[shown], y[shown], s=120, facecolor='none', edgecolor='black',
               lw=0.8, zorder=4, label='ROIs shown above')

    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_aspect('equal')
    ax.set_xlabel(f'Pentachrome {axis_label}')
    ax.set_ylabel(f'H&E {axis_label}')
    ax.spines[['top', 'right']].set_visible(False)
    if show_legend:
        ax.legend(frameon=False, fontsize=6.5, loc='lower right', handletextpad=0.3)


# ----------------------------------------------------------------------------
print('Finding files...')
files = {
    ('HE', 'normal', 'img'): find_file(CHECK_DIR, NORMAL_ID, ('_he.png',)),
    ('HE', 'wound', 'img'): find_file(CHECK_DIR, WOUND_ID, ('_he.png',)),
    ('Penta', 'normal', 'img'): find_file(CHECK_DIR, NORMAL_ID, ('_penta.png',)),
    ('Penta', 'wound', 'img'): find_file(CHECK_DIR, WOUND_ID, ('_penta.png',)),
    ('HE', 'normal', 'mask'): find_file(HE_MASK_DIR, NORMAL_ID, ('.tif', '.tiff', '.png')),
    ('HE', 'wound', 'mask'): find_file(HE_MASK_DIR, WOUND_ID, ('.tif', '.tiff', '.png')),
    ('Penta', 'normal', 'mask'): find_file(PENTA_MASK_DIR, NORMAL_ID, ('.tif', '.tiff', '.png')),
    ('Penta', 'wound', 'mask'): find_file(PENTA_MASK_DIR, WOUND_ID, ('.tif', '.tiff', '.png')),
}
df = pd.read_csv(B1_CSV)

fig = plt.figure(figsize=(7.2, 8.0))

# Top: image/mask grid. Columns: H&E image, H&E mask, gap, Penta image, Penta mask
top = GridSpec(2, 5, figure=fig, width_ratios=[1, 1, 0.12, 1, 1],
               left=0.07, right=0.98, top=0.93, bottom=0.47, wspace=0.04, hspace=0.06)
columns = [('HE', 'img', 0), ('HE', 'mask', 1), ('Penta', 'img', 3), ('Penta', 'mask', 4)]
rows = [('normal', NORMAL_ID, 'Normal skin'), ('wound', WOUND_ID, 'Wound')]

axes_top = {}
for r, (tissue, mask_id, row_label) in enumerate(rows):
    for stain, kind, c in columns:
        ax = fig.add_subplot(top[r, c])
        axes_top[(r, c)] = ax
        path = files[(stain, tissue, kind)]
        if kind == 'img':
            image = load_rgb(path)
            ax.imshow(image)
            add_scale_bar(ax, image.shape[1])
        else:
            ax.imshow(load_mask(path), cmap='gray', vmin=0, vmax=1, interpolation='nearest')
            label_mask(ax, metric_values(df, mask_id, stain))
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)
        if r == 0:
            ax.set_title('Image' if kind == 'img' else 'Collagen mask', fontsize=8, pad=3)
        if c == 0:
            ax.text(-0.08, 0.5, row_label, transform=ax.transAxes, rotation=90,
                    va='center', ha='center', fontsize=9, fontweight='bold')

# Block titles and panel letters
for (c0, c1), title, letter in [((0, 1), 'H&E', 'A'), ((3, 4), 'Pentachrome', 'B')]:
    p0, p1 = axes_top[(0, c0)].get_position(), axes_top[(0, c1)].get_position()
    fig.text((p0.x0 + p1.x1) / 2, p0.y1 + 0.03, title, ha='center', va='bottom',
             fontsize=10, fontweight='bold')
    fig.text(p0.x0 - (0.05 if letter == 'A' else 0.02), p0.y1 + 0.03, letter,
             fontsize=12, fontweight='bold', va='bottom')

# Bottom: scatter plots
bottom = GridSpec(1, 2, figure=fig, left=0.09, right=0.97, top=0.39, bottom=0.06, wspace=0.35)
ax_fd = fig.add_subplot(bottom[0, 0])
ax_lac = fig.add_subplot(bottom[0, 1])
scatter_panel(ax_fd, df, 'Fractal_Dimension', 'fractal dimension', show_legend=True)
scatter_panel(ax_lac, df, 'Lacunarity', 'lacunarity', show_legend=False)
for ax, letter in [(ax_fd, 'C'), (ax_lac, 'D')]:
    ax.text(-0.22, 1.05, letter, transform=ax.transAxes, fontsize=12,
            fontweight='bold', va='bottom')

os.makedirs(OUT_DIR, exist_ok=True)
for ext in ('png', 'tiff', 'pdf'):
    out = os.path.join(OUT_DIR, f'{OUT_NAME}.{ext}')
    kwargs = {'pil_kwargs': {'compression': 'tiff_lzw'}} if ext == 'tiff' else {}
    fig.savefig(out, dpi=600, **kwargs)
    print(f'Saved: {out}')
plt.close(fig)
