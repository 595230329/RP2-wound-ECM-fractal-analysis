"""Build the analysis mask folders (Penta/, H&E/) from the QuPath exports.

Each ROI was drawn on one stain (the reference) and propagated to the other, where
it may appear rotated. For each ROI this script reads its angle on both stains from
<QuPath masks>/roi_angles.csv (written by the two *_aligned*.groovy scripts) and picks:

  both stains axis-aligned  -> Penta/ and H&E/ from masks/Penta and masks/H&E
  pentachrome is reference  -> Penta/ from masks/Penta,         H&E/ from masks/H&E_aligned
  H&E is reference          -> Penta/ from masks/Penta_aligned, H&E/ from masks/H&E

An ROI is axis-aligned on a stain if its angle is within AXIS_ALIGNED_TOLERANCE_DEG of a
multiple of 90 degrees; otherwise the stain where it is closer to axis-aligned is the
reference. ROIs without a pentachrome angle (normal skin, which the pentachrome script
skips) are treated as drawn on pentachrome.

Every mask is written in one format (8-bit greyscale TIFF, collagen = 255, background = 0),
because QuPath saves masks from export_masks.groovy and from the *_aligned scripts in two
TIFF formats that FracLac reads with opposite polarity. Pixel values are not changed.

The check images for each ROI are copied to check_images/ from check_aligned/ (pentachrome
or both as reference) or check_aligned_penta/ (H&E as reference).

The same masks are also copied to "FracLac input/Penta" and "FracLac input/H&E": run FracLac
there, because it writes its own image files into the folder it reads.

Existing Penta/, H&E/, FracLac input and check folders are moved to
<name>_previous_<timestamp>/ (nothing is deleted). The choice for every ROI is written to outputs/.../mask_sources.csv.

Usage:
    py scripts/he_vs_pentachrome/analysis/assemble_masks.py
Set QUPATH_MASKS_DIR (environment variable or config.py) to the project's masks/ folder.
"""
import csv
import os
import shutil
import sys
import time

import numpy as np
from PIL import Image

from config import (AXIS_ALIGNED_TOLERANCE_DEG, CHECK_DIRS, DATA_DIR, FRACLAC_INPUT_DIRS,
                    HE_MASK_DIR, OUT_DIR, PENTA_MASK_DIR, QUPATH_MASKS_DIR)
from naming import index_mask_folder

# Orientation check (check images, averaged per image): flag if a 90-degree turn fits clearly better
ORIENTATION_MIN_R = 0.1
ORIENTATION_MARGIN = 0.1


def axis_deviation(angle):
    """Degrees away from the nearest multiple of 90 (0-45)."""
    a = abs(float(angle)) % 90.0
    return min(a, 90.0 - a)


def read_angles(path):
    """(stain, ROI name) -> angle and -> image name, keeping the last row for each
    (the log is appended)."""
    angles, images = {}, {}
    with open(path, newline='', encoding='utf-8') as handle:
        for row in csv.DictReader(handle):
            key = (row['Stain'], row['ROI'].strip())
            angles[key] = float(row['Angle_deg'])
            images[key] = row['Image'].strip()
    return angles, images


def roi_name(path):
    name = os.path.basename(path)
    for suffix in ('_mask.tif', '_mask.tiff'):
        if name.lower().endswith(suffix):
            return name[:-len(suffix)]
    return os.path.splitext(name)[0]


def write_standard_mask(src, dst):
    """Copy a mask as 8-bit greyscale, collagen = 255; return its collagen fraction."""
    with Image.open(src) as image:
        values = np.array(image)
    if values.ndim != 2:
        raise ValueError(f'{src}: expected a single-channel mask, got shape {values.shape}')
    unique = np.unique(values)
    if not set(unique.tolist()) <= {0, 255}:
        raise ValueError(f'{src}: pixel values {unique[:10]} are not only 0 and 255')
    mask = np.where(values > 0, 255, 0).astype(np.uint8)
    Image.fromarray(mask, mode='L').save(dst, format='TIFF')
    return float((mask == 255).mean())


def lumen_map(path, size=256, block=16):
    """Block map of white, unstained space (vessel lumens, gaps) in an RGB check image.
    These spaces are empty on both stains, so their layout can be compared across stains."""
    with Image.open(path) as image:
        rgb = np.asarray(image.convert('RGB').resize((size, size)), dtype=float)
    grey, spread = rgb.mean(axis=2), rgb.max(axis=2) - rgb.min(axis=2)
    white = ((grey > 215) & (spread < 40)).astype(float)
    n = size // block
    return white.reshape(n, block, n, block).mean(axis=(1, 3))


def orientation_scores(he_path, penta_path):
    """Correlation of the H&E and pentachrome lumen maps with the H&E image rotated by
    0, 90, 180 and 270 degrees (counter-clockwise)."""
    he, penta = lumen_map(he_path), lumen_map(penta_path)
    scores = []
    for k in range(4):
        a, b = np.rot90(he, k).ravel(), penta.ravel()
        scores.append(float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and b.std() > 0 else float('nan'))
    return scores


def move_aside(folder, stamp):
    if os.path.isdir(folder) and os.listdir(folder):
        target = f'{folder}_previous_{stamp}'
        shutil.move(folder, target)
        print(f'  Moved existing {os.path.basename(folder)}/ to {os.path.basename(target)}/')
    os.makedirs(folder, exist_ok=True)


def main():
    if not QUPATH_MASKS_DIR or not os.path.isdir(QUPATH_MASKS_DIR):
        sys.exit('Set QUPATH_MASKS_DIR to the masks/ folder of the QuPath project '
                 f'(now: "{QUPATH_MASKS_DIR}").')
    src = {name: os.path.join(QUPATH_MASKS_DIR, name)
           for name in ('Penta', 'H&E', 'Penta_aligned', 'H&E_aligned',
                        'check_aligned', 'check_aligned_penta')}
    angle_path = os.path.join(QUPATH_MASKS_DIR, 'roi_angles.csv')
    for label, path in list(src.items()) + [('roi_angles.csv', angle_path)]:
        if not os.path.exists(path):
            sys.exit(f'Missing {path}. Run the three QuPath scripts with Run for project first.')

    angles, images = read_angles(angle_path)
    penta = {roi_name(p): p for p, *_ in index_mask_folder(src['Penta']).values()}
    he = {roi_name(p): p for p, *_ in index_mask_folder(src['H&E']).values()}
    penta_al = {roi_name(p): p for p, *_ in index_mask_folder(src['Penta_aligned']).values()}
    he_al = {roi_name(p): p for p, *_ in index_mask_folder(src['H&E_aligned']).values()}

    names = sorted(set(penta) | set(he))
    missing = sorted(set(penta) ^ set(he))
    if missing:
        sys.exit(f'ROIs exported for one stain only: {missing}')

    # Decide the source of every mask before writing anything
    plan = []
    for name in names:
        if ('HE', name) not in angles:
            sys.exit(f'{name}: no H&E angle in roi_angles.csv. Rerun export_masks_aligned.groovy.')
        dev_he = axis_deviation(angles[('HE', name)])
        dev_penta = axis_deviation(angles.get(('Penta', name), 0.0))
        if dev_he <= AXIS_ALIGNED_TOLERANCE_DEG and dev_penta <= AXIS_ALIGNED_TOLERANCE_DEG:
            reference, penta_src, he_src = 'both', penta[name], he[name]
        elif dev_penta <= dev_he:
            if name not in he_al:
                sys.exit(f'{name}: pentachrome is the reference but there is no H&E_aligned mask.')
            reference, penta_src, he_src = 'Penta', penta[name], he_al[name]
        else:
            if name not in penta_al:
                sys.exit(f'{name}: H&E is the reference but there is no Penta_aligned mask.')
            reference, penta_src, he_src = 'H&E', penta_al[name], he[name]
        plan.append((name, angles[('HE', name)], angles.get(('Penta', name)),
                     reference, penta_src, he_src))

    stamp = time.strftime('%Y%m%d_%H%M%S')
    print(f'Writing {len(plan)} mask pairs to {DATA_DIR}')
    for folder in (PENTA_MASK_DIR, HE_MASK_DIR, *FRACLAC_INPUT_DIRS.values()):
        move_aside(folder, stamp)

    os.makedirs(OUT_DIR, exist_ok=True)
    log_path = os.path.join(OUT_DIR, 'mask_sources.csv')
    with open(log_path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(['ROI', 'Angle_HE_deg', 'Angle_Penta_deg', 'Reference',
                         'Penta_Source', 'HE_Source', 'Penta_Collagen_Fraction',
                         'HE_Collagen_Fraction'])
        for name, angle_he, angle_penta, reference, penta_src, he_src in plan:
            out_name = f'{name}_mask.tif'
            fg_penta = write_standard_mask(penta_src, os.path.join(PENTA_MASK_DIR, out_name))
            fg_he = write_standard_mask(he_src, os.path.join(HE_MASK_DIR, out_name))
            shutil.copy2(os.path.join(PENTA_MASK_DIR, out_name), FRACLAC_INPUT_DIRS['Penta'])
            shutil.copy2(os.path.join(HE_MASK_DIR, out_name), FRACLAC_INPUT_DIRS['H&E'])
            rel = lambda p: os.path.relpath(p, QUPATH_MASKS_DIR)
            writer.writerow([name, f'{angle_he:.3f}',
                             '' if angle_penta is None else f'{angle_penta:.3f}',
                             reference, rel(penta_src), rel(he_src),
                             f'{fg_penta:.4f}', f'{fg_he:.4f}'])
            print(f'  {name:32s} reference {reference:5s} | Penta <- {rel(penta_src)} | '
                  f'H&E <- {rel(he_src)}')

    # Normal skin is mostly collagen on both stains; a near-empty normal mask means the
    # export went wrong (e.g. the wrong classifier was applied), not a biological result.
    with open(log_path, newline='', encoding='utf-8') as handle:
        for row in csv.DictReader(handle):
            if 'NORMAL' in row['ROI'].upper():
                for stain in ('Penta', 'HE'):
                    if float(row[f'{stain}_Collagen_Fraction']) < 0.05:
                        print(f'  WARNING {row["ROI"]} ({stain}): collagen fraction '
                              f'{float(row[f"{stain}_Collagen_Fraction"]):.4f} in normal skin. '
                              'Check the QuPath export before running FracLac.')

    # Check images: the pair from the export that matches each ROI's reference stain
    check_dir = CHECK_DIRS[0]
    move_aside(check_dir, stamp)
    for name, _, _, reference, _, _ in plan:
        source = src['check_aligned_penta'] if reference == 'H&E' else src['check_aligned']
        for suffix in ('_HE.png', '_penta.png'):
            path = os.path.join(source, name + suffix)
            if os.path.exists(path):
                shutil.copy2(path, check_dir)
            else:
                print(f'  Note: no check image {os.path.relpath(path, QUPATH_MASKS_DIR)}')

    # Orientation check, averaged over the ROIs of each image: the H&E and pentachrome
    # check images should line up without rotation. A clearly better match after a
    # 90-degree turn means the start corner is wrong for that image.
    orient_path = os.path.join(OUT_DIR, 'orientation_check.csv')
    by_image = {}
    with open(orient_path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(['ROI', 'Reference', 'Rotated_Image', 'r_0', 'r_90', 'r_180', 'r_270'])
        for name, _, _, reference, _, _ in plan:
            if reference == 'both':
                continue   # no rotated export, nothing to check
            he_png = os.path.join(check_dir, name + '_HE.png')
            penta_png = os.path.join(check_dir, name + '_penta.png')
            if not (os.path.exists(he_png) and os.path.exists(penta_png)):
                continue
            # the rotated export comes from the non-reference stain
            image = images.get(('HE', name) if reference == 'Penta' else ('Penta', name), '?')
            scores = orientation_scores(he_png, penta_png)
            by_image.setdefault((image, reference), []).append(scores)
            writer.writerow([name, reference, image] + [f'{x:.3f}' for x in scores])
    flagged = False
    for (image, reference), rows in sorted(by_image.items()):
        mean = np.nanmean(np.array(rows), axis=0)
        best = int(np.nanargmax(mean))
        script = ('export_masks_aligned.groovy' if reference == 'Penta'
                  else 'export_masks_aligned_penta.groovy')
        line = (f'{image}: mean r over {len(rows)} ROIs = ' +
                ', '.join(f'{k * 90} deg {m:.2f}' for k, m in enumerate(mean)))
        if best != 0 and mean[best] >= ORIENTATION_MIN_R and mean[best] - mean[0] >= ORIENTATION_MARGIN:
            flagged = True
            print(f'  WARNING orientation {line}. The check images match better after turning '
                  f'the H&E image {best * 90} deg counter-clockwise. The sections may be '
                  f'mounted rotated: set START_OVERRIDE for this image in {script}, rerun '
                  f'that script and this one.')
        else:
            print(f'  Orientation OK  {line}')
    if not flagged:
        print(f'  Orientation check passed ({orient_path}).')

    print(f'\nDone. Sources of every mask: {log_path}')
    print('Next: run FracLac on "FracLac input/Penta" and "FracLac input/H&E" '
          'with "Lock black background" ticked.')


if __name__ == '__main__':
    main()
