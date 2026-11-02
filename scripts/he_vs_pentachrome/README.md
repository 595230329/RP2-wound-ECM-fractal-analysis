# H&E vs Movat pentachrome: fractal analysis of collagen masks

Code for testing whether box-counting fractal analysis of collagen masks from
H&E sections gives the same results as Movat pentachrome on adjacent sections
of human skin (unwounded skin and acute wounds).

All paths below are relative to the repository root.

## Contents

```
scripts/he_vs_pentachrome/
├── qupath/
│   ├── export_masks.groovy                 unrotated square ROIs, both stains
│   ├── export_masks_aligned.groovy         rotated ROIs on H&E (pentachrome = reference)
│   ├── export_masks_aligned_penta.groovy   rotated wound ROIs on pentachrome (H&E = reference)
│   └── classifiers/                        HE_Collagen_v1.json, Pentachrome_Collagen_v1.json
├── imagej/
│   └── invert_for_fraclac.ijm              inverted copies of masks for FracLac
└── analysis/
    ├── config.py                           paths and pixel size (edit only this file)
    ├── naming.py                           file name -> Mask_ID
    ├── compare_masks.py                    pixel overlap (QC) + paired FracLac outputs
    ├── make_figure.py                      paper figure
    └── coarse_check.py                     optional registration QC
```

## Requirements

- QuPath 0.7.0
- ImageJ 1.54g with FracLac 2015Sep090313a9330
- Python ≥ 3.9 and the packages in the root `requirements.txt`:
  ```bash
  pip install -r requirements.txt
  ```

## Pipeline

### 1. Export masks (QuPath)

Copy the two classifiers into `<QuPath project>/classifiers/pixel_classifiers/`.
Run each script with *Automate > Run for project*.

| Script | ROIs | Reference frame | Output |
|---|---|---|---|
| `export_masks.groovy` | Unrotated squares, both stains | — | `masks/Penta`, `masks/H&E` |
| `export_masks_aligned.groovy` | Drawn on pentachrome, rotated on H&E | Pentachrome | `masks/H&E_aligned`, `masks/check_aligned` |
| `export_masks_aligned_penta.groovy` | Wound ROIs redrawn on H&E, rotated on pentachrome | H&E | `masks/Penta_aligned`, `masks/check_aligned_penta` |

- ROIs are 500 µm squares. Their annotation names contain the sample code (see *File naming*).
- Masks are 8-bit at 0.2431 µm/px (about 2057 × 2057 px), with collagen = 255 and background = 0.
- The classifiers decide at 1.94 µm/px, so each classified pixel covers an 8 × 8 block in the mask.
- The pentachrome classifier takes yellow only as collagen. Green areas (collagen with mucin) are excluded.
- For rotated ROIs, the start corner (`START_BY_IMAGE` / `START_CORNER`) sets the orientation.
  Confirm it for each image by comparing the check images. Resampling is nearest neighbour.

### 2. Invert masks for FracLac (ImageJ)

FracLac treats black as foreground, whatever the ImageJ *Black background* setting.
The log line `changed inverted LUT on ...` shows this has happened. Run
`imagej/invert_for_fraclac.ijm` to write **inverted copies** to a separate folder,
and run FracLac on those copies. Keep the original masks (collagen = white) for step 4.

### 3. Fractal analysis (FracLac)

Settings: box counting; binary image; no filter; 12 grid positions; square boxes;
100 box sizes; smallest box 5 px; largest box 45 % of image width.

Save the summary tables (`.xlsx` or `*summary*.txt`) into the BOXCOUNT folders listed under *Data layout*.

### 4. Compare H&E and pentachrome

```bash
py scripts/he_vs_pentachrome/analysis/compare_masks.py
```

- **Part A** (QC only): pixel-level overlap between paired masks.
  It reports Dice, tolerant F at 2, 5 and 10 µm, chance-level Dice, and residual-shift correction.
  It writes `A1_Pixel_Agreement_Per_Mask.csv`, `A2_Pixel_Agreement_By_Region.csv`, plots and overlays.
- **Part B**: pairs the FracLac outputs by Mask_ID and writes
  `B1_FracLac_Paired_Table.csv` and `B2_FracLac_Agreement_Statistics.csv`, plus Bland–Altman plots.

Check the console output for the number of paired masks, and for duplicate or skipped messages.
If no pairs are found, the B1 table is **not** rewritten.

### 5. Figure

```bash
py scripts/he_vs_pentachrome/analysis/make_figure.py
```

Always run this after step 4, because the figure reads `B1_FracLac_Paired_Table.csv`.
Set the ROIs shown in panels A and B with `NORMAL_ID` and `WOUND_ID` at the top of the script.

### 6. Optional: coarse registration check

```bash
py scripts/he_vs_pentachrome/analysis/coarse_check.py
```

This prints the H&E vs pentachrome correlation after averaging the masks into 10–100 µm blocks.
The null value comes from the H&E map rotated by 90° steps, so it has no spatial correspondence.

## Data layout

Data are not in the repository. Put these folders in `inputs/he_vs_pentachrome/`,
or point the environment variable `FRACTAL_DATA_DIR` to a folder that contains them:

```
Penta/                   pentachrome masks (original, collagen = white)
H&E/                     H&E masks (original, collagen = white)
NORMAL SKIN BOXCOUNT/    pentachrome FracLac summaries, normal skin
WOUND SKIN BOXCOUNT/     pentachrome FracLac summaries, wounds
HE BOXCOUNT/             H&E FracLac summaries
check_aligned/           *_HE.png and *_penta.png check images (figure only)
```

All subfolders are read, so do not keep `_old` folders inside the BOXCOUNT folders.
Results are written to `outputs/he_vs_pentachrome/`.

## File naming

Masks and FracLac rows are paired by a Mask_ID built from the file name or row label:

- Sample code: `9M05`, `9M05-01`, `9M05-02`, `9M06`, `9M08`, `9M09` (pentachrome).
  The paired H&E sections are 9M10, 9M11, 9M13 and 9M14.
- Region: `NORMAL SKIN <n>`, `WOUND CENTRE` / `WOUND CENTER`, `PROXIMAL EDGE` or `DISTANT EDGE`.
  Spaces are optional.

Examples: `9M06 NORMAL SKIN 1` → `Normal_06_ROI_1`; `9M05-02 WOUND CENTRE` → `Wound_05_2_Wound_Centre`.

## FracLac columns read

`compare_masks.py` reads these 0-based columns from each FracLac summary row:

| Column | Value |
|---|---|
| 5 | Fractal dimension |
| 10 | R² of the fit |
| 86 | Lacunarity |
| 117 | Foreground pixels |
| 119 | Total pixels |

Foreground fraction = column 117 / column 119. Tables with fewer than 120 columns are skipped.
A row with a non-numeric R² is kept, with R² recorded as NaN.


## Licence

See `LICENSE` in the repository root.
