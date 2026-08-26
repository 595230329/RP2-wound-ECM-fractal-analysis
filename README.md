# Quantitative Spatial ECM Remodelling in Human Wound Histology

Reproducible Python workflows supporting the dissertation: **Quantitative Spatial Analysis of Extracellular Matrix Remodelling in Wound Healing Using Digital Pathology and Fractal Analysis**.

The repository contains three linked analyses:

1. **ECM coverage analysis** — aggregates QuPath tile measurements to area-weighted ROI coverage, fits linear mixed-effects models, applies Holm-adjusted contrasts, and creates coverage/correlation figures.
2. **FracLac analysis** — extracts prespecified fractal dimension, lacunarity and quality-control fields from FracLac Box Count Summary workbooks and fits equivalent mixed-effects models.
3. **HWHI integration** — joins ECM and fractal outputs by `Sample_ID` and anatomical region, standardises five components, and calculates an exploratory Histological Wound Healing Index (HWHI).

## Study design

The workflow was developed for ten digitised human histological sections: four normal-skin sections and six wound sections. Normal sections contribute three spatially separated ROIs per section; wound sections contribute one ROI each from Wound Centre, Proximal Edge and Distant Edge. The resulting analytical design contains 30 ROIs/masks.

Tile-level image measurements are **not** treated as independent biological replicates. ECM tile measurements are aggregated to area-weighted ROI-level values, and mixed-effects models use section `Sample_ID` as the grouping factor.

## Repository structure

```text
RP2-wound-ECM-fractal-analysis/
├── scripts/
│   ├── HistologyData_Pipeline_MixedEffects.py
│   ├── comprehensive_spatial_dynamics_plotter_mixed_effects.py
│   ├── FracLac_Spatial_Analysis_MixedEffects_HWHI_Ready.py
│   └── HWHI_Integrated_Analysis.py
├── inputs/
│   ├── README.md
│   ├── ecm_coverage/                         # private QuPath workbooks
│   └── fraclac/
│       ├── normal_skin_boxcount/              # private FracLac workbooks
│       └── wound_skin_boxcount/
├── outputs/                                   # created automatically; ignored by Git
├── requirements.txt
├── LICENSE
└── CITATION.cff
```

## Installation

This project was developed with Python 3.10+ on Windows. In PowerShell, from the repository root:

```powershell
py -m venv .venv
.\.venv\Scripts\Activate.ps1
py -m pip install --upgrade pip
py -m pip install -r requirements.txt
```

If PowerShell blocks activation, use Command Prompt or consult your institution's Python setup guidance. The scripts can also be executed from PyCharm after configuring the same interpreter and packages.

## Inputs

Raw histology images and original quantitative workbooks are not included. They may contain restricted human-tissue data and are excluded intentionally. See [`inputs/README.md`](inputs/README.md) for the required folder names, file naming conventions and columns.

## Run the complete workflow

Copy authorised input workbooks into the folders specified in `inputs/README.md`, then run these commands in order from the repository root:

```powershell
py scripts\HistologyData_Pipeline_MixedEffects.py
py scripts\comprehensive_spatial_dynamics_plotter_mixed_effects.py
py scripts\FracLac_Spatial_Analysis_MixedEffects_HWHI_Ready.py
py scripts\HWHI_Integrated_Analysis.py
```

Outputs are written to `outputs/ecm`, `outputs/fraclac`, and `outputs/hwhi`. The HWHI script requires the ECM and FracLac scripts to have run successfully because it reads their generated summary tables.

## HWHI definition

The exploratory equal-weight composite index is:

```text
HWHI = [z(Collagen) + z(Fractal Dimension) − z(Mucin)
        − z(Necrosis) − z(Lacunarity)] / 5
```

Each `z` value is calculated across all matched section-region observations as `(value − overall mean) / overall standard deviation`. Positive directions are collagen coverage and fractal dimension; negative directions are mucin coverage, necrosis coverage and lacunarity. HWHI is an exploratory research composite and is **not** a clinically validated diagnostic or prognostic score.

## Reproducibility and limitations

The scripts are designed to be auditable and portable, but reproduce results only when provided with data exported using the required QuPath and FracLac formats. The FracLac parser relies on fixed columns from the Box Count Summary layout; this is documented in the script. Results should be interpreted with the study's small number of independent sections in mind.

## Data availability

Raw human histology images and source workbooks are not deposited in this public code repository. Access, if appropriate, must be governed by the originating study's data-governance and ethical procedures.

## Citation

If you use or adapt this workflow, please cite the associated dissertation/publication and this repository. Update `CITATION.cff` with the final DOI or repository URL after publication.
