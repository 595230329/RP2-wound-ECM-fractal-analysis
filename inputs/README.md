# Input files

Do not commit human histology images, patient-identifiable information, or source workbooks unless you have written authorisation to make them public. The input folders are intentionally ignored by Git.

## 1. QuPath ROI-annotated ECM coverage workbooks

Place the 15 `.xlsx` workbooks in:

```text
inputs/ecm_coverage/
```

The pipeline identifies each metric from the file name. File names must contain one of the following phrases:

- `collagen coverage`
- `mucin coverage`
- `necrosis coverage`

Each workbook must include these columns:

```text
Sample_ID | region | ROI_ID
```

It must also include QuPath measurement columns for the relevant component area, background/non-component area and percentage coverage. The expected column suffixes are documented in `scripts/HistologyData_Pipeline_MixedEffects.py` under `METRICS_INFO`.

Recognised region labels include `Normal Skin`, `Wound Centre`, `Distant Edge`, `Proximal Edge`, and legacy equivalents documented in the script.

## 2. FracLac Box Count Summary workbooks

Place the 12 normal-skin exports in:

```text
inputs/fraclac/normal_skin_boxcount/
```

Place the 18 wound exports in:

```text
inputs/fraclac/wound_skin_boxcount/
```

The script expects 30 FracLac Box Count Summary `.xlsx` workbooks. It extracts the prespecified summary fields from the second row of the FracLac layout. Do not alter column ordering in the export.

### Required naming convention

Normal-skin filenames must contain the section code (`9M06`, `9M07`, `9M08` or `9M09`) and a normal ROI identifier such as `NORMAL SKIN01`.

Wound filenames must contain one of these unambiguous section codes:

```text
9M0501 -> Wound_05_1
9M0502 -> Wound_05_2
9M06   -> Wound_06
9M07   -> Wound_07
9M08   -> Wound_08
9M09   -> Wound_09
```

They must also contain one of `DISTANT EDGE`, `PROXIMAL EDGE`, or `WOUND CENTER`/`WOUND CENTRE`.

The Wound 05 files must use `9M0501` and `9M0502`; the ambiguous code `9M05` is deliberately not accepted because it cannot identify which of the two wound sections supplied the mask.

## Outputs used as HWHI inputs

Do not manually copy these files. After the first and third scripts run, the HWHI script automatically reads:

```text
outputs/ecm/0_SectionRegion_Summary.csv
outputs/fraclac/0_FracLac_QC_Passed_Data.csv
```
