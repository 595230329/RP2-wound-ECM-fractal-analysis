"""File-name parsing shared by all scripts.

Names must contain the sample code (e.g. '9M05', '9M05-01') and either
'NORMAL SKIN <n>' or one of 'WOUND CENTRE/CENTER', 'PROXIMAL EDGE', 'DISTANT EDGE'.
Spaces are optional.
"""
import glob
import os
import re

MASK_EXTENSIONS = ('.tif', '.tiff', '.png', '.bmp', '.jpg', '.jpeg')

WOUND_SAMPLE_MAP = {
    '0501': 'Wound_05_1', '0502': 'Wound_05_2', '05': 'Wound_05',
    '06': 'Wound_06', '07': 'Wound_07', '08': 'Wound_08', '09': 'Wound_09',
}


def parse_file_metadata(path_or_label, is_normal=None):
    """Return (Sample_ID, Region, Mask_ID); raise ValueError if the name cannot be parsed."""
    name = os.path.basename(str(path_or_label))
    upper_name = name.upper()
    sample_match = re.search(r'9M(05-?01|05-?02|05|06|07|08|09)', upper_name)
    if not sample_match:
        raise ValueError(f'Could not identify the sample code in {name}.')
    sample_code = sample_match.group(1).replace('-', '')
    # remove the sample code so its digits are not mistaken for an ROI number
    rest = upper_name[:sample_match.start()] + ' ' + upper_name[sample_match.end():]

    if is_normal is None:
        is_normal = re.search(r'NORMAL\s*SKIN', rest) is not None

    if is_normal:
        sample_id = f'Normal_{sample_code}'
        roi_match = re.search(r'NORMAL\s*SKIN\s*0?(\d{1,2})(?!\d)', rest)
        if not roi_match:
            raise ValueError(f'Could not identify the normal ROI number in {name}.')
        return sample_id, 'Normal Skin', f'{sample_id}_ROI_{int(roi_match.group(1))}'

    sample_id = WOUND_SAMPLE_MAP[sample_code]
    if re.search(r'DISTANT\s*EDGE', rest):
        region = 'Distant Edge'
    elif re.search(r'PROXIMAL\s*EDGE', rest):
        region = 'Proximal Edge'
    elif re.search(r'WOUND\s*CENT(ER|RE)', rest):
        region = 'Wound Centre'
    else:
        raise ValueError(f'Could not identify the wound region in {name}.')
    return sample_id, region, f"{sample_id}_{region.replace(' ', '_')}"


def index_mask_folder(folder):
    """Map Mask_ID -> (path, Sample_ID, Region) for every mask image in a folder tree."""
    index = {}
    for path in sorted(glob.glob(os.path.join(folder, '**', '*'), recursive=True)):
        if not path.lower().endswith(MASK_EXTENSIONS):
            continue
        is_normal = 'NORMAL' in os.path.basename(path).upper()
        try:
            sample_id, region, mask_id = parse_file_metadata(path, is_normal)
        except ValueError as error:
            print(f'  Skipped: {error}')
            continue
        if mask_id in index:
            raise ValueError(f'Duplicate Mask_ID {mask_id} in {folder}.')
        index[mask_id] = (path, sample_id, region)
    return index
