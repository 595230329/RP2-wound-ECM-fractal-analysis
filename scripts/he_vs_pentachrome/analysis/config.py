"""Paths and constants shared by all analysis scripts. Edit only this file."""
import os

# analysis/ -> he_vs_pentachrome/ -> scripts/ -> repository root
REPO_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))

# Data are not in the repository. Either set the environment variable
# FRACTAL_DATA_DIR, or put the folders listed in the README in
# <repo>/inputs/he_vs_pentachrome.
DATA_DIR = os.environ.get('FRACTAL_DATA_DIR',
                          os.path.join(REPO_DIR, 'inputs', 'he_vs_pentachrome'))

PENTA_MASK_DIR = os.path.join(DATA_DIR, 'Penta')
HE_MASK_DIR = os.path.join(DATA_DIR, 'H&E')
PENTA_BOX_DIRS = [os.path.join(DATA_DIR, 'NORMAL SKIN BOXCOUNT'),
                  os.path.join(DATA_DIR, 'WOUND SKIN BOXCOUNT')]
HE_BOX_DIRS = [os.path.join(DATA_DIR, 'HE BOXCOUNT')]
CHECK_DIR = os.path.join(DATA_DIR, 'check_aligned')   # *_HE.png / *_penta.png from QuPath

OUT_DIR = os.path.join(REPO_DIR, 'outputs', 'he_vs_pentachrome')

PIXEL_SIZE_UM = 0.2431   # mask export resolution (classifiers work at 1.94 µm/px)
