"""Paths and constants shared by all analysis scripts. Edit only this file."""
import os

# analysis/ -> he_vs_pentachrome/ -> scripts/ -> repository root
REPO_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))

# Data are not in the repository. Either set the environment variable
# FRACTAL_DATA_DIR, or put the folders listed in the README in
# <repo>/inputs/he_vs_pentachrome.
DATA_DIR = os.environ.get('FRACTAL_DATA_DIR',
                          os.path.join(REPO_DIR, 'inputs', 'he_vs_pentachrome'))

# masks/ folder of the QuPath project (written by the three QuPath scripts).
# Either set the environment variable QUPATH_MASKS_DIR, or put the path here.
# assemble_masks.py reads it and writes Penta/, H&E/ and the check folders below.
QUPATH_MASKS_DIR = os.environ.get('QUPATH_MASKS_DIR', '')

PENTA_MASK_DIR = os.path.join(DATA_DIR, 'Penta')
HE_MASK_DIR = os.path.join(DATA_DIR, 'H&E')
# Copies of the masks for FracLac, which writes its own image files into the folder it reads.
FRACLAC_INPUT_DIRS = {'Penta': os.path.join(DATA_DIR, 'FracLac input', 'Penta'),
                      'H&E': os.path.join(DATA_DIR, 'FracLac input', 'H&E')}
PENTA_BOX_DIRS = [os.path.join(DATA_DIR, 'NORMAL SKIN BOXCOUNT'),
                  os.path.join(DATA_DIR, 'WOUND SKIN BOXCOUNT')]
HE_BOX_DIRS = [os.path.join(DATA_DIR, 'HE BOXCOUNT')]
# *_HE.png / *_penta.png check images (figure only). assemble_masks.py copies, for each
# ROI, the pair from the same export as its masks (check_aligned or check_aligned_penta).
CHECK_DIRS = [os.path.join(DATA_DIR, 'check_images')]

OUT_DIR = os.path.join(REPO_DIR, 'outputs', 'he_vs_pentachrome')

# assemble_masks.py: an ROI counts as axis-aligned on a stain if its angle is within
# this many degrees of a multiple of 90. The stain where it is closer to axis-aligned
# is the reference; the other stain uses its *_aligned mask.
AXIS_ALIGNED_TOLERANCE_DEG = 0.25

PIXEL_SIZE_UM = 0.2431   # mask export resolution (classifiers work at 1.94 µm/px)
