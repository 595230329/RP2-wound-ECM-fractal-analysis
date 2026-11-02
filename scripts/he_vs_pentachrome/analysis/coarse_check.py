"""QC only: correlation between H&E and pentachrome masks after block-averaging.

Checks registration at 10-100 µm scale; not part of the fractal analysis.
'null' = mean correlation with the H&E map rotated (no spatial correspondence).
"""
import numpy as np
from PIL import Image

from config import HE_MASK_DIR, PENTA_MASK_DIR, PIXEL_SIZE_UM
from naming import index_mask_folder

BLOCKS_UM = [10, 25, 50, 100]


def load(path):
    with Image.open(path) as im:
        return (np.asarray(im.convert('L')) > 127).astype(np.float32)


def block_map(m, b):
    h, w = (m.shape[0] // b) * b, (m.shape[1] // b) * b
    return m[:h, :w].reshape(h // b, b, w // b, b).mean(axis=(1, 3))


def corr(a, b):
    a, b = a.ravel(), b.ravel()
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def null_corr(pm, hm):
    ks = (1, 2, 3) if pm.shape[0] == pm.shape[1] else (2,)   # 90° rotations need a square map
    return np.nanmean([corr(pm, np.rot90(hm, k)) for k in ks])


penta_index = index_mask_folder(PENTA_MASK_DIR)
he_index = index_mask_folder(HE_MASK_DIR)

for mask_id in sorted(set(penta_index) & set(he_index)):
    p, h = load(penta_index[mask_id][0]), load(he_index[mask_id][0])
    n0, n1 = min(p.shape[0], h.shape[0]), min(p.shape[1], h.shape[1])
    p, h = p[:n0, :n1], h[:n0, :n1]
    parts = [mask_id]
    for um in BLOCKS_UM:
        b = max(1, round(um / PIXEL_SIZE_UM))
        pm, hm = block_map(p, b), block_map(h, b)
        parts.append(f'{um}um r={corr(pm, hm):+.2f} (null {null_corr(pm, hm):+.2f})')
    print(' | '.join(parts))
