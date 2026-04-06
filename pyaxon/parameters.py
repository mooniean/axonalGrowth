"""
parameters.py – fixed numerical constants for the axonal growth model.

User-tunable parameters are loaded at runtime by calling parser() from
pyaxon.functions directly in main.py.
"""

import numpy as np

# ── Fixed constants ───────────────────────────────────────────────────────────

# Plotting defaults
kwargs = {"origin": "lower", "interpolation": "sinc", "cmap": "seismic"}

# Allen-Cahn regularisation / confinement strength (prevents division by zero)
alpha_ = 10e-4

# Boundary mode for scipy.ndimage.convolve
boundary_condition = "constant"

# Physical cell width (all length scales are in units of dL)
dL = 1.0


# 4th-order isotropic Laplacian stencil (5×5)
def _make_stencil(dl: float = 1.0) -> np.ndarray:
    return (1.0 / (12.0 * dl * dl)) * np.array(
        [
            [0, 0, -1, 0, 0],
            [0, 0, 16, 0, 0],
            [-1, 16, -60, 16, -1],
            [0, 0, 16, 0, 0],
            [0, 0, -1, 0, 0],
        ]
    )


stencil = _make_stencil(dL)
