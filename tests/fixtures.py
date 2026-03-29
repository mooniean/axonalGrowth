"""
tests/fixtures.py – shared helpers for building minimal simulation .npz files.

Import these from any test module that needs real or synthetic snapshot data.
All helpers write to a caller-supplied ``tmp_path`` (a ``pathlib.Path`` to a
temporary directory) so tests stay isolated and leave no artefacts behind.
"""

from pathlib import Path
from typing import List, Optional

import numpy as np

# Default field shape used in all synthetic snapshots.
SHAPE = (10, 12)

# RNG with a fixed seed so synthetic data is fully deterministic.
RNG = np.random.default_rng(7)


def make_snapshot(
    path: Path,
    shape: tuple = SHAPE,
    include_vm: bool = True,
    extra_fields: Optional[dict] = None,
) -> Path:
    """Write a minimal synthetic simulation .npz snapshot and return the path.

    Parameters
    ----------
    path : Path
        Destination file path (must have the ``.npz`` suffix).
    shape : tuple
        2-D grid shape ``(Nx, Ny)`` for all scalar fields.
    include_vm : bool
        Whether to include the ``v_m`` motor-velocity field ``(2, Nx, Ny)``.
    extra_fields : dict, optional
        Additional arrays to include (key → ndarray).

    Returns
    -------
    Path
        Same as ``path``, for convenience.
    """
    data = {
        "phi": RNG.random(shape).astype(np.float64),
        "psi": RNG.random(shape).astype(np.float64),
        "ngf": RNG.random(shape).astype(np.float64),
        "mf":  RNG.random(shape).astype(np.float64),
        "ml":  RNG.random(shape).astype(np.float64),
        "mtb": RNG.random(shape).astype(np.float64),
    }
    if include_vm:
        data["v_m"] = RNG.random((2,) + shape).astype(np.float64)
    if extra_fields:
        data.update(extra_fields)

    np.savez(path, **data)
    return path


def make_snapshot_dir(
    directory: Path,
    n_steps: int = 3,
    step_size: int = 10000,
    prefix: str = "sim_",
    shape: tuple = SHAPE,
) -> List[Path]:
    """Write *n_steps* sequentially numbered snapshots into *directory*.

    File names follow the pattern ``{prefix}{step}.npz`` so that
    ``_step_sort_key`` orders them correctly.

    Returns the list of written paths in chronological order.
    """
    directory.mkdir(parents=True, exist_ok=True)
    paths = []
    for i in range(n_steps):
        step = i * step_size
        p = make_snapshot(directory / f"{prefix}{step}.npz", shape=shape)
        paths.append(p)
    return paths


def make_snapshot_no_vm(path: Path, shape: tuple = SHAPE) -> Path:
    """Write a snapshot that deliberately omits the ``v_m`` field."""
    return make_snapshot(path, shape=shape, include_vm=False)


def make_snapshot_partial(path: Path, fields: List[str], shape: tuple = SHAPE) -> Path:
    """Write a snapshot that contains only the listed field names."""
    data = {f: RNG.random(shape).astype(np.float64) for f in fields}
    np.savez(path, **data)
    return path


