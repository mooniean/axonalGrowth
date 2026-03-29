"""
functions.py – helper functions for the axonal growth / mRNA transport model.

Contents
--------
  Phase-field helpers
      heaviside        – vectorised Heaviside step function
      h                – smoothed indicator  h(x) = x²(3-2x)
      smoothing        – Allen-Cahn relaxation to a diffuse interface
      init_cell        – initialise a filled disc phase field
      container        – confinement potential (keeps a field inside a carrier)

  Geometric helpers
      com              – centre of mass of a phase field
      interface_pos    – 1-D interface position (legacy, 1-D only)
      unit_vector      – vectorised unit vector

  Reaction helpers
      coeff_beta_      – spatially varying release rate prefactor  1/mtb

  Microtubule track helpers
      rasterized_line  – integer Bresenham-style grid path between two points
      build_transport_field – directed motor-velocity field along an ordered track

  mRNA transport
      conservative_upwind_advection – first-order upwind -div(field * velocity)
"""

import numpy as np
from scipy.ndimage import convolve

# ── Phase-field indicator functions ──────────────────────────────────────────


def heaviside(x):
    """Heaviside step function H(x): returns 1 for x ≥ 0, 0 otherwise."""
    if x >= 0:
        return 1
    else:
        return 0


heaviside = np.vectorize(heaviside)


# h(x) = x²(3 - 2x) is the standard smoothed-step (Hermite) interpolant.
# It maps [0, 1] → [0, 1] with h(0)=0, h(1)=1, h'(0)=h'(1)=0.
# Used to compute a smooth volume proxy for the Lagrange multiplier.
def h(x):
    return x**2 * (3.0 - 2.0 * x)


def smoothing(phi, tstep, stencil):
    """Relax a phase field to a diffuse interface via volume-conserving Allen-Cahn.

    Integrates the Allen-Cahn equation
        ∂_t φ = ∇²φ + φ(1-φ)(φ - 0.5 + λ)
    with a Lagrange multiplier λ = V_target - v(t) that pins the total volume
    h(φ).sum() to its initial value.

    Parameters
    ----------
    phi    : 2-D ndarray
        Phase-field to smooth (modified in-place and returned).
    tstep  : int
        Number of Allen-Cahn steps to perform.
    stencil : 2-D ndarray
        Discrete Laplacian kernel (see parameters.py).

    Returns
    -------
    phi : 2-D ndarray
        Smoothed phase field.
    """
    nstep = 0
    V_target = h(phi).sum()  # conserved volume target
    dt = 0.001
    while nstep <= tstep:
        v = h(phi).sum()
        Lagrange = np.ones(phi.shape) * (V_target - v)
        phi = phi + dt * (
            convolve(phi, stencil, mode="wrap") + phi * (1.0 - phi) * (phi - 0.5 + Lagrange)
        )
        nstep += 1
    return phi


def com(phi, volume, ndim):
    """Compute the centre of mass (centroid) of a 2-D phase field.

    Parameters
    ----------
    phi    : 2-D ndarray
        Phase field whose centroid is sought.
    volume : float
        Unused legacy parameter (recomputed internally as phi.sum()).
    ndim   : int
        Number of spatial dimensions (2 for this model).

    Returns
    -------
    r_cm : 1-D ndarray of shape (ndim,)
        Coordinates of the centroid.
    """
    r_cm = np.zeros(ndim)
    volume = phi.sum()  # recompute to avoid passing stale value
    for i in range(phi.size):
        s = np.asarray(np.unravel_index(i, phi.shape))
        r_cm = r_cm + s * phi.item(i)
    return r_cm / volume


def interface_pos(phi):
    """Return the 1-D centre-of-mass position of a 1-D phase field (legacy).

    Only meaningful for 1-D or when phi has a single active row/column.
    """
    pos = 0
    for i in range(phi.size):
        pos = pos + i * phi.item(i)
    return pos / phi.sum()


def init_cell(phi, cell_position, cell_radius):
    """Initialise a circular (disc) region of a 2-D phase field to 1.

    Sets phi[r] = 1 for all grid points r whose Euclidean distance from
    cell_position is ≤ cell_radius; leaves the rest at 0.

    Parameters
    ----------
    phi           : 2-D ndarray
        Phase field to modify (in-place).
    cell_position : array-like of shape (2,)
        Centre of the disc in grid coordinates.
    cell_radius   : float
        Radius of the disc in grid units.

    Returns
    -------
    phi : 2-D ndarray
        Modified phase field.
    """
    for i in range(phi.size):
        s = np.asarray(np.unravel_index(i, phi.shape))
        distance = np.sqrt(np.sum(np.power(s - cell_position, 2)))
        if distance <= cell_radius:
            phi.flat[i] = 1.0
    return phi


def init_growth_factor(ngf, neuron_position, neuron_radius, rng=None):
    """Initialise the NGF field as a random background outside the neuron soma.

    Sets ngf[r] ~ Uniform(0, 10) for grid points more than (neuron_radius + 2)
    away from neuron_position, creating a noisy chemo-attractive landscape.

    Parameters
    ----------
    ngf             : 2-D ndarray
        NGF field to modify (in-place).
    neuron_position : array-like of shape (2,)
        Centre of the neuron soma.
    neuron_radius   : float
        Radius of the soma exclusion zone.
    rng             : numpy.random.Generator or None, optional
        Random number generator to use.  Pass ``np.random.default_rng(seed)``
        for reproducible results.  Defaults to ``np.random.default_rng()``.

    Returns
    -------
    ngf : 2-D ndarray
        Initialised NGF field.
    """
    # Use the modern Generator API (NPY002); create a fresh one if none supplied.
    if rng is None:
        rng = np.random.default_rng()
    for i in range(ngf.size):
        s = np.asarray(np.unravel_index(i, ngf.shape))
        distance = np.sqrt(np.sum(np.power(s - neuron_position, 2)))
        if distance > neuron_radius + 2:
            ngf.flat[i] = rng.random() * 10.0
    return ngf


def container(field_1, field_2, alpha):
    """Confinement potential that keeps field_1 inside the support of field_2.

    Returns a normalised repulsion term
        field_1 * μ / sqrt(μ(field_1² + α))
    where  μ = (1 - field_2)² + α.

    When field_2 ≈ 1 (inside the axon), μ ≈ α (small), so the term
    is approximately field_1/sqrt(field_1²+α) ≈ ±1, and the Laplacian of
    this term acts like a confining force.  When field_2 ≈ 0 (outside),
    μ ≈ 1 + α, and the term is stronger, repelling field_1 back inward.

    Parameters
    ----------
    field_1 : ndarray
        Chemical field to be confined (mf, ml, or mtb).
    field_2 : ndarray
        Carrier phase field (psi – the axon body).
    alpha   : float
        Regularisation constant preventing division by zero (alpha_ ≈ 1e-3).

    Returns
    -------
    ndarray
        Confinement term, same shape as field_1.
    """
    mu = np.power(1.0 - field_2, 2) + alpha
    return field_1 * mu * np.power(mu * (np.power(field_1, 2) + alpha), -0.5)


def unit_vector(vector, norm):
    """Return vector / norm if norm > 0, else 0 (vectorised element-wise).

    Used to normalise the NGF gradient to a unit direction for the chemotactic
    forcing term.  Vectorised so it can operate on arrays.
    """
    if norm > 0:
        return vector / norm
    else:
        return 0


unit_vector = np.vectorize(unit_vector)


def coeff_beta_(mtb):
    """Compute the microtubule-density-dependent release prefactor 1/mtb.

    Returns 1/mtb when mtb > 0.1 (threshold avoids division near zero),
    and 0 otherwise.  Vectorised so it operates element-wise on arrays.

    This factor appears in the ml → mf release term:
        coeff_beta = (1/mtb) * beta_m * (1 - mtb) * ml
    The combination (1/mtb)*(1-mtb) = (1-mtb)/mtb decreases monotonically
    from ∞ at mtb→0 to 0 at mtb=1, so release is effectively zero where the
    MT rail is fully saturated.
    """
    if mtb > 10e-2:
        return 1.0 / mtb
    else:
        return 0.0


coeff_beta_ = np.vectorize(coeff_beta_)


# ── Microtubule track helpers ─────────────────────────────────────────────────


def rasterized_line(start, end):
    """Return all integer grid points along the segment from start to end.

    Uses linear interpolation with spacing 1 grid cell, then removes
    duplicate points (can arise along axis-aligned segments).  This is
    essentially a continuous-space version of Bresenham's line algorithm.

    Parameters
    ----------
    start : array-like of shape (2,)
        Start point in grid coordinates (rounded to nearest integer).
    end   : array-like of shape (2,)
        End point in grid coordinates (rounded to nearest integer).

    Returns
    -------
    line : 2-D ndarray of shape (N, 2), dtype int
        Ordered array of (row, col) grid coordinates along the segment.
    """
    start = np.asarray(np.rint(start), dtype=int)
    end = np.asarray(np.rint(end), dtype=int)
    # Number of steps = longest axis span, ensuring each step is ≤ 1 cell.
    nsteps = int(np.max(np.abs(end - start))) + 1

    if nsteps <= 1:
        return np.asarray([start], dtype=int)

    line = np.rint(np.linspace(start, end, nsteps)).astype(int)
    # Remove consecutive duplicates that can appear in axis-aligned segments.
    keep = np.ones(line.shape[0], dtype=bool)
    keep[1:] = np.any(np.diff(line, axis=0) != 0, axis=1)
    return line[keep]


def build_transport_field(track_points, shape, support_radius=1):
    """Build a directed motor-velocity field aligned with an ordered MT track.

    For each consecutive pair of track points, the unit tangent vector is
    accumulated into all grid cells within support_radius.  After processing
    all segments, the velocity at each cell is normalised by the accumulated
    weight so that |v_m| ≤ 1 everywhere.

    Parameters
    ----------
    track_points    : list of (row, col) tuples
        Ordered sequence of grid points defining the microtubule track,
        running from soma to current growth-cone position.
    shape           : tuple (Nx, Ny)
        Grid shape (same as L in parameters.py).
    support_radius  : int, optional
        Number of cells around each track point that receive the velocity.
        Default 1 means a 3×3 neighbourhood; 0 means single-cell width.

    Returns
    -------
    velocity : ndarray of shape (2, Nx, Ny)
        velocity[0] – row component of the motor field (soma→tip direction)
        velocity[1] – col component of the motor field
    """
    velocity = np.zeros((2,) + tuple(shape), dtype=float)
    weights = np.zeros(shape, dtype=float)

    if len(track_points) < 2:
        return velocity  # no track yet; return zero field

    track_points = np.asarray(track_points, dtype=int)
    max_i, max_j = shape[0] - 1, shape[1] - 1

    for current, nxt in zip(track_points[:-1], track_points[1:]):
        delta = nxt - current
        norm = np.linalg.norm(delta)
        if norm == 0:
            continue  # degenerate segment; skip

        direction = delta / norm  # unit tangent toward the tip
        for di in range(-support_radius, support_radius + 1):
            for dj in range(-support_radius, support_radius + 1):
                i = current[0] + di
                j = current[1] + dj
                if 0 <= i <= max_i and 0 <= j <= max_j:
                    velocity[:, i, j] += direction
                    weights[i, j] += 1.0

    # Normalise by accumulated weight to get a proper unit-vector field.
    mask = weights > 0
    velocity[0, mask] /= weights[mask]
    velocity[1, mask] /= weights[mask]
    return velocity


# ── Conservative mRNA transport ───────────────────────────────────────────────


def conservative_upwind_advection(field, velocity, dL=1.0):
    """First-order upwind approximation of -div(field · velocity).

    Implements the finite-volume update
        (d/dt) field ≈ -∇·(field · v)
    using upwind-biased face fluxes:
        F_{i+1/2} = v_{i+1/2} * field_i   if v_{i+1/2} ≥ 0
                  = v_{i+1/2} * field_{i+1} otherwise

    The face velocity v_{i+1/2} is the arithmetic average of adjacent cell
    velocities.  The resulting stencil is mass-conservative (sum of advection
    over all cells is zero to machine precision) and first-order accurate in
    space.

    Stability requires the Courant number  max|v| * dt / dL < 1.
    For the linked-mRNA equation: chi_ml * |v_m| * dt / dL < 1.

    Parameters
    ----------
    field    : 2-D ndarray
        Transported scalar field (e.g. coeff_ml * ml * mtb_effective).
    velocity : ndarray of shape (2, Nx, Ny)
        Motor-velocity field v_m scaled by chi_ml.
    dL       : float, optional
        Grid spacing (default 1.0).

    Returns
    -------
    advection : 2-D ndarray
        Discrete divergence term -div(field * velocity), same shape as field.
        Add this directly to the ml update:  ml += dt * advection.
    """
    advection = np.zeros_like(field, dtype=float)
    v0, v1 = velocity

    # ── Row (axis-0) fluxes ───────────────────────────────────────────────────
    # Face velocity: average of neighbouring cell velocities.
    face_v0 = 0.5 * (v0[:-1, :] + v0[1:, :])
    # Upwind flux: pick upstream cell depending on sign of face velocity.
    flux0 = np.where(face_v0 >= 0.0, face_v0 * field[:-1, :], face_v0 * field[1:, :])
    advection[:-1, :] -= flux0 / dL  # outflow from cell i
    advection[1:, :] += flux0 / dL  # inflow  to  cell i+1

    # ── Column (axis-1) fluxes ────────────────────────────────────────────────
    face_v1 = 0.5 * (v1[:, :-1] + v1[:, 1:])
    flux1 = np.where(face_v1 >= 0.0, face_v1 * field[:, :-1], face_v1 * field[:, 1:])
    advection[:, :-1] -= flux1 / dL
    advection[:, 1:] += flux1 / dL

    return advection
