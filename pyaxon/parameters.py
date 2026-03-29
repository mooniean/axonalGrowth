"""
parameters.py – all physical and numerical parameters for the axonal growth model.

Every symbol used in main.py is defined here and imported explicitly.
Changing a value here is the only place you need to edit to run a new parameter study.

Physical model overview
-----------------------
  Grid
      L       = [Nx, Ny]   – domain size in grid cells
      dL      = 1.0        – cell width (physical length per cell)
      dt      = 0.001      – time step (explicit Euler; decrease if unstable)

  Geometry
      neuron_position   – soma centre (row, col)
      neuron_radius     – soma disc radius
      gc_position       – initial growth-cone centre (on the soma boundary)
      gc_radius         – growth-cone disc radius
      small_box_size    – half-width of the local computation window around the GC

  Phase-field
      alpha_           – Allen-Cahn regularisation / confinement strength
      boundary_condition – SciPy convolution boundary mode ('constant' = zero-flux)
      stencil          – 4th-order isotropic finite-difference Laplacian kernel

  NGF (Nerve Growth Factor)
      D_ngf              – diffusion coefficient
      ngf_source_position – point source location (row, col)
      ngf_source_value   – Dirichlet injection value

  Microtubules
      D_mtb         – diffusion coefficient
      lambda_mtb    – linear decay rate
      mtb_source    – value pinned at tracked rail cells each step
      mtb_positions – flat-index list of the MT rail (populated at run time)

  Free mRNA (mf)
      kf              – diffusion coefficient
      lambda_f        – linear decay rate
      mf_source_position – Dirichlet source location (soma centre)
      mf_source_value    – Dirichlet injection value
      gamma           – binding rate  mf → ml  (requires mtb and free rail capacity)

  Linked mRNA (ml)
      kl              – diffusion coefficient
      lambda_l        – linear decay rate
      Mm              – maximum MT rail capacity (saturation concentration)
      chi_ml          – motor transport coefficient; scales v_m in the upwind flux
                        CFL stability: chi_ml * dt / dL < 1
      beta_everywhere – baseline release rate  ml → mf  along the axon shaft
      beta_growth_cone – enhanced release rate  ml → mf  at the growth-cone tip
      beta_m          – spatially varying release field (set each step in main.py)

  Growth-cone chemotaxis and axon growth
      chi             – chemotactic speed coefficient  (v = chi * mf/(1+mf) * ∇NGF̂)
                        reduced from 10 → 3 to prevent the cone from outrunning the shaft
      alpha_p         – proliferative coupling: rate at which phi/ngf drive psi growth
                        increased from 10 → 15 to help the shaft keep up with the cone

  Time integration
      nstep       – current step counter (starts at 0; advanced in main.py)
      tstep       – total number of steps to run
      dt          – time step size
      print_period – save a snapshot every this many steps
      nprint      – counter for the next save (initialised to print_period)
"""

import numpy as np

# Functions imported from pyaxon.functions and re-exported so that main.py
# can obtain both parameters and helper functions from a single import block.
# The __all__ list below documents every public name this module exposes.
from pyaxon.functions import (
    build_transport_field,
    coeff_beta_,
    com,
    conservative_upwind_advection,
    container,
    h,
    init_cell,
    init_growth_factor,
    rasterized_line,
    smoothing,
    unit_vector,
)

__all__ = [
    # re-exported helpers
    "build_transport_field",
    "coeff_beta_",
    "com",
    "conservative_upwind_advection",
    "container",
    "h",
    "init_cell",
    "init_growth_factor",
    "rasterized_line",
    "smoothing",
    "unit_vector",
    # grid / numerics
    "L",
    "dL",
    "dt",
    "stencil",
    "boundary_condition",
    "alpha_",
    # geometry
    "neuron_position",
    "neuron_radius",
    "gc_position",
    "gc_radius",
    "small_box_size",
    # NGF
    "D_ngf",
    "ngf_source_position",
    "ngf_source_value",
    # microtubules
    "D_mtb",
    "lambda_mtb",
    "mtb_source",
    "mtb_positions",
    # free mRNA
    "kf",
    "lambda_f",
    "gamma",
    "mf_source_position",
    "mf_source_value",
    # linked mRNA
    "kl",
    "lambda_l",
    "Mm",
    "chi_ml",
    "beta_everywhere",
    "beta_growth_cone",
    "beta_m",
    # growth cone / axon
    "chi",
    "alpha_p",
    # time integration
    "nstep",
    "tstep",
    "nprint",
    "print_period",
]

# ── Plotting defaults ─────────────────────────────────────────────────────────
kwargs = {"origin": "lower", "interpolation": "sinc", "cmap": "seismic"}

# ── Spatial discretisation ────────────────────────────────────────────────────
# alpha_ is used as the regularisation parameter in the container() potential;
# a small positive value prevents division by zero.
alpha_ = 10e-4

# Boundary mode for scipy.ndimage.convolve: 'constant' pads with zeros,
# which is equivalent to a Neumann (no-flux) boundary for these stencils.
boundary_condition = "constant"

# Physical cell width (grid units).  All length scales are in units of dL.
dL = 1.0

# 4th-order isotropic Laplacian stencil (5×5, cross-shaped).
# Prefactor 1/(12 dL²) gives O(dL⁴) accuracy on a square grid.
# The central coefficient is -60; the four axis-aligned neighbours are +16;
# the four diagonal neighbours at distance 2 are -1.
stencil = (1.0 / (12.0 * dL * dL)) * np.array(
    [
        [0, 0, -1, 0, 0],
        [0, 0, 16, 0, 0],
        [-1, 16, -60, 16, -1],
        [0, 0, 16, 0, 0],
        [0, 0, -1, 0, 0],
    ]
)

# ── Domain geometry ───────────────────────────────────────────────────────────
# L = [Nx, Ny]: number of grid cells in the row and column directions.
# The axon is expected to grow along the column (y) axis.
L = np.array([40, 100])
L0 = np.array([0, 0])  # lower-left corner (unused legacy variable)

# ── Neuron geometry ───────────────────────────────────────────────────────────
# Soma is placed at the left of the domain.  The growth cone starts at the
# soma surface (neuron_position + neuron_radius in the column direction).
neuron_position = np.array([20, 20])
neuron_radius = 15  # soma disc radius (grid cells)
gc_radius = 3  # growth-cone disc radius (grid cells)
gc_position = np.array([neuron_position[0], neuron_position[1] + neuron_radius])

# Half-width of the local computation box used around the growth cone for
# gradient / chemotaxis calculations.  Should be > gc_radius.
small_box_size = 15

# ── NGF source ────────────────────────────────────────────────────────────────
# The NGF point source is placed near the right boundary to act as a distant
# chemoattractant target for the growth cone.
ngf_source_position = (30, 90)  # (row, col) – target the cone navigates toward
ngf_source_value = 10.0  # Dirichlet value injected at the source each step
D_ngf = 100.0  # fast diffusion so NGF gradient spans the domain

# ── Free mRNA source ──────────────────────────────────────────────────────────
# Free mRNA is continuously produced at the soma centre, mimicking ribosomal
# synthesis in the cell body.
mf_source_position = (neuron_position[0], neuron_position[1])
mf_source_value = 10.0  # Dirichlet injection value

# ── Microtubules ──────────────────────────────────────────────────────────────
# The MT rail is maintained by pinning all tracked cells to mtb_source each step
# (see main.py: np.put(mtb, mtb_positions, mtb_source)).
# mtb_source = 1.0 corresponds to full saturation of the rail.
# Values < 1 leave some capacity for additional polymerisation.
mtb_source = 1.0  # value pinned to all rail cells each timestep
mtb_positions = []  # flat indices of rail cells (populated at runtime)
lambda_mtb = 0.1  # linear MT decay rate (dt * lambda_mtb ≪ 1 for stability)
D_mtb = 0.1  # MT diffusion coefficient (small: MTs spread slowly)

# ── mRNA parameters ───────────────────────────────────────────────────────────

# --- Free mRNA (mf) ---
kf = 1.0  # diffusion coefficient  (larger → faster soma-to-tip spread)
lambda_f = 0.0001  # linear decay rate      (very slow; mf is long-lived)
gamma = 10.0  # binding rate  mf + MT → ml
# larger γ → more of the pool shifts to linked form

# --- Linked mRNA (ml) ---
kl = 0.5  # diffusion coefficient  (smaller than kf; ml is less mobile)
lambda_l = 0.0001  # linear decay rate

# Maximum MT rail capacity.  When ml → Mm the rail is saturated and further
# binding (coeff_gamma) and motor transport (coeff_ml) both saturate.
Mm = 1.0

# Motor transport coefficient for the active flux  -∇·(chi_ml * coeff_ml * ml * mtb * v_m).
# CFL stability requires:  chi_ml * dt / dL < 1
# With dt=0.001 and dL=1.0, chi_ml < 1000.  In practice chi_ml=100 is safe.
chi_ml = 100.0

# Spatially varying release rate  ml → mf:
#   beta_everywhere  – baseline release along the whole axon shaft
#   beta_growth_cone – enhanced release at the growth cone tip
# The large contrast (0.001 vs 10) ensures that free mRNA accumulates
# preferentially at the tip rather than being released prematurely along the shaft.
beta_everywhere = 0.001
beta_growth_cone = 10.0
beta_m = np.zeros(L)  # full release-rate field; updated each step in main.py

# ── Growth-cone chemotaxis ────────────────────────────────────────────────────
# Chemotactic speed law:  v = chi * mf_GC/(1+mf_GC) * ∇NGF̂
# The Hill-type saturation factor mf_GC/(1+mf_GC) bounds the speed at chi.
# chi was reduced from 10 → 3 to prevent the cone from pulling ahead of the
# shaft faster than alpha_p can refill it.
chi = 3.0

# ── Axon shaft proliferation ──────────────────────────────────────────────────
# The proliferative source in the psi (Cahn-Hilliard) equation:
#   alpha_p * phi * ngf * psi * (1 - psi)
# drives new axon material to grow at the interface between phi and psi.
# alpha_p was increased from 10 → 15 so the shaft keeps up with the cone.
alpha_p = 15.0

# ── Time integration ──────────────────────────────────────────────────────────
nstep = 0  # current step (advanced inside the while loop in main.py)
tstep = 200000  # total steps to run  (physical time T = tstep * dt = 200)
dt = 0.001  # time step size
print_period = 10000  # save a snapshot every this many steps
nprint = print_period  # counter; snapshot saved when nprint >= print_period
