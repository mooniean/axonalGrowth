"""
main.py – time integration loop for the axonal growth / mRNA transport model.

Overview
--------
Three interacting subsystems are evolved on a 2-D rectangular grid:

  1. Phase-field morphology
       phi  – growth cone  (Allen-Cahn + volume-constraint + chemotaxis)
       psi  – neuron body / axon shaft  (Cahn-Hilliard + NGF-driven proliferation)

  2. Chemical guidance and microtubule scaffold
       ngf  – Nerve Growth Factor  (diffusion + uptake by psi)
       mtb  – microtubule density  (diffusion + decay; seeded along the axon track)

  3. mRNA transport (two phenotypes)
       mf   – free mRNA  (diffusion + binding → ml; released from ml at the GC)
       ml   – MT-linked mRNA  (diffusion + active motor transport toward tip;
                                binding from mf; release → mf at the growth cone)

All fields are updated with an explicit forward-Euler scheme using time step dt.
Spatial derivatives use a 4th-order isotropic Laplacian stencil (see parameters.py).

Motor transport of ml uses a conservative first-order upwind flux on an ordered
microtubule track that is extended step by step as the growth cone advances.

Once the growth cone reaches the NGF source the simulation enters a dormant phase:
chemotaxis, axon growth (alpha_p), and NGF re-injection are switched off so that
the cone stops moving and the shaft stabilises.
"""

import os

import numpy as np
from scipy.ndimage import convolve

from pyaxon.parameters import (
    # ── Microtubules ──────────────────────────────────────────────────────────
    D_mtb,
    # ── NGF ───────────────────────────────────────────────────────────────────
    D_ngf,
    # ── Grid / numerics ───────────────────────────────────────────────────────
    L,
    Mm,
    alpha_,
    alpha_p,
    beta_everywhere,
    beta_growth_cone,
    beta_m,
    boundary_condition,
    # ── Functions (re-exported from pyaxon.functions via parameters) ──────────
    build_transport_field,
    # ── Growth cone / axon ────────────────────────────────────────────────────
    chi,
    chi_ml,
    coeff_beta_,
    com,
    conservative_upwind_advection,
    container,
    dL,
    dt,
    gamma,
    gc_position,
    gc_radius,
    h,
    init_cell,
    init_growth_factor,
    # ── Free mRNA ─────────────────────────────────────────────────────────────
    kf,
    # ── Linked mRNA ───────────────────────────────────────────────────────────
    kl,
    lambda_f,
    lambda_l,
    lambda_mtb,
    mf_source_position,
    mf_source_value,
    mtb_positions,
    mtb_source,
    # ── Geometry ──────────────────────────────────────────────────────────────
    neuron_position,
    neuron_radius,
    ngf_source_position,
    ngf_source_value,
    nprint,
    # ── Time integration ─────────────────────────────────────────────────────
    nstep,
    print_period,
    rasterized_line,
    small_box_size,
    smoothing,
    stencil,
    tstep,
    unit_vector,
)

# ── Output directory ──────────────────────────────────────────────────────────
prefix_ = "deb_5_"
os.makedirs(prefix_, exist_ok=True)

# ── Field allocation ──────────────────────────────────────────────────────────
# All fields share the same 2-D grid of shape L = [Nx, Ny].
ngf = np.zeros(L)  # Nerve Growth Factor concentration
phi = np.zeros(L)  # growth-cone phase field  (0 = outside, 1 = inside)
psi = np.zeros(L)  # neuron-body/axon phase field
mtb = np.zeros(L)  # microtubule density
mf = np.zeros(L)  # free mRNA concentration
ml = np.zeros(L)  # MT-linked mRNA concentration

# ── Initial conditions ────────────────────────────────────────────────────────
# Free mRNA starts uniformly inside the neuron soma (slightly smaller radius
# than psi so that mf is strictly interior).
mf = 1.0 * init_cell(mf, neuron_position, neuron_radius - 5)

# NGF is initialised as a random field outside the neuron body; it provides
# the chemo-attractive gradient that steers the growth cone.
ngf = init_growth_factor(ngf, neuron_position, neuron_radius)

# Growth cone and neuron body are initialised as filled discs.
phi = init_cell(phi, gc_position, gc_radius)
psi = init_cell(psi, neuron_position, neuron_radius)

# ── Microtubule track initialisation ─────────────────────────────────────────
# Seed an ordered list of grid points from the soma centre to the initial
# growth-cone position.  This gives ml a transport rail from the very start
# rather than waiting for the cone to advance before any delivery occurs.
track_points = [tuple(point) for point in rasterized_line(neuron_position, gc_position)]
# Convert track points to flat indices so np.put() can reinforce the track
# each timestep without an explicit loop.
mtb_positions[:] = [np.ravel_multi_index(point, mtb.shape) for point in track_points]
# Build the directed motor-velocity field v_m along the seeded track.
# support_radius=0 keeps the field strictly on the track skeleton (one cell wide).
v_m = build_transport_field(track_points, L, support_radius=0)

# ── Visualisation coordinate arrays (for streamplot etc.) ────────────────────
y = np.linspace(0, 40, 40)
x = np.linspace(0, 100, 100)

# ── Allen-Cahn smoothing of initial interfaces ────────────────────────────────
# Relax phi and psi to diffuse interfaces before the main loop.
# 2000 smoothing steps are sufficient to establish a well-resolved interface.
phi = smoothing(phi, 2000, stencil)  # growth cone
psi = smoothing(psi, 2000, stencil)  # neuron body / axon

# ── Pre-equilibration of mf and ngf ──────────────────────────────────────────
# Run a short pre-loop so that the free-mRNA and NGF fields are already spread
# inside/outside the neuron before the coupled dynamics start.
for _i in range(0, 100):
    # mf diffuses with coefficient 10 and is confined by the container potential.
    mf = mf + dt * convolve(10 * mf + container(mf, psi, 0.0001), stencil, mode=boundary_condition)
    # NGF diffuses freely during pre-equilibration (no cone/psi present yet).
    ngf = ngf + dt * D_ngf * convolve(ngf, stencil, mode=boundary_condition)

# ── Volume conservation target ────────────────────────────────────────────────
# V_target is the smoothed volume of the initial growth cone.  The Lagrange
# multiplier in the phi equation enforces phi.sum() ≈ V_target throughout.
V_target = h(phi).sum()

# ── Loop control variables ────────────────────────────────────────────────────
r_cm = com(phi, V_target, ndim=2)  # centre of mass of the growth cone
icalc = 0
ncalc = 100  # recompute r_cm and extend MT track every ncalc steps
# (previously 250; smaller value keeps the track more up-to-date)

# Arrival flag – latched to True once the cone overlaps the NGF source.
# Once set, chemotaxis, shaft growth, and NGF injection are all turned off.
source_reached = False
stop_radius = gc_radius  # arrival distance threshold (grid units)

# ─────────────────────────────────────────────────────────────────────────────
#  MAIN TIME-INTEGRATION LOOP
# ─────────────────────────────────────────────────────────────────────────────
while nstep <= tstep:
    # ── Volume / Lagrange multiplier ──────────────────────────────────────────
    v = h(phi).sum()
    # The Lagrange multiplier corrects the phase field to conserve volume:
    #   L(x) = V_target - v(t)   (uniform field, same sign everywhere)
    Lagrange = np.ones(L) * (V_target - v)

    # ── Arrival detection (latched) ───────────────────────────────────────────
    # Check whether the growth cone has reached the NGF point source.
    # Two equivalent criteria are combined (OR):
    #   1. The phase-field value at the source location exceeds 0.5
    #      (i.e., the cone interior overlaps the source cell).
    #   2. The Euclidean distance from the cone centre of mass to the source
    #      is within stop_radius (≈ one cone radius).
    if not source_reached:
        reached_now = (
            phi[ngf_source_position] > 0.5
            or np.linalg.norm(
                np.asarray(r_cm, dtype=float) - np.asarray(ngf_source_position, dtype=float)
            )
            <= stop_radius
        )
        if reached_now:
            source_reached = True
            print(f"Growth cone reached NGF source at step {nstep}")

    # ── Local window around the growth cone ───────────────────────────────────
    # All chemotactic gradients are computed on a small sub-domain centred on
    # the cone.  ic, jc are the clipped integer centre coordinates; clipping
    # to [1, L-2] guarantees the gradient window is always at least 2 cells wide.
    ic = int(np.clip(np.rint(r_cm[0]), 1, L[0] - 2))
    jc = int(np.clip(np.rint(r_cm[1]), 1, L[1] - 2))

    # Large box [il:ih, jl:jh] – used for gradient computation and phi update.
    il, icl = max(0, ic - small_box_size), max(0, ic - gc_radius)
    ih, ich = min(L[0], ic + small_box_size), min(L[0], ic + gc_radius)
    # Small box [icl:ich, jcl:jch] – used for the beta_m pattern and mf_GC mean.
    jl, jcl = max(0, jc - small_box_size), max(0, jc - gc_radius)
    jh, jch = min(L[1], jc + small_box_size), min(L[1], jc + gc_radius)

    # Safety guard: if the window collapses (e.g., cone leaves the domain or
    # r_cm becomes NaN after a numerical instability), raise an informative error
    # instead of letting np.gradient produce a cryptic exception.
    if (ih - il) < 2 or (jh - jl) < 2:
        raise RuntimeError(
            f"Gradient window too small at step {nstep}: r_cm={r_cm}, window=({il}:{ih}, {jl}:{jh})"
        )

    # ── Chemotactic gradients ─────────────────────────────────────────────────
    # Gradients are computed only on the local window to reduce cost and focus
    # the chemotactic signal on the cone region.
    grad_phi = np.gradient(phi[il:ih, jl:jh])
    grad_ngf = np.gradient(ngf[il:ih, jl:jh])
    # Normalise the NGF gradient to a unit vector so that chi is the only
    # parameter controlling speed magnitude.
    grad_ngf = unit_vector(grad_ngf, norm=np.linalg.norm(grad_ngf, axis=0))

    # ── mRNA reaction / transport coefficients ────────────────────────────────
    # Clamp mtb to [0, 1] to prevent saturated values from destabilising the
    # reaction terms (mtb_source = 1 by default; diffusion can create overshoots).
    mtb_effective = np.clip(mtb, 0.0, 1.0)

    # coeff_ml = (Mm - ml) / Mm  is the fractional free capacity of the MT rail.
    # When ml → Mm the rail is fully occupied and binding/transport saturates.
    coeff_ml = np.clip((Mm - ml) / Mm, 0.0, 1.0)

    # beta_m is the spatially varying mRNA release rate:
    #   beta_everywhere  – small baseline release along the axon shaft
    #   beta_growth_cone – large release specifically at the growth cone tip,
    #                      producing the local spike of free mRNA observed there.
    beta_m[:, :] = beta_everywhere
    beta_m[icl:ich, jcl:jch] = beta_growth_cone

    # coeff_beta: effective release rate  ml → mf
    #   formula: (1/mtb) * beta_m * (1 - mtb) * ml
    #   The factor (1/mtb)*(1-mtb) vanishes when mtb → 1 (fully occupied rail)
    #   so release is suppressed where microtubules are saturated.
    coeff_beta = coeff_beta_(mtb_effective) * beta_m * (1 - mtb_effective) * ml

    # coeff_gamma: effective binding rate  mf → ml
    #   formula: gamma * mf * mtb * (Mm - ml) / Mm
    #   Binding is proportional to free mRNA, local MT density, and available
    #   capacity on the rail.
    coeff_gamma = gamma * mf * mtb_effective * coeff_ml

    # Motor transport of linked mRNA toward the growth cone tip.
    # The transported quantity is:  chi_ml * (Mm-ml)/Mm * ml * mtb
    # using a conservative first-order upwind flux -div(flux * v_m).
    # The upwind scheme is stable for |chi_ml * v_m * dt / dL| < 1
    # (Courant–Friedrichs–Lewy condition) and exactly conserves total ml.
    transport_ml = conservative_upwind_advection(coeff_ml * ml * mtb_effective, chi_ml * v_m, dL=dL)

    # ── Chemotactic velocity of the growth cone ───────────────────────────────
    grad_ngf_norm = np.linalg.norm(grad_ngf, axis=0, keepdims=True)
    unit_grad_ngf = unit_vector(grad_ngf, grad_ngf_norm)
    # chem: projection of grad(phi) onto the NGF gradient direction.
    # Positive values pull the cone toward increasing NGF.
    chem = grad_phi[0] * unit_grad_ngf[0] + grad_phi[1] * unit_grad_ngf[1]

    # Mean free-mRNA concentration inside the growth-cone core box.
    # Using mean (not sum) makes the chemotactic speed depend on local
    # concentration rather than total mRNA mass in the box, preventing
    # runaway acceleration as more mRNA accumulates.
    mf_GC = mf[icl:ich, jcl:jch].mean()

    # Gate chemotaxis and axon proliferation after the cone reaches its target.
    chi_eff = 0.0 if source_reached else chi
    alpha_p_eff = 0.0 if source_reached else alpha_p

    # Saturating speed law:  v = chi_eff * mf_GC/(1+mf_GC) * chem
    # The Hill-type factor mf_GC/(1+mf_GC) saturates at chi_eff, so the cone
    # speed is bounded even when tip mRNA is large.
    chi_ = chi_eff * (mf_GC / (1.0 + mf_GC)) * chem

    # ── Source / boundary conditions ─────────────────────────────────────────
    # Constant Dirichlet injection of free mRNA at the soma centre.
    mf[mf_source_position] = mf_source_value

    # NGF point source: injected only while the cone is still growing.
    # After arrival the source is turned off so that NGF diffuses away and
    # the residual gradient can no longer drive migration.
    if not source_reached:
        ngf[ngf_source_position] = ngf_source_value
    # Uncomment the line below to use a line source instead of a point source:
    # ngf[:, 80] = 10.0

    # Reinforce the MT track: set all tracked grid cells to mtb_source each step.
    # This keeps the microtubule rail from decaying away due to the lambda_mtb term.
    np.put(mtb, mtb_positions, mtb_source)

    # ── Forward-Euler field updates ───────────────────────────────────────────
    # All six fields are updated simultaneously from the values at step n.
    # The tuple unpacking on the left ensures no field reads its own new value.

    phi[il:ih, jl:jh], psi, mtb, ngf, mf, ml = (
        # phi – growth-cone Allen-Cahn equation with volume constraint and
        #        chemotactic forcing (-chi_ drives the cone toward NGF).
        # Only the local window is updated for efficiency; the rest of phi is
        # zero and not physically active.
        phi[il:ih, jl:jh]
        + dt
        * (
            -chi_
            + convolve(phi[il:ih, jl:jh], stencil, mode="constant")
            + phi[il:ih, jl:jh]
            * (1.0 - phi[il:ih, jl:jh])
            * ((phi[il:ih, jl:jh] - 0.5) + Lagrange[il:ih, jl:jh])
        ),
        # psi – neuron/axon Cahn-Hilliard equation.
        # First term: standard CH free-energy relaxation (maintains sharp interface).
        # Second term: NGF-dependent proliferative source at the phi/psi interface,
        #              which extends the axon shaft behind the advancing cone.
        #              Turned off (alpha_p_eff = 0) once the target is reached.
        psi
        + dt
        * (
            -convolve(
                8 * (psi * (1.0 - psi) * (psi - 0.5))
                + convolve(psi, stencil, mode=boundary_condition),
                stencil,
                mode=boundary_condition,
            )
            + alpha_p_eff * phi * ngf * psi * (1.0 - psi)
        ),
        # mtb – microtubule density: diffusion + decay.
        # The container term confines mtb within the axon body (psi ≈ 1).
        mtb
        + dt
        * (
            D_mtb
            * convolve(mtb + 2 * container(mtb, psi, alpha_), stencil, mode=boundary_condition)
            - lambda_mtb * mtb
        ),
        # ngf – Nerve Growth Factor: diffusion + consumption by the axon (psi * ngf).
        ngf + dt * (D_ngf * convolve(ngf, stencil, mode=boundary_condition) - psi * ngf),
        # mf – free mRNA:
        #   + diffusion  (kf∇²mf)
        #   + container confinement inside the axon
        #   - decay  (lambda_f * mf)
        #   - binding to microtubules  (coeff_gamma)
        #   + release from linked mRNA (coeff_beta)
        mf
        + dt
        * (
            convolve(
                kf * mf + 2 * container(mf, psi, alpha_),
                stencil,
                mode=boundary_condition,
            )
            - lambda_f * mf
            - coeff_gamma
            + coeff_beta
        ),
        # ml – MT-linked mRNA:
        #   + diffusion  (kl∇²ml)
        #   + container confinement
        #   - decay  (lambda_l * ml)
        #   + binding from free mRNA (coeff_gamma)
        #   - release to free mRNA   (coeff_beta)
        #   + active motor transport along the MT track (transport_ml)
        #     implements: -∇·[ chi_ml * (Mm-ml)/Mm * ml * mtb * v_m ]
        ml
        + dt
        * (
            convolve(
                kl * ml + 2 * container(ml, psi, alpha_),
                stencil,
                mode=boundary_condition,
            )
            - lambda_l * ml
            + coeff_gamma
            - coeff_beta
            + transport_ml
        ),
    )

    # ── Periodic centre-of-mass and MT track update ───────────────────────────
    if icalc >= ncalc:
        icalc = 0

        new_r_cm = com(phi, v, 2)
        new_r_cm = (np.int64(new_r_cm[0]), np.int64(new_r_cm[1]))

        # Extend the MT track from the last known position to the new cone centre.
        # rasterized_line returns all integer grid points along the segment;
        # we append only the new points (segment[1:]) to avoid duplicates.
        segment = rasterized_line(r_cm, new_r_cm)
        if (not source_reached) and segment.shape[0] > 1:
            track_points.extend(tuple(point) for point in segment[1:])
            # Rebuild flat-index list for np.put() and recompute the full v_m field.
            mtb_positions[:] = [np.ravel_multi_index(point, mtb.shape) for point in track_points]
            v_m = build_transport_field(track_points, L, support_radius=0)

        r_cm = new_r_cm

    # ── I/O ───────────────────────────────────────────────────────────────────
    if nprint >= print_period:
        print(nstep)
        np.savez(
            prefix_ + "/" + prefix_ + str(nstep) + ".npz",
            psi=psi,
            mtb=mtb,
            v_m=v_m,
            ml=ml,
            mf=mf,
            phi=phi,
            ngf=ngf,
        )
        nprint = 0

    icalc += 1
    nstep += 1
    nprint += 1
