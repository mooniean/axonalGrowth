Introduction
============

``pyaxon`` is a computational implementation of a multi-phase-field model for
axonal growth and mRNA transport.

The model describes three coupled subsystems on a 2-D rectangular grid:

1. **Phase-field morphology** — a growth cone (φ) migrates chemotactically
   toward a Nerve Growth Factor (NGF) source while a Cahn-Hilliard equation
   extends the axon shaft (ψ) behind it.

2. **Microtubule scaffold** — a directed microtubule rail is seeded from the
   soma to the current growth-cone position and is maintained throughout the
   simulation.

3. **mRNA transport** — free mRNA (m\ :sub:`f`) produced at the soma binds to
   the microtubule rail, is actively transported as linked mRNA (m\ :sub:`l`)
   toward the tip by molecular motors, and is released as free mRNA at the
   growth cone.  The local concentration of free mRNA at the tip controls the
   chemotactic speed of the cone.

See :doc:`tutorial` for a step-by-step guide to running simulations, plotting
snapshots, and generating animated GIFs.
