Tutorial
========

This page walks through a complete workflow: editing parameters, running a
simulation, visualising a snapshot, and generating an animated GIF.


Project layout
--------------

.. code-block:: text

   axonalGrowth/
   ├── pyaxon/
   │   ├── parameters.py   ← edit this to change the simulation
   │   ├── functions.py    ← helper functions (do not edit unless extending the model)
   │   ├── main.py         ← time-integration loop (run this)
   │   └── plot_npz.py     ← CLI for plots and GIFs
   ├── tests/
   │   └── general.py
   ├── docs/
   ├── pyproject.toml
   ├── Makefile
   └── README.md


Step 1 – Choose your parameters
--------------------------------

Open ``pyaxon/parameters.py``.  Every physical and numerical knob is defined
there with an inline comment explaining its role.  The most commonly adjusted
parameters are:

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Parameter
     - Default
     - Description
   * - ``tstep``
     - ``200000``
     - Total number of time steps.
   * - ``print_period``
     - ``10000``
     - Save a ``.npz`` snapshot every this many steps.
   * - ``chi``
     - ``3.0``
     - Chemotactic speed coefficient.  Increase to make the cone move faster.
   * - ``alpha_p``
     - ``15.0``
     - Axon shaft proliferation rate.  Increase if the shaft lags behind the cone.
   * - ``chi_ml``
     - ``100.0``
     - Motor transport strength for linked mRNA.
       **CFL stability:** ``chi_ml * dt / dL < 1``.
   * - ``gamma``
     - ``10.0``
     - Binding rate free mRNA → linked mRNA.
   * - ``beta_growth_cone``
     - ``10.0``
     - mRNA release rate at the growth cone tip.
   * - ``ngf_source_position``
     - ``(30, 90)``
     - ``(row, col)`` location of the NGF point source — the cone navigates toward it.
   * - ``dt``
     - ``0.001``
     - Time step.  Halve this if the simulation becomes numerically unstable.

.. note::

   ``prefix_`` at the top of ``main.py`` sets the output directory name
   (default ``deb_4_/``).  Change it before each run to avoid overwriting
   previous results.


Step 2 – Run the simulation
----------------------------

From the repository root:

.. code-block:: bash

   uv run python pyaxon/main.py

Or, with the virtual environment activated:

.. code-block:: bash

   python pyaxon/main.py

The terminal prints the current step number every ``print_period`` steps::

   0
   10000
   20000
   ...
   Growth cone reached NGF source at step 87500
   ...

Compressed NumPy archives are written to the output directory:

.. code-block:: text

   deb_4_/
   ├── deb_4_0.npz
   ├── deb_4_10000.npz
   ├── deb_4_20000.npz
   └── …

Each archive contains the arrays:

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Key
     - Description
   * - ``phi``
     - Growth-cone phase field (0 = outside, 1 = inside).
   * - ``psi``
     - Neuron body / axon shaft phase field.
   * - ``ngf``
     - Nerve Growth Factor concentration.
   * - ``mf``
     - Free mRNA concentration.
   * - ``ml``
     - MT-linked mRNA concentration.
   * - ``mtb``
     - Microtubule density along the axon track.
   * - ``v_m``
     - Motor-velocity field, shape ``(2, Nx, Ny)``.

**Typical run time:** 10–30 minutes on a modern laptop CPU at the default
settings (``tstep = 200 000``, grid ``40 × 100``).


Step 3 – Plot a snapshot
-------------------------

The ``plot-npz`` entry-point is installed automatically by ``uv``.

**Interactive plot**

.. code-block:: bash

   uv run plot-npz deb_4_/deb_4_10000.npz

**Save to a PNG without opening a window**

.. code-block:: bash

   uv run plot-npz deb_4_/deb_4_10000.npz --save snapshot_10000.png --no-show

**Use a different colormap**

.. code-block:: bash

   uv run plot-npz deb_4_/deb_4_10000.npz --cmap plasma

The output figure is a 2 × 4 grid showing φ, ψ, NGF, m_f, m_l, mtb, and the
magnitude of v_m.

You can also use the Python API directly:

.. code-block:: python

   from pyaxon.plot_npz import plot_npz

   plot_npz("deb_4_/deb_4_10000.npz", save_path="snapshot.png", show=False)


Step 4 – Generate an animated GIF
-----------------------------------

Pass the output directory (or a glob pattern) together with ``--gif``.

**From the whole output directory (all snapshots in chronological order)**

.. code-block:: bash

   uv run plot-npz deb_4_ --gif growth.gif

**Control the frame rate** (default 6 fps)

.. code-block:: bash

   uv run plot-npz deb_4_ --gif growth.gif --fps 10

**Choose a colormap**

.. code-block:: bash

   uv run plot-npz deb_4_ --gif growth.gif --fps 8 --cmap inferno

**Select a subset of frames with a glob pattern**

.. code-block:: bash

   uv run plot-npz "deb_4_/deb_4_*.npz" --gif growth.gif --fps 6

Frames are sorted automatically by the trailing step number, so the animation
always plays in chronological order.

You can also call the Python API directly:

.. code-block:: python

   from pyaxon.plot_npz import collect_npz_files, make_gif_from_npz_files

   files = collect_npz_files("deb_4_")
   make_gif_from_npz_files(files, "growth.gif", fps=8, cmap="viridis")


Step 5 – Run the tests
-----------------------

.. code-block:: bash

   make test          # full test suite
   make cov           # tests + branch coverage report
   make cov-report    # open htmlcov/index.html in the browser

Or directly with pytest:

.. code-block:: bash

   uv run pytest tests/ -v


Step 6 – Build the Sphinx documentation
-----------------------------------------

.. code-block:: bash

   make docs

The HTML output is written to ``docs/_build/html/``.  Open
``docs/_build/html/index.html`` in your browser:

.. code-block:: bash

   open docs/_build/html/index.html       # macOS
   xdg-open docs/_build/html/index.html  # Linux


Makefile reference
------------------

.. code-block:: bash

   make help          # print all available targets

The most useful targets:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Target
     - Action
   * - ``make install``
     - Create ``.venv`` and install runtime + dev dependencies via ``uv sync``.
   * - ``make install-docs``
     - Also install Sphinx and nbsphinx.
   * - ``make check``
     - Run ``ruff check`` (read-only; exits non-zero on violations).
   * - ``make fix``
     - Auto-fix lint violations, then auto-format with ``ruff``.
   * - ``make fmt``
     - Auto-format source files with ``ruff format``.
   * - ``make test``
     - Run the full pytest test suite.
   * - ``make cov``
     - Run tests with branch coverage measurement.
   * - ``make build``
     - Build sdist + wheel into ``dist/`` using ``uv build``.
   * - ``make docs``
     - Build the Sphinx HTML documentation.
   * - ``make clean``
     - Remove caches and build artefacts (keeps ``.venv``).
   * - ``make all``
     - ``check`` + ``test`` — safe for CI, no auto-modification.

