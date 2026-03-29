![pyaxon](https://github.com/mooniean/axonalGrowth/workflows/pyaxon/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/axonalgrowth/badge/?version=latest)](https://axonalgrowth.readthedocs.io/en/latest/?badge=latest)
[![Python 3](https://pyup.io/repos/github/mooniean/axonalGrowth/python-3-shield.svg)](https://pyup.io/repos/github/mooniean/axonalGrowth)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# pyaxon – Phase-field model for axonal growth and mRNA transport

`pyaxon` simulates axonal growth using a multi-phase-field model.
A growth cone navigates a Nerve Growth Factor (NGF) gradient while free and
microtubule-linked mRNA are transported along the growing axon shaft.

---

## Table of contents

1. [Requirements](#requirements)
2. [Installation](#installation)
3. [Project structure](#project-structure)
4. [Running a simulation](#running-a-simulation)
5. [Key parameters](#key-parameters)
6. [Plotting a snapshot](#plotting-a-snapshot)
7. [Generating a GIF](#generating-a-gif)
8. [Running the tests](#running-the-tests)
9. [Building the documentation](#building-the-documentation)
10. [Model equations](#model-equations)

---

## Requirements

| Tool | Minimum version | Purpose |
|---|---|---|
| Python | 3.9 | runtime |
| [uv](https://github.com/astral-sh/uv) | 0.4 | dependency & environment manager |
| numpy | 1.14.3 | numerical arrays |
| scipy | 1.3.0 | convolutions |
| matplotlib | 3.5 | plotting |
| Pillow | 9.0 | GIF encoding |

Development extras (`ruff`, `pytest`, `coverage`) and documentation extras
(`sphinx`) are managed through `pyproject.toml` dependency groups.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/mooniean/axonalGrowth.git
cd axonalGrowth
```

### 2. Create the virtual environment and install all dependencies

```bash
make install          # runtime + dev tools (ruff, pytest, …)
make install-docs     # additionally installs Sphinx and nbsphinx
```

Under the hood this runs `uv sync --group dev`, which creates `.venv/` and
installs every dependency from the lock file.

### 3. Activate the environment (optional – uv handles it automatically)

```bash
source .venv/bin/activate
```

---

## Project structure

```
axonalGrowth/
├── pyaxon/
│   ├── parameters.py   # all physical & numerical parameters (edit here)
│   ├── functions.py    # helper functions (phase-field, transport, geometry)
│   ├── main.py         # time-integration loop – run this to simulate
│   └── plot_npz.py     # CLI tool for plots and GIFs
├── tests/
│   └── general.py      # unit tests
├── docs/               # Sphinx documentation source
├── pyproject.toml      # project metadata, dependencies, ruff & pytest config
├── Makefile            # developer shortcuts
└── README.md
```

---

## Running a simulation

All parameters are in `pyaxon/parameters.py`.
Edit that file, then run the integration loop directly:

```bash
uv run python pyaxon/main.py
```

Or with the virtual environment activated:

```bash
python pyaxon/main.py
```

The simulation prints the current step number every `print_period` steps and
writes compressed NumPy archives (`.npz`) to the output directory defined by
`prefix_` at the top of `main.py` (default: `deb_4_/`).

Each archive contains the fields:

| Key | Description |
|---|---|
| `phi` | growth-cone phase field |
| `psi` | neuron body / axon shaft phase field |
| `ngf` | Nerve Growth Factor concentration |
| `mf` | free mRNA concentration |
| `ml` | MT-linked mRNA concentration |
| `mtb` | microtubule density |
| `v_m` | motor-velocity field (2 × Nx × Ny) |

### Typical run time

With the default settings (`tstep = 200 000`, `dt = 0.001`, grid `40 × 100`)
a full simulation takes roughly **10–30 minutes** on a modern laptop CPU.

---

## Key parameters

All of these live in `pyaxon/parameters.py` and are the only file you need to
change for a parameter study.

| Parameter | Default | Description |
|---|---|---|
| `L` | `[40, 100]` | Grid size `[Nx, Ny]` |
| `dt` | `0.001` | Time step. Decrease if the simulation goes unstable. |
| `tstep` | `200000` | Total number of time steps |
| `print_period` | `10000` | Save a snapshot every this many steps |
| `chi` | `3.0` | Chemotactic speed coefficient |
| `alpha_p` | `15.0` | Axon shaft proliferation rate |
| `chi_ml` | `100.0` | MT motor transport strength. CFL: `chi_ml * dt / dL < 1` |
| `gamma` | `10.0` | Binding rate free mRNA → linked mRNA |
| `beta_growth_cone` | `10.0` | mRNA release rate at the growth cone tip |
| `ngf_source_position` | `(30, 90)` | Location of the NGF point source |
| `mtb_source` | `1.0` | Microtubule density pinned to the rail each step |
| `prefix_` *(in main.py)* | `'deb_4_'` | Output directory for `.npz` files |

---

## Plotting a snapshot

Use the `plot-npz` CLI (installed as an entry-point) or call the module
directly.

### Plot a single `.npz` file interactively

```bash
uv run plot-npz deb_4_/deb_4_10000.npz
```

### Save the plot to an image file without opening a window

```bash
uv run plot-npz deb_4_/deb_4_10000.npz --save snapshot_10000.png --no-show
```

### Choose a different colormap

```bash
uv run plot-npz deb_4_/deb_4_10000.npz --cmap plasma
```

The figure shows all saved fields (φ, ψ, NGF, m_f, m_l, mtb) plus the
magnitude of the motor-velocity field v_m in an 2 × 4 grid.

---

## Generating a GIF

Pass an output directory (or a glob pattern) together with `--gif`:

### From an entire output directory

```bash
uv run plot-npz deb_4_ --gif growth.gif
```

### Control the frame rate (default 6 fps)

```bash
uv run plot-npz deb_4_ --gif growth.gif --fps 10
```

### Choose a colormap

```bash
uv run plot-npz deb_4_ --gif growth.gif --fps 8 --cmap inferno
```

### From a glob pattern

```bash
uv run plot-npz "deb_4_/deb_4_*.npz" --gif growth.gif --fps 6
```

Frames are sorted by the trailing step number in the filename, so the GIF
always plays in chronological order.

---

## Running the tests

```bash
make test          # run the full test suite
make cov           # run tests with branch coverage
make cov-report    # open the HTML coverage report in the browser
```

Or directly:

```bash
uv run pytest tests/ -v
```

---

## Building the documentation

```bash
make docs
```

The HTML output is written to `docs/_build/html/index.html`.
Open it in your browser:

```bash
open docs/_build/html/index.html          # macOS
xdg-open docs/_build/html/index.html     # Linux
```

---

## Model equations

### Growth cone (Allen-Cahn + chemotaxis)

The growth cone is described by the order parameter φ.
Its velocity toward the NGF source depends on the local free-mRNA
concentration $m_f$ at the cone tip:

$$\vec{v}_\mathrm{GC} = \chi(m_f)\,\frac{\nabla[\mathrm{NGF}]}{||\nabla[\mathrm{NGF}]||}\,,\qquad \chi(m_f) = \chi_0\,\frac{m_f}{1+m_f}$$

### Axon shaft (Cahn-Hilliard + proliferation)

The neuron body and axon shaft are described by the order parameter ψ,
governed by the Cahn-Hilliard equation with an NGF-dependent proliferative
source that extends the shaft behind the advancing cone:

![equation](https://latex.codecogs.com/gif.latex?%5Cpartial_t%20%5Cpsi%20%3D%20-M%20%5Cnabla%5E2%20%5Cleft%28%208%20%5Cpsi%20%281-%5Cpsi%29%5Cleft%28%5Cpsi%20-%20%5Cfrac%7B1%7D%7B2%7D%5Cright%29%20&plus;%20%5Cvarepsilon%5E2%20%5Cnabla%5E2%20%5Cpsi%20%5Cright%29%20&plus;%20%5Calpha_p%20%5Cpsi%20%5Cphi%20%281-%5Cphi%29%20%5B%5Cmathrm%7BNGF%7D%5D%20%5C%2C%20%2C)

### mRNA transport

Two mRNA phenotypes are coupled by binding/release reactions and
active motor transport along the microtubule (MT) rail:

**Free mRNA** $m_f$ (diffuses, binds to MT rail, released at the tip):

![equation](https://latex.codecogs.com/gif.latex?%5Cpartial_t%20m_f%20%3D%20k_f%5Cnabla%5E2m_f%20-%20%5Cgamma%20m_f%20m_t%20%5Cleft%28%5Cfrac%7BM_m-m_l%7D%7BM_m%7D%5Cright%29%20-%5Clambda_f%20m_f%20%5C%5C&plus;%20%5Cfrac%7B%5Cbeta_m%20%5Cleft%281-m_t%5Cright%29m_l%7D%7Bm_t%7D%20&plus;%20%5Cnabla%5E2%20%5CXi%5E%7B%5Cnatural%7D_C%5Cleft%5Bm_f%2C%20%5Cpsi%20%5Cright%5D%20%5C%2C%20%2C)

**Linked mRNA** $m_l$ (diffuses, actively transported soma → tip by molecular motors):

![equation](https://latex.codecogs.com/gif.latex?%5Cpartial_t%20m_l%20%3D%20k_l%20%5Cnabla%5E2%20m_l%20-%5Cnabla%5Ccdot%5Cleft%28%5Cleft%28%5Cfrac%7BM_m-m_l%7D%7BM_m%7D%5Cright%29m_l%20m_t%20%5Cmathbf%7Bv%7D_m%5Cright%29-%5Clambda_l%20m_l%20%5C%5C%20&plus;%20%5Cgamma%20m_f%20m_t%20%5Cleft%28%5Cfrac%7BM_m-m_l%7D%7BM_m%7D%5Cright%29%20-%20%5Cfrac%7B%5Cbeta_m%20%5Cleft%281-m_t%5Cright%29m_l%7D%7Bm_t%7D%20&plus;%20%5Cnabla%5E2%20%5CXi%5E%7B%5Cnatural%7D_C%5Cleft%5Bm_l%2C%20%5Cpsi%20%5Cright%5D%20%5C%2C%20%2C)

The divergence term $-\nabla\cdot(\cdots)$ is discretised with a
**conservative first-order upwind scheme** along an ordered microtubule track
that is extended step by step as the growth cone advances.

### Arrival condition

Once the growth cone reaches the NGF point source the simulation switches to
a dormant phase: chemotaxis (χ), shaft proliferation (α_p), and NGF
re-injection are all set to zero so the cone stops and the shaft stabilises.
