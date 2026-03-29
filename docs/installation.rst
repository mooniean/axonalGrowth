Installation
============

Requirements
------------

============  ================  ================================
Tool          Minimum version   Purpose
============  ================  ================================
Python        3.9               runtime
uv_           0.4               dependency & environment manager
numpy         1.14.3            numerical arrays
scipy         1.3.0             convolutions
matplotlib    3.5               plotting
Pillow        9.0               GIF encoding
============  ================  ================================

.. _uv: https://github.com/astral-sh/uv

Development extras (``ruff``, ``pytest``, ``coverage``) and documentation
extras (``sphinx``, ``nbsphinx``) are declared as dependency groups in
``pyproject.toml`` and are installed automatically by ``uv``.


Setting up the environment
--------------------------

1. **Clone the repository**

   .. code-block:: bash

      git clone https://github.com/mooniean/axonalGrowth.git
      cd axonalGrowth

2. **Create the virtual environment and install all dependencies**

   .. code-block:: bash

      make install           # runtime + dev tools (ruff, pytest, …)
      make install-docs      # additionally installs Sphinx and nbsphinx

   Under the hood this runs ``uv sync --group dev``, which creates ``.venv/``
   and pins every dependency from ``uv.lock``.

3. **Activate the environment** *(optional — uv handles activation automatically)*

   .. code-block:: bash

      source .venv/bin/activate   # macOS / Linux

Verifying the installation
--------------------------

Run the test suite to confirm everything works:

.. code-block:: bash

   make test

All four tests should pass with no errors.

