# xiaoyu-cytoflow

`xiaoyu-cytoflow` is a self-contained flow-cytometry analysis package for
96-well plate FCS workflows. It packages Xiaoyu's analysis helpers together
with the small Cytoflow runtime subset they need, so users do not have to
install the unmaintained upstream `cytoflow` package separately.

The package focuses on the workflow used in `xiaoyu_CF.py`: importing plate
FCS files, attaching sample and concentration metadata, removing outliers,
gating cell populations and doublets, plotting channel distributions, subsetting
by well/row/column, and exporting 96-well median signal tables.

## Why This Exists

The original analysis code depended on a locally patched Cytoflow install that
was difficult to recreate on modern Python. This repository vendors the working
pieces of that local runtime and updates them for Python 3.12, making the
analysis easier to share, test, and publish.

The vendored Cytoflow subset is intentionally small. It keeps the pieces used by
the current workflow:

- FCS import and experiment storage
- range, polygon, and density gates
- autofluorescence and linear bleedthrough correction
- histogram, scatter, and density views
- `linear` and `log` scales

The compiled Cytoflow logicle extension was removed, so this package is pure
Python.

## Installation

The recommended development setup uses mamba:

```bash
mamba env create -f environment.yml
mamba activate xiaoyu-cytoflow
```

The environment installs this repository in editable mode, so local changes are
immediately visible.

In an existing Python 3.12 environment, install with:

```bash
pip install -e .
```

After the package is published to PyPI, users will be able to install it with:

```bash
pip install xiaoyu-cytoflow
```

## Quick Start

```python
import matplotlib.pyplot as plt
import xiaoyu_CF as flow

expr = flow.xiaoyu_Expr("path/to/fcs_directory")
expr.add_conc_condition(1e-5, dilu_dir="down")

without_outliers = flow.gate_outliers(expr.expr)
cell_population = flow.auto_gate_FSC_SSC(without_outliers, if_plot=True)
single_cells = flow.auto_gate_FSC_A_H(cell_population, if_plot=True, keep=0.7)

flow.channel_histogram(flow.subset_by_num(single_cells, "12"), huefacet="conc")
expr.median_96well(single_cells, channel="Alexa 647-A")

plt.show()
```

## Documentation

See [PROJECT_GUIDE.md](PROJECT_GUIDE.md) for the full project documentation,
including data layout assumptions, API reference, example workflows, testing,
packaging, and maintenance notes.

## Test Data

The `test/` directory may contain local experiment data and scratch analysis
files. The included smoke test exercises the packaged workflow against the
provided example data:

```bash
MPLCONFIGDIR=/tmp/mplconfig python test/package_smoke_test.py
```

## Licensing

This project includes a vendored subset of Cytoflow, which is GPL-licensed.
The package is therefore distributed under GPL-compatible terms. See
[LICENSE.txt](LICENSE.txt) and [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
