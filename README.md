# xiaoyu-cytoflow

`xiaoyu-cytoflow` is a self-contained flow-cytometry analysis package for
96-well plate FCS workflows. It packages Xiaoyu's analysis helpers together
with the slim Cytoflow runtime subset they need — no separate upstream
`cytoflow` installation required.

Includes FCS import, metadata assignment, outlier removal, cell/doublet
gating, histograms, density plots, well/row/column subsetting, and 96-well
median export.

The vendored Cytoflow subset includes:

- FCS import and experiment storage
- range, polygon, and density gates
- autofluorescence and linear bleedthrough correction
- histogram, scatter, and density views
- `linear` and `log` scales

The compiled logicle extension was removed, so this package is pure Python.

## Installation

```bash
mamba env create -f environment.yml
mamba activate xiaoyu-cytoflow
```

Or in an existing Python >=3.9 environment:

```bash
pip install -e .
```

## Quick Start

```python
import matplotlib.pyplot as plt
import xiaoyu_CF as flow

expr = flow.xiaoyu_Expr("path/to/fcs_directory")
expr.add_conc_condition(1e-5, dilu_dir="down")

cells = flow.auto_gate_FSC_SSC(expr.expr, keep=0.6)
single_cells = flow.auto_gate_FSC_A_H(cells, keep=0.7)

flow.channel_histogram(single_cells, channel="Alexa 647-A", huefacet="conc")
expr.median_96well(single_cells, channel="Alexa 647-A")

plt.show()
```

## Documentation

- [PROJECT_GUIDE.md](PROJECT_GUIDE.md) — full API reference, workflows, testing, publishing
- [AGENTS.md](AGENTS.md) — quick guide for AI coding agents
- [llms.txt](llms.txt) — compact machine-readable project map

## Testing

```bash
MPLCONFIGDIR=/tmp/mplconfig python test/package_smoke_test.py
```

## License

GPL-compatible. See [LICENSE.txt](LICENSE.txt) and [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
