# Agent Guide for xiaoyu-cytoflow

This file is for AI coding agents, notebook assistants, and other tools helping
users write analysis code with `xiaoyu-cytoflow`.

## Package Identity

- PyPI package name: `xiaoyu-cytoflow`
- Main import: `import xiaoyu_CF as flow`
- The former `cytoflow` and `fcsparser` packages are now internal modules (prefixed with `_`)
- Python target: `>=3.9`
- Primary domain: 96-well plate flow-cytometry FCS analysis

## Internal Module Map

For codebase exploration, the vendored code lives under `src/xiaoyu_CF/`:

| Module | Contents |
|--------|----------|
| `__init__.py` | Public API: `xiaoyu_Expr`, all gating/plotting helpers |
| `_experiment.py` | `Experiment` class (core data structure) |
| `_operations.py` | `ImportOp`, `Tube`, `RangeOp`, `PolygonOp`, `DensityGateOp`, `AutofluorescenceOp`, `BleedthroughLinearOp` |
| `_views.py` | `HistogramView`, `ScatterplotView`, `DensityView`, base view classes |
| `_utility.py` | Scales, errors, algorithms, custom traits, docstring helpers |
| `_fcsparser.py` | FCS file parser (`parse` function)

## Fast Mental Model

Users typically have one directory of `.fcs` files named with well identifiers:

```text
A1 sample.fcs
A2 sample.fcs
...
H12 sample.fcs
```

The usual code pattern is:

```python
import matplotlib.pyplot as plt
import xiaoyu_CF as flow

expr = flow.xiaoyu_Expr("path/to/fcs_directory")
expr.add_conc_condition(1e-5, dilu_dir="down")

cells = flow.auto_gate_FSC_SSC(expr.expr, keep=0.6)
single_cells = flow.auto_gate_FSC_A_H(cells, keep=0.7)

flow.density_plot(single_cells, xchannel="FSC 488/10-A",
                  ychannel="Alexa 647-A", xscale="log", yscale="log")
expr.median_96well(single_cells, channel="Alexa 647-A")

plt.show()
```

`expr.expr` is a Cytoflow `Experiment`. Event data lives in
`expr.expr.data`, a pandas `DataFrame`.

## Public Helpers

Use these `xiaoyu_CF` helpers first:

- `xiaoyu_Expr(working_dir, merge_dirs=[])`
- `xiaoyu_Expr.add_conc_condition(ori_conc, dilu_dir="down", dilu_factor=5)`
- `xiaoyu_Expr.add_condition_with_sample_array(condition_name, dtype, conditions)`
- `xiaoyu_Expr.sample_index_to_well_index_output()`
- `xiaoyu_Expr.median_96well(experiment, channel="Alexa 647-A", interval="\t")`
- `gate_outliers(experiment, channel="Alexa 647-A", low=2000, high=2147483647)`
- `assist_gate(experiment, channel, low, high)`
- `auto_gate_FSC_SSC(experiment, keep=0.6)`
- `auto_gate_FSC_A_H(experiment, keep=0.9)`
- `gate_FSC_SSC(experiment, vertices)`
- `gate_FSC_A_H(experiment, vertices)`
- `polygon_gate(experiment, xchannel, ychannel, vertices=...)`
- `interactive_gate_preview(experiment, xchannel, ychannel)` — marimo-compatible polygon picker
- `apply_drawn_gate(experiment, name)` — apply the polygon drawn in interactive_gate_preview
- `density_plot(experiment, xchannel, ychannel)`
- `subset_by_char(experiment, char)`
- `subset_by_num(experiment, num)`
- `subset_by_well(experiment, well)`
- `cell_count(experiment)`

See `PROJECT_GUIDE.md` for fuller descriptions.

## Channel Name Defaults

Many helpers assume channel names from Xiaoyu's instrument:

- `FSC 488/10-A`
- `SSC 488/10-A`
- `FSC 488/10-H`
- `Alexa 647-A`
- `mCherry-A`
- `TagBFP-A`
- `549/15-488 nm citrine-A`

If a user has different channel names, inspect `experiment.data.columns` and
pass explicit `channel`, `xchannel`, or `ychannel` arguments.

## Important Constraints

- Available scales are `linear` and `log`.
- Logicle support was removed.
- Use `interactive_gate_preview()` + `apply_drawn_gate()` for polygon gating;
  the old `DraggablePolygon`-based interactive mode is removed.
- Headless scripts should use `matplotlib.use("Agg")` before importing pyplot.
- `subset_by_num` expects plate column values as strings, for example `"12"`.
- `median_96well` writes `median_values.txt` into the FCS working directory.

## Safe Code Generation Rules

When generating user scripts:

- Start with `import xiaoyu_CF as flow`.
- Use `xiaoyu_Expr(...)` for loading FCS directories.
- Print `expr.expr.data.columns` when channel names are uncertain.
- Prefer explicit channel arguments instead of hard-coding new names.
- Use `auto_gate_FSC_SSC` before `auto_gate_FSC_A_H` for the standard cleanup.
- Use `plt.show()` only for interactive scripts.
- Avoid using removed Cytoflow features such as logicle scales or unused
  advanced views.

## Marimo Notebook Support

This package works with [marimo](https://marimo.io) notebooks. Install with:

```bash
pip install "xiaoyu-cytoflow[marimo]"
```

See `examples/marimo_template.py` and `examples/flow_analysis.py` for a working
template.

**Module-style pattern** — Put frequently-changed plotting/analysis logic in a
separate `.py` module (like `examples/flow_analysis.py`), import it into the
notebook, and enable **"On module change → autorun"** in marimo settings. This
lets you edit the module in any editor (opencode, VS Code, etc.) and marimo
auto-refreshes affected cells in the browser.

Key marimo conventions:
- Every cell-level variable name must be unique across the notebook
- Prefix temporary variables with `_` to keep them cell-local
- Use `mo.mpl.interactive(fig)` to display matplotlib figures
- Use `mo.ui.table(df)` for interactive data tables
- For interactive polygon gating, use `flow.interactive_gate_preview()` +
  `flow.apply_drawn_gate()` — ``polygon_gate`` with ``vertices=[]`` does not
  work in marimo (it launches a blocking ``plt.show()`` window).

## Useful Local Commands

Development environment:

```bash
mamba env create -f environment.yml
mamba activate xiaoyu-cytoflow
```

Smoke test:

```bash
MPLCONFIGDIR=/tmp/mplconfig python test/package_smoke_test.py
```

Build:

```bash
poetry check
poetry build
```

## Documentation Map

- `README.md`: high-level overview and quick start.
- `PROJECT_GUIDE.md`: detailed user and maintainer documentation.
- `THIRD_PARTY_NOTICES.md`: vendored dependency and license notes.
- `test/package_smoke_test.py`: executable example of a non-interactive
  analysis workflow.
