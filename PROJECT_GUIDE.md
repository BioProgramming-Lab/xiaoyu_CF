# xiaoyu-cytoflow Project Guide

This guide documents the structure, assumptions, workflow, and maintenance
practice for `xiaoyu-cytoflow`.

The short version: this package wraps Xiaoyu's flow-cytometry analysis workflow
and vendors the Cytoflow pieces required to run it on Python 3.12 without a
separate upstream Cytoflow installation.

## Project Goals

`xiaoyu-cytoflow` is designed for reproducible analysis of plate-based flow
cytometry experiments. The current workflow is centered on 96-well FCS data with
well names such as `A1`, `A2`, ..., `H12`.

The package aims to:

- load FCS files from one experiment directory;
- assign plate metadata automatically from file names;
- add concentration or custom sample conditions;
- perform common cleanup gates, including outlier removal, FSC/SSC gating, and
  FSC-A/FSC-H doublet removal;
- make quick diagnostic plots;
- subset experiments by row, column, or well;
- export median channel values in a 96-well table format;
- run without installing the original upstream `cytoflow` package.

This package is intentionally practical rather than a complete replacement for
all Cytoflow features.

## Repository Layout

| Path | Purpose |
| --- | --- |
| `xiaoyu_CF.py` | Main user-facing analysis helpers. Import this as `xiaoyu_CF` or `flow`. |
| `analysis.py` | Minimal example script using the original workflow style. |
| `cytoflow/` | Slim vendored Cytoflow runtime subset. |
| `fcsparser/` | Vendored FCS parser used by Cytoflow import operations. |
| `environment.yml` | Reproducible mamba environment for Python 3.12 development and testing. |
| `pyproject.toml` | Package metadata and Poetry build configuration. |
| `test/package_smoke_test.py` | Non-interactive smoke test for the self-contained package. |
| `THIRD_PARTY_NOTICES.md` | Third-party license and vendoring notes. |
| `LICENSE.txt` | Project license text. |

The `test/` directory can also contain local experiment data and exploratory
analysis scripts. Those files are useful for validating real workflows, but they
are not part of the public package API.

## Installation

### Recommended Local Environment

From the repository root:

```bash
mamba env create -f environment.yml
mamba activate xiaoyu-cytoflow
```

The environment installs this repository in editable mode with `pip -e .`.
That means edits to `xiaoyu_CF.py` or the vendored runtime are picked up without
reinstalling.

### Existing Python 3.12 Environment

If the dependencies are already available:

```bash
pip install -e .
```

The package metadata currently supports Python `>=3.12,<3.13`. It may work on
nearby Python versions, but Python 3.12 is the tested target.

### Optional Table Output

Some helper methods print prettier tables when `tabulate` is installed:

```bash
pip install "xiaoyu-cytoflow[tables]"
```

For editable development:

```bash
pip install -e ".[tables]"
```

## Data Layout Assumptions

The loader expects FCS files in a single directory:

```text
experiment_directory/
  A1 sample name.fcs
  A2 sample name.fcs
  ...
  H12 sample name.fcs
```

Only the first whitespace-separated part of each file name is used as the well
identifier. For example:

```text
A1 EGFR.fcs -> A1
B12 condition.fcs -> B12
```

The well identifier is split into:

- `char`: the plate row, such as `A`, `B`, or `H`;
- `num`: the plate column, stored as a string such as `"1"` or `"12"`;
- `well`: the combined well name, such as `"A1"`;
- `sample`: a numeric sample index assigned according to sorted file order.

The workflow assumes a standard 8 by 12 plate when using concentration mapping
and 96-well median export.

## Main Workflow

The common workflow has five phases:

1. Import FCS files.
2. Attach experimental metadata.
3. Gate unwanted events.
4. Plot or subset the gated experiment.
5. Export summary values.

Example:

```python
import matplotlib.pyplot as plt
import xiaoyu_CF as flow

expr = flow.xiaoyu_Expr("test/26Mar2025/HuaiDan-New Experiment-20250326-1722")
expr.add_conc_condition(1e-5, dilu_dir="down")

cell_population = flow.auto_gate_FSC_SSC(expr.expr, keep=0.6)
single_cells = flow.auto_gate_FSC_A_H(cell_population, keep=0.7)

row_a = flow.subset_by_char(single_cells, "A")
flow.channel_histogram(row_a, channel="Alexa 647-A", huefacet="conc")

expr.median_96well(single_cells, channel="Alexa 647-A")
plt.show()
```

The object returned by `xiaoyu_Expr(...).expr` is a vendored Cytoflow
`Experiment`. Its event table is available as `experiment.data`, a pandas
`DataFrame`.

## Public API

### `xiaoyu_Expr`

```python
expr = flow.xiaoyu_Expr(working_dir, merge_dirs=[])
```

Loads all `.fcs` files in `working_dir`, sorts them by file path, extracts well
IDs from file names, builds Cytoflow `Tube` objects, and imports them into one
experiment.

Important attributes:

- `expr.working_dir`: directory used for loading and output;
- `expr.filelist`: loaded FCS file paths;
- `expr.name_locus_mapping`: sample index to `[row, zero_based_column]`;
- `expr.expr`: imported Cytoflow experiment.

`merge_dirs` can be used to include additional FCS files from other directories.
The current metadata mapping is based on the primary directory file names, so
merged use should be checked carefully for matching order and plate layout.

### `add_conc_condition`

```python
expr.add_conc_condition(ori_conc, dilu_dir="down", dilu_factor=5)
```

Adds a floating-point `conc` condition to the experiment.

Supported dilution directions:

- `"down"`: concentration decreases from row A toward row H;
- `"up"`: concentration decreases from row H toward row A;
- `"right"`: concentration decreases from column 1 toward column 12;
- `"left"`: concentration decreases from column 12 toward column 1.

This method is designed for 96-well plates.

### `add_condition_with_sample_array`

```python
expr.add_condition_with_sample_array("cell_type", "category", values)
```

Adds a custom condition by sample index. `values[0]` is assigned to sample 1,
`values[1]` to sample 2, and so on.

The `dtype` argument is passed to Cytoflow's `add_condition`; common values are
`"category"`, `"str"`, `"float"`, and `"int"`.

### `sample_index_to_well_index_output`

```python
expr.sample_index_to_well_index_output()
```

Prints the mapping from numeric sample index to well name. If `tabulate` is
installed, the output is formatted as a table.

### `median_96well`

```python
median_values = expr.median_96well(
    experiment,
    channel="Alexa 647-A",
    interval="\t",
)
```

Computes the median of `channel` for each sample and writes a plate-shaped table
to:

```text
<working_dir>/median_values.txt
```

It also returns the values as an 8 by 12 Python list.

## Gating Helpers

### `gate_outliers`

```python
filtered = flow.gate_outliers(
    experiment,
    channel="Alexa 647-A",
    if_plot=False,
    low=2000,
    high=2147483647,
)
```

Applies a one-dimensional range gate and keeps events where `channel` is between
`low` and `high`.

### `assist_gate`

```python
positive = flow.assist_gate(
    experiment,
    channel="mCherry-A",
    low=5e6,
    high=1e9,
)
```

Applies a named range gate for quick manual thresholding. This function stores
created gate objects in a module-level list so names stay unique during a
session.

### `auto_gate_FSC_SSC`

```python
cells = flow.auto_gate_FSC_SSC(experiment, keep=0.6, if_plot=True)
```

Uses a Cytoflow density gate on:

- `FSC 488/10-A`
- `SSC 488/10-A`

The `keep` value controls the fraction of events retained by the density gate.

### `auto_gate_FSC_A_H`

```python
single_cells = flow.auto_gate_FSC_A_H(experiment, keep=0.9)
```

Uses a Cytoflow density gate on:

- `FSC 488/10-A`
- `FSC 488/10-H`

This is typically used after FSC/SSC gating to reduce doublets or abnormal
events.

### `gate_FSC_SSC` and `gate_FSC_A_H`

```python
cells = flow.gate_FSC_SSC(experiment, vertices)
single_cells = flow.gate_FSC_A_H(experiment, vertices)
```

Apply explicit polygon gates using user-provided vertices:

```python
vertices = [
    [100000, 50000],
    [250000, 70000],
    [260000, 180000],
    [120000, 160000],
]
```

### `polygon_gate`

```python
gated = flow.polygon_gate(
    experiment,
    xchannel="FSC 488/10-A",
    ychannel="SSC 488/10-A",
    vertices=vertices,
)
```

Applies a general polygon gate for any two channels.

If `vertices=[]`, the function opens an interactive matplotlib plot where points
can be added and dragged. This requires a graphical matplotlib backend, so it is
not suitable for headless test runs.

## Plotting Helpers

### `FSC_SSC_ploting`

```python
flow.FSC_SSC_ploting(experiment, type="scatter")
flow.FSC_SSC_ploting(experiment, type="density")
```

Plots FSC-A versus SSC-A using either Cytoflow scatter or density views.

The function name keeps the original spelling for compatibility.

### `FSC_A_H_ploting`

```python
flow.FSC_A_H_ploting(experiment)
```

Plots FSC-A versus FSC-H.

### `channel_histogram`

```python
flow.channel_histogram(
    experiment,
    channel="Alexa 647-A",
    huefacet="sample",
    bins=200,
)
```

Plots a log-scaled histogram for a channel. If `huefacet="conc"`, the hue scale
is also log-scaled.

### `density_plot`

```python
fig, ax = flow.density_plot(
    experiment,
    xchannel="FSC 488/10-A",
    ychannel="FSC 488/10-H",
    xscale="log",
    yscale="log",
)
```

Creates a `mpl-scatter-density` plot and returns the matplotlib figure and axes.

### `log_density_plot`

```python
flow.log_density_plot(
    experiment,
    xchannel="Alexa 647-A",
    ychannel="mCherry-A",
    huefacet="sample",
)
```

Uses Cytoflow's `DensityView` with log-scaled x and y channels.

## Subsetting and Counting

### `subset_by_char`

```python
row_a = flow.subset_by_char(experiment, "A")
```

Returns events from one plate row.

### `subset_by_num`

```python
column_12 = flow.subset_by_num(experiment, "12")
```

Returns events from one plate column. Column values are stored as strings.

### `subset_by_well`

```python
a1 = flow.subset_by_well(experiment, "A1")
selected = flow.subset_by_well(experiment, ["A1", "A2", "B1"])
```

Returns events from one well or a list of wells.

### `cell_count`

```python
n = flow.cell_count(experiment)
```

Returns the number of events in an experiment.

## Vendored Cytoflow Runtime

The package includes a reduced `cytoflow/` tree copied from a working local
environment and updated for modern dependencies. This was done because upstream
Cytoflow was difficult to install and needed local fixes.

Included Cytoflow pieces:

- `Experiment`
- `Tube`
- `ImportOp`
- `RangeOp`
- `PolygonOp`
- `DensityGateOp`
- `AutofluorescenceOp`
- `BleedthroughLinearOp`
- `HistogramView`
- `ScatterplotView`
- `DensityView`

Removed Cytoflow pieces include many statistical views, advanced operations,
script entry points, and the compiled logicle extension.

Available scales:

- `linear`
- `log`

The old compiled logicle extension is not included. If future work needs
logicle-like visualization, the cleanest options are either implementing a pure
Python scale or adding a maintained external transform dependency.

## Testing

Run the smoke test from the repository root:

```bash
MPLCONFIGDIR=/tmp/mplconfig python test/package_smoke_test.py
```

With mamba:

```bash
MPLCONFIGDIR=/tmp/mplconfig mamba run -n xiaoyu-cytoflow python -u test/package_smoke_test.py
```

The smoke test checks that:

- no compiled `.so` files are present in the repository;
- only `linear` and `log` scales are registered;
- test FCS files can be imported;
- automatic FSC/SSC and FSC-A/FSC-H gates run;
- expected fluorescence channels exist;
- representative threshold gates run;
- a JSON summary is written to `test/package_smoke_summary.json`.

The smoke test uses the non-interactive `Agg` matplotlib backend, so it can run
without opening plot windows.

## Building and Publishing

The project uses Poetry's build backend through `poetry-core`. Poetry does not
need to be installed inside the runtime mamba environment.

Recommended publishing setup:

```bash
mamba install -n base -c conda-forge pipx
pipx ensurepath
pipx install poetry
poetry config pypi-token.pypi <your-token>
```

Build and check:

```bash
poetry check
poetry build
python -m twine check dist/*
```

Publish:

```bash
poetry publish
```

Or rebuild and publish in one step:

```bash
poetry publish --build
```

Before publishing a new release, update the version in `pyproject.toml`.

## Maintenance Notes

### Keep the Public API Stable

Existing analysis scripts import `xiaoyu_CF` directly and call functions such as
`auto_gate_FSC_SSC`, `subset_by_num`, and `median_96well`. Avoid renaming these
without adding compatibility aliases.

### Be Careful With Channel Names

The helpers currently use channel names from the source instrument, including:

- `FSC 488/10-A`
- `SSC 488/10-A`
- `FSC 488/10-H`
- `Alexa 647-A`
- `mCherry-A`
- `TagBFP-A`
- `549/15-488 nm citrine-A`

If another instrument exports different names, users should pass explicit
channel arguments where supported. A future improvement could add a channel name
configuration layer.

### Interactive Plotting Needs a GUI Backend

Most plotting helpers work in normal matplotlib sessions. The interactive
polygon picker requires a GUI-capable backend and will not work in headless
CI-style runs.

### Generated Outputs

Common generated files include:

- `median_values.txt` in the experiment directory;
- `test/package_smoke_summary.json`;
- build artifacts under `dist/`;
- Python cache directories such as `__pycache__/`.

Do not treat generated files as source unless they are intentionally checked in.

## Troubleshooting

### Import Fails With Missing Dependencies

Use the mamba environment:

```bash
mamba env create -f environment.yml
mamba activate xiaoyu-cytoflow
```

This is the most reproducible path because scientific Python packages can be
more reliable from conda-forge than from source builds.

### No FCS Files Loaded

Check that the directory path points directly to the folder containing `.fcs`
files. The loader does not recursively search subdirectories.

### Well Mapping Looks Wrong

Check file names. The loader expects the first whitespace-separated token to be
the well name, such as `A1` or `H12`.

### Plots Do Not Appear

Check the matplotlib backend. In notebooks or scripts, make sure `plt.show()` is
called after plotting.

### Logicle Scale Is Missing

That is expected. The compiled logicle extension was removed to keep this
package pure Python. Use `linear` or `log` scales with the current workflow.

## Example End-to-End Script

```python
import matplotlib.pyplot as plt
import xiaoyu_CF as flow

DATA_DIR = "test/26Mar2025/HuaiDan-New Experiment-20250326-1722"

expr = flow.xiaoyu_Expr(DATA_DIR)
expr.add_conc_condition(1e-5, dilu_dir="down")

print("Imported events:", flow.cell_count(expr.expr))

without_outliers = flow.gate_outliers(
    expr.expr,
    channel="Alexa 647-A",
    low=2000,
)

cells = flow.auto_gate_FSC_SSC(without_outliers, keep=0.6, if_plot=True)
single_cells = flow.auto_gate_FSC_A_H(cells, keep=0.7, if_plot=True)

print("Single-cell events:", flow.cell_count(single_cells))

column_12 = flow.subset_by_num(single_cells, "12")
flow.channel_histogram(column_12, channel="Alexa 647-A", huefacet="conc")

expr.median_96well(single_cells, channel="Alexa 647-A")

plt.show()
```
