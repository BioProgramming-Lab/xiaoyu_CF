import mpl_scatter_density
import matplotlib.pyplot as plt
from matplotlib.widgets import PolygonSelector
import glob
import os.path
import numpy as np
import pandas as pd
import scipy.ndimage as ndi
import warnings

from ._experiment import Experiment
from ._operations import Tube, ImportOp, RangeOp, PolygonOp, DensityGateOp
from ._views import HistogramView
from ._utility import *
try:
    from .flow_style import figure as _style_fig, density as _style_den
except ImportError:
    _style_fig = {}
    _style_den = {}

# ===========================================================================
# xiaoyu_Expr - experiment loader and metadata manager
# ===========================================================================

class xiaoyu_Expr:

    # load fcs files in the working directory
    # usage: xiaoyu_Expr("working_directory", ["merge_directory1", "merge_directory2", ...])
    def __init__(self, working_dir : str, merge_dirs = []):
        self.working_dir = ""
        self.filelist = []
        self.name_locus_mapping = {}
        self.expr = None
        if working_dir != "" and working_dir[-1] != "/":
            working_dir += "/"
        self.working_dir = working_dir
        filelist = glob.glob(working_dir + "*.fcs")
        filelist = sorted(filelist)
        self.filelist = filelist
        print (f"{len(filelist)} fsc files loaded")

        # load other fcs files in the merge_dirs, and merge them to the current experiment by well index.
        # this is useful when you have some additive data to the current experiment.
        if len(merge_dirs) != 0:
            for merge_dir in merge_dirs:
                if merge_dir[-1] != "/":
                    merge_dir += "/"
                filelist_merge = glob.glob(merge_dir + "*.fcs")
                filelist_merge = sorted(filelist_merge)
                self.filelist += filelist_merge
                print (f"other {len(filelist_merge)} fsc files loaded from {merge_dir}. will be merged to the current experiment by well index.")

        # classify samples with name (A1, A2, A3, ..., B1, B2, B3, ...)
        file_names = []
        for i in filelist:
            file_names.append(os.path.splitext(os.path.split(i)[-1])[0])
            file_names[-1] = file_names[-1].split(" ")[0]   # filter out the well name. e.g. "A1 EGFR.fcs" -> "A1"
        for i in range(1, len(file_names)+1):
            self.name_locus_mapping[i] = [file_names[i-1][0], int(file_names[i-1][1:])-1]
        #print (f"mapping of sample name to char and num: {name_locus_mapping}")

        # generate experiment instance
        tubes = []
        count = 1
        # for merged data, the sample index (num, char, well) should be the same with the original data, but the sample index (sample) should be different for clarification.
        for i in filelist:
            tube = Tube(file = i, conditions = {"sample" : count, "char" : self.name_locus_mapping[count][0], "num" : str(self.name_locus_mapping[count][1]+1), "well" : self.name_locus_mapping[count][0]+str(self.name_locus_mapping[count][1]+1)})
            tubes.append(tube)
            count += 1
        experiment = ImportOp(conditions = {"sample" : "int", "char" : "category", "num" : "category", "well" : "str"}, tubes = tubes).apply()
        print (f"properties of the experiment: {experiment.data.columns}")
        self.expr = experiment

    # add concetration condition to the experiment
    # usage: add_conc_condition(ori_conc, dilu_dir = "down", dilu_factor = 5)
    # this function will add a new condition to the experiment with the name "conc" and the dtype "float".
    # the concentration of the samples will be assigned according to the dilution direction and dilution factor.
    # the concentration of the first sample will be the ori_conc, and the concentration of the other samples will be calculated by the dilution factor.
    # this sample only works for 96-well plate.
    def add_conc_condition(self, ori_conc, dilu_dir = "down", dilu_factor = 5):
        conc = []
        condition_list = []
        if dilu_dir == "left" or dilu_dir == "right":
            for i in range(12):
                if dilu_dir == "left":
                    conc.append(ori_conc / (dilu_factor ** (12-i)))
                else:
                    conc.append(ori_conc / (dilu_factor ** i))
            for i in self.expr.data["num"]:
                condition_list.append( round(conc[int(i)-1], 18) )
            condition_list = np.array(condition_list)
            self.expr.add_condition(name = "conc", dtype = "float", data = condition_list)
        if dilu_dir == "down" or dilu_dir == "up":
            for i in range(8):
                if dilu_dir == "up":
                    conc.append(ori_conc / (dilu_factor ** (8-i)))
                else:
                    conc.append(ori_conc / (dilu_factor ** i))
            for i in self.expr.data["char"]:
                condition_list.append( round(conc[ord(i)-ord('A')], 18) )
            condition_list = np.array(condition_list)
            self.expr.add_condition(name = "conc", dtype = "float", data = condition_list)

    # add other conditions
    # usage: add_condition_with_sample_array("condition_name", "dtype", [condition1, condition2, condition3, ...])
    # this function will add a new condition to the experiment with the name "condition_name" and the dtype "dtype".
    # by assigning the list of conditions to their corresponding samples.
    # the length of the conditions should be the same with the number of samples in the experiment.
    # if module tabulate is installed, the mapping of the sample index to the new condition value will be shown.
    def add_condition_with_sample_array(self, condition_name, dtype, conditions):
        condition_list = []
        for i in self.expr.data["sample"]:
            condition_list.append( conditions[i-1] )
        condition_list = np.array(condition_list)
        self.expr.add_condition(name = condition_name, dtype = dtype, data = condition_list)
        print (f"condition {condition_name} added to the experiment")
        print ("mapping of the sample index to the new condition value is shown below:")
        try:
            from tabulate import tabulate
            print (tabulate([[i, f"{self.name_locus_mapping[i][0]}{self.name_locus_mapping[i][1]+1}", conditions[i-1]] for i in range(1, len(self.name_locus_mapping)+1)], headers = ["sample", "well", condition_name], tablefmt='fancy_grid'))
        except:
            print ("tabulate module not found. please install tabulate module to show the mapping of the sample index to the new condition value.")

    def sample_index_to_well_index_output(self):
        try:
            from tabulate import tabulate
            print (tabulate([[i, f"{self.name_locus_mapping[i][0]}{self.name_locus_mapping[i][1]+1}"] for i in range(1, len(self.name_locus_mapping)+1)], headers = ["sample", "well"], tablefmt='fancy_grid'))
        except:
            print ("tabulate module not found. please install tabulate module to show the mapping of the sample index to the well index.")

    # median value output
    # usage: median_96well(experiment, channel = "Alexa 647-A", interval = "\t")
    # this function will output the median values of the channel in the 96-well plate format.
    # the output will be saved to the working directory with the name "median_values.txt".
    # the interval between the values can be set by the parameter "interval".
    # the median values will also be returned as a 2D list.
    # this function only works for 96-well plate.
    def median_96well(self, experiment, channel = "Alexa 647-A", interval = "\t"):
        median_values_raw = experiment.data.groupby('sample')[channel].median()
        median_values = [[0 for i in range(12)] for j in range(8)]
        for key in median_values_raw.index:
            median_values[ord(self.name_locus_mapping[key][0])-ord('A')][self.name_locus_mapping[key][1]] = median_values_raw[key]
        file = open(self.working_dir + "median_values.txt", "w")
        for i in median_values:
            for j in i:
                file.write(str(j) + interval)
            file.write("\n")
        file.close()
        print (f"median values of {channel} output to {self.working_dir}median_values.txt")
        return median_values

# ===========================================================================
# Gating helpers
# ===========================================================================

# gate out outliers
def gate_outliers(experiment, channel = "Alexa 647-A", if_plot = False, low = 2000, high = 2147483647):
    gate = RangeOp(name = "de_outliers", channel = channel, low = low, high = high)
    if if_plot:
        HistogramView(channel = channel, scale = "log", huefacet = "well").plot(experiment, num_bins = 200, density = True, title = f"Histogram of {channel}")
    return gate.apply(experiment).subset("de_outliers", True)

# assist auto gate
assistants = []
def assist_gate(experiment, channel = "FSC 488/10-A", if_plot = False, low = 2000, high = 2147483647):
    global assistants
    gate = RangeOp(name = f"assist{len(assistants)+1}", channel = channel, low = low, high = high)
    assistants.append(gate)
    if if_plot:
        gate.default_view(scale = "log").plot(experiment, title = channel)
    return gate.apply(experiment).subset(f"assist{len(assistants)}", True)  # without +1 because already conducted appending

def plot_density_gate(experiment, xchannel, ychannel, gate_op,
                      xscale=None, yscale=None, contour_color="black"):
    if xscale is None:
        xscale, _, _ = channel_info(experiment, xchannel)
    if yscale is None:
        yscale, _, _ = channel_info(experiment, ychannel)

    fig, ax = density_plot(experiment, xchannel=xchannel, ychannel=ychannel,
                            xscale=xscale, yscale=yscale)

    hist = gate_op._histogram.get(True)
    if hist is not None and hist.size > 0:
        level = hist[gate_op._keep_xbins[True][-1], gate_op._keep_ybins[True][-1]]
        xbins = gate_op._xbins[:-1]
        ybins = gate_op._ybins[:-1]
        ax.contour(xbins, ybins, hist.T, [level],
                   colors=contour_color, linewidths=2)
    return fig, ax


def _warn_saturation(experiment):
    skip = {"TIME", "TLSW", "TMSW", "Event Info"}
    for ch in experiment.channels:
        if ch in skip:
            continue
        rng = experiment.metadata[ch].get("range", experiment[ch].max())
        near_max = (experiment[ch] > rng * 0.999).sum()
        if near_max > len(experiment) * 0.01:
            print(f"[WARN] {ch}: {near_max:,} events near detector max "
                  f"({near_max/len(experiment)*100:.1f}%). "
                  "Run gate_saturation() first.")
            return


def plot_density_gate(experiment, xchannel, ychannel, gate_op,
                      xscale=None, yscale=None, contour_color="black"):
    """Plot a density gate: density of full experiment + contour at gate boundary.

    Parameters
    ----------
    experiment : Experiment
        The experiment BEFORE gating (original data).
    gate_op : DensityGateOp
        The gate operation (after estimate).
    """
    if xscale is None:
        xscale, _, _ = channel_info(experiment, xchannel)
    if yscale is None:
        yscale, _, _ = channel_info(experiment, ychannel)

    fig, ax = density_plot(experiment, xchannel=xchannel, ychannel=ychannel,
                            xscale=xscale, yscale=yscale)

    hist = gate_op._histogram.get(True)
    if hist is not None and hist.size > 0:
        level = hist[gate_op._keep_xbins[True][-1], gate_op._keep_ybins[True][-1]]
        xbins = gate_op._xbins[:-1]
        ybins = gate_op._ybins[:-1]
        ax.contour(xbins, ybins, hist.T, [level],
                   colors=contour_color, linewidths=2)
    return fig, ax


def _gate_with_polygon(experiment, gate_op, name="cells"):
    """Apply a DensityGateOp as a polygon gate for smooth boundaries."""
    import warnings
    poly = _polygon_from_gate(gate_op)
    if poly and len(poly) >= 3:
        gate = PolygonOp(name=name,
                         xchannel=gate_op.xchannel,
                         ychannel=gate_op.ychannel,
                         vertices=poly)
        gated = gate.apply(experiment)
        return gated.subset(name, True)
    else:
        # fallback: use the original density gate
        warnings.warn("Could not extract polygon; falling back to bin gate")
        gated = gate_op.apply(experiment)
        return gated.subset(name, True)


def gate_saturation(experiment):
    """Remove events where any channel is at the detector saturation ceiling.

    This cleans up events whose signal hit the maximum detection range,
    a common artifact in flow cytometry data.
    """
    skip = {"TIME", "TLSW", "TMSW", "Event Info"}
    keep = pd.Series(True, index=experiment.data.index)

    for ch in experiment.channels:
        if ch in skip:
            continue
        rng = experiment.metadata[ch].get("range", experiment.data[ch].max())
        keep &= experiment.data[ch] < rng * 0.999

    result = experiment.clone(deep=False)
    result.data = experiment.data[keep].reset_index(drop=True)
    return result


def watershed_gate_FSC_SSC(experiment, sigma=12, grid=512,
                            h_min=0.005, if_plot=False):
    """Gate cells via histogram blur + watershed on FSC-A vs SSC-A.

    Returns the largest basin (cell population) as the gated experiment.

    Parameters
    ----------
    experiment : Experiment
    sigma : float or int
        Gaussian blur sigma for the KDE smoothing (default 8).
    grid : int
        Grid resolution (default 512).
    h_min : float
        Density threshold fraction (default 0.005 = 0.5% of max).
    if_plot : bool
    """
    fsc = _find_channel(experiment, "FSC", "-A")
    ssc = _find_channel(experiment, "SSC", "-A")
    if not fsc or not ssc:
        raise ValueError("Cannot find scatter channels")

    data = experiment.data[[fsc, ssc]].values
    data_log = np.log10(data.clip(1))

    # 2D histogram
    xmin, xmax = data_log[:, 0].min(), data_log[:, 0].max()
    ymin, ymax = data_log[:, 1].min(), data_log[:, 1].max()
    H, xe, ye = np.histogram2d(data_log[:, 0], data_log[:, 1], bins=grid,
                                range=[[xmin, xmax], [ymin, ymax]])

    # gaussian blur = KDE equivalent
    Z = ndi.gaussian_filter(H, sigma=sigma)

    # density threshold → connected components = cell populations
    from scipy.ndimage import label
    mask = Z > Z.max() * h_min
    if not mask.any():
        raise ValueError("Density threshold too high — lower h_min")
    labels = label(mask)[0]

    # optional: fill holes in each label
    labels = ndi.binary_fill_holes(labels > 0).astype(np.int32)
    labels = label(labels)[0]

    # map each event to its basin
    bx = np.digitize(data_log[:, 0], xe) - 1
    by = np.digitize(data_log[:, 1], ye) - 1
    bx = bx.clip(0, grid - 1)
    by = by.clip(0, grid - 1)
    event_labels = labels[by, bx]

    # keep largest basin
    unique, counts = np.unique(event_labels[event_labels > 0], return_counts=True)
    if len(unique) == 0:
        raise ValueError("Watershed found no basins — lower h_min or increase sigma")
    best = unique[np.argmax(counts)]

    result = experiment.clone(deep=False)
    result.data = experiment.data[event_labels == best].reset_index(drop=True)

    if if_plot:
        fig, axes = plt.subplots(1, 2, figsize=(16, 7))

        ax = axes[0]
        ax.imshow(Z.T, origin='lower', extent=[xe[0], xe[-1], ye[0], ye[-1]],
                  cmap=_cytoflow_cmap(), aspect='auto')
        from matplotlib.colors import ListedColormap
        _unq = np.unique(labels)
        _unq = _unq[_unq > 0]
        masked = np.ma.masked_where(labels == 0, labels)
        ax.imshow(masked.T, origin='lower', extent=[xe[0], xe[-1], ye[0], ye[-1]],
                  cmap=ListedColormap(plt.cm.tab10(np.linspace(0, 1, len(_unq)))),
                  alpha=0.3, aspect='auto')
        ax.set_xlabel(f"log10 {fsc}"); ax.set_ylabel(f"log10 {ssc}")
        ax.set_title(f"Watershed basins ({len(_unq)})")

        _fig, _ax = density_plot(result, xchannel=fsc, ychannel=ssc,
                                  xscale="linear", yscale="linear")
        _ax.set_title(f"Kept basin ({len(result):,})")

        plt.tight_layout()
        plt.show()

    return result


def _find_channel(experiment, *patterns):
    """Return the first channel name containing all given patterns."""
    for ch in experiment.channels:
        if all(p in ch for p in patterns):
            return ch
    return None


def auto_gate_FSC_SSC(experiment, keep=0.6, return_fig=False):
    fsc = _find_channel(experiment, "FSC", "-A")
    ssc = _find_channel(experiment, "SSC", "-A")
    if not fsc or not ssc:
        raise ValueError(f"Cannot find FSC-A ({fsc}) or SSC-A ({ssc}) channels")

    _warn_saturation(experiment)
    dens_gate = DensityGateOp(name = "auto_FSCA_SSCA", xchannel = fsc, ychannel = ssc,
                              xscale = "linear", yscale = "linear", keep = keep)
    dens_gate.estimate(experiment)
    gated = dens_gate.apply(experiment)
    result = gated.subset("auto_FSCA_SSCA", True)

    if return_fig:
        fig, _ = plot_density_gate(experiment, fsc, ssc, dens_gate,
                                    xscale="linear", yscale="linear")
        return result, fig
    return result


def auto_gate_FSC_A_H(experiment, keep=0.9, return_fig=False):
    fsc_a = _find_channel(experiment, "FSC", "-A")
    fsc_h = _find_channel(experiment, "FSC", "-H")
    if not fsc_a or not fsc_h:
        raise ValueError(f"Cannot find FSC-A ({fsc_a}) or FSC-H ({fsc_h}) channels")
    dens_gate = DensityGateOp(name = "auto_FSCA_FSCH", xchannel = fsc_a, ychannel = fsc_h,
                              xscale = "linear", yscale = "linear", keep = keep)
    dens_gate.estimate(experiment)
    gated = dens_gate.apply(experiment)
    result = gated.subset("auto_FSCA_FSCH", True)

    if return_fig:
        fig, _ = plot_density_gate(experiment, fsc_a, fsc_h, dens_gate,
                                    xscale="linear", yscale="linear")
        return result, fig
    return result

# get sub set by char
def subset_by_char(exper, char):
    subset = exper.subset("char", char)
    return subset

# get sub set by num
def subset_by_num(exper, num):
    subset = exper.subset('num', num)
    return subset

# get sub set by well index
def subset_by_well(exper, well):
    if (type(well) == str):
        subset = exper.subset('well', well)
        return subset
    elif (type(well) == list):
        subset = exper.clone(deep = False)
        subset.data = exper.data[exper.data['well'].isin(well)]
        subset.data.reset_index(drop = True, inplace = True)
        return subset
    else:
        raise TypeError("well should be a string or a list of strings, but got " + str(type(well)))

# get count of cells in the experiment
def cell_count(experiment):
    return len(experiment.data)

def channel_info(experiment, channel):
    """Return (scale, lo, hi) for plotting a channel.

    Default: linear for scatter (FSC/SSC), log for fluorescence.
    """
    scale = "linear" if "FSC" in channel or "SSC" in channel else "log"
    rng = experiment.metadata[channel].get("range", experiment.data[channel].max())
    if scale == "log":
        lo = float(max(experiment.data[channel][experiment.data[channel] > 0].min(), 1))
    else:
        lo = 0
    return scale, lo, float(rng)


# 2D Density plot
def density_plot(experiment, xchannel = "FSC 488/10-A", ychannel = "FSC 488/10-H",
                  xscale = None, yscale = None, cmap = None, colorbar = False,
                  **kwargs):
    if xscale is None:
        xscale, _xlo, _xhi = channel_info(experiment, xchannel)
    else:
        _, _xlo, _xhi = channel_info(experiment, xchannel)
    if yscale is None:
        yscale, _ylo, _yhi = channel_info(experiment, ychannel)
    else:
        _, _ylo, _yhi = channel_info(experiment, ychannel)

    if cmap is None:
        cmap = _cytoflow_cmap()

    fig = plt.figure(**_style_fig)
    fig.patch.set_alpha(0)
    ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
    ax.patch.set_alpha(0)
    DensityPlot = ax.scatter_density(experiment.data[xchannel], experiment.data[ychannel], cmap = cmap, **kwargs)
    ax.set_box_aspect(1)
    ax.set_xlabel(xchannel)
    ax.set_ylabel(ychannel)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)

    def _format_ticks(axis, scale):
        if scale == "linear":
            axis.set_major_formatter(
                plt.FuncFormatter(lambda v, _: f"{v/1000:.0f}k" if v >= 1000 else f"{v:.0f}")
            )

    _format_ticks(ax.xaxis, xscale)
    _format_ticks(ax.yaxis, yscale)

    ax.set_xlim(_xlo, _xhi)
    ax.set_ylim(_ylo, _yhi)

    if colorbar:
        fig.colorbar(DensityPlot, ax=ax, label='Density')
    ax.grid(False)
    return fig, ax


def _cytoflow_cmap():
    """Blue-through-transparent to red colormap, like professional cytometry software."""
    import matplotlib.colors as mcolors
    return mcolors.LinearSegmentedColormap.from_list(
        "cytoflow_density",
        [(0, 0, 0, 0),       # transparent (zero density)
         (0, 0.3, 0.8, 0.6), # blue
         (0.2, 0.7, 1.0, 1.0), # light blue
         (0.5, 0.9, 0.5, 1.0), # green
         (1.0, 0.9, 0.2, 1.0), # yellow
         (1.0, 0.3, 0.0, 1.0), # orange-red
         (0.7, 0.0, 0.0, 1.0), # red
        ],
        N=256,
    )

# polygon gate (pre-defined vertices)
polygon_gate_number = 0
def polygon_gate(experiment, xchannel = "FSC 488/10-A", ychannel = "SSC 488/10-A", xscale = "linear", yscale = "linear", vertices = None, if_plot = False):
    global polygon_gate_number
    if vertices is None or len(vertices) < 3:
        raise ValueError("vertices must be a list of at least 3 [x, y] pairs")
    gate = PolygonOp(name = f"polygon_gate_{polygon_gate_number}", xchannel = xchannel, ychannel = ychannel, vertices = vertices)
    gatedDATA = gate.apply(experiment)
    polygon_gate_number += 1
    if if_plot:
        gate.default_view(xscale = xscale, yscale = yscale).plot(gatedDATA, s = 10, alpha = 0.1)
    return gatedDATA.subset(f"polygon_gate_{polygon_gate_number-1}", True)

# interactive polygon gate (marimo / notebook compatible)
_gate_vertices = []
_gate_xchannel = "FSC 488/10-A"
_gate_ychannel = "SSC 488/10-A"


def interactive_gate_preview(experiment, xchannel="FSC 488/10-A",
                              ychannel="SSC 488/10-A",
                              xscale=None, yscale=None,
                              cmap=None):
    """Interactive polygon gate preview for marimo notebooks.

    Displays a density plot with a polygon selector. After drawing
    a polygon (click points, then re-click the first point to close),
    call ``apply_drawn_gate()`` to apply the gate.

    Parameters
    ----------
    experiment : Experiment
        The experiment to gate.
    xchannel, ychannel : str
        Channel names for X and Y axes.
    xscale, yscale : str
        Scale for axes ('linear' or 'log').
    cmap : str or Colormap, optional
        Colormap override. Default: transparent-blue-red.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Display with ``mo.mpl.interactive(fig)`` in marimo.
    """
    global _gate_vertices, _gate_xchannel, _gate_ychannel
    _gate_vertices.clear()
    _gate_xchannel = xchannel
    _gate_ychannel = ychannel

    fig, ax = density_plot(experiment, xchannel=xchannel, ychannel=ychannel,
                            xscale=xscale, yscale=yscale, cmap=cmap,
                            colorbar=False)

    def _onselect(verts):
        _gate_vertices.clear()
        _gate_vertices.extend([[float(x), float(y)] for x, y in verts])

    selector = PolygonSelector(ax, _onselect, useblit=True)
    fig._gate_selector = selector

    return fig


def apply_drawn_gate(experiment, name="marimo_gate"):
    """Apply the polygon gate drawn in ``interactive_gate_preview()``.

    Call this after the polygon has been drawn in the preview plot.
    Works both in marimo and in standard Python scripts.

    Parameters
    ----------
    experiment : Experiment
        The experiment to gate.
    name : str
        Name for the gate condition column.

    Returns
    -------
    Experiment
        Gated experiment subset containing only events inside the polygon.
    """
    global _gate_vertices, _gate_xchannel, _gate_ychannel
    if len(_gate_vertices) < 3:
        raise ValueError(
            "No polygon drawn yet. Draw a polygon in the interactive "
            "preview cell first, then re-run this cell."
        )
    gate = PolygonOp(name=name,
                     xchannel=_gate_xchannel,
                     ychannel=_gate_ychannel,
                     vertices=_gate_vertices)
    gated = gate.apply(experiment)
    return gated.subset(name, True)
