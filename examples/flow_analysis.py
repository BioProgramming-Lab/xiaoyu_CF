"""Shared plotting logic for marimo notebooks.

Edit this file in any editor (VS Code, opencode, etc.) and marimo
auto-re-runs affected cells if "On module change → autorun" is
enabled in marimo settings.

All plotting functions accept keyword overrides that take priority
over the defaults in DEFAULT_PLOT_STYLE.toml.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, FuncFormatter
from scipy.stats import gaussian_kde

import xiaoyu_CF as flow
import xiaoyu_CF.flow_style as flow_style

flow_style.apply()


def _ticks_k(axis):
    axis.set_major_formatter(
        FuncFormatter(lambda v, _: f"{v/1000:.0f}k" if v >= 1000 else f"{v:.0f}")
    )


def plot_histogram(experiment, channel="Alexa 647-A", huefacet="sample",
                   **kwargs):
    """Log-scale histogram with KDE overlay, faceted by huefacet.

    KDE lines are hidden when peak density < kde_threshold (avoids
    flat noisy lines for low-signal samples).

    Parameters
    ----------
    experiment : xiaoyu_CF.Experiment
        Gated experiment to plot.
    channel : str
        Channel name.
    huefacet : str
        Condition to color by.
    kde_threshold : float (default: 0.01)
        Don't draw KDE line if max density below this.
    **kwargs : dict
        Overrides for DEFAULT_PLOT_STYLE.histogram.
    """
    data = experiment.data
    groups = data[huefacet].unique()
    palette = flow_style.PALETTE_10
    kde_threshold = kwargs.pop("kde_threshold", flow_style.histogram.get("kde_threshold", 0.01))

    opts = flow_style.merge("histogram", **kwargs)

    fig, ax = plt.subplots(**flow_style.figure)
    _scale, _lo, _hi = flow.channel_info(experiment, channel)
    _max_val = float(max(_hi, data[channel].max()))

    if _scale == "log":
        _bins = np.logspace(np.log10(_lo), np.log10(_max_val), opts["bins"])
        _log_bin_width = np.log10(_bins[1]) - np.log10(_bins[0])
        _kde_x = np.logspace(np.log10(_lo), np.log10(_max_val), 500)
    else:
        _bins = np.linspace(_lo, _max_val, opts["bins"])
        _lin_bin_width = _bins[1] - _bins[0]
        _kde_x = np.linspace(_lo, _max_val, 500)

    for i, group in enumerate(groups):
        subset = data[data[huefacet] == group][channel]
        color = palette[i % len(palette)]

        _counts, _edges = np.histogram(subset, bins=_bins)
        if len(subset) > 0:
            _density = _counts / (len(subset) * (_log_bin_width if _scale == "log" else _lin_bin_width))
            ax.stairs(_density, _edges, fill=True,
                      alpha=opts["alpha"], color=color, linewidth=0)

            if opts.get("show_kde", True):
                _vals = np.log10(subset[subset > 0]) if _scale == "log" else subset[subset > 0]
                _grid = np.log10(_kde_x) if _scale == "log" else _kde_x
                if len(_vals) > 3:
                    _kde = gaussian_kde(_vals)
                    _kde_density = _kde(_grid)
                    _masked = np.where(_kde_density >= kde_threshold,
                                       _kde_density, np.nan)
                    if not np.all(np.isnan(_masked)):
                        ax.plot(_kde_x, _masked, color=color,
                                linewidth=opts["linewidth"], label=str(group))

    ax.set_xscale(_scale)
    ax.set_xlabel(channel)
    ax.set_ylabel("Density")
    ax.set_title(f"{channel} by {huefacet}")
    if _scale == "linear":
        _ticks_k(ax.xaxis)
    ax.legend()
    return fig


def plot_ridgeline(experiment, channel="Alexa 647-A", huefacet="sample",
                            **kwargs):
    """Ridgeline histogram — samples offset vertically with KDE overlay.

    Parameters
    ----------
    experiment : xiaoyu_CF.Experiment
        Gated experiment to plot.
    channel : str
        Channel name.
    huefacet : str
        Condition to color by.
    overlap : float (default: 0.6)
        How much consecutive traces overlap (0 = no overlap, 1 = full).
    **kwargs : dict
        Overrides for DEFAULT_PLOT_STYLE.histogram.
    """
    data = experiment.data
    groups = list(data[huefacet].unique())
    palette = flow_style.PALETTE_10
    overlap = kwargs.pop("overlap", 0.6)

    opts = flow_style.merge("histogram", **kwargs)

    fig, ax = plt.subplots(**flow_style.figure)
    _scale, _lo, _hi = flow.channel_info(experiment, channel)
    _max_val = float(max(_hi, data[channel].max()))

    if _scale == "log":
        _bins = np.logspace(np.log10(_lo), np.log10(_max_val), opts["bins"])
        _bin_width = np.log10(_bins[1]) - np.log10(_bins[0])
        _kde_x = np.logspace(np.log10(_lo), np.log10(_max_val), 500)
    else:
        _bins = np.linspace(_lo, _max_val, opts["bins"])
        _bin_width = _bins[1] - _bins[0]
        _kde_x = np.linspace(_lo, _max_val, 500)

    _all_densities = []
    _all_subsets = []
    for group in groups:
        subset = data[data[huefacet] == group][channel]
        _counts, _edges = np.histogram(subset, bins=_bins)
        _all_subsets.append(subset)
        if len(subset) > 0:
            _all_densities.append(_counts / (len(subset) * _bin_width))
        else:
            _all_densities.append(np.zeros(len(_bins) - 1))

    _max_y = max(d.max() for d in _all_densities) if _all_densities else 1
    _step = _max_y * (1 - overlap)
    _kde_threshold = kwargs.pop("kde_threshold", flow_style.histogram.get("kde_threshold", 0.01))

    for i, group in enumerate(reversed(groups)):
        color = palette[(len(groups) - 1 - i) % len(palette)]
        offset = i * _step
        _den = _all_densities[len(groups) - 1 - i]
        ax.fill_between(_bins[:-1], offset, _den + offset,
                         color=color, alpha=opts["alpha"], linewidth=0)

        if opts.get("show_kde", True):
            _subset = _all_subsets[len(groups) - 1 - i]
            _vals = np.log10(_subset[_subset > 0]) if _scale == "log" else _subset[_subset > 0]
            _grid = np.log10(_kde_x) if _scale == "log" else _kde_x
            if len(_vals) > 3:
                _kde = gaussian_kde(_vals)
                _kde_density = _kde(_grid)
                _masked = np.where(_kde_density >= _kde_threshold,
                                   _kde_density, np.nan)
                if not np.all(np.isnan(_masked)):
                    ax.plot(_kde_x, _masked + offset, color=color,
                            linewidth=opts["linewidth"])

    ax.set_xscale(_scale)
    ax.set_xlabel(channel)
    ax.set_yticks([])
    ax.set_ylabel("")
    ax.set_title(f"{channel} by {huefacet}")
    if _scale == "linear":
        _ticks_k(ax.xaxis)

    # label each trace on the left
    for i, group in enumerate(reversed(groups)):
        offset = i * _step
        ax.text(_lo * 0.9, offset + _step * 0.5, str(group),
                color=palette[(len(groups) - 1 - i) % len(palette)],
                ha="right", va="center", fontsize=9, fontweight="bold")

    return fig
