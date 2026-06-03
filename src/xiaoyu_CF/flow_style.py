"""Load DEFAULT_PLOT_STYLE.toml and expose defaults + merge utilities.

Usage:
    import flow_style

    flow_style.apply()                          # apply global rcParams
    fig, ax = plt.subplots(**flow_style.figure)  # figure presets

    # merge: user kwargs override defaults
    ax.hist(data, **flow_style.merge("histogram", bins=100, alpha=0.5))
"""
import os
import tomllib
import matplotlib as mpl

_TOML_PATH = os.path.join(os.path.dirname(__file__), "DEFAULT_PLOT_STYLE.toml")

with open(_TOML_PATH, "rb") as _f:
    _cfg = tomllib.load(_f)

# ===========================================================================
# Exposed defaults — access as flow_style.figure, flow_style.histogram, etc.
# ===========================================================================

figure = _cfg["figure"]
histogram = _cfg["histogram"]
scatter = _cfg["scatter"]
density = _cfg["density"]
palette = _cfg["palette"]

PALETTE_10 = _cfg["palette"]["palette_10"]

# ===========================================================================
# Global rcParams
# ===========================================================================

_GLOBAL_RC = {
    "figure.figsize": tuple(_cfg["figure"]["figsize"]),
    "figure.dpi": _cfg["figure"]["dpi"],
    "font.family": _cfg["global"]["font_family"],
    "font.size": _cfg["global"]["font_size"],
    "axes.titlesize": _cfg["global"]["axes_titlesize"],
    "axes.titleweight": _cfg["global"]["axes_titleweight"],
    "axes.labelsize": _cfg["global"]["axes_labelsize"],
    "axes.spines.top": _cfg["global"]["axes_spines_top"],
    "axes.spines.right": _cfg["global"]["axes_spines_right"],
    "axes.grid": _cfg["global"]["axes_grid"],
    "legend.fontsize": _cfg["global"]["legend_fontsize"],
    "legend.frameon": _cfg["global"]["legend_frameon"],
    "legend.fancybox": _cfg["global"]["legend_fancybox"],
    "legend.framealpha": _cfg["global"]["legend_framealpha"],
    "xtick.labelsize": _cfg["global"]["xtick_labelsize"],
    "ytick.labelsize": _cfg["global"]["ytick_labelsize"],
    "savefig.bbox": _cfg["global"]["savefig_bbox"],
    "savefig.dpi": _cfg["global"]["savefig_dpi"],
}


def apply():
    """Apply global rcParams from DEFAULT_PLOT_STYLE.toml."""
    mpl.rcParams.update(_GLOBAL_RC)


def merge(section, **kwargs):
    """Return default values from `section` merged with user kwargs.

    User kwargs take priority over defaults.

    Parameters
    ----------
    section : str
        One of "figure", "histogram", "scatter", "density".
    **kwargs : dict
        User overrides (e.g., bins=100, alpha=0.5).

    Returns
    -------
    dict
        Merged keyword-argument dict.
    """
    defaults = dict(_cfg.get(section, {}))
    defaults.update(kwargs)
    return defaults


def reset():
    """Reset rcParams to matplotlib defaults."""
    mpl.rcdefaults()
