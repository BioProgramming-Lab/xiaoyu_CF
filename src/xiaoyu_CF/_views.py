# _views.py - merged views module for xiaoyu_CF
# Vendored from cytoflow.views.*

from traits.api import Interface, Constant, Bool, HasStrictTraits, Str, Tuple, List, Dict, provides
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from natsort import natsorted
from warnings import warn
import numpy as np
import math
import bottleneck
import scipy.ndimage.filters
import copy

from . import _utility as util


class IView(Interface):
    """
    An interface for a visualization of flow data.

    Could be a histogram, a density plot, a scatter plot, a statistical
    visualization like a bar chart of population means; even a textual
    representation like a table.

    Attributes
    ----------
    id : Constant
        A unique id for this view.  Prefix: "edu.mit.cytoflow.views"

    friendly_id : Constant
        The human-readable id of this view: eg, "Histogram"
    """

    id = Constant("FIXME")
    friendly_id = Constant("FIXME")

    def plot(self, experiment, **kwargs):
        """Plot a visualization of flow data using the pyplot stateful interface

        Parameters
        ----------
        experiment : Experiment
            the Experiment containing the data to plot
        kwargs : dict
            additional arguments to pass to the underlying plotting function.
        """


class ISelectionView(IView):
    """A decorator that lets you add (possibly interactive) selections to an IView.

    Note that this is a Decorator *design pattern*, not a Python ``@decorator``.

    Attributes
    ----------
    interactive : Bool
        Is this view's interactivity turned on?

    """

    interactive = Bool(False, transient=True)


class BaseView(HasStrictTraits):
    """
    The base class for facetted plotting.

    Attributes
    ----------
    xfacet : String
        Set to one of the `Experiment.conditions` in the `Experiment`, and
        a new column of subplots will be added for every unique value
        of that condition.

    yfacet : String
        Set to one of the `Experiment.conditions` in the `Experiment`, and
        a new row of subplots will be added for every unique value
        of that condition.

    huefacet : String
        Set to one of the `Experiment.conditions` in the in the `Experiment`,
        and a new color will be added to the plot for every unique value of
        that condition.

    huescale : {'linear', 'log'}
        How should the color scale for `huefacet` be scaled?
    """

    xfacet = Str
    yfacet = Str
    huefacet = Str
    huescale = util.ScaleEnum

    def plot(self, experiment, data, **kwargs):
        """
        Base function for facetted plotting

        Parameters
        ----------
        experiment: Experiment
            The `Experiment` to plot using this view.

        title : str
            Set the plot title

        xlabel : str
            Set the X axis label

        ylabel : str
            Set the Y axis label

        huelabel : str
            Set the label for the hue facet (in the legend)

        legend : bool
            Plot a legend for the color or hue facet?  Defaults to `True`.

        sharex : bool
            If there are multiple subplots, should they share X axes?  Defaults
            to `True`.

        sharey : bool
            If there are multiple subplots, should they share Y axes?  Defaults
            to `True`.

        row_order : list
            Override the row facet value order with the given list.
            If a value is not given in the ordering, it is not plotted.
            Defaults to a "natural ordering" of all the values.

        col_order : list
            Override the column facet value order with the given list.
            If a value is not given in the ordering, it is not plotted.
            Defaults to a "natural ordering" of all the values.

        hue_order : list
            Override the hue facet value order with the given list.
            If a value is not given in the ordering, it is not plotted.
            Defaults to a "natural ordering" of all the values.

        height : float
            The height of *each row* in inches.  Default = 3.0

        aspect : float
            The aspect ratio *of each subplot*.  Default = 1.5

        col_wrap : int
            If `xfacet` is set and `yfacet` is not set, you can "wrap" the
            subplots around so that they form a multi-row grid by setting
            this to the number of columns you want.

        sns_style : {"darkgrid", "whitegrid", "dark", "white", "ticks"}
            Which ``seaborn`` style to apply to the plot?  Default is ``whitegrid``.

        sns_context : {"paper", "notebook", "talk", "poster"}
            Which ``seaborn`` context to use?  Controls the scaling of plot
            elements such as tick labels and the legend.  Default is ``talk``.

        palette : palette name, list, or dict
            Colors to use for the different levels of the hue variable.
            Should be something that can be interpreted by
            `seaborn.color_palette`, or a dictionary mapping hue levels to
            matplotlib colors.

        despine : Bool
            Remove the top and right axes from the plot?  Default is ``True``.

        Other Parameters
        ----------------
        cmap : matplotlib colormap
            If plotting a huefacet with many values, use this color map instead
            of the default.

        norm : matplotlib.colors.Normalize
            If plotting a huefacet with many values, use this object for color
            scale normalization.

        """

        if experiment is None:
            raise util.CytoflowViewError('experiment',
                                         "No experiment specified")

        col_wrap = kwargs.pop('col_wrap', None)

        if col_wrap is not None and self.yfacet:
            raise util.CytoflowViewError('yfacet',
                                         "Can't set yfacet and col_wrap at the same time.")

        if col_wrap is not None and not self.xfacet:
            raise util.CytoflowViewError('xfacet',
                                         "Must set xfacet to use col_wrap.")

        if col_wrap is not None and col_wrap < 2:
            raise util.CytoflowViewError(None,
                                         "col_wrap must be None or > 1")

        title = kwargs.pop("title", None)
        xlabel = kwargs.pop("xlabel", None)
        ylabel = kwargs.pop("ylabel", None)
        huelabel = kwargs.pop("huelabel", self.huefacet)
        if huelabel == "":
            huelabel = self.huefacet

        sharex = kwargs.pop("sharex", True)
        sharey = kwargs.pop("sharey", True)

        height = kwargs.pop("height", 3)
        aspect = kwargs.pop("aspect", 1.5)

        legend = kwargs.pop('legend', True)

        despine = kwargs.pop('despine', False)
        palette = kwargs.pop('palette', None)

        if 'sns_style' in kwargs:
            kwargs.pop('sns_style')
        if 'sns_context' in kwargs:
            kwargs.pop('sns_context')

        col_order = kwargs.pop("col_order", (natsorted(data[self.xfacet].unique()) if self.xfacet else None))
        row_order = kwargs.pop("row_order", (natsorted(data[self.yfacet].unique()) if self.yfacet else None))
        hue_order = kwargs.pop("hue_order", (natsorted(data[self.huefacet].unique()) if self.huefacet else None))
        g = sns.FacetGrid(data,
                          height=height,
                          aspect=aspect,
                          col=(self.xfacet if self.xfacet else None),
                          row=(self.yfacet if self.yfacet else None),
                          hue=(self.huefacet if self.huefacet else None),
                          col_order=col_order,
                          row_order=row_order,
                          hue_order=hue_order,
                          col_wrap=col_wrap,
                          legend_out=False,
                          sharex=sharex,
                          sharey=sharey,
                          despine=despine,
                          palette=palette)

        plot_ret = self._grid_plot(experiment=experiment, grid=g, **kwargs)

        kwargs.update(plot_ret)

        xscale = kwargs.pop("xscale", None)
        yscale = kwargs.pop("yscale", None)
        xlim = kwargs.pop("xlim", None)
        ylim = kwargs.pop("ylim", None)

        for ax in g.axes.flatten():
            if xscale:
                ax.set_xscale(xscale.name, **xscale.get_mpl_params(ax.get_xaxis()))
            if yscale:
                ax.set_yscale(yscale.name, **yscale.get_mpl_params(ax.get_yaxis()))
            if xlim:
                ax.set_xlim(xlim)
            if ylim:
                ax.set_ylim(ylim)

        if sharex:
            fig = plt.gcf()
            fig_x_min = float("inf")
            fig_x_max = float("-inf")

            for ax in fig.get_axes():
                ax_x_min, ax_x_max = ax.get_xlim()
                if ax_x_min < fig_x_min:
                    fig_x_min = ax_x_min
                if ax_x_max > fig_x_max:
                    fig_x_max = ax_x_max

            for ax in fig.get_axes():
                ax.set_xlim(fig_x_min, fig_x_max)

        if sharey:
            fig = plt.gcf()
            fig_y_max = float("-inf")

            for ax in fig.get_axes():
                _, ax_y_max = ax.get_ylim()
                if ax_y_max > fig_y_max:
                    fig_y_max = ax_y_max

            for ax in fig.get_axes():
                ax.set_ylim(None, fig_y_max)

        cmap = kwargs.pop('cmap', None)
        norm = kwargs.pop('norm', None)
        legend_data = kwargs.pop('legend_data', None)

        if legend:
            if cmap and norm:
                plot_ax = plt.gca()
                cax, _ = mpl.colorbar.make_axes(plt.gcf().get_axes())
                mpl.colorbar.ColorbarBase(cax,
                                          cmap=cmap,
                                          norm=norm)
                plt.sca(plot_ax)
            elif self.huefacet:

                current_palette = mpl.rcParams['axes.prop_cycle']

                if util.is_numeric(data[self.huefacet]) and \
                   len(g.hue_names) > len(current_palette):

                    cmap = mpl.colors.ListedColormap(sns.color_palette("husl",
                                                                       n_colors=len(g.hue_names)))
                    hue_scale = util.scale_factory(self.huescale,
                                                   experiment,
                                                   data=data[self.huefacet].values)

                    plot_ax = plt.gca()

                    cax, _ = mpl.colorbar.make_axes(plt.gcf().get_axes())

                    mpl.colorbar.ColorbarBase(cax,
                                              cmap=cmap,
                                              norm=hue_scale.norm(),
                                              label=huelabel)
                    plt.sca(plot_ax)
                else:
                    g.add_legend(title=huelabel, legend_data=legend_data)
                    ax = g.axes.flat[0]
                    legend = ax.legend_
                    self._update_legend(legend)

        if title:
            plt.suptitle(title)

        if xlabel == "":
            xlabel = None

        if ylabel == "":
            ylabel = None

        g.set_axis_labels(xlabel, ylabel)

    def _grid_plot(self, experiment, grid, xlim, ylim, xscale, yscale, **kwargs):
        raise NotImplementedError("You must override _grid_plot in a derived class")

    def _update_legend(self, legend):
        pass


class BaseDataView(BaseView):
    """
    The base class for data views (as opposed to statistics views).

    Attributes
    ----------
    subset : str
        An expression that specifies the subset of the statistic to plot.
        Passed unmodified to `pandas.DataFrame.query`.
    """

    subset = Str

    def plot(self, experiment, **kwargs):
        """
        Plot some data from an experiment.  This function takes care of
        checking for facet name validity and subsetting, then passes the
        underlying dataframe to `BaseView.plot`

        Parameters
        ----------
        min_quantile : float (>0.0 and <1.0, default = 0.001)
            Clip data that is less than this quantile.

        max_quantile : float (>0.0 and <1.0, default = 1.00)
            Clip data that is greater than this quantile.

        Other Parameters
        ----------------
        lim : Dict(Str : (float, float))
            Set the range of each channel's axis.  If unspecified, assume
            that the limits are the minimum and maximum of the clipped data.
            Required.

        scale : Dict(Str : IScale)
            Scale the data on each axis.  Required.

        """
        if experiment is None:
            raise util.CytoflowViewError('experiment',
                                         "No experiment specified")

        if self.xfacet and self.xfacet not in experiment.conditions:
            raise util.CytoflowViewError('xfacet',
                                         "X facet {0} not in the experiment"
                                         .format(self.xfacet))

        if self.yfacet and self.yfacet not in experiment.conditions:
            raise util.CytoflowViewError('yfacet',
                                         "Y facet {0} not in the experiment"
                                         .format(self.yfacet))

        if self.huefacet and self.huefacet not in experiment.conditions:
            raise util.CytoflowViewError('huefacet',
                                         "Hue facet {0} not in the experiment"
                                         .format(self.huefacet))

        min_quantile = kwargs.pop("min_quantile", 0.001)
        max_quantile = kwargs.pop("max_quantile", 1.0)

        if min_quantile < 0.0 or min_quantile > 1:
            raise util.CytoflowViewError('min_quantile',
                                         "min_quantile must be between 0 and 1")

        if max_quantile < 0.0 or max_quantile > 1:
            raise util.CytoflowViewError('max_quantile',
                                         "max_quantile must be between 0 and 1")

        if min_quantile >= max_quantile:
            raise util.CytoflowViewError('min_quantile',
                                         "min_quantile must be less than max_quantile")

        lim = kwargs.get('lim')
        scale = kwargs.get('scale')

        for c in lim.keys():
            if lim[c] is None:
                lim[c] = (experiment[c].quantile(min_quantile),
                          experiment[c].quantile(max_quantile))
            elif isinstance(lim[c], list) or isinstance(lim[c], tuple):
                if len(lim[c]) != 2:
                    raise util.CytoflowError('lim',
                                             'Length of lim[{}] must be 2'
                                             .format(c))
                if lim[c][0] is None:
                    lim[c] = (experiment[c].quantile(min_quantile),
                              lim[c][1])

                if lim[c][1] is None:
                    lim[c] = (lim[c][0],
                              experiment[c].quantile(max_quantile))

            else:
                raise util.CytoflowError('lim',
                                         "lim[{}] is an unknown data type"
                                         .format(c))

            lim[c] = [scale[c].clip(x) for x in lim[c]]

        facets = [x for x in [self.xfacet, self.yfacet, self.huefacet] if x]

        if len(facets) != len(set(facets)):
            raise util.CytoflowViewError(None,
                                         "Can't reuse facets")

        if self.subset:
            try:
                experiment = experiment.query(self.subset)
            except util.CytoflowError as e:
                raise util.CytoflowViewError('subset', str(e)) from e
            except Exception as e:
                raise util.CytoflowViewError('subset',
                                             "Subset string '{0}' isn't valid"
                                             .format(self.subset)) from e

            if len(experiment) == 0:
                raise util.CytoflowViewError('subset',
                                             "Subset string '{0}' returned no events"
                                             .format(self.subset))

        super().plot(experiment,
                     experiment.data,
                     **kwargs)


class Base1DView(BaseDataView):
    """
    A data view that plots data from a single channel.

    Attributes
    ----------
    channel : Str
        The channel to view

    scale : {'linear', 'log'}
        The scale applied to the data before plotting it.
    """

    channel = Str
    scale = util.ScaleEnum

    def plot(self, experiment, **kwargs):
        """
        Parameters
        ----------
        lim : (float, float)
            Set the range of the plot's data axis.

        orientation : {'vertical', 'horizontal'}
        """

        if experiment is None:
            raise util.CytoflowViewError('experiment',
                                         "No experiment specified")

        if not self.channel:
            raise util.CytoflowViewError('channel',
                                         "Must specify a channel")

        if self.channel not in experiment.data:
            raise util.CytoflowViewError('channel',
                                         "Channel {0} not in the experiment"
                                         .format(self.channel))

        scale = kwargs.pop('scale', None)
        if scale is None:
            scale = util.scale_factory(self.scale, experiment, channel=self.channel)

        lim = kwargs.pop("lim", None)

        super().plot(experiment,
                     lim={self.channel: lim},
                     scale={self.channel: scale},
                     **kwargs)


class Base2DView(BaseDataView):
    """
    A data view that plots data from two channels.

    Attributes
    ----------
    xchannel : Str
        The channel to view on the X axis

    ychannel : Str
        The channel to view on the Y axis

    xscale : {'linear', 'log'} (default = 'linear')
        The scales applied to the `xchannel` data before plotting it.

    yscale : {'linear', 'log'} (default = 'linear')
        The scales applied to the `ychannel` data before plotting it.
    """

    xchannel = Str
    xscale = util.ScaleEnum
    ychannel = Str
    yscale = util.ScaleEnum

    def plot(self, experiment, **kwargs):
        """
        Parameters
        ----------
        xlim : (float, float)
            Set the range of the plot's X axis.

        ylim : (float, float)
            Set the range of the plot's Y axis.
        """

        if experiment is None:
            raise util.CytoflowViewError('experiment',
                                         "No experiment specified")

        if not self.xchannel:
            raise util.CytoflowViewError('xchannel',
                                         "Must specify an xchannel")

        if self.xchannel not in experiment.data:
            raise util.CytoflowViewError('xchannel',
                                         "Channel {} not in the experiment"
                                         .format(self.xchannel))

        if not self.ychannel:
            raise util.CytoflowViewError('ychannel',
                                         "Must specify a ychannel")

        if self.ychannel not in experiment.data:
            raise util.CytoflowViewError('ychannel',
                                         "Channel {} not in the experiment"
                                         .format(self.ychannel))

        xscale = kwargs.pop('xscale', None)
        if xscale is None:
            xscale = util.scale_factory(self.xscale, experiment, channel=self.xchannel)

        yscale = kwargs.pop('yscale', None)
        if yscale is None:
            yscale = util.scale_factory(self.yscale, experiment, channel=self.ychannel)

        xlim = kwargs.pop('xlim', None)
        ylim = kwargs.pop('ylim', None)

        super().plot(experiment,
                     lim={self.xchannel: xlim,
                          self.ychannel: ylim},
                     scale={self.xchannel: xscale,
                            self.ychannel: yscale},
                     **kwargs)


class BaseNDView(BaseDataView):
    """
    A data view that plots data from one or more channels.

    Attributes
    ----------
    channels : List(Str)
        The channels to view

    scale : Dict(Str : {"linear", "log"})
        Re-scale the data in the specified channels before plotting.  If a
        channel isn't specified, assume that the scale is linear.
    """

    channels = List(Str)
    scale = Dict(Str, util.ScaleEnum)

    def plot(self, experiment, **kwargs):
        """
        Parameters
        ----------
        lim : Dict(Str : (float, float))
            Set the range of each channel's axis.  If unspecified, assume
            that the limits are the minimum and maximum of the clipped data
        """

        if experiment is None:
            raise util.CytoflowViewError('experiment',
                                         "No experiment specified")

        if len(self.channels) == 0:
            raise util.CytoflowOpError('channels',
                                       "Must set at least one channel")

        if len(self.channels) != len(set(self.channels)):
            raise util.CytoflowOpError('channels',
                                       "Must not duplicate channels")

        for c in self.channels:
            if c not in experiment.data:
                raise util.CytoflowOpError('channels',
                                           "Channel {0} not found in the experiment"
                                           .format(c))

        for c in self.scale:
            if c not in self.channels:
                raise util.CytoflowOpError('scale',
                                           "Scale set for channel {0}, but it isn't "
                                           "in 'channels'"
                                           .format(c))

        scale = {}
        for c in self.channels:
            if c in self.scale:
                scale[c] = util.scale_factory(self.scale[c], experiment, channel=c)
            else:
                scale[c] = util.scale_factory(util.get_default_scale(), experiment, channel=c)

        lim = kwargs.pop("lim", {})
        for c in self.channels:
            if c not in lim:
                lim[c] = None

        super().plot(experiment,
                     lim=lim,
                     scale=scale,
                     **kwargs)


@provides(IView)
class BaseStatisticsView(BaseView):
    """
    The base class for statistics views (as opposed to data views).

    Attributes
    ----------
    variable : str
        The condition that varies when plotting this statistic: used for the
        x axis of line plots, the bar groups in bar plots, etc.

    subset : str
        An expression that specifies the subset of the statistic to plot.
        Passed unmodified to `pandas.DataFrame.query`.

    """

    variable = Str
    subset = Str

    def enum_plots(self, experiment, data):
        """
        Enumerate the named plots we can make from this set of statistics.

        Returns
        -------
        iterator
            An iterator across the possible plot names. The iterator ALSO has an instance
            attribute called ``by``, which holds a list of the facets that are
            not yet set (and thus need to be specified in the plot name.)
        """

        if experiment is None:
            raise util.CytoflowViewError('experiment',
                                         "No experiment specified")

        if not self.variable:
            raise util.CytoflowViewError('variable',
                                         "variable not set")

        if self.variable not in experiment.conditions:
            raise util.CytoflowViewError('variable',
                                         "variable {0} not in the experiment"
                                         .format(self.variable))

        data, facets, names = self._subset_data(data)

        by = list(set(names) - set(facets))

        class plot_enum(object):

            def __init__(self, data, by):
                self.by = by
                self._iter = None
                self._returned = False

                if by:
                    self._iter = data.groupby(by).__iter__()

            def __iter__(self):
                return self

            def __next__(self):
                if self._iter:
                    return next(self._iter)[0]
                else:
                    if self._returned:
                        raise StopIteration
                    else:
                        self._returned = True
                        return None

        return plot_enum(data.reset_index(), by)

    def plot(self, experiment, data, plot_name=None, **kwargs):
        """
        Plot some data from a statistic.

        This function takes care of checking for facet name validity and
        subsetting, then passes the dataframe to `BaseView.plot`

        Parameters
        ----------
        plot_name : str
            If this `IView` can make multiple plots, ``plot_name`` is
            the name of the plot to make.  Must be one of the values retrieved
            from `enum_plots`.

        """

        if experiment is None:
            raise util.CytoflowViewError('experiment',
                                         "No experiment specified")

        if not self.variable:
            raise util.CytoflowViewError('variable',
                                         "variable not set")

        if self.variable not in experiment.conditions:
            raise util.CytoflowViewError('variable',
                                         "variable {0} not in the experiment"
                                         .format(self.variable))

        data, facets, names = self._subset_data(data)

        unused_names = list(set(names) - set(facets))

        if plot_name is not None and not unused_names:
            raise util.CytoflowViewError('plot_name',
                                         "You specified a plot name, but all "
                                         "the facets are already used")

        if unused_names:
            groupby = data.groupby(unused_names)

            if plot_name is None:
                raise util.CytoflowViewError('plot_name',
                                             "You must use facets {} in either the "
                                             "plot variables or the plot name. "
                                             "Possible plot names: {}"
                                             .format(unused_names, list(groupby.groups.keys())))

            if plot_name not in set(groupby.groups.keys()):
                raise util.CytoflowViewError('plot_name',
                                             "Plot {} not from plot_enum; must "
                                             "be one of {}"
                                             .format(plot_name, list(groupby.groups.keys())))

            data = groupby.get_group(plot_name)

        data.reset_index(inplace=True)
        super().plot(experiment, data, **kwargs)

    def _subset_data(self, data):

        if self.subset:
            try:
                data = data.query(self.subset)
            except Exception as e:
                raise util.CytoflowViewError('subset',
                                             "Subset string '{0}' isn't valid"
                                             .format(self.subset)) from e

            if len(data) == 0:
                raise util.CytoflowViewError('subset',
                                             "Subset string '{0}' returned no values"
                                             .format(self.subset))

        names = list(data.index.names)

        for name in names:
            unique_values = data.index.get_level_values(name).unique()
            if len(unique_values) == 1:
                warn("Only one value for level {}; dropping it.".format(name),
                     util.CytoflowViewWarning)
                try:
                    data.index = data.index.droplevel(name)
                except AttributeError as e:
                    raise util.CytoflowViewError(None,
                                                 "Must have more than one "
                                                 "value to plot.") from e

        names = list(data.index.names)

        if self.xfacet and self.xfacet not in data.index.names:
            raise util.CytoflowViewError('xfacet',
                                         "X facet {} not in statistics; must be one of {}"
                                         .format(self.xfacet, data.index.names))

        if self.yfacet and self.yfacet not in data.index.names:
            raise util.CytoflowViewError('yfacet',
                                         "Y facet {} not in statistics; must be one of {}"
                                         .format(self.yfacet, data.index.names))

        if self.huefacet and self.huefacet not in data.index.names:
            raise util.CytoflowViewError('huefacet',
                                         "Hue facet {} not in statistics; must be one of {}"
                                         .format(self.huefacet, data.index.names))

        facets = [x for x in [self.variable, self.xfacet, self.yfacet, self.huefacet] if x]
        if len(facets) != len(set(facets)):
            raise util.CytoflowViewError(None, "Can't reuse facets")

        return data, facets, names


class Base1DStatisticsView(BaseStatisticsView):
    """
    The base class for 1-dimensional statistic views -- ie, the `variable`
    attribute is on the x axis, and the statistic value is on the y axis.

    Attributes
    ----------
    statistic : (str, str)
        The name of the statistic to plot.  Must be a key in the
        `Experiment.statistics` attribute of the `Experiment`
        being plotted.

    error_statistic : (str, str)
        The name of the statistic used to plot error bars.  Must be a key in the
        `Experiment.statistics` attribute of the `Experiment`
        being plotted.

    scale : {'linear', 'log'}
        The scale applied to the data before plotting it.
    """

    statistic = Tuple(Str, Str)
    error_statistic = Tuple(Str, Str)

    scale = util.ScaleEnum

    def enum_plots(self, experiment):
        if experiment is None:
            raise util.CytoflowViewError('experiment',
                                         "No experiment specified")
        data = self._make_data(experiment)
        return super().enum_plots(experiment, data)

    def plot(self, experiment, plot_name=None, **kwargs):
        """
        Parameters
        ----------
        orientation : {'vertical', 'horizontal'}

        lim : (float, float)
            Set the range of the plot's axis.
        """

        if experiment is None:
            raise util.CytoflowViewError('experiment',
                                         "No experiment specified")

        data = self._make_data(experiment)

        if not self.variable:
            raise util.CytoflowViewError('variable',
                                         "variable not set")

        if self.variable not in experiment.conditions:
            raise util.CytoflowViewError('variable',
                                         "variable {0} not in the experiment"
                                         .format(self.variable))

        scale = util.scale_factory(self.scale,
                                   experiment,
                                   statistic=self.statistic,
                                   error_statistic=self.error_statistic)

        super().plot(experiment,
                     data,
                     plot_name=plot_name,
                     scale=scale,
                     **kwargs)

    def _make_data(self, experiment):
        if experiment is None:
            raise util.CytoflowViewError('experiment', "No experiment specified")

        if not self.statistic:
            raise util.CytoflowViewError('statistic', "Statistic not set")

        if self.statistic not in experiment.statistics:
            raise util.CytoflowViewError('statistic',
                                         "Can't find the statistic {} in the experiment"
                                         .format(self.statistic))
        else:
            stat = experiment.statistics[self.statistic]

        if not util.is_numeric(stat):
            raise util.CytoflowViewError('statistic',
                                         "Statistic must be numeric")

        if self.error_statistic[0]:
            if self.error_statistic not in experiment.statistics:
                raise util.CytoflowViewError('error_statistic',
                                             "Can't find the error statistic in the experiment")
            else:
                error_stat = experiment.statistics[self.error_statistic]
        else:
            error_stat = None

        if error_stat is not None:

            if set(stat.index.names) != set(error_stat.index.names):
                raise util.CytoflowViewError('error_statistic',
                                             "Data statistic and error statistic "
                                             "don't have the same index.")

            try:
                error_stat.index = error_stat.index.reorder_levels(stat.index.names)
                error_stat.sort_index(inplace=True)
            except AttributeError:
                pass

            if not stat.index.equals(error_stat.index):
                raise util.CytoflowViewError('error_statistic',
                                             "Data statistic and error statistic "
                                             " don't have the same index.")

            if stat.name == error_stat.name:
                raise util.CytoflowViewError('error_statistic',
                                             "Data statistic and error statistic can "
                                             "not have the same name.")

        data = pd.DataFrame(index=stat.index)
        data[stat.name] = stat

        if error_stat is not None:
            data[error_stat.name] = error_stat

        return data


class Base2DStatisticsView(BaseStatisticsView):
    """
    The base class for 2-dimensional statistic views -- ie, the `variable`
    attribute varies independently, and the corresponding values from the x and
    y statistics are plotted on the x and y axes.

    Attributes
    ----------
    xstatistic : (str, str)
        The name of the statistic to plot on the X axis.  Must be a key in the
        `Experiment.statistics` attribute of the `Experiment`
        being plotted.

    ystatistic : (str, str)
        The name of the statistic to plot on the Y axis.  Must be a key in the
        `Experiment.statistics` attribute of the `Experiment`
        being plotted.

    x_error_statistic : (str, str)
        The name of the statistic used to plot error bars on the X axis.
        Must be a key in the `Experiment.statistics` attribute of the `Experiment`
        being plotted.

    y_error_statistic : (str, str)
        The name of the statistic used to plot error bars on the Y axis.
        Must be a key in the `Experiment.statistics` attribute of the `Experiment`
        being plotted.

    xscale : {'linear', 'log'}
        The scale applied to `xstatistic` before plotting it.

    yscale : {'linear', 'log'}
        The scale applied to `ystatistic` before plotting it.
    """

    xstatistic = Tuple(Str, Str)
    ystatistic = Tuple(Str, Str)
    x_error_statistic = Tuple(Str, Str)
    y_error_statistic = Tuple(Str, Str)

    xscale = util.ScaleEnum
    yscale = util.ScaleEnum

    def enum_plots(self, experiment):
        if experiment is None:
            raise util.CytoflowViewError('experiment',
                                         "No experiment specified")
        data = self._make_data(experiment)
        return super().enum_plots(experiment, data)

    def plot(self, experiment, plot_name=None, **kwargs):
        """
        Parameters
        ----------
        xlim : (float, float)
            Set the range of the plot's X axis.

        ylim : (float, float)
            Set the range of the plot's Y axis.
        """
        if experiment is None:
            raise util.CytoflowViewError('experiment',
                                         "No experiment specified")

        data = self._make_data(experiment)

        xscale = util.scale_factory(self.xscale,
                                    experiment,
                                    statistic=self.xstatistic,
                                    error_statistic=self.x_error_statistic)

        yscale = util.scale_factory(self.yscale,
                                    experiment,
                                    statistic=self.ystatistic,
                                    error_statistic=self.y_error_statistic)

        super().plot(experiment,
                     data,
                     plot_name,
                     xscale=xscale,
                     yscale=yscale,
                     **kwargs)

    def _make_data(self, experiment):
        if experiment is None:
            raise util.CytoflowViewError('experiment',
                                         "No experiment specified")

        if not self.xstatistic:
            raise util.CytoflowViewError('xstatistic',
                                         "X Statistic not set")

        if self.xstatistic not in experiment.statistics:
            raise util.CytoflowViewError('xstatistic',
                                         "Can't find the statistic {} in the experiment"
                                         .format(self.xstatistic))
        else:
            xstat = experiment.statistics[self.xstatistic]

        if not util.is_numeric(xstat):
            raise util.CytoflowViewError('xstatistic',
                                         "X statistic must be numeric")

        if self.x_error_statistic[0]:
            if self.x_error_statistic not in experiment.statistics:
                raise util.CytoflowViewError('x_error_statistic',
                                             "Can't find the X error statistic in the experiment")
            else:
                x_error_stat = experiment.statistics[self.x_error_statistic]
        else:
            x_error_stat = None

        if x_error_stat is not None:

            if set(xstat.index.names) != set(x_error_stat.index.names):
                raise util.CytoflowViewError('x_error_statistic',
                                             "X data statistic and error statistic "
                                             "don't have the same index.")

            try:
                x_error_stat.index = x_error_stat.index.reorder_levels(xstat.index.names)
                x_error_stat.sort_index(inplace=True)
            except AttributeError:
                pass

            if not xstat.index.equals(x_error_stat.index):
                raise util.CytoflowViewError('x_error_statistic',
                                             "X data statistic and error statistic "
                                             " don't have the same index.")

            if xstat.name == x_error_stat.name:
                raise util.CytoflowViewError('x_error_statistic',
                                             "X data statistic and error statistic can "
                                             "not have the same name.")

        if not self.ystatistic:
            raise util.CytoflowViewError('ystatistic',
                                         "Y statistic not set")

        if self.ystatistic not in experiment.statistics:
            raise util.CytoflowViewError('ystatistic',
                                         "Can't find the Y statistic {} in the experiment"
                                         .format(self.ystatistic))
        else:
            ystat = experiment.statistics[self.ystatistic]

        if not util.is_numeric(ystat):
            raise util.CytoflowViewError('ystatistic',
                                         "Y statistic must be numeric")

        if self.y_error_statistic[0]:
            if self.y_error_statistic not in experiment.statistics:
                raise util.CytoflowViewError('y_error_statistic',
                                             "Can't find the Y error statistic in the experiment")
            else:
                y_error_stat = experiment.statistics[self.y_error_statistic]
        else:
            y_error_stat = None

        if y_error_stat is not None:

            if set(ystat.index.names) != set(y_error_stat.index.names):
                raise util.CytoflowViewError('y_error_statistic',
                                             "Y data statistic and error statistic "
                                             "don't have the same index.")

            try:
                y_error_stat.index = y_error_stat.index.reorder_levels(ystat.index.names)
                y_error_stat.sort_index(inplace=True)
            except AttributeError:
                pass

            if not ystat.index.equals(y_error_stat.index):
                raise util.CytoflowViewError('y_error_statistic',
                                             "Y data statistic and error statistic "
                                             " don't have the same index.")

            if ystat.name == y_error_stat.name:
                raise util.CytoflowViewError('y_error_statistic',
                                             "Data statistic and error statistic can "
                                             "not have the same name.")

        if xstat.name == ystat.name:
            raise util.CytoflowViewError('ystatistic',
                                         "X and Y statistics can "
                                         "not have the same name.")

        if set(xstat.index.names) != set(ystat.index.names):
            raise util.CytoflowViewError('ystatistic',
                                         "X and Y data statistics "
                                         "don't have the same index.")

        try:
            ystat.index = ystat.index.reorder_levels(xstat.index.names)
            ystat.sort_index(inplace=True)
        except AttributeError:
            pass

        intersect_idx = xstat.index.intersection(ystat.index)
        xstat = xstat.reindex(intersect_idx)
        xstat.sort_index(inplace=True)
        ystat = ystat.reindex(intersect_idx)
        ystat.sort_index(inplace=True)

        data = pd.DataFrame(index=xstat.index)
        data[xstat.name] = xstat
        data[ystat.name] = ystat

        if x_error_stat is not None:
            data[x_error_stat.name] = x_error_stat

        if y_error_stat is not None:
            data[y_error_stat.name] = y_error_stat

        return data


@provides(IView)
class HistogramView(Base1DView):
    """
    Plots a one-channel histogram

    Attributes
    ----------

    Examples
    --------

    Make a little data set.

    .. plot::
        :context: close-figs

        >>> import cytoflow as flow
        >>> import_op = flow.ImportOp()
        >>> import_op.tubes = [flow.Tube(file = "Plate01/RFP_Well_A3.fcs",
        ...                              conditions = {'Dox' : 10.0}),
        ...                    flow.Tube(file = "Plate01/CFP_Well_A4.fcs",
        ...                              conditions = {'Dox' : 1.0})]
        >>> import_op.conditions = {'Dox' : 'float'}
        >>> ex = import_op.apply()

    Plot a histogram

    .. plot::
        :context: close-figs

        >>> flow.HistogramView(channel = 'Y2-A',
        ...                    scale = 'log',
        ...                    huefacet = 'Dox').plot(ex)
    """

    id = Constant("edu.mit.synbio.cytoflow.view.histogram")
    friendly_id = Constant("Histogram")

    def plot(self, experiment, **kwargs):
        """
        Plot a faceted histogram view of a channel

        Parameters
        ----------
        num_bins : int
            The number of bins to plot in the histogram.  Clipped to [100, 1000]

        histtype : {'stepfilled', 'step', 'bar'}
            The type of histogram to draw.  ``stepfilled`` is the default, which
            is a line plot with a color filled under the curve.

        density: bool
            If `True`, re-scale the histogram to form a probability density
            function, so the area under the histogram is 1.

        orientation : {'horizontal', 'vertical'}
            The orientation of the histogram.  ``horizontal`` gives a histogram
            with the intensity on the Y axis and the count on the X axis;
            default is ``vertical``.

        linewidth : float
            The width of the histogram line (in points)

        linestyle : ['-' | '--' | '-.' | ':' | "None"]
            The style of the line to plot

        alpha : float (default = 0.5)
            The alpha blending value, between 0 (transparent) and 1 (opaque).


        Notes
        -----
        Other ``kwargs`` are passed to `matplotlib.pyplot.hist <https://matplotlib.org/devdocs/api/_as_gen/matplotlib.pyplot.hist.html>`_


        """
        ylabel = 'Density' if kwargs.get('density', False) else 'Count'
        if kwargs.get('orientation', 'vertical') == 'vertical':
            kwargs.setdefault('xlabel', self.channel)
            kwargs.setdefault('ylabel', ylabel)
        else:
            kwargs.setdefault('xlabel', ylabel)
            kwargs.setdefault('ylabel', self.channel)
        super().plot(experiment, **kwargs)

    def _grid_plot(self, experiment, grid, **kwargs):

        kwargs.setdefault('histtype', 'stepfilled')
        kwargs.setdefault('alpha', 0.5)
        kwargs.setdefault('antialiased', True)

        scale = kwargs.pop('scale')[self.channel]
        lim = kwargs.pop('lim')[self.channel]

        scaled_data = scale(experiment[self.channel])
        num_bins = kwargs.pop('num_bins', util.num_hist_bins(scaled_data))
        num_bins = util.num_hist_bins(scaled_data) if num_bins is None else num_bins

        num_bins = max(min(num_bins, 1000), 100)

        if (self.huefacet
            and "bins" in experiment.metadata[self.huefacet]
            and experiment.metadata[self.huefacet]["bin_scale"] == self.scale):

            bins = experiment.metadata[self.huefacet]["bins"]
            scaled_bins = scale(bins)

            num_hues = len(experiment[self.huefacet].unique())
            bins_per_hue = math.floor(num_bins / num_hues)

            if bins_per_hue == 1:
                new_bins = scaled_bins
            else:
                new_bins = []
                for idx in range(1, len(scaled_bins)):
                    new_bins = np.append(new_bins,
                                         np.linspace(scaled_bins[idx - 1],
                                                     scaled_bins[idx],
                                                     bins_per_hue + 1,
                                                     endpoint=False))

            bins = scale.inverse(new_bins)
        else:
            xmin = bottleneck.nanmin(scaled_data)
            xmax = bottleneck.nanmax(scaled_data)
            bins = scale.inverse(np.linspace(xmin, xmax, num=int(num_bins), endpoint=True))

        kwargs.setdefault('bins', bins)
        kwargs.setdefault('orientation', 'vertical')

        if ('linewidth' not in kwargs) or ('linewidth' in kwargs and kwargs['linewidth'] is None):
            kwargs['linewidth'] = 0 if kwargs['histtype'] == "stepfilled" else 2

        count_max = []

        def hist_lims(*args, **kwargs):
            bins = kwargs.get('bins')
            new_args = []
            for x in args:
                x = x[x > bins[0]]
                x = x[x < bins[-1]]
                new_args.append(x)

            if scale.name != "linear" and kwargs.get("density"):
                kwargs["density"] = False
                counts, _ = np.histogram(new_args, bins=kwargs["bins"])
                kwargs["weights"] = counts / np.sum(counts)
                n, _, _ = plt.hist(kwargs["bins"][:-1], **kwargs)
            else:
                n, _, _ = plt.hist(*new_args, **kwargs)

            count_max.append(max(n))

        grid.map(hist_lims, self.channel, **kwargs)

        ret = {}
        if kwargs['orientation'] == 'vertical':
            ret['xscale'] = scale
            ret['xlim'] = lim
            ret['ylim'] = (0, 1.05 * max(count_max))
        else:
            ret['yscale'] = scale
            ret['ylim'] = lim
            ret['xlim'] = (0, 1.05 * max(count_max))

        return ret


@provides(IView)
class ScatterplotView(Base2DView):
    """
    Plots a 2-d scatterplot.

    Attributes
    ----------

    Examples
    --------

    Make a little data set.

    .. plot::
        :context: close-figs

        >>> import cytoflow as flow
        >>> import_op = flow.ImportOp()
        >>> import_op.tubes = [flow.Tube(file = "Plate01/RFP_Well_A3.fcs",
        ...                              conditions = {'Dox' : 10.0}),
        ...                    flow.Tube(file = "Plate01/CFP_Well_A4.fcs",
        ...                              conditions = {'Dox' : 1.0})]
        >>> import_op.conditions = {'Dox' : 'float'}
        >>> ex = import_op.apply()

    Plot a density plot

    .. plot::
        :context: close-figs

        >>> flow.ScatterplotView(xchannel = 'V2-A',
        ...                      xscale = 'log',
        ...                      ychannel = 'Y2-A',
        ...                      yscale = 'log',
        ...                      huefacet = 'Dox').plot(ex)

    """

    id = Constant('edu.mit.synbio.cytoflow.view.scatterplot')
    friend_id = Constant("Scatter Plot")

    def plot(self, experiment, **kwargs):
        """
        Plot a faceted scatter plot view of a channel

        Parameters
        ----------

        alpha : float (default = 0.25)
            The alpha blending value, between 0 (transparent) and 1 (opaque).

        s : int (default = 2)
            The size in points^2.

        marker : a matplotlib marker style, usually a string
            Specfies the glyph to draw for each point on the scatterplot.
            See `matplotlib.markers <http://matplotlib.org/api/markers_api.html#module-matplotlib.markers>`_ for examples.  Default: 'o'


        Notes
        -----
        Other ``kwargs`` are passed to `matplotlib.pyplot.scatter <https://matplotlib.org/devdocs/api/_as_gen/matplotlib.pyplot.scatter.html>`_

        """

        super().plot(experiment, **kwargs)

    def _grid_plot(self, experiment, grid, **kwargs):

        kwargs.setdefault('alpha', 0.25)
        kwargs.setdefault('s', 2)
        kwargs.setdefault('marker', 'o')
        kwargs.setdefault('antialiased', True)

        lim = kwargs.pop('lim')
        xlim = lim[self.xchannel]
        ylim = lim[self.ychannel]

        scale = kwargs.pop('scale')
        xscale = scale[self.xchannel]
        yscale = scale[self.ychannel]

        grid.map(plt.scatter, self.xchannel, self.ychannel, **kwargs)

        return dict(xlim=xlim,
                    xscale=xscale,
                    ylim=ylim,
                    yscale=yscale)

    def _update_legend(self, legend):
        for lh in legend.legendHandles:
            lh.set_alpha(0.5)
            lh.set_sizes([10.0])


@provides(IView)
class DensityView(Base2DView):
    """
    Plots a 2-d density plot.

    Attributes
    ----------

    huefacet : None
        You must leave the hue facet unset!

    Examples
    --------

    Make a little data set.

    .. plot::
        :context: close-figs

        >>> import cytoflow as flow
        >>> import_op = flow.ImportOp()
        >>> import_op.tubes = [flow.Tube(file = "Plate01/RFP_Well_A3.fcs",
        ...                              conditions = {'Dox' : 10.0}),
        ...                    flow.Tube(file = "Plate01/CFP_Well_A4.fcs",
        ...                              conditions = {'Dox' : 1.0})]
        >>> import_op.conditions = {'Dox' : 'float'}
        >>> ex = import_op.apply()

    Plot a density plot

    .. plot::
        :context: close-figs

        >>> flow.DensityView(xchannel = 'V2-A',
        ...                 xscale = 'log',
        ...                 ychannel = 'Y2-A',
        ...                 yscale = 'log').plot(ex)

    The same plot, smoothed, with a log color scale.  *Note - you can change*
    *the hue scale, even if you don't have control over the hue facet!*

    .. plot::
        :context: close-figs

        >>> flow.DensityView(xchannel = 'V2-A',
        ...                  xscale = 'log',
        ...                  ychannel = 'Y2-A',
        ...                  yscale = 'log',
        ...                  huescale = 'log').plot(ex, smoothed = True)
    """

    id = Constant('edu.mit.synbio.cytoflow.view.density')
    friend_id = Constant("Density Plot")

    huefacet = Constant(None)

    def plot(self, experiment, **kwargs):
        """
        Plot a faceted density plot view of a channel

        Parameters
        ----------
        gridsize : int
            The size of the grid on each axis.  Default = 50

        smoothed : bool
            Should the resulting mesh be smoothed?

        smoothed_sigma : int
            The standard deviation of the smoothing kernel.  default = 1.

        cmap : cmap
            An instance of matplotlib.colors.Colormap.  By default, the
            ``viridis`` colormap is used

        under_color : matplotlib color
            Sets the color to be used for low out-of-range values.

        bad_color : matplotlib color
            Set the color to be used for masked values.

        Notes
        -----
        Other ``kwargs`` are passed to `matplotlib.axes.Axes.pcolormesh`

        """

        super().plot(experiment, **kwargs)

    def _grid_plot(self, experiment, grid, **kwargs):

        kwargs.setdefault('antialiased', False)
        kwargs.setdefault('linewidth', 0)
        kwargs.setdefault('edgecolors', 'face')
        kwargs.setdefault('cmap', plt.get_cmap('viridis'))

        lim = kwargs.pop('lim')
        xlim = lim[self.xchannel]
        ylim = lim[self.ychannel]

        scale = kwargs.pop('scale')
        xscale = scale[self.xchannel]
        yscale = scale[self.ychannel]

        cmap = copy.copy(kwargs['cmap'])

        under_color = kwargs.pop('under_color', None)
        if under_color is not None:
            cmap.set_under(color=under_color)
        else:
            cmap.set_under(cmap(0.0))

        bad_color = kwargs.pop('bad_color', None)
        if bad_color is not None:
            cmap.set_bad(color=bad_color)
        else:
            cmap.set_bad(color=cmap(0.0))

        gridsize = kwargs.pop('gridsize', 50)

        xbins = xscale.inverse(np.linspace(xscale(xlim[0]), xscale(xlim[1]), gridsize))
        ybins = yscale.inverse(np.linspace(yscale(ylim[0]), yscale(ylim[1]), gridsize))

        if 'norm' not in kwargs:
            data_max = 0
            for _, data_ijk in grid.facet_data():
                x = data_ijk[self.xchannel]
                y = data_ijk[self.ychannel]
                h, _, _ = np.histogram2d(x, y, bins=[xbins, ybins])
                data_max = max(data_max, h.max())

            hue_scale = util.scale_factory(self.huescale,
                                           experiment,
                                           data=np.array([1, data_max]))
            kwargs['norm'] = hue_scale.norm()

        grid.map(_densityplot, self.xchannel, self.ychannel, xbins=xbins, ybins=ybins, **kwargs)

        return dict(xlim=xlim,
                    xscale=xscale,
                    ylim=ylim,
                    yscale=yscale,
                    cmap=kwargs['cmap'],
                    norm=kwargs['norm'])


def _densityplot(x, y, xbins, ybins, **kwargs):

    h, X, Y = np.histogram2d(x, y, bins=[xbins, ybins])

    smoothed = kwargs.pop('smoothed', False)
    smoothed_sigma = kwargs.pop('smoothed_sigma', 1)

    if smoothed:
        h = scipy.ndimage.filters.gaussian_filter(h, sigma=smoothed_sigma)

    ax = plt.gca()
    ax.pcolormesh(X, Y, h.T, **kwargs)


util.expand_class_attributes(HistogramView)
util.expand_method_parameters(HistogramView, HistogramView.plot)

util.expand_class_attributes(ScatterplotView)
util.expand_method_parameters(ScatterplotView, ScatterplotView.plot)

util.expand_class_attributes(DensityView)
util.expand_method_parameters(DensityView, DensityView.plot)
