"""
Xiaoyu-CF operations module.

This module contains all operation classes implementing `IOperation`
and their associated views, merged from cytoflow/operations/*.py.
"""

import collections
import math
import os
import warnings
from pathlib import Path
from warnings import warn

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.widgets import PolygonSelector, SpanSelector
from natsort import natsorted
import numpy as np
import pandas as pd
import scipy.ndimage.filters
import scipy.optimize
import scipy.stats
from traits.api import (HasTraits, HasStrictTraits, Str, List, Any, Dict,
                        File, Constant, Enum, Int, Float, Instance, Bool,
                        provides, observe, Tuple, Array, DelegatesTo,
                        Property, Interface)

from . import _utility as util
from ._experiment import Experiment
from ._fcsparser import parse
from ._views import (IView, ISelectionView, HistogramView, ScatterplotView,
                      DensityView, BaseDataView, Base1DView, Base2DView)


# ---------------------------------------------------------------------------
# IOperation interface
# ---------------------------------------------------------------------------

class IOperation(Interface):
    """
    The basic interface for an operation on cytometry data.

    Attributes
    ----------
    id : Str
        a unique identifier for this class. prefix: ``edu.mit.synbio.cytoflow.operations``

    friendly_id : Str
        The operation's human-readable id (like ``Threshold`` or ``K-means``).
        Used for UI implementations.

    name : Str
        The name of this IOperation instance (like ``Debris_Filter``).  Useful
        for UI implementations; sometimes used for naming gates' metadata
    """

    # interface traits
    id = Constant('FIXME')
    friendly_id = Constant('FIXME')

    def estimate(self, experiment, subset = None):
        """
        Estimate this operation's parameters from some data.

        For operations that are data-driven (for example, a mixture model),
        estimate the operation's parameters from an experiment.

        Parameters
        ----------
        experiment : `Experiment`
            the `Experiment` to use in the estimation.

        subset : Str (optional)
            a string passed to `pandas.DataFrame.query` to select the subset
            of data on which to run the parameter estimation.

        Raises
        ------
        `CytoflowOpException`
            If the operation can't be be completed because of bad op
            parameters.
        """

    def apply(self, experiment):
        """
        Apply an operation to an experiment.

        Parameters
        ----------
        experiment : `Experiment`
            the `Experiment` to apply this op to

        Returns
        -------
        `Experiment`
            the old `Experiment` with this operation applied

        Raises
        ------
        `CytoflowOpException`
            If the operation can't be be completed because of bad op
            parameters.
        """

    def default_view(self, **kwargs):
        """
        Many operations have a "default" view.  This can either be a diagnostic
        for the operation's `estimate` method, an interactive for setting
        gates, etc.  Frequently it makes sense to link the properties of the
        view to the properties of the `IOperation`; sometimes,
        `default_view` is the only way to get the view (ie, it's not
        useful when it doesn't reference an `IOperation` instance.)

        Parameters
        ----------
        **kwargs : Dict
            The keyword args passed to the view's constructor

        Returns
        -------
        `IView`
            the `IView` instance
        """


# ---------------------------------------------------------------------------
# Helper functions (from import_op.py)
# ---------------------------------------------------------------------------

def autodetect_name_metadata(filename, data_set = 0):
    """
    Tries to determine whether the channel names should come from
    $PnN or $PnS.

    Parameters
    ----------
    filename : string
        The name of the FCS file to operate on

    data_set : int (optional, default = 0)
        Which data set in the FCS file to operate on

    Returns
    -------
    The name of the parameter to parse channel names from,
    either "$PnN" or "$PnS"

    """
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            metadata = parse(filename,
                             data_set = data_set,
                             meta_data_only = True,
                             reformat_meta = True)
    except Exception as e:
        warnings.warn("Trouble getting metadata from {}: {}".format(filename, str(e)),
                      util.CytoflowWarning)
        return '$PnS'

    meta_channels = metadata["_channels_"]

    if "$PnN" in meta_channels and not "$PnS" in meta_channels:
        name_metadata = "$PnN"
    elif "$PnN" not in meta_channels and "$PnS" in meta_channels:
        name_metadata = "$PnS"
    else:
        PnN = meta_channels["$PnN"]
        PnS = meta_channels["$PnS"]

        if None in PnS:
            name_metadata = "$PnN"

        if (len(set(PnN)) == len(PnN) and
            len(set(PnS)) != len(PnS)):
            name_metadata = "$PnN"
        elif (len(set(PnN)) != len(PnN) and
              len(set(PnS)) == len(PnS)):
            name_metadata = "$PnS"
        else:
            name_metadata = "$PnS"

    return name_metadata


def check_tube(filename, experiment, data_set = 0):
    """
    Check to see if an FCS file can be parsed, and that the tube's parameters
    are the same as those already in the `Experiment`.  If not, raises `CytoflowError`.
    At the moment, only checks $PnV, the detector voltages.

    Parameters
    ----------
    filename : string
        An FCS filename

    experiment : `Experiment`
        The `Experiment` to check ``filename`` against.

    data_set : int (optional, default = 0)
        The FCS standard allows for multiple data sets; ``data_set``
        specifies which one to check.

    Raises
    ------
    `CytoflowError`
        If the FCS file can't be read, or if the voltages in ``filename`` are
        different than those in ``experiment``.

    """

    if experiment is None:
        raise util.CytoflowError("No experiment specified")

    ignore_v = experiment.metadata['ignore_v']

    try:
        tube_meta = parse( filename,
                           channel_naming = experiment.metadata["name_metadata"],
                           data_set = data_set,
                           meta_data_only = True,
                           reformat_meta = True)
    except Exception as e:
        raise util.CytoflowError("FCS reader threw an error reading metadata "
                                 "for tube {0}"
                                 .format(filename)) from e

    if not set([experiment.metadata[c]["fcs_name"] for c in experiment.channels]) <= set(tube_meta["_channel_names_"]):
        raise util.CytoflowError("Tube {0} doesn't have the same channels as the rest of the experiment"
                                 .format(filename))

    tube_channels = tube_meta["_channels_"]
    tube_channels.set_index(experiment.metadata["name_metadata"],
                            inplace = True)

    for channel in experiment.channels:
        fcs_name = experiment.metadata[channel]["fcs_name"]
        if "voltage" in experiment.metadata[channel]:
            if not "$PnV" in tube_channels.loc[fcs_name]:
                raise util.CytoflowError("Didn't find a voltage for channel {0}" \
                                   "in tube {1}".format(channel, filename))

            old_v = experiment.metadata[channel]["voltage"]
            new_v = tube_channels.loc[fcs_name]['$PnV']

            if old_v != new_v and not channel in ignore_v:
                raise util.CytoflowError("Tube {0}, channel {1} doesn't have the same voltages as the rest of the experiment"
                                    .format(filename, channel))


def parse_tube(filename, experiment = None, data_set = 0, metadata_only = False):
    """
    Parses an FCS file.  A thin wrapper over `parse`.

    Parameters
    ----------
    filename : string
        The file to parse.

    experiment : `Experiment` (optional, default: None)
        If provided, check the tube's parameters against this experiment
        first.

    data_set : int (optional, default: 0)
        Which data set in the FCS file to parse?

    metadata_only : bool (optional, default: False)
        If ``True``, only parse the metadata.  Because this is at the beginning
        of the FCS file, this happens much faster than parsing the entire file.

    Returns
    -------
    tube_metadata : dict
        The metadata from the FCS file

    tube_data : `pandas.DataFrame`
        The actual tabular data from the FCS file.  Each row is an event, and
        each column is a channel.

    """

    if experiment:
        check_tube(filename, experiment)
        name_metadata = experiment.metadata["name_metadata"]
    else:
        name_metadata = '$PnS'


    try:
        if metadata_only:
            tube_data = None
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                tube_meta = parse(
                                filename,
                                meta_data_only = True,
                                data_set = data_set,
                                channel_naming = name_metadata)
        else:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                tube_meta, tube_data = parse(
                                        filename,
                                        meta_data_only = metadata_only,
                                        data_set = data_set,
                                        channel_naming = name_metadata)
    except Exception as e:
        raise util.CytoflowError("FCS reader threw an error reading data for tube {}"
                                 .format(filename)) from e

    del tube_meta['__header__']

    return tube_meta, tube_data


# ---------------------------------------------------------------------------
# Tube class (from import_op.py)
# ---------------------------------------------------------------------------

class Tube(HasTraits):
    """
    Represents a tube or plate well we want to import.

    Attributes
    ----------
    file : File
        The file name of the FCS file to import

    conditions : Dict(Str, Any)
        A dictionary containing this tube's experimental conditions.  Keys
        are condition names, values are condition values.

    Examples
    --------
    >>> tube1 = flow.Tube(file = 'RFP_Well_A3.fcs', conditions = {"Dox" : 10.0})
    >>> tube2 = flow.Tube(file='CFP_Well_A4.fcs', conditions = {"Dox" : 1.0})
    """

    file = File

    conditions = Dict(Str, Any)

    def conditions_equal(self, other):
        return other and len(set(self.conditions.items()) ^
                   set(other.conditions.items())) == 0

    def __eq__(self, other):
        return other and (self.file == other.file and
                self.conditions == other.conditions)

    def __hash__(self):
        return hash((self.file,
                     (frozenset(self.conditions.keys()),
                       frozenset(self.conditions.values()))))


# ---------------------------------------------------------------------------
# ImportOp (from import_op.py)
# ---------------------------------------------------------------------------

@provides(IOperation)
class ImportOp(HasStrictTraits):
    """
    An operation for importing data and making an `Experiment`.

    To use, set the `conditions` dict to a mapping between condition name
    and NumPy ``dtype``.  Useful dtypes include ``category``, ``float``,
    ``int``, ``bool``.

    Next, set `tubes` to a list of `Tube` containing FCS filenames
    and the corresponding conditions.

    If you would rather not analyze every single event in every FCS file,
    set `events` to the number of events from each FCS file you want to
    load.

    Call `apply` to load the data.  The usual ``experiment`` parameter
    can be ``None``.

    Attributes
    ----------
    conditions : Dict(Str, Str)
        A dictionary mapping condition names (keys) to NumPy ``dtype``s (values).
        Useful ``dtype``s include ``category``, ``float``, ``int``, and ``bool``.

    tubes : List(Tube)
        A list of ``Tube`` instances, which map FCS files to their corresponding
        experimental conditions.  Each ``Tube`` must have a
        ``Tube.conditions`` dict whose keys match those of
        `conditions`.

    channels : Dict(Str, Str)
        If you only need a subset of the channels available in the data set,
        specify them here.  Each ``(key, value)`` pair specifies a channel to
        include in the output experiment.  The key is the channel name in the
        FCS file, and the value is the name of the channel in the Experiment.
        You can use this to rename channels as you import data (because flow
        channel names are frequently not terribly informative.)  New channel
        names must be valid Python identifiers: start with a letter or ``_``, and
        all characters must be letters, numbers or ``_``.  If `channels` is
        empty, load all channels in the FCS files.

    events : Int
        If not None, import only a random subset of events of size `events`.
        Presumably the analysis will go faster but less precisely; good for
        interactive data exploration.  Then, unset `events` and re-run
        the analysis non-interactively.

    name_metadata : {None, "$PnN", "$PnS"} (default = None)
        Which FCS metadata is the channel name?  If ``None``, attempt to
        autodetect.

    data_set : Int (default = 0)
        The FCS standard allows you to encode multiple data sets in a single
        FCS file.  Some software (such as the Beckman-Coulter software)
        also encode the same data in two different formats -- for example,
        FCS2.0 and FCS3.0.  To access a data set other than the first one,
        set `data_set` to the 0-based index of the data set you
        would like to use.  This will be used for *all FCS files imported by
        this operation.*

    ignore_v : List(Str)
        `cytoflow` is designed to operate on an `Experiment` containing
        tubes that were all collected under the same instrument settings.
        In particular, the same PMT voltages ensure that data can be
        compared across samples.

        *Very rarely*, you may need to set up an `Experiment` with
        different voltage settings on different `Tube` instances.  This is likely
        only to be the case when you are trying to figure out which voltages
        should be used in future experiments.  If so, set `ignore_v` to a
        list of channel names to ignore particular channels.

        .. warning::

            THIS WILL BREAK REAL EXPERIMENTS

    Examples
    --------
    >>> tube1 = flow.Tube(file = 'RFP_Well_A3.fcs', conditions = {"Dox" : 10.0})
    >>> tube2 = flow.Tube(file='CFP_Well_A4.fcs', conditions = {"Dox" : 1.0})
    >>> import_op = flow.ImportOp(conditions = {"Dox" : "float"},
    ...                           tubes = [tube1, tube2])
    >>> ex = import_op.apply()
    """

    id = Constant("edu.mit.synbio.cytoflow.operations.import")
    friendly_id = Constant("Import")
    name = Constant("Import Data")

    conditions = Dict(Str, Str)

    tubes = List(Tube)

    channels = Dict(Str, Str)

    name_metadata = Enum(None, "$PnN", "$PnS")

    data_set = Int(0)

    events = Int(None)

    ignore_v = List(Str)

    def apply(self, experiment = None, metadata_only = False):
        """
        Load a new `Experiment`.

        Parameters
        ----------
        experiment : Experiment
            Ignored

        metadata_only : bool (default = False)
            Only "import" the metadata, creating an Experiment with all the
            expected metadata and structure but 0 events.

        Returns
        -------
        Experiment
            The new `Experiment`.  New channels have the following
            metadata:

            - **voltage** - int
                The voltage that this channel was collected at.  Determined
                by the ``$PnV`` field from the first FCS file.

            - **range** - int
                The maximum range of this channel.  Determined by the ``$PnR``
                field from the first FCS file.

            New experimental conditions do not have **voltage** or **range**
            metadata, obviously.  Instead, they have **experiment** set to
            ``True``, to distinguish the experimental variables from the
            conditions that were added by gates, etc.

            If `ignore_v` is set, it is added as a key to the
            `Experiment`-wide metadata.

        """

        if not self.tubes or len(self.tubes) == 0:
            raise util.CytoflowOpError('tubes',
                                       "Must specify some tubes!")

        if self.channels:
            for old_name, new_name in self.channels.items():
                if not new_name:
                    raise util.CytoflowOpError('channels',
                                               'Can\'t leave the new name for {} empty'
                                               .format(old_name))
                if old_name != new_name and new_name != util.sanitize_identifier(new_name):
                    raise util.CytoflowOpError('channels',
                                               "Channel name {} must be a "
                                               "valid Python identifier."
                                               .format(new_name))

        tube0_conditions = set(self.tubes[0].conditions)
        for tube in self.tubes:
            tube_conditions = set(tube.conditions)
            if len(tube0_conditions ^ tube_conditions) > 0:
                raise util.CytoflowOpError('tubes',
                                           "Tube {0} didn't have the same "
                                           "conditions as tube {1}"
                                           .format(tube.file, self.tubes[0].file))

        for idx, i in enumerate(self.tubes[0:-1]):
            for j in self.tubes[idx+1:]:
                if i.conditions_equal(j):
                    raise util.CytoflowOpError('tubes',
                                               "The same conditions specified for "
                                               "tube {0} and tube {1}"
                                               .format(i.file, j.file))

        experiment = Experiment()

        experiment.metadata["ignore_v"] = self.ignore_v

        for condition, dtype in list(self.conditions.items()):
            experiment.add_condition(condition, dtype)
            experiment.metadata[condition]['experiment'] = True

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                tube0_meta = parse(self.tubes[0].file,
                                   data_set = self.data_set,
                                   meta_data_only = True,
                                   reformat_meta = True)
        except Exception as e:
            raise util.CytoflowOpError('tubes',
                                       "FCS reader threw an error reading metadata "
                                       "for tube {}: {}"
                                       .format(self.tubes[0].file, str(e))) from e

        # ZE5 cytometer stores data at 8192x BD-scale; normalize
        is_ze5 = tube0_meta.get('$CYT') == 'YETI'
        experiment.metadata['is_ze5'] = is_ze5

        meta_channels = tube0_meta["_channels_"]

        if self.name_metadata:
            experiment.metadata["name_metadata"] = self.name_metadata
        else:
            experiment.metadata["name_metadata"] = autodetect_name_metadata(self.tubes[0].file,
                                                                             data_set = self.data_set)

        meta_channels['Index'] = meta_channels.index
        meta_channels.set_index(experiment.metadata["name_metadata"],
                                inplace = True)

        channels = list(self.channels.keys()) if self.channels \
                   else list(meta_channels.index.values)

        for channel in channels:
            if channel not in meta_channels.index:
                raise util.CytoflowOpError('channels',
                                           "Channel {0} not in tube {1}"
                                           .format(channel, self.tubes[0].file))

        for channel in channels:
            experiment.add_channel(channel)

            experiment.metadata[channel]["fcs_name"] = channel

            if("$PnV" in meta_channels.loc[channel]):
                v = meta_channels.loc[channel]['$PnV']
                if v: experiment.metadata[channel]["voltage"] = v

            data_range = meta_channels.loc[channel]['$PnR']
            data_range = float(data_range)
            experiment.metadata[channel]['range'] = data_range


        experiment.metadata['fcs_metadata'] = {}
        for tube in self.tubes:
            if metadata_only:
                tube_meta, tube_data = parse_tube(tube.file,
                                                   experiment,
                                                   data_set = self.data_set,
                                                   metadata_only = True)
            else:
                tube_meta, tube_data = parse_tube(tube.file,
                                                   experiment,
                                                   data_set = self.data_set)

                if self.events is not None:
                    if self.events <= len(tube_data):
                        tube_data = tube_data.loc[np.random.choice(tube_data.index,
                                                                   self.events,
                                                                   replace = False)]
                    else:
                        warnings.warn("Only {0} events in tube {1}"
                                      .format(len(tube_data), tube.file),
                                      util.CytoflowWarning)

                experiment.add_events(tube_data[channels], tube.conditions)

            if 'WELL ID' in tube_meta:
                pos = tube_meta['WELL ID']
                tube_meta['CF_Row'] = pos[0]
                tube_meta['CF_Col'] = int(pos[1:3])

            for i, channel in enumerate(channels):

                if '$P{}V'.format(i+1) in tube_meta:
                    del tube_meta['$P{}V'.format(i+1)]

                pnr = '$P{}R'.format(i+1)
                if pnr in tube_meta and float(tube_meta[pnr]) > experiment.metadata[channel]['range']:
                    experiment.metadata[channel]['range'] = float(tube_meta[pnr])


            tube_meta['CF_File'] = Path(tube.file).stem

            experiment.metadata['fcs_metadata'][tube.file] = tube_meta

        for channel in channels:
            if tube0_meta['$DATATYPE'] == 'I':
                data_bits  = int(meta_channels.loc[channel]['$PnB'])
                data_range = float(meta_channels.loc[channel]['$PnR'])
                range_bits = int(math.ceil(math.log(data_range, 2)))

                if range_bits < data_bits:
                    warnings.warn('The data range $PnR doesn\'t match the data bits $PnB for channel {}, masking out {} bits'
                                  .format(channel, data_bits - range_bits),
                                  util.CytoflowWarning)
                    mask = 1
                    for _ in range(1, range_bits):
                        mask = mask << 1 | 1

                    experiment.data[channel] = experiment.data[channel].values.astype('int') & mask

            data_range = float(meta_channels.loc[channel]['$PnR'])
            f1 = float(meta_channels.loc[channel]['$PnE'][0])
            f2 = float(meta_channels.loc[channel]['$PnE'][1])

            # store per-channel amplifier scale: log if $PnE indicates log amp
            if f1 > 0.0:
                experiment.metadata[channel]["amp_scale"] = "log"
            else:
                experiment.metadata[channel]["amp_scale"] = "linear"

            if f1 > 0.0 and f2 == 0.0:
                warnings.warn('Invalid $PnE = {},{} for channel {}, changing it to {},1.0'
                              .format(f1, f2, channel, f1),
                              util.CytoflowWarning)
                f2 = 1.0

            if f1 > 0.0 and f2 > 0.0 and tube0_meta['$DATATYPE'] == 'I':
                warnings.warn('Converting channel {} from logarithmic to linear'
                              .format(channel),
                              util.CytoflowWarning)
                experiment.data[channel] = 10 ** (f1 * experiment.data[channel] / data_range) * f2

        if is_ze5:
            experiment.data[channels] = experiment.data[channels].astype("float64") / 8192
            for channel in channels:
                experiment.metadata[channel]["range"] = float(experiment.metadata[channel]["range"]) / 8192

        for channel in channels:
            if self.channels and channel in self.channels:
                new_name = self.channels[channel]
                if channel == new_name:
                    continue
                experiment.data.rename(columns = {channel : new_name}, inplace = True)
                experiment.metadata[new_name] = experiment.metadata[channel]
                experiment.metadata[new_name]["fcs_name"] = channel
                del experiment.metadata[channel]

        return experiment


# ---------------------------------------------------------------------------
# RangeOp (from range.py)
# ---------------------------------------------------------------------------

@provides(IOperation)
class RangeOp(HasStrictTraits):
    """
    Apply a range gate to a cytometry experiment.

    Attributes
    ----------
    name : Str
        The operation name.  Used to name the new metadata field in the
        experiment that's created by `apply`

    channel : Str
        The name of the channel to apply the range gate.

    low : Float
        The lowest value to include in this gate.

    high : Float
        The highest value to include in this gate.

    Examples
    --------

    .. plot::
        :context: close-figs

        Make a little data set.

        >>> import cytoflow as flow
        >>> import_op = flow.ImportOp()
        >>> import_op.tubes = [flow.Tube(file = "Plate01/RFP_Well_A3.fcs",
        ...                              conditions = {'Dox' : 10.0}),
        ...                    flow.Tube(file = "Plate01/CFP_Well_A4.fcs",
        ...                              conditions = {'Dox' : 1.0})]
        >>> import_op.conditions = {'Dox' : 'float'}
        >>> ex = import_op.apply()

    Create and parameterize the operation.

    .. plot::
        :context: close-figs

        >>> range_op = flow.RangeOp(name = 'Range',
        ...                         channel = 'Y2-A',
        ...                         low = 2000,
        ...                         high = 10000)


    Plot a diagnostic view

    .. plot::
        :context: close-figs

        >>> rv = range_op.default_view(scale = 'log')
        >>> rv.plot(ex)


    .. note::
       If you want to use the interactive default view in a Jupyter notebook,
       make sure you say ``%matplotlib notebook`` in the first cell
       (instead of ``%matplotlib inline`` or similar).  Then call
       ``default_view()`` with ``interactive = True``::

            rv = range_op.default_view(scale = 'log',
                                       interactive = True)
            rv.plot(ex)

    Apply the gate, and show the result

    .. plot::
        :context: close-figs

        >>> ex2 = range_op.apply(ex)
        >>> ex2.data.groupby('Range').size()
        Range
        False    16042
        True      3958
        dtype: int64
    """

    # traits
    id = Constant('edu.mit.synbio.cytoflow.operations.range')
    friendly_id = Constant('Range')

    name = Str
    channel = Str
    low = Float(None)
    high = Float(None)

    _selection_view = Instance('RangeSelection', transient = True)

    def apply(self, experiment):
        """
        Applies the range gate to an experiment.

        Parameters
        ----------
        experiment : `Experiment`
            the `Experiment` to which this op is applied

        Returns
        -------
        Experiment
            a new    `Experiment`, the same as old `Experiment` but with a new
            column of type ``bool`` with the same as the operation name.  The
            bool is ``True`` if the event's measurement in `channel` is
            greater than `low` and less than `high`; it is ``False``
            otherwise.
        """

        if experiment is None:
            raise util.CytoflowOpError('experiment', "No experiment specified")

        if not self.name:
            raise util.CytoflowOpError('name',
                                       "You have to set the gate's name "
                                       "before applying it!")

        if self.name != util.sanitize_identifier(self.name):
            raise util.CytoflowOpError('name',
                                       "Name can only contain letters, numbers and underscores."
                                       .format(self.name))

        if self.name in experiment.data.columns:
            raise util.CytoflowOpError('name',
                                       "Experiment already has a column named {0}"
                                       .format(self.name))

        if not self.channel:
            raise util.CytoflowOpError('channel', "Channel not specified")

        if not self.channel in experiment.channels:
            raise util.CytoflowOpError('channel',
                                       "Channel {0} not in the experiment"
                                       .format(self.channel))

        if self.low is None:
            raise util.CytoflowOpError('low',
                                       "must set 'low'")

        if self.high is None:
            raise util.CytoflowOpError('high',
                                       "must set 'high'")

        if self.high <= self.low:
            raise util.CytoflowOpError('high',
                                       "range high must be > range low")

        if self.high <= experiment[self.channel].min():
            raise util.CytoflowOpError('high',
                                       "range high must be > {0}"
                                       .format(experiment[self.channel].min()))
        if self.low >= experiment[self.channel].max():
            raise util.CytoflowOpError('low',
                                       "range low must be < {0}"
                                       .format(experiment[self.channel].max()))

        gate = experiment[self.channel].between(self.low, self.high)
        new_experiment = experiment.clone(deep = False)
        new_experiment.add_condition(self.name, "bool", gate)
        new_experiment.history.append(self.clone_traits(transient = lambda _: True))

        return new_experiment

    def default_view(self, **kwargs):
        self._selection_view = RangeSelection(op = self)
        self._selection_view.trait_set(**kwargs)
        return self._selection_view


# ---------------------------------------------------------------------------
# PolygonOp (from polygon.py)
# ---------------------------------------------------------------------------

@provides(IOperation)
class PolygonOp(HasStrictTraits):
    """
    Apply a polygon gate to a cytometry experiment.

    Attributes
    ----------
    name : Str
        The operation name.  Used to name the new metadata field in the
        experiment that's created by `apply`

    xchannel, ychannel : Str
        The names of the x and y channels to apply the gate.

    xscale, yscale : {'linear', 'log'} (default = 'linear')
        The scales applied to the data before drawing the polygon.

    vertices : List((Float, Float))
        The polygon verticies.  An ordered list of 2-tuples, representing
        the x and y coordinates of the vertices.

    Notes
    -----
    You can set the verticies by hand, I suppose, but it's much easier to use
    the interactive view you get from `default_view` to do so.
    Set `ScatterplotPolygonSelectionView.interactive` to `True`, then
    single-click to set vertices. Click the first vertex a second time to
    close the polygon.  You'll need to do this in a Jupyter notebook with
    ``%matplotlib notebook`` -- see the ``Interactive Plots`` demo for an example.


    Examples
    --------

    .. plot::
        :context: close-figs

        Make a little data set.

        >>> import cytoflow as flow
        >>> import_op = flow.ImportOp()
        >>> import_op.tubes = [flow.Tube(file = "Plate01/RFP_Well_A3.fcs",
        ...                              conditions = {'Dox' : 10.0}),
        ...                    flow.Tube(file = "Plate01/CFP_Well_A4.fcs",
        ...                              conditions = {'Dox' : 1.0})]
        >>> import_op.conditions = {'Dox' : 'float'}
        >>> ex = import_op.apply()

    Create and parameterize the operation.

    .. plot::
        :context: close-figs

        >>> p = flow.PolygonOp(name = "Polygon",
        ...                    xchannel = "V2-A",
        ...                    ychannel = "Y2-A")
        >>> p.vertices = [(23.411982294776319, 5158.7027015021222),
        ...               (102.22182270573683, 23124.058843387455),
        ...               (510.94519955277201, 23124.058843387455),
        ...               (1089.5215641232173, 3800.3424832180476),
        ...               (340.56382570202402, 801.98947404942271),
        ...               (65.42597937575897, 1119.3133482602157)]


    Show the default view.

    .. plot::
        :context: close-figs

        >>> df = p.default_view(huefacet = "Dox",
        ...                     xscale = 'log',
        ...                     yscale = 'log',
        ...                     density = True)

        >>> df.plot(ex)


    .. note::
       If you want to use the interactive default view in a Jupyter notebook,
       make sure you say ``%matplotlib notebook`` in the first cell
       (instead of ``%matplotlib inline`` or similar).  Then call
       `default_view` with ``interactive = True``::

            df = p.default_view(huefacet = "Dox",
                                xscale = 'log',
                                yscale = 'log',
                                interactive = True)
            df.plot(ex)

    Apply the gate, and show the result

    .. plot::
        :context: close-figs

        >>> ex2 = p.apply(ex)
        >>> ex2.data.groupby('Polygon').size()
        Polygon
        False    15875
        True      4125
        dtype: int64

    You can also get (or draw) the polygon on a density plot instead of a
    scatterplot:

    .. plot::
        :context: close-figs

        >>> df = p.default_view(huefacet = "Dox",
        ...                     xscale = 'log',
        ...                     yscale = 'log')

        >>> df.plot(ex)
    """

    # traits
    id = Constant('edu.mit.synbio.cytoflow.operations.polygon')
    friendly_id = Constant("Polygon")

    name = Str
    xchannel = Str
    ychannel = Str
    vertices = List((Float, Float))

    xscale = util.ScaleEnum()
    yscale = util.ScaleEnum()

    _selection_view = Instance('_PolygonSelection', transient = True)

    def apply(self, experiment):
        """Applies the threshold to an experiment.

        Parameters
        ----------
        experiment : Experiment
            the old `Experiment` to which this op is applied

        Returns
        -------
        Experiment
            a new 'Experiment`, the same as ``experiment`` but with
            a new column of type `bool` with the same as the operation name.
            The bool is ``True`` if the event's measurement is within the
            polygon, and ``False`` otherwise.

        Raises
        ------
        CytoflowOpError
            if for some reason the operation can't be applied to this
            experiment. The reason is in the ``args`` attribute.
        """

        if experiment is None:
            raise util.CytoflowOpError('experiment',
                                       "No experiment specified")

        if not self.name:
            raise util.CytoflowOpError('name',
                                       "You have to set the Polygon gate's name "
                                       "before applying it!")

        if self.name in experiment.data.columns:
            raise util.CytoflowOpError('name',
                                       "{} is in the experiment already!"
                                       .format(self.name))

        if self.name != util.sanitize_identifier(self.name):
            raise util.CytoflowOpError('name',
                                       "Name can only contain letters, numbers and underscores."
                                       .format(self.name))

        if not self.xchannel:
            raise util.CytoflowOpError('xchannel',
                                       "Must specify an x channel")

        if not self.ychannel:
            raise util.CytoflowOpError('ychannel',
                                       "Must specify a y channel")

        if not self.xchannel in experiment.channels:
            raise util.CytoflowOpError('xchannel',
                                       "xchannel {0} is not in the experiment"
                                       .format(self.xchannel))

        if not self.ychannel in experiment.channels:
            raise util.CytoflowOpError('ychannel',
                                       "ychannel {0} is not in the experiment"
                                       .format(self.ychannel))

        if len(self.vertices) < 3:
            raise util.CytoflowOpError('vertices',
                                       "Must have at least 3 vertices")

        if any([len(x) != 2 for x in self.vertices]):
            return util.CytoflowOpError('vertices',
                                        "All vertices must be lists or tuples "
                                        "of length = 2")

        xscale = util.scale_factory(self.xscale, experiment, channel = self.xchannel)
        yscale = util.scale_factory(self.yscale, experiment, channel = self.ychannel)

        vertices = [(xscale(x), yscale(y)) for (x, y) in self.vertices]
        vertices.append(vertices[0])
        data = experiment.data[[self.xchannel, self.ychannel]].copy()
        data[self.xchannel] = xscale(data[self.xchannel])
        data[self.ychannel] = yscale(data[self.ychannel])

        xy_data = data[[self.xchannel, self.ychannel]].values
        in_polygon = util.polygon_contains(xy_data, np.array(vertices))

        new_experiment = experiment.clone(deep = False)
        new_experiment.add_condition(self.name, "bool", in_polygon)
        new_experiment.history.append(self.clone_traits(transient = lambda _: True))

        return new_experiment

    def default_view(self, **kwargs):
        """
        Returns an `IView` that allows a user to view the polygon or interactively draw it.

        Parameters
        ----------

        density : bool, default = False
            If `True`, return a density plot instead of a scatterplot.
        """

        density = kwargs.pop('density', False)
        if density:
            self._selection_view = DensityPolygonSelectionView(op = self)
        else:
            self._selection_view = ScatterplotPolygonSelectionView(op = self)

        self._selection_view.trait_set(**kwargs)
        return self._selection_view


# ---------------------------------------------------------------------------
# DensityGateOp (from density.py)
# ---------------------------------------------------------------------------

@provides(IOperation)
class DensityGateOp(HasStrictTraits):
    """
    This module computes a gate based on a 2D density plot.  The user chooses
    what proportion of events to keep, and the module creates a gate that selects
    that proportion of events in the highest-density bins of the 2D density
    histogram.

    Attributes
    ----------
    name : Str
        The operation name; determines the name of the new metadata column

    xchannel : Str
        The X channel to apply the binning to.

    ychannel : Str
        The Y channel to apply the binning to.

    xscale : {"linear", "log"} (default = "linear")
        Re-scale the data on the X acis before fitting the data?

    yscale : {"linear", "log"} (default = "linear")
        Re-scale the data on the Y axis before fitting the data?

    keep : Float (default = 0.9)
        What proportion of events to keep?  Must be ``>0`` and ``<1``

    bins : Int (default = 100)
        How many bins should there be on each axis?  Must be positive.

    min_quantile : Float (default = 0.001)
        Clip values below this quantile

    max_quantile : Float (default = 1.0)
        Clip values above this quantile

    sigma : Float (default = 1.0)
        What standard deviation to use for the gaussian blur?

    by : List(Str)
        A list of metadata attributes to aggregate the data before estimating
        the gate.  For example, if the experiment has two pieces of metadata,
        ``Time`` and ``Dox``, setting ``by = ["Time", "Dox"]`` will fit a
        separate gate to each subset of the data with a unique combination of
        ``Time`` and ``Dox``.

    Notes
    -----
    This gating method was developed by John Sexton, in Jeff Tabor's lab at
    Rice University.

    From http://taborlab.github.io/FlowCal/fundamentals/density_gate.html,
    the method is as follows:

    1. Determines the number of events to keep, based on the user specified
       gating fraction and the total number of events of the input sample.

    2. Divides the 2D channel space into a rectangular grid, and counts the
       number of events falling within each bin of the grid. The number of
       counts per bin across all bins comprises a 2D histogram, which is a
       coarse approximation of the underlying probability density function.

    3. Smoothes the histogram generated in Step 2 by applying a Gaussian Blur.
       Theoretically, the proper amount of smoothing results in a better
       estimate of the probability density function. Practically, smoothing
       eliminates isolated bins with high counts, most likely corresponding to
       noise, and smoothes the contour of the gated region.

    4. Selects the bins with the greatest number of events in the smoothed
       histogram, starting with the highest and proceeding downward until the
       desired number of events to keep, calculated in step 1, is achieved.

    Examples
    --------

    .. plot::
        :context: close-figs

        Make a little data set.

        >>> import cytoflow as flow
        >>> import_op = flow.ImportOp()
        >>> import_op.tubes = [flow.Tube(file = "Plate01/RFP_Well_A3.fcs",
        ...                              conditions = {'Dox' : 10.0}),
        ...                    flow.Tube(file = "Plate01/CFP_Well_A4.fcs",
        ...                              conditions = {'Dox' : 1.0})]
        >>> import_op.conditions = {'Dox' : 'float'}
        >>> ex = import_op.apply()

    Create and parameterize the operation.

    .. plot::
        :context: close-figs

        >>> dens_op = flow.DensityGateOp(name = 'Density',
        ...                              xchannel = 'FSC-A',
        ...                              xscale = 'log',
        ...                              ychannel = 'SSC-A',
        ...                              yscale = 'log',
        ...                              keep = 0.5)

    Find the bins to keep

    .. plot::
        :context: close-figs

        >>> dens_op.estimate(ex)

    Plot a diagnostic view

    .. plot::
        :context: close-figs

        >>> dens_op.default_view().plot(ex)

    Apply the gate

    .. plot::
        :context: close-figs

        >>> ex2 = dens_op.apply(ex)

    """

    id = Constant('edu.mit.synbio.cytoflow.operations.density')
    friendly_id = Constant("Density Gate")

    name = Str
    xchannel = Str
    ychannel = Str
    xscale = util.ScaleEnum
    yscale = util.ScaleEnum
    keep = util.PositiveFloat(0.9, allow_zero = False)
    bins = util.PositiveInt(100, allow_zero = False)
    min_quantile = util.PositiveFloat(0.001, allow_zero = True)
    max_quantile = util.PositiveFloat(1.0, allow_zero = False)
    sigma = util.PositiveFloat(1.0, allow_zero = False)
    by = List(Str)

    _xscale = Instance(util.IScale, transient = True)
    _yscale = Instance(util.IScale, transient = True)

    _xbins = Array(transient = True)
    _ybins = Array(transient = True)

    _keep_xbins = Dict(Any, Array, transient = True)
    _keep_ybins = Dict(Any, Array, transient = True)
    _histogram = Dict(Any, Array, transient = True)

    def estimate(self, experiment, subset = None):
        """
        Split the data set into bins and determine which ones to keep.

        Parameters
        ----------
        experiment : `Experiment`
            The `Experiment` to use to estimate the gate parameters.

        subset : Str (default = None)
            If set, determine the gate parameters on only a subset of the
            ``experiment`` parameter.
        """

        if experiment is None:
            raise util.CytoflowOpError('experiment',
                                       "No experiment specified")

        if self.xchannel not in experiment.data:
            raise util.CytoflowOpError('xchannel',
                                       "Column {0} not found in the experiment"
                                       .format(self.xchannel))

        if self.ychannel not in experiment.data:
            raise util.CytoflowOpError('ychannel',
                                       "Column {0} not found in the experiment"
                                       .format(self.ychannel))

        if self.min_quantile > 1.0:
            raise util.CytoflowOpError('min_quantile',
                                       "min_quantile must be <= 1.0")

        if self.max_quantile > 1.0:
            raise util.CytoflowOpError('max_quantile',
                                       "max_quantile must be <= 1.0")

        if not (self.max_quantile > self.min_quantile):
            raise util.CytoflowOpError('max_quantile',
                                       "max_quantile must be > min_quantile")

        if self.keep > 1.0:
            raise util.CytoflowOpError('keep',
                                       "keep must be <= 1.0")

        for b in self.by:
            if b not in experiment.conditions:
                raise util.CytoflowOpError('by',
                                           "Aggregation metadata {} not found, "
                                           "must be one of {}"
                                           .format(b, experiment.conditions))

        if subset:
            try:
                experiment = experiment.query(subset)
            except:
                raise util.CytoflowOpError('subset',
                                            "Subset string '{0}' isn't valid"
                                            .format(subset))

            if len(experiment) == 0:
                raise util.CytoflowOpError('subset',
                                           "Subset string '{0}' returned no events"
                                           .format(subset))

        if self.by:
            groupby = experiment.data.groupby(self.by)
        else:
            groupby = experiment.data.groupby(lambda _: True)

        self._xscale = xscale = util.scale_factory(self.xscale, experiment, channel = self.xchannel)
        self._yscale = yscale = util.scale_factory(self.yscale, experiment, channel = self.ychannel)


        xlim = (xscale.clip(experiment[self.xchannel].quantile(self.min_quantile)),
                xscale.clip(experiment[self.xchannel].quantile(self.max_quantile)))

        ylim = (yscale.clip(experiment[self.ychannel].quantile(self.min_quantile)),
                yscale.clip(experiment[self.ychannel].quantile(self.max_quantile)))

        self._xbins = xbins = xscale.inverse(np.linspace(xscale(xlim[0]),
                                                         xscale(xlim[1]),
                                                         self.bins))
        self._ybins = ybins = yscale.inverse(np.linspace(yscale(ylim[0]),
                                                         yscale(ylim[1]),
                                                         self.bins))

        histogram = {}
        for group, group_data in groupby:
            if len(group_data) == 0:
                raise util.CytoflowOpError('by',
                                           "Group {} had no data"
                                           .format(group))

            h, _, _ = np.histogram2d(group_data[self.xchannel],
                                     group_data[self.ychannel],
                                     bins=[xbins, ybins])

            h = scipy.ndimage.filters.gaussian_filter(h, sigma = self.sigma)

            i = scipy.stats.rankdata(h, method = "ordinal") - 1
            i = np.unravel_index(np.argsort(-i), h.shape)

            goal_count = self.keep * len(group_data)
            curr_count = 0
            num_bins = 0

            while(curr_count < goal_count and num_bins < i[0].size):
                curr_count += h[i[0][num_bins], i[1][num_bins]]
                num_bins += 1

            self._keep_xbins[group] = i[0][0:num_bins]
            self._keep_ybins[group] = i[1][0:num_bins]
            histogram[group] = h

        self._histogram = histogram


    def apply(self, experiment):
        """
        Creates a new condition based on membership in the gate that was
        parameterized with `estimate`.

        Parameters
        ----------
        experiment : `Experiment`
            the `Experiment` to apply the gate to.

        Returns
        -------
        `Experiment`
            a new `Experiment` with the new gate applied.
        """

        if experiment is None:
            raise util.CytoflowOpError('experiment',
                                       "No experiment specified")

        if not self.xchannel:
            raise util.CytoflowOpError('xchannel',
                                       "Must set X channel")

        if not self.ychannel:
            raise util.CytoflowOpError('ychannel',
                                       "Must set Y channel")

        if not self.name:
            raise util.CytoflowOpError('name',
                                       "You have to set the gate's name "
                                       "before applying it!")

        if self.name != util.sanitize_identifier(self.name):
            raise util.CytoflowOpError('name',
                                       "Name can only contain letters, numbers and underscores."
                                       .format(self.name))

        if self.name in experiment.data.columns:
            raise util.CytoflowOpError('name',
                                       "Experiment already has a column named {0}"
                                       .format(self.name))

        if not (self._xbins.size and self._ybins.size and self._keep_xbins):
            raise util.CytoflowOpError(None,
                                       "No gate estimate found.  Did you forget to "
                                       "call estimate()?")

        if not self._xscale:
            raise util.CytoflowOpError(None,
                                       "Couldn't find _xscale.  What happened??")

        if not self._yscale:
            raise util.CytoflowOpError(None,
                                       "Couldn't find _yscale.  What happened??")

        if self.xchannel not in experiment.data:
            raise util.CytoflowOpError('xchannel',
                                       "Column {0} not found in the experiment"
                                       .format(self.xchannel))

        if self.ychannel not in experiment.data:
            raise util.CytoflowOpError('ychannel',
                                       "Column {0} not found in the experiment"
                                       .format(self.ychannel))

        for b in self.by:
            if b not in experiment.conditions:
                raise util.CytoflowOpError('by',
                                           "Aggregation metadata {} not found, "
                                           "must be one of {}"
                                           .format(b, experiment.conditions))

        if self.by:
            groupby = experiment.data.groupby(self.by)
        else:
            groupby = experiment.data.groupby(lambda _: True)

        new_experiment = experiment.clone(deep=False)

        from ._utility import density_gate_to_polygon
        poly = density_gate_to_polygon(self)
        if poly and len(poly) >= 3:
            poly_gate = PolygonOp(name=self.name,
                                  xchannel=self.xchannel,
                                  ychannel=self.ychannel,
                                  vertices=poly)
            gated = poly_gate.apply(experiment)
            new_experiment.add_condition(self.name, "bool",
                                          gated.data[self.name])
        else:
            event_assignments = pd.Series([False] * len(experiment), dtype="bool")
            for group, group_data in groupby:
                if group not in self._keep_xbins:
                    continue
                group_idx = groupby.groups[group]
                cX = pd.cut(group_data[self.xchannel], self._xbins, include_lowest=True, labels=False).reset_index(drop=True)
                cY = pd.cut(group_data[self.ychannel], self._ybins, include_lowest=True, labels=False).reset_index(drop=True)
                group_keep = pd.Series([False] * len(group_data))
                keep_x = self._keep_xbins[group]
                keep_y = self._keep_ybins[group]
                for (xbin, ybin) in zip(keep_x, keep_y):
                    group_keep = group_keep | ((cX == xbin) & (cY == ybin))
                event_assignments.iloc[group_idx] = group_keep
            new_experiment.add_condition(self.name, "bool", event_assignments)

        new_experiment.history.append(self.clone_traits(transient = lambda _: True))
        return new_experiment

    def default_view(self, **kwargs):
        """
        Returns a diagnostic plot of the Gaussian mixture model.

        Returns
        -------
        `IView`
            a diagnostic view, call `DensityGateView.plot` to see the
            diagnostic plot.
        """
        v = DensityGateView(op = self)
        v.trait_set(**kwargs)
        return v


# ---------------------------------------------------------------------------
# AutofluorescenceOp (from autofluorescence.py)
# ---------------------------------------------------------------------------

@provides(IOperation)
class AutofluorescenceOp(HasStrictTraits):
    """
    Apply autofluorescence correction to a set of fluorescence channels.

    The `estimate` function loads a separate FCS file (not part of the input
    `Experiment`) and computes the untransformed median and standard deviation
    of the blank cells.  Then, `apply` subtracts the median from the
    experiment data.

    To use, set the `blank_file` property to point to an FCS file with
    unstained or nonfluorescing cells in it; set the `channels`
    property to a  list of channels to correct.

    `apply` also adds the ``af_median`` and ``af_stdev`` metadata to the
    corrected channels, representing the median and standard deviation of the
    measured blank distributions.

    Attributes
    ----------
    channels : List(Str)
        The channels to correct.

    blank_file : File
        The filename of a file with "blank" cells (not fluorescent).  Used
        to `estimate` the autofluorescence.

    blank_file_conditions : Dict
        Occasionally, you'll need to specify the experimental conditions that
        the blank tube was collected under (to apply the operations in the
        history.)  Specify them here.


    Examples
    --------
    Create a small experiment:

    .. plot::
        :context: close-figs

        >>> import cytoflow as flow
        >>> import_op = flow.ImportOp()
        >>> import_op.tubes = [flow.Tube(file = "tasbe/rby.fcs")]
        >>> ex = import_op.apply()

    Create and parameterize the operation

    .. plot::
        :context: close-figs

        >>> af_op = flow.AutofluorescenceOp()
        >>> af_op.channels = ["Pacific Blue-A", "FITC-A", "PE-Tx-Red-YG-A"]
        >>> af_op.blank_file = "tasbe/blank.fcs"

    Estimate the model parameters

    .. plot::
        :context: close-figs

        >>> af_op.estimate(ex)

    Plot the diagnostic plot

    .. plot::
        :context: close-figs

        >>> af_op.default_view().plot(ex)

    Apply the operation to the experiment

    .. plot::
        :context: close-figs

        >>> ex2 = af_op.apply(ex)

    """

    # traits
    id = Constant('edu.mit.synbio.cytoflow.operations.autofluorescence')
    friendly_id = Constant("Autofluorescence correction")

    name = Constant("Autofluorescence")
    channels = List(Str)
    blank_file = File(exists = True)
    blank_file_conditions = Dict({})

    _af_median = Dict(Str, Float, transient = True)
    _af_stdev = Dict(Str, Float, transient = True)
    _af_histogram = Dict(Str, Tuple(Array, Array), transient = True)


    def estimate(self, experiment, subset = None):
        """
        Estimate the autofluorescence from `blank_file` in channels
        specified in `channels`.

        Parameters
        ----------
        experiment : `Experiment`
            The experiment to which this operation is applied

        subset : Str (default = "")
            An expression that specifies the events used to compute the
            autofluorescence

        """
        if experiment is None:
            raise util.CytoflowOpError('experiment',
                                       "No experiment specified")

        if not self.blank_file:
            raise util.CytoflowOpError('blank_file',
                                       "No blank tube specified")

        if not self.channels:
            raise util.CytoflowOpError('channels', "No channels specified")

        if not set(self.channels) <= set(experiment.channels):
            raise util.CytoflowOpError('channels',
                                       "Specified channels that weren't found "
                                       "in the experiment.")

        self._af_median.clear()
        self._af_stdev.clear()
        self._af_histogram.clear()

        check_tube(self.blank_file, experiment)
        exp_conditions = {k: experiment.data[k].dtype.name for k in self.blank_file_conditions.keys()}

        blank_exp = ImportOp(tubes = [Tube(file = self.blank_file,
                                           conditions = self.blank_file_conditions)],
                             conditions = exp_conditions,
                             channels = {experiment.metadata[c]["fcs_name"] : c for c in experiment.channels},
                             name_metadata = experiment.metadata['name_metadata']).apply()

        for op in experiment.history:
            if hasattr(op, 'by'):
                for by in op.by:
                    if 'experiment' in experiment.metadata[by]:
                        warn("Found experimental metadata {} in experiment history; "
                             "make sure it's specified in blank_file_conditions."
                             .format(by),
                             util.CytoflowOpWarning)


            blank_exp = op.apply(blank_exp)

        if subset:
            try:
                blank_exp = blank_exp.query(subset)
            except Exception as exc:
                raise util.CytoflowOpError('subset',
                                           "Subset string '{0}' isn't valid"
                                           .format(subset)) from exc

            if len(blank_exp.data) == 0:
                raise util.CytoflowOpError('subset',
                                           "Subset string '{0}' returned no events"
                                      .format(subset))

        af_median = {}
        af_stdev = {}
        for channel in self.channels:
            self._af_histogram[channel] = np.histogram(blank_exp[channel], bins = 250)

            channel_min = blank_exp[channel].quantile(0.025)
            channel_max = blank_exp[channel].quantile(0.975)

            blank_exp[channel] = blank_exp[channel].clip(channel_min,
                                                         channel_max)

            af_median[channel] = np.median(blank_exp[channel])
            af_stdev[channel] = np.std(blank_exp[channel])

        self._af_stdev = af_stdev
        self._af_median = af_median

    def apply(self, experiment):
        """
        Applies the autofluorescence correction to channels in an experiment.

        Parameters
        ----------
        experiment : `Experiment`
            the experiment to which this op is applied

        Returns
        -------
        Experiment
            a new experiment with the autofluorescence median subtracted.  The
            corrected channels have the following metadata added to them:

            - **af_median** : Float
              The median of the non-fluorescent distribution

            - **af_stdev** : Float
              The standard deviation of the non-fluorescent distribution
        """

        if experiment is None:
            raise util.CytoflowOpError('experiment', "No experiment specified")

        if not self.channels:
            raise util.CytoflowOpError('channels', "No channels specified")

        if not self._af_median:
            raise util.CytoflowOpError(None, "Autofluorescence values aren't set. Did "
                                       "you forget to run estimate()?")

        if not set(self._af_median.keys()) <= set(experiment.channels) or \
           not set(self._af_stdev.keys()) <= set(experiment.channels):
            raise util.CytoflowOpError(None, "Autofluorescence estimates aren't set, or are "
                               "different than those in the experiment "
                               "parameter. Did you forget to run estimate()?")

        if not set(self._af_median.keys()) == set(self._af_stdev.keys()):
            raise util.CytoflowOpError(None, "Median and stdev keys are different! "
                                  "What the hell happened?!")

        if not set(self.channels) == set(self._af_median.keys()):
            raise util.CytoflowOpError('channels', "Estimated channels differ from the channels "
                               "parameter.  Did you forget to (re)run estimate()?")

        new_experiment = experiment.clone(deep = True)

        for channel in self.channels:
            new_experiment[channel] = \
                experiment[channel] - self._af_median[channel]

            new_experiment.metadata[channel]['af_median'] = self._af_median[channel]
            new_experiment.metadata[channel]['af_stdev'] = self._af_stdev[channel]

        new_experiment.history.append(self.clone_traits(transient = lambda t: True))

        return new_experiment

    def default_view(self, **kwargs):
        """
        Returns a diagnostic plot to see if the autofluorescence estimation
        is working.

        Returns
        -------
        IView
            An diagnostic view, call `AutofluorescenceDiagnosticView.plot`
            to see the diagnostic plots
        """
        v = AutofluorescenceDiagnosticView(op = self)
        v.trait_set(**kwargs)
        return v


# ---------------------------------------------------------------------------
# BleedthroughLinearOp (from bleedthrough_linear.py)
# ---------------------------------------------------------------------------

@provides(IOperation)
class BleedthroughLinearOp(HasStrictTraits):
    """
    Apply matrix-based bleedthrough correction to a set of fluorescence channels.

    This is a traditional matrix-based compensation for bleedthrough.  For each
    pair of channels, the user specifies the proportion of the first channel
    that bleeds through into the second; then, the module performs a matrix
    multiplication to compensate the raw data.

    The module can also estimate the bleedthrough matrix using one
    single-color control per channel.

    This works best on data that has had autofluorescence removed first;
    if that is the case, then the autofluorescence will be subtracted from
    the single-color controls too.

    To use, set up the `controls` dict with the single color controls;
    call `estimate` to parameterize the operation; check that the bleedthrough
    plots look good by calling `BleedthroughLinearDiagnostic.plot` on
    the `BleedthroughLinearDiagnostic` instance returned by
    `default_view`; and then `apply` on an `Experiment`.

    Attributes
    ----------
    controls : Dict(Str, File)
        The channel names to correct, and corresponding single-color control
        FCS files to estimate the correction splines with.  Must be set to
        use `estimate`.

    spillover : Dict(Tuple(Str, Str), Float)
        The spillover "matrix" to use to correct the data.  The keys are pairs
        of channels, and the values are proportions of spectral overlap.  If
        ``("channel1", "channel2")`` is present as a key,
        ``("channel2", "channel1")`` must also be present.  The module does not
        assume that the matrix is symmetric.

    control_conditions : Dict(Str, Dict(Str, Any))
        Occasionally, you'll need to specify the experimental conditions that
        the bleedthrough tubes were collected under (to apply the operations in the
        history.)  Specify them here.  The key is the channel name; they value
        is a dictionary of the conditions (same as you would specify for a
        `Tube` )

    Examples
    --------

    Create a small experiment:

    .. plot::
        :context: close-figs

        >>> import cytoflow as flow
        >>> import_op = flow.ImportOp()
        >>> import_op.tubes = [flow.Tube(file = "tasbe/rby.fcs")]
        >>> ex = import_op.apply()

    Correct for autofluorescence

    .. plot::
        :context: close-figs

        >>> af_op = flow.AutofluorescenceOp()
        >>> af_op.channels = ["Pacific Blue-A", "FITC-A", "PE-Tx-Red-YG-A"]
        >>> af_op.blank_file = "tasbe/blank.fcs"

        >>> af_op.estimate(ex)
        >>> af_op.default_view().plot(ex)

        >>> ex2 = af_op.apply(ex)

    Create and parameterize the operation

    .. plot::
        :context: close-figs

        >>> bl_op = flow.BleedthroughLinearOp()
        >>> bl_op.controls = {'Pacific Blue-A' : 'tasbe/ebfp.fcs',
        ...                   'FITC-A' : 'tasbe/eyfp.fcs',
        ...                   'PE-Tx-Red-YG-A' : 'tasbe/mkate.fcs'}

    Estimate the model parameters

    .. plot::
        :context: close-figs

        >>> bl_op.estimate(ex2)

    Plot the diagnostic plot

    .. note::
       The diagnostic plots look really bad in the online documentation.
       They're better in a real-world example, I promise!

    .. plot::
        :context: close-figs

        >>> bl_op.default_view().plot(ex2)

    Apply the operation to the experiment

    .. plot::
        :context: close-figs

        >>> ex2 = bl_op.apply(ex2)

    """

    # traits
    id = Constant('edu.mit.synbio.cytoflow.operations.bleedthrough_linear')
    friendly_id = Constant("Linear Bleedthrough Correction")

    name = Constant("Bleedthrough")

    controls = Dict(Str, File)
    spillover = Dict(Tuple(Str, Str), Float)
    control_conditions = Dict(Str, Dict(Str, Any), {})

    _sample = Dict(Str, Any, transient = True)

    def estimate(self, experiment, subset = None):
        """
        Estimate the bleedthrough from simgle-channel controls in `controls`
        """
        if experiment is None:
            raise util.CytoflowOpError('experiment', "No experiment specified")

        channels = list(self.controls.keys())

        if len(channels) < 2:
            raise util.CytoflowOpError('channels',
                                       "Need at least two channels to correct bleedthrough.")

        for channel in channels:
            if not os.path.isfile(self.controls[channel]):
                raise util.CytoflowOpError('channels',
                                           "Can't find file {0} for channel {1}."
                                           .format(self.controls[channel], channel))

        self.spillover = {}
        self._sample.clear()

        spillover = {}
        for channel in channels:

            check_tube(self.controls[channel], experiment)
            tube_conditions = self.control_conditions[channel] if channel in self.control_conditions else {}
            exp_conditions = {k: experiment.data[k].dtype.name for k in tube_conditions.keys()}

            tube_exp = ImportOp(tubes = [Tube(file = self.controls[channel],
                                              conditions = tube_conditions)],
                                conditions = exp_conditions,
                                channels = {experiment.metadata[c]["fcs_name"] : c for c in experiment.channels},
                                name_metadata = experiment.metadata['name_metadata']).apply()

            for op in experiment.history:
                if hasattr(op, 'by'):
                    for by in op.by:
                        if 'experiment' in experiment.metadata[by]:
                            raise util.CytoflowOpError('experiment',
                                                       "Prior to applying this operation, "
                                                       "you must not apply any operation with 'by' "
                                                       "set to an experimental condition.")
                tube_exp = op.apply(tube_exp)

            if subset:
                try:
                    tube_exp = tube_exp.query(subset)
                except Exception as exc:
                    raise util.CytoflowOpError('subset',
                                               "Subset string '{0}' isn't valid"
                                               .format(self.subset)) from exc

                if len(tube_exp.data) == 0:
                    raise util.CytoflowOpError('subset',
                                               "Subset string '{0}' returned no events"
                                               .format(self.subset))

            tube_data = tube_exp.data

            tube_data.sort_values(channel, inplace = True)

            self._sample[channel] = tube_data.sample(n = 1000)

            from_channel = channel

            from_min = np.min(tube_data[from_channel]) * 1.025
            from_max = np.max(tube_data[from_channel]) * 0.975
            tube_data[from_channel] = \
                tube_data[from_channel].clip(from_min, from_max)
            for to_channel in channels:

                if from_channel == to_channel:
                    continue

                to_min = np.min(tube_data[to_channel]) * 1.025
                to_max = np.max(tube_data[to_channel]) * 0.975
                tube_data[to_channel] = \
                    tube_data[to_channel].clip(to_min, to_max)

                tube_data.reset_index(drop = True, inplace = True)

                f = lambda x, k: x * k

                popt, _ = scipy.optimize.curve_fit(f,
                                                   tube_data[from_channel],
                                                   tube_data[to_channel],
                                                   0)

                spillover[(from_channel, to_channel)] = popt[0]

        self.spillover = spillover

    def apply(self, experiment):
        """
        Applies the bleedthrough correction to an experiment.

        Parameters
        ----------
        experiment : `Experiment`
            The experiment to which this operation is applied

        Returns
        -------
        Experiment
            A new `Experiment` with the bleedthrough subtracted out.
            The corrected channels have the following metadata added:

            - **linear_bleedthrough** : Dict(Str : Float)
              The values for spillover from other channels into this channel.

            - **bleedthrough_channels** : List(Str)
              The channels that were used to correct this one.

            - **bleedthrough_fn** : Callable (Tuple(Float) --> Float)
              The function that will correct one event in this channel.  Pass it
              the values specified in `controls` and it will return
              the corrected value for this channel.
        """
        if experiment is None:
            raise util.CytoflowOpError('experiment', "No experiment specified")

        if not self.spillover:
            raise util.CytoflowOpError('spillover',
                                       "Spillover matrix isn't set. "
                                       "Did you forget to run estimate()?")

        for (from_channel, to_channel) in self.spillover:
            if not from_channel in experiment.data:
                raise util.CytoflowOpError('spillover',
                                           "Can't find channel {0} in experiment"
                                           .format(from_channel))
            if not to_channel in experiment.data:
                raise util.CytoflowOpError('spillover',
                                           "Can't find channel {0} in experiment"
                                           .format(to_channel))

            if not (to_channel, from_channel) in self.spillover:
                raise util.CytoflowOpError('spillover',
                                           "Must have both (from, to) and "
                                           "(to, from) keys in self.spillover")

        new_experiment = experiment.clone(deep = True)

        channels = list(set([x for (x, _) in list(self.spillover.keys())]))

        a = [  [self.spillover[(y, x)] if x != y else 1.0 for x in channels]
               for y in channels]

        a_inv = np.linalg.pinv(a)

        new_channels = np.dot(experiment.data[channels], a_inv)

        for i, c in enumerate(channels):
            new_experiment[c] = pd.Series(new_channels[:, i])

        new_experiment.data = new_experiment.data[list(experiment.data.columns)]

        for channel in channels:
            new_experiment.metadata[channel]['linear_bleedthrough'] = \
                {x : self.spillover[(x, channel)]
                     for x in channels if x != channel}
            new_experiment.metadata[channel]['bleedthrough_channels'] = list(channels)
            new_experiment.metadata[channel]['bleedthrough_fn'] = lambda x, a_inv = a_inv: np.dot(x, a_inv)

        new_experiment.history.append(self.clone_traits(transient = lambda _: True))
        return new_experiment

    def default_view(self, **kwargs):
        """
        Returns a diagnostic plot to make sure spillover estimation is working.

        Returns
        -------
        `IView`
            An `IView`, call `BleedthroughLinearDiagnostic.plot` to see the diagnostic plots
        """

        channels = list(set([x for (x, _) in list(self.spillover.keys())]))

        if set(self.controls.keys()) != set(channels):
            raise util.CytoflowOpError('controls',
                                       "Must have both the controls and bleedthrough to plot")

        v = BleedthroughLinearDiagnostic(op = self)
        v.trait_set(**kwargs)
        return v


# ---------------------------------------------------------------------------
# Base operation views (from base_op_views.py)
# ---------------------------------------------------------------------------

@provides(IView)
class OpView(BaseDataView):
    """
    Attributes
    ----------
    op : Instance(`IOperation`)
        The `IOperation` that this view is associated with.  If you
        created the view using `default_view`, this is already set.
    """

    op = Instance(IOperation)

@provides(IView)
class Op1DView(OpView, Base1DView):
    """
    Attributes
    ----------
    channel : Str
        The channel this view is viewing.  If you created the view using
        `default_view`, this is already set.

    scale : {'linear', 'log'}
        The way to scale the x axes.  If you created the view using
        `default_view`, this may be already set.
    """

    channel = DelegatesTo('op')
    scale = DelegatesTo('op')

@provides(IView)
class Op2DView(OpView, Base2DView):
    """
    Attributes
    ----------
    xchannel : Str
        The channels to use for this view's X axis.  If you created the
        view using `default_view`, this is already set.

    ychannel : Str
        The channels to use for this view's Y axis.  If you created the
        view using `default_view`, this is already set.

    xscale : {'linear', 'log'}
        The way to scale the x axis.  If you created the view using
        `default_view`, this may be already set.

    yscale : {'linear', 'log'}
        The way to scale the y axis.  If you created the view using
        `default_view`, this may be already set.
    """
    xchannel = DelegatesTo('op')
    xscale = DelegatesTo('op')
    ychannel = DelegatesTo('op')
    yscale = DelegatesTo('op')

@provides(IView)
class ByView(OpView):
    """
    A view that can plot various plots based on the ``plot_name`` parameter
    of `plot`.

    Attributes
    ----------
    facets : List(Str)
        A read-only list of the conditions used to facet this view.

    by : List(Str)
        A read-only list of the conditions used to group this view's data before
        plotting.
    """

    facets = Property(List)
    by = Property(List)

    def _get_facets(self):
        return natsorted([x for x in [self.xfacet, self.yfacet, self.huefacet] if x])

    def _get_by(self):
        if self.op.by:
            return self.op.by
        else:
            return []

    def enum_plots(self, experiment):
        """
        Returns an iterator over the possible plots that this View can
        produce.  The values returned can be passed to the ``plot_name``
        keyword of `plot`.

        Parameters
        ----------
        experiment : `Experiment`
            The `Experiment` that will be producing the plots.
        """

        if experiment is None:
            raise util.CytoflowViewError('experiment',
                                         "No experiment specified")

        if len(self.by) == 0 and len(self.facets) > 1:
            raise util.CytoflowViewError('facets',
                                         "You can only facet this view if you "
                                         "specify some variables in `by`")

        for facet in self.facets:
            if facet not in experiment.conditions:
                raise util.CytoflowViewError('facets',
                                             "Facet {} not in the experiment"
                                            .format(facet))

        if len(self.facets) != len(set(self.facets)):
            raise util.CytoflowViewError('facets',
                                         "You can't reuse facets!")

        for b in self.by:
            if b not in experiment.conditions:
                raise util.CytoflowOpError('by',
                                           "Aggregation metadata {} not found, "
                                           "must be one of {}"
                                           .format(b, experiment.conditions))

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

        by = [x for x in self.by if x not in self.facets]

        class plot_enum(object):

            def __init__(self, by, experiment):
                self.by = by
                self._iter = None
                self._returned = False

                if by:
                    self._iter = experiment.data.groupby(by).__iter__()

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

        return plot_enum(by, experiment)

    def plot(self, experiment, **kwargs):
        """
        Make the plot.

        Parameters
        ----------
        plot_name : Str
            If this `IView` can make multiple plots, ``plot_name`` is
            the name of the plot to make.  Must be one of the values retrieved
            from `enum_plots`.
        """
        if experiment is None:
            raise util.CytoflowViewError('experiment',
                                         "No experiment specified")

        if len(self.by) == 0 and len(self.facets) > 1:
            raise util.CytoflowViewError('facets',
                                         "You can only facet this view if you "
                                         "specify some variables in `by`")

        for facet in self.facets:
            if facet not in experiment.conditions:
                raise util.CytoflowViewError('facets',
                                             "Facet {} not in the experiment"
                                             .format(facet))

            if facet not in self.by:
                raise util.CytoflowViewError('facets',
                                             "Facet {} must be one of {}"
                                             .format(facet, self.by))

        if len(self.facets) != len(set(self.facets)):
            raise util.CytoflowViewError('facets',
                                         "You can't reuse facets!")

        for b in self.by:
            if b not in experiment.conditions:
                raise util.CytoflowOpError('by',
                                           "Aggregation metadata {} not found, "
                                           "must be one of {}"
                                           .format(b, experiment.conditions))

        if self.subset:
            try:
                experiment = experiment.query(self.subset)
                experiment.data.reset_index(drop = True, inplace = True)
            except Exception as e:
                raise util.CytoflowViewError('subset',
                                             "Subset string '{0}' isn't valid"
                                             .format(self.subset)) from e

            if len(experiment) == 0:
                raise util.CytoflowViewError('subset',
                                             "Subset string '{0}' returned no events"
                                             .format(self.subset))

        by = [x for x in self.by if x not in self.facets]

        plot_name = kwargs.get('plot_name', None)

        if by and plot_name is None:
            raise util.CytoflowViewError('plot_name',
                                         "You must use facets {} in either the "
                                         "plot facets or the plot name. "
                                         "Possible plot names: {}"
                                         .format(by, [x for x in self.enum_plots(experiment)]))

        if plot_name is not None:
            if plot_name is not None and not by:
                raise util.CytoflowViewError('plot_name',
                                             "Don't set view.plot_name if you don't also set operation.by"
                                             .format(plot_name))

            groupby = experiment.data.groupby(by)

            if plot_name not in groupby.groups.keys():
                raise util.CytoflowViewError('plot_name',
                                             "Plot {} must be one of the values "
                                             "returned by enum_plots(). "
                                             "Possible plot names: {} "
                                             "(DEBUG: groupby keys: {}"
                                             .format(plot_name,
                                                     [x for x in self.enum_plots(experiment)],
                                                     groupby.groups.keys()))

            experiment = experiment.clone()
            experiment.data = groupby.get_group(plot_name)
            experiment.data.reset_index(drop = True, inplace = True)

        super().plot(experiment, **kwargs)

@provides(IView)
class By1DView(ByView, Op1DView):
    pass

@provides(IView)
class By2DView(ByView, Op2DView):
    pass

@provides(IView)
class NullView(BaseDataView):
    """
    An `IView` that doesn't actually do any plotting.
    """

    def _grid_plot(self, experiment, grid, **kwargs):
        return {}


@provides(IView)
class AnnotatingView(BaseDataView):
    """
    A `IView` that plots an underlying data plot, then plots some
    annotations on top of it.  See `gaussian.GaussianMixture1DView` for an
    example.  By default, it assumes that the annotations are to be plotted
    in the same color as the view's `huefacet`, and sets `huefacet`
    accordingly if the annotation isn't already set to a different facet.

    .. note::

        The ``annotation_facet`` and ``annotation_plot`` parameters that the
        `plot` method consumes are only for internal use, which is why
        they're not documented in the `plot` docstring.
    """

    def plot(self, experiment, **kwargs):
        """
        Parameters
        ----------
        color : matplotlib color
            The color to plot the annotations.  Overrides the default color
            cycle.
        """
        if experiment is None:
            raise util.CytoflowViewError('experiment',
                                         "No experiment specified")

        annotation_facet = kwargs.pop('annotation_facet', None)
        annotation_trait = kwargs.pop('annotation_trait', None)

        if annotation_facet is not None and annotation_facet in experiment.data:
            if annotation_trait:
                self.trait_set(**{annotation_trait : annotation_facet})
            elif not self.huefacet:
                warn("Setting 'huefacet' to '{}'".format(annotation_facet),
                     util.CytoflowViewWarning)
                annotation_trait = 'huefacet'
                self.trait_set(**{'huefacet' : annotation_facet})

        super().plot(experiment,
                     annotation_facet = annotation_facet,
                     **kwargs)

    def _grid_plot(self, experiment, grid, **kwargs):

        annotation_facet = kwargs.pop('annotation_facet', None)
        annotations = kwargs.pop('annotations', None)
        plot_name = kwargs.pop('plot_name', None)
        color = kwargs.get('color', None)

        plot_ret = super()._grid_plot(experiment, grid, **kwargs)
        kwargs.update(plot_ret)

        for (i, j, k), _ in grid.facet_data():
            ax = grid.facet_axis(i, j)

            row_name = grid.row_names[i] if grid.row_names and grid._row_var is not annotation_facet else None
            col_name = grid.col_names[j] if grid.col_names and grid._col_var is not annotation_facet else None
            hue_name = grid.hue_names[k] if grid.hue_names and grid._hue_var is not annotation_facet else None

            facets = [x for x in [row_name, col_name, hue_name] if x is not None]

            if plot_name is not None:
                if isinstance(plot_name, collections.abc.Iterable) and not isinstance(plot_name, str):
                    plot_name = list(plot_name)
                else:
                    plot_name = [plot_name]

                annotation_name = plot_name + facets
            else:
                annotation_name = facets

            annotation = None
            for group, a in annotations.items():
                if isinstance(group, collections.abc.Iterable) and not isinstance(group, str):
                    g_set = set(group)
                else:
                    g_set = set([group])

                if g_set == set(annotation_name):
                    annotation = a

            if (annotation is None
                and len(annotations.keys()) == 1
                and list(annotations.keys())[0] is True):
                annotation = annotations[True]

            if annotation is None:
                continue

            if annotation_facet is not None:
                if annotation_facet == grid._row_var:
                    annotation_value = grid.row_names[i]
                elif annotation_facet == grid._col_var:
                    annotation_value = grid.col_names[j]
                elif annotation_facet == grid._hue_var:
                    annotation_value = grid.hue_names[k]
                else:
                    annotation_value = None
            else:
                annotation_value = None

            annotation_color = grid._facet_color(k, color)

            self._annotation_plot(ax,
                                  annotation,
                                  annotation_facet,
                                  annotation_value,
                                  annotation_color,
                                  **kwargs)

        return plot_ret

    def _strip_trait(self, val):
        if val:
            trait_name = self._find_trait_name(val)
            if trait_name is not None:
                view = self.clone_traits('all')
                view.trait_set(**{trait_name : ""})
                return view, trait_name
        return self, None

    def _find_trait_name(self, val):
        traits = self.trait_get()
        for n, v in traits.items():
            if v == val:
                return n


# ---------------------------------------------------------------------------
# RangeSelection (from range.py)
# ---------------------------------------------------------------------------

@provides(ISelectionView)
class RangeSelection(Op1DView, HistogramView):
    """
    Plots, and lets the user interact with, a selection on the X axis.

    Attributes
    ----------

    interactive : Bool
        is this view interactive?  Ie, can the user set min and max
        with a mouse drag?

    Examples
    --------

    In an IPython notebook with ``%matplotlib notebook``

    >>> r = RangeOp(name = "RangeGate",
    ...             channel = 'Y2-A')
    >>> rv = r.default_view()
    >>> rv.interactive = True
    >>> rv.plot(ex2)
    >>> ### draw a range on the plot ###
    >>> print r.low, r.high
    """

    id = Constant('edu.mit.synbio.cytoflow.views.range')
    friendly_id = Constant("Range Selection")

    xfacet = Constant(None)
    yfacet = Constant(None)

    scale = util.ScaleEnum

    interactive = Bool(False, transient = True)

    # internal state.
    _ax = Any(transient = True)
    _span = Instance(SpanSelector, transient = True)
    _low_line = Instance(Line2D, transient = True)
    _high_line = Instance(Line2D, transient = True)
    _hline = Instance(Line2D, transient = True)
    _line_props = Dict()

    def plot(self, experiment, **kwargs):
        """
        Plot the underlying histogram and then plot the selection on top of it.

        Parameters
        ----------

        line_props : Dict
           The properties of the `matplotlib.lines.Line2D` that are drawn
           on top of the histogram.  They're passed directly to the
           `matplotlib.lines.Line2D` constructor.
           Default: ``{color : 'black', linewidth : 2}``
        """

        if experiment is None:
            raise util.CytoflowViewError('experiment',
                                         "No experiment specified")

        self._line_props = kwargs.pop('line_props',
                                        {'color' : 'black',
                                         'linewidth' : 2})

        super(RangeSelection, self).plot(experiment, **kwargs)
        self._ax = plt.gca()
        self._draw_span(None)
        self._interactive(None)

    @observe('[op.low,op.high]', post_init = True)
    def _draw_span(self, _):
        if not (self._ax and self.op.low and self.op.high):
            return

        if self._low_line and self._low_line in self._ax.lines:
            self._low_line.remove()

        if self._high_line and self._high_line in self._ax.lines:
            self._high_line.remove()

        if self._hline and self._hline in self._ax.lines:
            self._hline.remove()

        self._low_line = plt.axvline(self.op.low, **self._line_props)
        self._high_line = plt.axvline(self.op.high, **self._line_props)

        ymin, ymax = plt.ylim()
        y = (ymin + ymax) / 2.0
        self._hline = plt.plot([self.op.low, self.op.high],
                               [y, y], **self._line_props)[0]

        plt.draw()

    @observe('interactive', post_init = True)
    def _interactive(self, _):
        if self._ax and self.interactive:
            self._span = SpanSelector(self._ax,
                                      onselect = self._onselect,
                                      direction = "horizontal",
                                      interactive = False,
                                      props = dict(facecolor = 'none',
                                                   edgecolor = 'blue',
                                                   linewidth = 2),
                                      useblit = True)
        else:
            self._span = None


    def _onselect(self, xmin, xmax):
        """Update selection traits"""
        self.op.low = xmin
        self.op.high = xmax

util.expand_class_attributes(RangeSelection)
util.expand_method_parameters(RangeSelection, RangeSelection.plot)


# ---------------------------------------------------------------------------
# _PolygonSelection (from polygon.py)
# ---------------------------------------------------------------------------

class _PolygonSelection(Op2DView):
    xfacet = Constant(None)
    yfacet = Constant(None)

    interactive = Bool(False, transient = True)

    # internal state.
    _ax = Any(transient = True)
    _widget = Instance(PolygonSelector, transient = True)
    _patch = Instance(mpl.patches.PathPatch, transient = True)
    _patch_props = Dict()

    def plot(self, experiment, **kwargs):

        if experiment is None:
            raise util.CytoflowViewError('experiment',
                                         "No experiment specified")

        self._patch_props = kwargs.pop('patch_props',
                                        {'edgecolor' : 'black',
                                         'linewidth' : 2,
                                         'fill' : False})

        super(_PolygonSelection, self).plot(experiment, **kwargs)
        self._ax = plt.gca()
        self._draw_poly(None)
        self._interactive(None)


    @observe('op.vertices', post_init = True)
    def _draw_poly(self, _):
        if not self._ax:
            return

        if self._patch and self._patch in self._ax.patches:
            self._patch.remove()

        if not self.op.vertices or len(self.op.vertices) < 3:
            return

        patch_vert = np.concatenate((np.array(self.op.vertices),
                                    np.array((0,0), ndmin = 2)))

        self._patch = \
            mpl.patches.PathPatch(mpl.path.Path(patch_vert, closed = True), **self._patch_props)

        self._ax.add_patch(self._patch)
        plt.draw()

    @observe('interactive', post_init = True)
    def _interactive(self, _):
        if self._ax and self.interactive:
            self._widget = PolygonSelector(self._ax,
                                           self._onselect,
                                           useblit = True,
                                           grab_range = 20)
        elif self._widget:
            self._widget.set_active(False)
            self._widget = None

    def _onselect(self, vertices):
        self.op.vertices = vertices
        self.interactive = False


# ---------------------------------------------------------------------------
# ScatterplotPolygonSelectionView (from polygon.py)
# ---------------------------------------------------------------------------

@provides(ISelectionView)
class ScatterplotPolygonSelectionView(_PolygonSelection, ScatterplotView):
    """
    Plots, and lets the user interact with, a 2D polygon selection on a scatterplot.

    Attributes
    ----------
    interactive : bool
        is this view interactive?  Ie, can the user set the polygon verticies
        with mouse clicks?

    Examples
    --------

    In a Jupyter notebook with ``%matplotlib notebook``

    >>> s = flow.PolygonOp(xchannel = "V2-A",
    ...                    ychannel = "Y2-A")
    >>> poly = s.default_view()
    >>> poly.plot(ex2)
    >>> poly.interactive = True
    """

    id = Constant('edu.mit.synbio.cytoflow.views.polygon')
    friendly_id = Constant("Polygon Selection")

    def plot(self, experiment, **kwargs):
        """
        Plot the default view, and then draw the selection on top of it.

        Parameters
        ----------

        patch_props : Dict
           The properties of the `matplotlib.patches.Patch` that are drawn
           on top of the scatterplot or density view.  They're passed
           directly to the `matplotlib.patches.Patch` constructor.
           Default: ``{edgecolor : 'black', linewidth : 2, fill : False}``

        """
        super().plot(experiment, **kwargs)

util.expand_class_attributes(ScatterplotPolygonSelectionView)
util.expand_method_parameters(ScatterplotPolygonSelectionView, ScatterplotPolygonSelectionView.plot)


# ---------------------------------------------------------------------------
# DensityPolygonSelectionView (from polygon.py)
# ---------------------------------------------------------------------------

@provides(ISelectionView)
class DensityPolygonSelectionView(_PolygonSelection, DensityView):
    """
    Plots, and lets the user interact with, a 2D polygon selection on a density plot.

    Attributes
    ----------
    interactive : bool
        is this view interactive?  Ie, can the user set the polygon verticies
        with mouse clicks?

    Examples
    --------

    In a Jupyter notebook with ``%matplotlib notebook``

    >>> s = flow.PolygonOp(xchannel = "V2-A",
    ...                    ychannel = "Y2-A")
    >>> poly = s.default_view(density = True)
    >>> poly.plot(ex2)
    >>> poly.interactive = True
    """

    id = Constant('edu.mit.synbio.cytoflow.views.polygon_density')
    friendly_id = Constant("Polygon Selection")

    def plot(self, experiment, **kwargs):
        """
        Plot the default view, and then draw the selection on top of it.

        Parameters
        ----------

        patch_props : Dict
           The properties of the `matplotlib.patches.Patch` that are drawn
           on top of the scatterplot or density view.  They're passed
           directly to the `matplotlib.patches.Patch` constructor.
           Default: {edgecolor : 'black', linewidth : 2, fill : False}

        """
        super().plot(experiment, **kwargs)

util.expand_class_attributes(DensityPolygonSelectionView)
util.expand_method_parameters(ScatterplotPolygonSelectionView, ScatterplotPolygonSelectionView.plot)


# ---------------------------------------------------------------------------
# DensityGateView (from density.py)
# ---------------------------------------------------------------------------

@provides(IView)
class DensityGateView(By2DView, AnnotatingView, DensityView):
    """
    A diagnostic view for `DensityGateOp`.  Draws a density plot,
    then outlines the selected bins in white.

    Attributes
    ----------

    """

    id = Constant('edu.mit.synbio.cytoflow.view.densitygateview')
    friendly_id = Constant("Density Gate Diagnostic Plot")

    huefacet = Constant(None)

    def plot(self, experiment, **kwargs):
        """
        Plot the plots.

        Parameters
        ----------

        contour_props : Dict
            The keyword arguments passed to the
            `matplotlib.axes.Axes.contour` constructor, which controls the
            visual properties of the contour that's plotted on top of the
            density plot. Default: ``{'colors' : 'w'}``

        """

        contour_props = kwargs.pop('contour_props',
                                   {'colors' : 'w'})

        annotations = {}
        for i in self.op._keep_xbins:
            annotations[i] = (self.op._keep_xbins[i],
                              self.op._keep_ybins[i],
                              self.op._histogram[i],
                              contour_props)

        super().plot(experiment,
                     annotations = annotations,
                     xscale = self.op._xscale,
                     yscale = self.op._yscale,
                     **kwargs)

    def _annotation_plot(self,
                         axes,
                         annotation,
                         annotation_facet,
                         annotation_value,
                         annotation_color,
                         **kwargs):


        keep_x = annotation[0]
        keep_y = annotation[1]
        h = annotation[2]
        contour_props = annotation[3]
        xbins = self.op._xbins[0:-1]
        ybins = self.op._ybins[0:-1]
        last_level = h[keep_x[-1], keep_y[-1]]

        axes.contour(xbins, ybins, h.T, [last_level], **contour_props)

util.expand_class_attributes(DensityGateView)
util.expand_method_parameters(DensityGateView, DensityGateView.plot)


# ---------------------------------------------------------------------------
# AutofluorescenceDiagnosticView (from autofluorescence.py)
# ---------------------------------------------------------------------------

@provides(IView)
class AutofluorescenceDiagnosticView(HasStrictTraits):
    """
    Plots a histogram of each channel, and its median in red.  Serves as a
    diagnostic for the autofluorescence correction.

    Attributes
    ----------
    op : Instance(`AutofluorescenceOp`)
        The `AutofluorescenceOp` whose parameters we're viewing. Set
        automatically if you created the instance using
        `AutofluorescenceOp.default_view`.

    """

    # traits
    id = Constant('edu.mit.synbio.cytoflow.view.autofluorescencediagnosticview')
    friendly_id = Constant("Autofluorescence Diagnostic")

    op = Instance(AutofluorescenceOp)

    def plot(self, experiment, **kwargs):
        """
        Plot a faceted histogram view of a channel
        """

        if experiment is None:
            raise util.CytoflowViewError('experiment', "No experiment specified")

        if not self.op.channels:
            raise util.CytoflowViewError('op', "No channels specified")

        if not self.op._af_median:
            raise util.CytoflowViewError('op',
                                         "Autofluorescence values aren't set. Did "
                                         "you forget to run estimate()?")

        if not set(self.op._af_median.keys()) <= set(experiment.channels) or \
           not set(self.op._af_stdev.keys()) <= set(experiment.channels) or \
           not set(self.op._af_histogram.keys()) <= set(experiment.channels):
            raise util.CytoflowViewError('op',
                                       "Autofluorescence estimates aren't set, or are "
                                       "different than those in the experiment "
                                       "parameter. Did you forget to run estimate()?")

        if not set(self.op._af_median.keys()) == set(self.op._af_stdev.keys()):
            raise util.CytoflowOpError('op',
                                       "Median and stdev keys are different! "
                                       "What the heck happened?!")

        if not set(self.op.channels) == set(self.op._af_median.keys()):
            raise util.CytoflowOpError('channels', "Estimated channels differ from the channels "
                               "parameter.  Did you forget to (re)run estimate()?")

        plt.figure()

        for idx, channel in enumerate(self.op.channels):
            hist, bin_edges = self.op._af_histogram[channel]
            hist = hist[1:-1]
            bin_edges = bin_edges[1:-1]
            plt.subplot(len(self.op.channels), 1, idx+1)
            plt.title(channel)
            plt.bar(bin_edges[:-1], hist, width = bin_edges[2] - bin_edges[1], linewidth = 0)
            plt.axvline(self.op._af_median[channel], color = 'r')

        plt.tight_layout(pad = 0.8)


# ---------------------------------------------------------------------------
# BleedthroughLinearDiagnostic (from bleedthrough_linear.py)
# ---------------------------------------------------------------------------

@provides(IView)
class BleedthroughLinearDiagnostic(HasStrictTraits):
    """
    Plots a scatterplot of each channel vs every other channel and the
    bleedthrough line

    Attributes
    ----------
    op : Instance(BleedthroughPiecewiseOp)
        The operation whose parameters we're viewing.  If you made the instance
        with `BleedthroughLinearOp.default_view`, this is set for you
        already.

    subset : str
        If set, only plot this subset of the underlying data.

    """

    # traits
    id = Constant("edu.mit.synbio.cytoflow.view.linearbleedthroughdiagnostic")
    friendly_id = Constant("Linear Bleedthrough Diagnostic")

    subset = Str

    op = Instance(IOperation)

    def plot(self, experiment = None, **kwargs):
        """
        Plot a diagnostic of the bleedthrough model computation.
        """

        if experiment is None:
            raise util.CytoflowViewError('experiment',
                                         "No experiment specified")

        if not self.op.controls:
            raise util.CytoflowViewError('op',
                                         "No controls specified")

        if not self.op.spillover:
            raise util.CytoflowViewError('op',
                                         "No spillover matrix specified")

        kwargs.setdefault('histtype', 'stepfilled')
        kwargs.setdefault('alpha', 0.5)
        kwargs.setdefault('antialiased', True)

        channels = natsorted(list(set([x for (x, _) in list(self.op.spillover.keys())])))
        num_channels = len(channels)
        _, axes2d = plt.subplots(nrows=num_channels, ncols=num_channels)

        for to_idx, row in enumerate(axes2d):
            for from_idx, ax in enumerate(row):

                plt.sca(ax)

                from_channel = channels[from_idx]
                to_channel = channels[to_idx]

                if to_idx == len(axes2d) - 1 or (to_idx == len(axes2d) - 2 and from_idx == len(row) - 1):
                    plt.xlabel(from_channel)
                if from_idx == 0 or (from_idx == 1 and to_idx == 0):
                    plt.ylabel(to_channel)

                if from_idx == to_idx:
                    ax.set_visible(False)
                    continue

                tube_data = self.op._sample[from_channel]

                scale_name = 'log'

                xscale = util.scale_factory(scale_name, experiment, channel = from_channel)
                yscale = util.scale_factory(scale_name, experiment, channel = to_channel)

                plt.xscale(scale_name, **xscale.get_mpl_params(ax.get_xaxis()))
                plt.yscale(scale_name, **yscale.get_mpl_params(ax.get_yaxis()))

                plt.scatter(tube_data[from_channel],
                            tube_data[to_channel],
                            alpha = 1,
                            s = 1,
                            marker = 'o')

                xs = np.logspace(-1, math.log(tube_data[from_channel].max(), 10))
                ys = xs * self.op.spillover[(from_channel, to_channel)]

                plt.plot(xs, ys, 'g-', lw=3)


        plt.tight_layout(pad = 0.8)
