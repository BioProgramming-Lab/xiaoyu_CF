# _utility.py - merged utility module for xiaoyu_CF
# Vendored from cytoflow.utility.*

from __future__ import print_function, unicode_literals, division

import random
import string
import numbers
import re
import logging
import warnings
from warnings import warn
import inspect
import struct

import numpy as np
import pandas as pd
from scipy import stats
import scipy.optimize

import matplotlib.colors
import matplotlib.scale
from matplotlib import transforms
from matplotlib.ticker import NullFormatter, LogFormatterMathtext
from matplotlib.ticker import Locator

from traits.api import (Interface, Str, Instance, Tuple, Array,
                        BaseInt, BaseCInt, BaseFloat, BaseCFloat, BaseEnum,
                        TraitType, HasTraits, Float, Property, Undefined,
                        provides, Constant, Enum, Dict)
from traits.has_traits import HasStrictTraits

from numba import jit, njit
import numba

# ===========================================================================
# cytoflow_errors
# ===========================================================================

import warnings
formatwarning_orig = warnings.formatwarning
warnings.formatwarning = lambda message, category, filename, lineno, line=None: \
    formatwarning_orig(message, category, filename, lineno, line='')

class CytoflowError(RuntimeError):
    """
    A general error
    """

class CytoflowOpError(CytoflowError):
    """
    An error raised by an operation.

    Parameters
    ----------
    args[0] : string
        The attribute or parameter whose bad value caused the error, or ``None``
        if there isn't one.

    args[1] : string
        A more verbose error message.
    """

class CytoflowViewError(CytoflowError):
    """
    An error raised by a view.

    Parameters
    ----------
    args[0] : string
        The attribute or parameter whose bad value caused the error, or ``None``
        if there isn't one.

    args[1] : string
        A more verbose error message.
    """

class CytoflowWarning(UserWarning):
    """
    A general warning.
    """

class CytoflowOpWarning(CytoflowWarning):
    """
    A warning raised by an operation.

    Parameters
    ----------
    args[0] : string
        A verbose warning message
    """

class CytoflowViewWarning(CytoflowWarning):
    """
    A warning raised by a view.

    Parameters
    ----------
    args[0] : string
        A verbose warning message
    """

warnings.simplefilter('always', CytoflowWarning)
warnings.simplefilter('always', CytoflowOpWarning)
warnings.simplefilter('always', CytoflowViewWarning)

# ===========================================================================
# util_functions
# ===========================================================================

def iqr(a):
    """
    Calculate the inter-quartile range for an array of numbers.

    Parameters
    ----------
    a : array_like
        The array of numbers to compute the IQR for.

    Returns
    -------
    float
        The IQR of the data.
    """
    a = np.asarray(a)
    q1 = np.nanpercentile(a, 25)
    q3 = np.nanpercentile(a, 75)
    return q3 - q1

def num_hist_bins(a):
    """
    Calculate number of histogram bins using Freedman-Diaconis rule.

    From http://stats.stackexchange.com/questions/798/

    Parameters
    ----------
    a : array_like
        The data to make a histogram of.

    Returns
    -------
    int
        The number of bins in the histogram
    """
    a = np.asarray(a)
    h = 2 * iqr(a) / (len(a) ** (1 / 3))

    # fall back to 10 bins if iqr is 0
    if h == 0:
        return 10.
    else:
        return np.ceil((np.nanpercentile(a, 99) -
                        np.nanpercentile(a, 1)) / h)

def geom_mean(a):
    """
    Compute the geometric mean for an "arbitrary" data set, ie one that
    contains zeros and negative numbers.

    Parameters
    ----------

    a : array-like
        A numpy.ndarray, or something that can be converted to an ndarray

    Returns
    -------
    The geometric mean of the input array

    Notes
    -----
    The traditional geometric mean can not be computed on a mixture of positive
    and negative numbers.  The approach here, validated rigorously in the
    cited paper[1], is to compute the geometric mean of the absolute value of
    the negative numbers separately, and then take a weighted arithmetic mean
    of that and the geometric mean of the positive numbers.  We're going to
    discard 0 values, operating under the assumption that in this context
    there are going to be few or no observations with a value of exactly 0.

    References
    ----------
    [1] Geometric mean for negative and zero values
        Elsayed A. E. Habib
        International Journal of Research and Reviews in Applied Sciences
        11:419 (2012)
        http://www.arpapress.com/Volumes/Vol11Issue3/IJRRAS_11_3_08.pdf

    """

    a = np.array(a)
    pos = a[a > 0]
    pos_mean = stats.gmean(pos)
    pos_prop = pos.size / a.size

    neg = a[a < 0]
    neg = np.abs(neg)
    neg_mean = stats.gmean(neg) if neg.size > 0 else 0
    neg_prop = neg.size / a.size

    return (pos_mean * pos_prop) - (neg_mean * neg_prop)

def geom_sd(a):
    """
    Compute the geometric standard deviation for an "abitrary" data set, ie one
    that contains zeros and negative numbers.  Since we're in log space, this
    gives a *dimensionless scaling factor*, not a measure.  If you want
    traditional "error bars", don't plot ``[geom_mean - geom_sd, geom_mean + geom_sd]``;
    rather, plot ``[geom_mean / geom_sd, geom_mean * geom_sd]``.

    Parameters
    ----------

    a : array-like
        A numpy.ndarray, or something that can be converted to an ndarray

    Returns
    -------
    The geometric standard deviation of the distribution.

    Notes
    -----
    As with `geom_mean`, non-positive numbers pose a problem.  The approach
    here, though less rigorously validated than the one above, is to replace
    negative numbers with their absolute value plus 2 * geometric mean, then
    go about our business as per the Wikipedia page for geometric sd[1].

    References
    ----------
    [1] https://en.wikipedia.org/wiki/Geometric_standard_deviation
    """

    a = np.array(a)
    u = geom_mean(a)
    a[a <= 0] = np.abs(a[a <= 0]) + 2 * u

    return np.exp(np.std(np.log(a)))

def geom_sd_range(a):
    """
    A convenience function to compute [geom_mean / geom_sd, geom_mean * geom_sd].

    Parameters
    ----------

    a : array-like
        A numpy.ndarray, or something that can be converted to an ndarray

    Returns
    -------
    A tuple, with ``(geom_mean / geom_sd, geom_mean * geom_sd)``
    """

    u = geom_mean(a)
    sd = geom_sd(a)

    return (u / sd, u * sd)

def geom_sem(a):
    """
    Compute the geometric standard error of the mean for an "arbirary" data set,
    ie one that contains zeros and negative numbers.

    Parameters
    ----------

    a : array-like
        A numpy.ndarray, or something that can be converted to an ndarray

    Returns
    -------
    The geometric mean of the distribution.

    Notes
    -----
    As with `geom_mean`, non-positive numbers pose a problem.  The approach
    here, though less rigorously validated than the one above, is to replace
    negative numbers with their absolute value plus 2 * geometric mean.  The
    geometric SEM is computed as in [1].

    References
    ----------
    [1] The Standard Errors of the Geometric and Harmonic Means and Their Application to Index Numbers
        Nilan Norris
        The Annals of Mathematical Statistics
        Vol. 11, No. 4 (Dec., 1940), pp. 445-448

        http://www.jstor.org/stable/2235723?seq=1#page_scan_tab_contents
    """

    a = np.array(a)
    u = geom_mean(a)
    a[a <= 0] = np.abs(a[a <= 0]) + 2 * u

    return u * np.std(np.log(a)) / np.sqrt(a.size)


def geom_sem_range(a):
    """
    A convenience function to compute [geom_mean / geom_sem, geom_mean * geom_sem].

    Parameters
    ----------

    a : array-like
        A numpy.ndarray, or something that can be converted to an ndarray

    Returns
    -------
    A tuple, with ``(geom_mean / geom_sem, geom_mean * geom_sem)``
    """

    u = geom_mean(a)
    sem = geom_sem(a)

    return (u / sem, u * sem)


def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    References
    ----------
    Originally from http://stackoverflow.com/a/1235363/4755587
    """

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n // arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in range(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

def sanitize_identifier(name):
    """
    Makes name a Python identifier by replacing all nonsafe characters with '_'
    """

    new_name = list(name)
    for i, c in enumerate(list(name)): # character by character
        if i == 0 and not (c.isalpha() or c == '_'):
            new_name[i] = '_'
        if i > 0 and not (c.isalnum() or c == '_'):
            new_name[i] = '_'

    return  "".join(new_name)


def random_string(n):
    """
    Makes a random string of ascii digits and lowercase letters of length ``n``

    from http://stackoverflow.com/questions/2257441/random-string-generation-with-upper-case-letters-and-digits-in-python
    """
    return ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(n))

def is_numeric(s):
    """
    Determine if a `pandas.Series` or `numpy.ndarray` is numeric
    from its dtype.
    """
    return s.dtype.kind in 'iufc'

def cov2corr(covariance):
    '''
    Compute the correlation matrix from the covariance matrix.

    From https://github.com/AndreaCensi/procgraph/blob/master/src/procgraph_statistics/cov2corr.py
    '''

    sigma = np.sqrt(covariance.diagonal())
    M = np.multiply.outer(sigma, sigma)
    correlation = covariance / M

    return sigma, correlation

# ===========================================================================
# algorithms
# ===========================================================================

def ci(data, func, which=95, boots=1000):
    """
    Determine the confidence interval of a function applied to a data set by
    bootstrapping.

    Parameters
    ----------
    data : pandas.DataFrame
        The data to resample.

    func : callable
        A function that is called on a resampled ``data``

    which : int
        The percentile to use for the confidence interval

    boots : int (default = 1000):
        How many times to bootstrap

    Returns
    -------
    (float, float)
        The confidence interval.

    """
    boots = bootstrap(data, func = func, n_boot = boots)
    p = 50 - which / 2, 50 + which / 2
    return tuple(percentiles(boots, p))

def percentiles(a, pcts, axis=None):
    """
    Like `scipy.stats.scoreatpercentile` but can take and return array of percentiles.

    from seaborn: https://github.com/mwaskom/seaborn/blob/master/seaborn/utils.py

    Parameters
    ----------
    a : array
        data

    pcts : sequence of percentile values
        percentile or percentiles to find score at

    axis : int or None
        if not None, computes scores over this axis

    Returns
    -------
    scores: array
        array of scores at requested percentiles
        first dimension is length of object passed to ``pcts``

    """
    scores = []
    try:
        n = len(pcts)
    except TypeError:
        pcts = [pcts]
        n = 0
    for p in pcts:
        if axis is None:
            score = stats.scoreatpercentile(a.ravel(), p)
        else:
            score = np.apply_along_axis(stats.scoreatpercentile, axis, a, p)
        scores.append(score)
    scores = np.asarray(scores)
    if not n:
        scores = scores.squeeze()
    return scores


def bootstrap(*args, **kwargs):
    """
    Resample one or more arrays with replacement and store aggregate values.
    Positional arguments are a sequence of arrays to bootstrap along the first
    axis and pass to a summary function.


    Parameters
    ----------
    n_boot : int, default 10000
        Number of iterations

    axis : int, default None
        Will pass axis to ``func`` as a keyword argument.

    units : array, default None
        Array of sampling unit IDs. When used the bootstrap resamples units
        and then observations within units instead of individual
        datapoints.

    smooth : bool, default False
        If True, performs a smoothed bootstrap (draws samples from a kernel
        destiny estimate); only works for one-dimensional inputs and cannot
        be used `units` is present.

    func : callable, default np.mean
        Function to call on the args that are passed in.

    random_seed : int | None, default None
        Seed for the random number generator; useful if you want
        reproducible resamples.

    Returns
    -------
    array
        array of bootstrapped statistic values

    from seaborn: https://github.com/mwaskom/seaborn/blob/master/seaborn/algorithms.py
    """
    # Ensure list of arrays are same length
    if len(np.unique(list(map(len, args)))) > 1:
        raise ValueError("All input arrays must have the same length")
    n = len(args[0])

    # Default keyword arguments
    n_boot = kwargs.get("n_boot", 10000)
    func = kwargs.get("func", np.mean)
    axis = kwargs.get("axis", None)
    units = kwargs.get("units", None)
    smooth = kwargs.get("smooth", False)
    random_seed = kwargs.get("random_seed", None)
    if axis is None:
        func_kwargs = dict()
    else:
        func_kwargs = dict(axis=axis)

    # Initialize the resampler
    rs = np.random.RandomState(random_seed)

    # Coerce to arrays
    args = list(map(np.asarray, args))
    if units is not None:
        units = np.asarray(units)

    # Do the bootstrap
    if smooth:
        return _smooth_bootstrap(args, n_boot, func, func_kwargs)

    if units is not None:
        return _structured_bootstrap(args, n_boot, units, func,
                                     func_kwargs, rs)

    boot_dist = []
    for _ in range(int(n_boot)):
        resampler = rs.randint(0, n, n)
        sample = [a.take(resampler, axis=0) for a in args]
        boot_dist.append(func(*sample, **func_kwargs))
    return np.array(boot_dist)


def _structured_bootstrap(args, n_boot, units, func, func_kwargs, rs):
    """Resample units instead of datapoints."""
    unique_units = np.unique(units)
    n_units = len(unique_units)

    args = [[a[units == unit] for unit in unique_units] for a in args]

    boot_dist = []
    for _ in range(int(n_boot)):
        resampler = rs.randint(0, n_units, n_units)
        sample = [np.take(a, resampler, axis=0) for a in args]
        lengths = list(map(len, sample[0]))
        resampler = [rs.randint(0, n, n) for n in lengths]
        sample = [[c.take(r, axis=0) for c, r in zip(a, resampler)]
                  for a in sample]
        sample = list(map(np.concatenate, sample))
        boot_dist.append(func(*sample, **func_kwargs))
    return np.array(boot_dist)


def _smooth_bootstrap(args, n_boot, func, func_kwargs):
    """Bootstrap by resampling from a kernel density estimate."""
    n = len(args[0])
    boot_dist = []
    kde = [stats.gaussian_kde(np.transpose(a)) for a in args]
    for _ in range(int(n_boot)):
        sample = [a.resample(n).T for a in kde]
        boot_dist.append(func(*sample, **func_kwargs))
    return np.array(boot_dist)


@jit(nopython=True)
def is_inside_sm(polygon, point):
    length = len(polygon)-1
    dy2 = point[1] - polygon[0][1]
    intersections = 0
    ii = 0
    jj = 1

    while ii<length:
        dy  = dy2
        dy2 = point[1] - polygon[jj][1]

        # consider only lines which are not completely above/bellow/right from the point
        if dy*dy2 <= 0.0 and (point[0] >= polygon[ii][0] or point[0] >= polygon[jj][0]):

            # non-horizontal line
            if dy<0 or dy2<0:
                F = dy*(polygon[jj][0] - polygon[ii][0])/(dy-dy2) + polygon[ii][0]

                if point[0] > F: # if line is left from the point - the ray moving towards left, will intersect it
                    intersections += 1
                elif point[0] == F: # point on line
                    return 2

            # point on upper peak (dy2=dx2=0) or horizontal line (dy=dy2=0 and dx*dx2<=0)
            elif dy2==0 and (point[0]==polygon[jj][0] or (dy==0 and (point[0]-polygon[ii][0])*(point[0]-polygon[jj][0])<=0)):
                return 2

        ii = jj
        jj += 1

    #print 'intersections =', intersections
    return intersections & 1


@njit(parallel=True)
def polygon_contains(points, polygon):
    ln = len(points)
    D = np.empty(ln, dtype=numba.boolean)
    for i in numba.prange(ln):
        D[i] = is_inside_sm(polygon,points[i])
    return D

# ===========================================================================
# logging
# ===========================================================================

class MplFilter(logging.Filter):
    """A `logging.Filter` that removes nuisance warnings"""
    def filter(self, record):
        if record.msg == "posx and posy should be finite values":
            return 0
        else:
            return 1

# ===========================================================================
# scale infrastructure (IScale, ScaleMixin, registry)
# ===========================================================================

class IScale(Interface):
    """
    An interface for various ways we could rescale flow data.

    Attributes
    ----------
    name : Str
        The name of this view (for serialization, UI, etc.)

    experiment : Instance("xiaoyu_CF._experiment.Experiment")
        The experiment this scale is to be applied to.  Needed because some
        scales have parameters estimated from data.

    channel : Str
        Which channel to scale.  Needed because some scales have parameters
        estimated from data.

    condition : Str
        What condition to scale.  Needed because some scales have parameters
        estimated from the a condition.  Must be a numeric condition; else
        instantiating the scale should fail.

    statistic : Tuple(Str, Str)
        What statistic to scale.  Needed because some scales have parameters
        estimated from a statistic.  The statistic must be numeric or an
        iterable of numerics; else instantiating the scale should fail.

    data : array_like
        What raw data to scale.
    """

    id = Str
    name = Str

    experiment = Instance("xiaoyu_CF._experiment.Experiment")

    # what are we using to parameterize the scale?  set one of these; if
    # multiple are set, the first is used.
    channel = Str
    condition = Str
    statistic = Tuple(Str, Str)
    error_statistic = Tuple(Str, Str)
    data = Array

    def __call__(self, data):
        """
        Transforms `data` using this scale.  Must know how to handle int, float,
        and lists, tuples, numpy.ndarrays and pandas.Series of int or float.
        Must return the same type passed.

        Careful!  May return `NaN` if the scale domain doesn't match the data
        (ie, applying a log10 scale to negative numbers.
        """

    def inverse(self, data):
        """
        Transforms 'data' using the inverse of this scale.  Must know how to
        handle int, float, and list, tuple, numpy.ndarray and pandas.Series of
        int or float.  Returns the same type as passed.
        """

    def clip(self, data):
        """
        Clips the data to the scale's domain.
        """

    def norm(self, vmin = None, vmax = None):
        """
        Return an instance of matplotlib.colors.Normalize, which normalizes
        this scale to [0, 1].  Among other things, this is used to correctly
        scale a color bar.
        """


class ScaleMixin(HasStrictTraits):
    """
    Provides useful functionality that scales can inherit.
    """

    def __init__(self, **kwargs):

        # run the traits constructor
        super().__init__(**kwargs)

        # check that the channel, condition, or statistic is either numeric or
        # an iterable of numerics
        if self.condition:
            if not is_numeric(self.experiment[self.condition]):
                raise CytoflowError("Tried to scale the non-numeric condition {}"
                                    .format(self.condition))

        elif self.statistic[0]:
            stat = self.experiment.statistics[self.statistic]
            if is_numeric(stat):
                return
            else:
                try:
                    for x in stat:
                        for y in x:
                            if not isinstance(y, numbers.Number):
                                raise CytoflowError("Tried to scale a non-numeric "
                                                    "statistic {}"
                                                    .format(self.statistic))
                except TypeError as e:
                    raise CytoflowError("Error scaling statistic {}"
                                        .format(self.statistic)) from e

# maps name -> scale object
_scale_mapping = {}
_scale_default = "linear"

def scale_factory(scale, experiment, **scale_params):
    """
    Make a new instance of a named scale.

    Parameters
    ----------
    scale : string
        The name of the scale to build

    experiment : Experiment
        The experiment to use to parameterize the new scale.

    **scale_params : kwargs
        Other parameters to pass to the scale constructor
    """

    scale = scale.lower()

    if scale not in _scale_mapping:
        raise CytoflowError("Unknown scale type {0}".format(scale))

    return _scale_mapping[scale](experiment = experiment, **scale_params)

def register_scale(scale_class):
    """
    Register a new scale for the `scale_factory` and `ScaleEnum`.
    """

    _scale_mapping[scale_class.name] = scale_class

def set_default_scale(scale):
    """
    Set a default scale for `ScaleEnum`
    """

    global _scale_default

    scale = scale.lower()
    if scale not in _scale_mapping:
        raise CytoflowError("Unknown scale type {0}".format(scale))

    _scale_default = scale

def get_default_scale():
    """Get the defaults scale set with `set_default_scale`"""
    return _scale_default


# ===========================================================================
# custom_traits
# ===========================================================================

class PositiveInt(BaseInt):
    """
    Defines a trait whose value must be a positive integer
    """

    info_text = 'a positive integer'

    def validate(self, obj, name, value):
        if self.allow_none and value == None:
            return None

        value = super().validate(obj, name, value)
        if (value > 0 or (self.allow_zero and value >= 0)):
            return value

        self.error(obj, name, value)

class PositiveCInt(BaseCInt):
    """
    Defines a trait whose converted value must be a positive integer
    """

    info_text = 'a positive integer'

    def validate(self, obj, name, value):
        if self.allow_none and (value == "" or value == None):
            return None

        value = super().validate(obj, name, value)
        if (value > 0 or (self.allow_zero and value >= 0)):
            return value

        self.error(obj, name, value)


class PositiveFloat(BaseFloat):
    """
    Defines a trait whose value must be a positive float
    """

    info_text = 'a positive float'

    def validate(self, obj, name, value):
        if self.allow_none and value == None:
            return None

        value = super().validate(obj, name, value)
        if (value > 0.0 or (self.allow_zero and value >= 0.0)):
            return value

        self.error(obj, name, value)

class PositiveCFloat(BaseCFloat):
    """
    Defines a trait whose converted value must be a positive float
    """

    info_text = 'a positive float'

    def validate(self, obj, name, value):
        if self.allow_none and (value == "" or value == None):
            return None

        value = super().validate(obj, name, value)
        if (value > 0.0 or (self.allow_zero and value >= 0.0)):
            return value

        self.error(obj, name, value)

class FloatOrNone(BaseFloat):
    """
    Defines a trait whose value must be a float or None
    """

    info_text = 'a float or None'

    def validate(self, obj, name, value):
        if value == "" or value == None:
            return None
        else:
            return super().validate(obj, name, value)

class CFloatOrNone(BaseCFloat):
    """
    Defines a trait whose converted value must be a float or None
    """

    info_text = 'a float or None'

    def validate(self, obj, name, value):
        if value == None or value == "":
            return None
        else:
            return super().validate(obj, name, value)

class IntOrNone(BaseInt):
    """
    Defines a trait whose value must be an integer or None
    """

    info_text = 'an int or None'

    def validate(self, obj, name, value):
        if value == None:
            return None
        else:
            return super().validate(obj, name, value)

class CIntOrNone(BaseCInt):
    """
    Defines a trait whose converted value must be an integer or None
    """

    info_text = 'an int or None'

    def validate(self, obj, name, value):
        if value == None or value == "":
            return None
        else:
            return super().validate(obj, name, value)


class ScaleEnum(BaseEnum):
    """
    Defines an enumeration that contains one of the registered data scales
    """

    info_text = 'an enum containing one of the registered scales'

    def __init__ ( self, *args, **metadata ):
        """ Returns an Enum trait with values from the registered scales
        """
        self.name = None
        self.values = list(_scale_mapping.keys())
        super(BaseEnum, self).__init__(_scale_default, **metadata )

    def get_default_value(self):
        # this is so silly.  get_default_value is ... called once?  as traits
        # sets up?  idunno.  anyways, instead of returning _scale_default, we
        # need to return a reference to a function that returns _scale_default.

        return (7, (self._get_default_value, (), None))

    def _get_default_value(self):
        return _scale_default

class Removed(TraitType):
    """
    Names a trait that was present in a previous version but was removed.

    Trait metadata:

        - **err_string** : the error string in the error

        - **warning** : if ``True``, raise a warning when the trait is referenced.
            Otherwise, raise an exception.

    """

    def __init__(self, **metadata):
        metadata.setdefault('err_string', 'Trait {} has been removed')
        metadata.setdefault('transient', True)
        super().__init__(**metadata)

    def get(self, obj, name):
        # TODO - this is quite slow.  come up with a better way.
        curframe = inspect.currentframe()
        calframe = inspect.getouterframes(curframe)
        if calframe[1][3] == "copy_traits" or calframe[1][3] == "trait_get":
            return

        if self.warning:
            warn(self.err_string.format(name), CytoflowWarning)
        else:
            raise CytoflowError(self.err_string.format(name))

    def set(self, obj, name, value):
        curframe = inspect.currentframe()
        calframe = inspect.getouterframes(curframe, 2)
        if calframe[1][3] == "copy_traits":
            return

        if self.warning:
            warn(self.err_string.format(name), CytoflowWarning)
        else:
            raise CytoflowError(self.err_string.format(name))

class Deprecated(TraitType):
    """
    Names a trait that was present in a previous version but was renamed.
    When the trait is accessed, a warning is raised, and the access
    is passed through to the new trait.

    Trait metadata:
        - **new** : the name of the new trait

        - **err_string** : the error string in the error

    """

    def __init__(self, **metadata):
        metadata.setdefault('err_string', 'Trait {} is deprecated; please use {}')
        metadata.setdefault('transient', True)
        super().__init__(**metadata)

    def get(self, obj, name):
        curframe = inspect.currentframe()
        calframe = inspect.getouterframes(curframe)
        if calframe[1][3] != "copy_traits" and calframe[1][3] != 'trait_get':
            warn(self.err_string.format(name, self.new), CytoflowWarning)

        return getattr(obj, self.new)

    def set(self, obj, name, value):
        curframe = inspect.currentframe()
        calframe = inspect.getouterframes(curframe)
        if calframe[1][3] != "copy_traits":
            warn(self.err_string.format(name, self.new), CytoflowWarning)
        setattr(obj, self.new, value)


# ===========================================================================
# docstring
# ===========================================================================

def expand_class_attributes(cls):
    """
    Takes entries in the ``Attributes`` section of a class's ancestors'
    docstrings and adds them to the ``Attributes`` section of this class's
    docstring.

    All the classes must have docstrings formatted with the using the ``numpy``
    docstring style.

    Parameters
    ----------
    cls : class
        The class whose docstring is to be expanded.
    """

    if cls.__doc__ is None:
        return

    lines = cls.__doc__.split("\n")

    first_attr_line, last_attr_line = find_section("Attributes", lines)
    if first_attr_line is None:
        warn("Couldn't find 'Attributes' section for {}".format(cls.__name__))
        return

    cls_attr = []

    for mro_cls in cls.__mro__:
        mro_attr = get_class_attributes(mro_cls)

        for mro_attr_name, mro_attr_type, mro_attr_body in mro_attr:
            do_add = True
            for cls_attr_name, _, _ in cls_attr:
                if cls_attr_name == mro_attr_name:
                    do_add = False
            if do_add:
                cls_attr.append((mro_attr_name, mro_attr_type, mro_attr_body))

    new_cls_attr = []
    for i in range(len(cls_attr)):
        for j in range(i+1, len(cls_attr)):
            if cls_attr[i] is not None and cls_attr[j] is not None and \
               cls_attr[i][1] == cls_attr[j][1] and \
               cls_attr[i][2] == cls_attr[j][2]:
                name = cls_attr[i][0] + ', ' + cls_attr[j][0]
                cls_attr[i] = (name, cls_attr[i][1], cls_attr[i][2])
                cls_attr[j] = None
        if cls_attr[i] is not None:
            new_cls_attr.append(cls_attr[i])

    attr_section = []
    for attr_name, attr_value, attr_body in new_cls_attr:
        attr_section.append("    " + attr_name + " : " + attr_value)
        attr_section.extend(attr_body)
        attr_section.append("")

    lines[first_attr_line:last_attr_line + 1] = attr_section

    cls.__doc__ = "\n".join(lines)

def expand_method_parameters(cls, method):
    """
    Expand the ``Parameters`` section of a method's docstring with
    ``Parameters`` from the overridden methods in the class's ancestors.

    All the methods must have docstrings formatted with the using the ``numpy``
    docstring style.

    Parameters
    ----------
    cls : class
        The class whose ancestors are to be parsed for more parameters.

    method : callable
        The method whose docstring is to be expanded.
    """
    if method.__doc__ is None:
        return

    lines = method.__doc__.split("\n")

    first_param_line, last_param_line = find_section("Parameters", lines)
    if first_param_line is None:
        warn("Couldn't find a 'Parameters' section for {} : {}".format(cls.__name__, method.__name__))
        return

    method_params = []

    for mro_cls in reversed(cls.__mro__):
        try:
            mro_method = getattr(mro_cls, method.__name__)
        except:
            continue

        mro_params = get_method_parameters(mro_method)

        for mro_param_name, mro_param_type, mro_param_body in mro_params:
            do_add = True
            for method_param_name, _, _ in method_params:
                if method_param_name == mro_param_name:
                    do_add = False
            if do_add:
                method_params.append((mro_param_name, mro_param_type, mro_param_body))

    new_method_params = []
    for i in range(len(method_params)):
        for j in range(i+1, len(method_params)):
            if method_params[i] is not None and method_params[j] is not None and \
               method_params[i][1] == method_params[j][1] and \
               method_params[i][2] == method_params[j][2]:
                name = method_params[i][0] + ', ' + method_params[j][0]
                method_params[i] = (name, method_params[i][1], method_params[i][2])
                method_params[j] = None
        if method_params[i] is not None:
            new_method_params.append(method_params[i])

    params_section = []
    for param_name, param_type, param_body in new_method_params:
        params_section.append("        " + param_name + " : " + param_type)
        params_section.extend(param_body)
        params_section.append("")

    lines[first_param_line:last_param_line + 1] = params_section

    method.__doc__ = "\n".join(lines)


def find_section(section, lines):
    """
    Find a named section in a ``numpy``-formatted docstring.

    Parameters
    ----------
    section : string
        The name of the section to find

    lines : array of string
        The docstring, split into lines

    Returns
    -------
    int, int
        The indices of the first and last lines of the section.
    """
    first_line = None
    for i, line in enumerate(lines):
        if "----" in line and i > 0 and section in lines[i - 1]:
            first_line = i + 1
            break

    if first_line is None:
        return None, None

    last_line = None

    for i in range(first_line, len(lines)):
        if "---" in lines[i]:
            last_line = i - 2
            break

    if last_line is None:
        last_line = len(lines) - 1

    return first_line, last_line


def get_class_attributes(cls):
    """
    Gets the entries from the ``Attributes`` section of a class's
    ``numpy``-formated docstring.

    Parameters
    ----------
    cls : class
        The class whose docstring to parse

    Returns
    ------
    array of ``(name, type, body)``

        - name : the attribute's name

        - type : the attribute's type

        - body : the attribute's description body
    """

    if not cls.__doc__:
        return []

    lines = cls.__doc__.split('\n')
    first_attr_line, last_attr_line = find_section("Attributes", lines)

    attributes = []

    # consume the attributes
    for  attr_name, attr_value, attr_body in _get_params_and_attrs(lines, first_attr_line, last_attr_line):
        attributes.append((attr_name, attr_value, attr_body))

    return attributes

def get_method_parameters(method):
    """
    Gets the entries from the ``Parameters`` section of a method's
    ``numpy``-formatted docstring.

    Parameters
    ----------
    method : callable
        The method whose docstring to parse.

    Returns
    ------
    array of ``(name, type, body)``

        - name : the attribute's name

        - type : the attribute's type

        - body : the attribute's description body

    """
    if not method.__doc__:
        return []

    lines = method.__doc__.split('\n')
    first_param_line, last_param_line = find_section("Parameters", lines)

    params = []

    # consume the attributes
    for param_name, param_value, param_body in _get_params_and_attrs(lines, first_param_line, last_param_line):
        params.append((param_name, param_value, param_body))

    return params


def _get_params_and_attrs(lines, first_attr_line, last_attr_line):
    if first_attr_line is None:
        return

    attr_names = None
    attr_type = None
    attr_body = []

    for i in range(first_attr_line, last_attr_line):
        if attr_names is None:
            colon_idx = lines[i].find(':')
            if colon_idx == -1:
                continue
            attr_names = lines[i][:colon_idx]
            attr_names = attr_names.split(',')
            attr_names = [x.strip() for x in attr_names]
            attr_type = lines[i][colon_idx + 1:].strip()

        if re.match(r'^\s*$', lines[i]):
            for attr_name in attr_names:
                yield attr_name, attr_type, attr_body[1:]
            attr_names = None
            attr_body = []
        else:
            attr_body.append(lines[i])

    if attr_names is not None:
        for attr_name in attr_names:
            yield attr_name, attr_type, attr_body[1:]


# ===========================================================================
# fcswrite
# ===========================================================================

def write_fcs(filename, chn_names, chn_ranges, data,
              compat_chn_names=True,
              compat_percent=True,
              compat_negative=True,
              compat_copy=True,
              verbose=0,
              **kws):
    """
    Write numpy data to an .fcs file (FCS3.0 file format)


    Parameters
    ----------
    filename : str
        Path to the output .fcs file

    chn_names : list of str, length C
        Names of the output channels

    chn_ranges : dictionary
        Keys: channel names.  Values: ranges

    data : 2d ndarray of shape (N,C)
        The numpy array data to store as .fcs file format.

    compat_chn_names : bool
        Compatibility mode for 3rd party flow analysis software:
        The characters " ", "?", and "_" are removed in the output
        channel names.

    compat_percent : bool
        Compatibliity mode for 3rd party flow analysis software:
        If a column in ``data`` contains values only between 0 and 1,
        they are multiplied by 100.

    compat_negative : bool
        Compatibliity mode for 3rd party flow analysis software:
        Flip the sign of ``data`` if its mean is smaller than zero.

    compat_copy : bool
        Do not override the input array ``data`` when modified in
        compatibility mode.

    kwargs : Str
        Additional keyword arguments are written as keyword/value pairs in
        the TEXT segment of the FCS file.

    Notes
    -----
    These commonly used unicode characters are replaced: "µ", "²"

    """
    if not isinstance(data, np.ndarray):
        data = np.array(data)

    msg="length of `chn_names` must match length of 2nd axis of `data`"
    assert len(chn_names) == data.shape[1], msg

    rpl = [["µ", "u"],
           ["²", "2"],
          ]

    if compat_chn_names:
        # Compatibility mode: Clean up headers.
        rpl += [[" ", ""],
                ["?", ""],
                ["_", ""],
                ]

    for i in range(len(chn_names)):
        for (a, b) in rpl:
            chn_names[i] = chn_names[i].replace(a, b)

    if compat_percent:
        # Compatibility mode: Scale values b/w 0 and 1 to percent
        toscale = []
        for ch in range(data.shape[1]):
            if data[:,ch].min() > 0 and data[:,ch].max() < 1:
                toscale.append(ch)
        if len(toscale):
            if compat_copy:
                # copy if requested
                data = data.copy()
            for ch in toscale:
                data[:,ch] *= 100

    if compat_negative:
        toflip = []
        for ch in range(data.shape[1]):
            if np.mean(data[:,ch]) < 0:
                toflip.append(ch)
        if len(toflip):
            if compat_copy:
                # copy if requested
                data = data.copy()
            for ch in toflip:
                data[:,ch] *= -1


    # DATA segment
    data1 = data.flatten().tolist()
    DATA = struct.pack('>%sf' % len(data1), *data1)

    # TEXT segment
    # fix length of TEXT to 4 kilo bytes
    ltxt = 4096
    ver='FCS3.0'
    textfirst= '{0: >8}'.format(256)
    datafirst= '{0: >8}'.format(256+ltxt)
    datalast = '{0: >8}'.format(256+ltxt+len(DATA)-1)
    anafirst = '{0: >8}'.format(0)
    analast  = '{0: >8}'.format(0)
    # use little endian
    #byteord = '1,2,3,4'
    # use big endian
    byteord = '4,3,2,1'
    TEXT ='/$BEGINANALYSIS/0/$ENDANALYSIS/0'
    TEXT+='/$BEGINSTEXT/0/$ENDSTEXT/0'
    TEXT+='/$BEGINDATA/{0}/$ENDDATA/{1}'.format(256+ltxt, 256+ltxt+len(DATA)-1)
    TEXT+='/$BYTEORD/{0}/$DATATYPE/F'.format(byteord)
    TEXT+='/$MODE/L/$NEXTDATA/0/$TOT/{0}'.format(data.shape[0])
    TEXT+='/$PAR/{0}'.format(data.shape[1])

    for i in range(data.shape[1]):
        pnrange = chn_ranges[chn_names[i]]
        # TODO:
        # - Set log/lin
        TEXT+='/$P{0}B/32/$P{0}E/0,0/$P{0}N/{1}/$P{0}R/{2}/$P{0}D/Linear'.format(i+1, chn_names[i], pnrange)

    for kw, val in kws.items():
        kw = kw.replace('/', '//')
        val = val.replace('/', '//')
        TEXT+='/{0}/{1}'.format(kw, val)

    TEXT += '/'

    if len(TEXT) > ltxt:
        raise RuntimeError("TEXT segment is too long; specify fewer keywords")

    textlast = '{0: >8}'.format(len(TEXT)+256-1)
    TEXT = TEXT.ljust(ltxt, ' ')

    # HEADER segment
    HEADER = '{0: <256}'.format(ver+'    '+
                                textfirst +
                                textlast  +
                                datafirst +
                                datalast  +
                                anafirst  +
                                analast)

    # Write data
    with open(filename, "wb") as fd:
        fd.write(HEADER.encode("ascii"))
        fd.write(TEXT.encode("ascii"))
        fd.write(DATA)
        fd.write(b'00000000')


# ===========================================================================
# linear_scale
# ===========================================================================

@provides(IScale)
class LinearScale(ScaleMixin):
    """
    A scale that doesn't transform the data at all.
    """

    id = Constant("edu.mit.synbio.cytoflow.utility.linear_scale")
    name = "linear"

    experiment = Instance("xiaoyu_CF._experiment.Experiment")

    # none of these are actually used
    channel = Str
    condition = Str
    statistic = Tuple(Str, Str)
    error_statistic = Tuple(Str, Str)
    data = Array

    def __call__(self, data):
        return data

    def inverse(self, data):
        return data

    def clip(self, data):
        return data

    def norm(self, vmin = None, vmax = None):
        if vmin is not None and vmax is not None:
            pass
        elif self.channel:
            vmin = self.experiment[self.channel].min()
            vmax = self.experiment[self.channel].max()
        elif self.condition:
            vmin = self.experiment[self.condition].min()
            vmax = self.experiment[self.condition].max()
        elif self.statistic in self.experiment.statistics:
            stat = self.experiment.statistics[self.statistic]
            try:
                vmin = min([min(x) for x in stat])
                vmax = max([max(x) for x in stat])
            except (TypeError, IndexError):
                vmin = stat.min()
                vmax = stat.max()

            if self.error_statistic in self.experiment.statistics:
                err_stat = self.experiment.statistics[self.error_statistic]
                try:
                    vmin = min([min(x) for x in err_stat])
                    vmax = max([max(x) for x in err_stat])
                except (TypeError, IndexError):
                    vmin = vmin - err_stat.min()
                    vmax = vmax + err_stat.max()
        elif self.data.size > 0:
            vmin = self.data.min()
            vmax = self.data.max()
        else:
            raise CytoflowError("Must set one of 'channel', 'condition' "
                                "or 'statistic'.")

        return matplotlib.colors.Normalize(vmin = vmin, vmax = vmax)

    def get_mpl_params(self, ax):
        return dict()


# ===========================================================================
# log_scale
# ===========================================================================

@provides(IScale)
class LogScale(ScaleMixin):
    id = Constant("edu.mit.synbio.cytoflow.utility.log_scale")
    name = "log"

    experiment = Instance("xiaoyu_CF._experiment.Experiment")

    # must set one of these.  they're considered in order.
    channel = Str
    condition = Str
    statistic = Tuple(Str, Str)
    error_statistic = Tuple(Str, Str)
    data = Array

    mode = Enum("mask", "clip")
    threshold = Property(Float, depends_on = "[experiment, condition, channel, statistic, error_statistic]")
    _channel_threshold = Float(0.1)

    def get_mpl_params(self, ax):
        """
        Returns a dict with the traits needed to initialize an instance of
        `matplotlib.scale.ScaleBase`
        """
        return {"nonpositive" : self.mode}

    def _set_threshold(self, threshold):
        self._channel_threshold = threshold

    def _get_threshold(self):
        if self.channel:
            return self._channel_threshold
        elif self.condition:
            cond = self.experiment[self.condition][self.experiment[self.condition] > 0]
            return cond.min()
        elif self.statistic in self.experiment.statistics \
            and not self.error_statistic in self.experiment.statistics:
            stat = self.experiment.statistics[self.statistic]
            assert is_numeric(stat)
            return stat[stat > 0].min()
        elif self.statistic in self.experiment.statistics \
            and self.error_statistic in self.experiment.statistics:
            stat = self.experiment.statistics[self.statistic]
            err_stat = self.experiment.statistics[self.error_statistic]
            stat_min = stat[stat > 0].min()

            try:
                err_min = min([x for x in [min(x) for x in err_stat] if x > 0])
                return err_min

            except (TypeError, IndexError):
                err_min = min([x for x in err_stat if stat_min - x > 0])
                return stat_min - err_min

        elif self.data.size > 0:
            return self.data[self.data > 0].min()


    def __call__(self, data):
        """
        Transforms `data` using this scale.

        Careful!  May return `NaN` if the scale domain doesn't match the data
        (ie, applying a log10 scale to negative numbers.)
        """

        # this function should work with: int, float, tuple, list, pd.Series,
        # np.ndframe.  it should return the same data type as it was passed.

        if isinstance(data, (int, float)):
            if self.mode == "mask":
                if data < self.threshold:
                    raise CytoflowError("data <= scale.threshold (currently: {})".format(self.threshold))
                else:
                    return np.log10(data)
            else:
                if data < self.threshold:
                    return np.log10(self.threshold)
                else:
                    return np.log10(data)
        elif isinstance(data, (list, tuple)):
            ret = [self.__call__(x) for x in data]
            if isinstance(data, tuple):
                return tuple(ret)
            else:
                return ret
        elif isinstance(data, (np.ndarray, pd.Series)):
            mask_value = np.nan if self.mode == "mask" else self.threshold
            x = pd.Series(data)
            x = x.mask(x < self.threshold, other = mask_value)
            ret = np.log10(x)

            if isinstance(data, pd.Series):
                return ret
            else:
                return ret.values
        else:
            raise CytoflowError("Unknown type {} passed to log_scale.__call__"
                                .format(type(data)))

    def inverse(self, data):
        """
        Transforms 'data' using the inverse of this scale.
        """

        # this function should work with: int, float, tuple, list, pd.Series,
        # np.ndframe
        if isinstance(data, (int, float)):
            return np.power(10, data)
        elif isinstance(data, (list, tuple)):
            ret = [np.power(10, x) for x in data]
            if isinstance(data, tuple):
                return tuple(ret)
            else:
                return ret
        elif isinstance(data, (np.ndarray, pd.Series)):
            return np.power(10, data)
        else:
            raise CytoflowError("Unknown type {} passed to log_scale.inverse"
                                .format(type(data)))

    def clip(self, data):
        """
        Clips data to the range of the scale function
        """

        if isinstance(data, pd.Series):
            return data.clip(lower = self.threshold)
        elif isinstance(data, np.ndarray):
            return data.clip(min = self.threshold)
        elif isinstance(data, float):
            return max(data, self.threshold)
        elif isinstance(data, int):
            data = float(data)
            return max(data, self.threshold)
        else:
            try:
                return [max(x, self.threshold) for x in data]
            except TypeError as e:
                raise CytoflowError("Unknown data type in LogScale.clip") from e

    def norm(self, vmin = None, vmax = None):
        """
        A factory function that returns `matplotlib.colors.Normalize` instance,
        which normalizes values for a `matplotlib` color palette.
        """

        if vmin is not None and vmax is not None:
            pass
        elif self.channel:
            vmin = self.experiment[self.channel].min()
            vmax = self.experiment[self.channel].max()

        elif self.condition:
            vmin = self.experiment[self.condition].min()
            vmax = self.experiment[self.condition].max()

        elif self.statistic in self.experiment.statistics:
            stat = self.experiment.statistics[self.statistic]
            try:
                vmin = min([min(x) for x in stat])
                vmax = max([max(x) for x in stat])
            except (TypeError, IndexError):
                vmin = stat.min()
                vmax = stat.max()

            if self.error_statistic in self.experiment.statistics:
                err_stat = self.experiment.statistics[self.error_statistic]
                try:
                    vmin = min([min(x) for x in err_stat])
                    vmax = max([max(x) for x in err_stat])
                except (TypeError, IndexError):
                    vmin = vmin - err_stat.min()
                    vmax = vmax + err_stat.max()
        elif self.data.size > 0:
            vmin = self.data.min()
            vmax = self.data.max()
        else:
            raise CytoflowError("Must set one of 'channel', 'condition' "
                                "or 'statistic'.")

        return matplotlib.colors.LogNorm(vmin = self.clip(vmin), vmax = self.clip(vmax))


# ===========================================================================
# hlog_scale
# ===========================================================================

@provides(IScale)
class HlogScale(ScaleMixin):
    """
    A scale that transforms the data using the ``hyperlog`` function.

    This scaling method implements a "linear-like" region around 0, and a
    "log-like" region for large values, with a smooth transition between
    them.

    The transformation has one parameter, `b`, which specifies the location of
    the transition from linear to log-like.  The default, ``500``, is good for
    18-bit scales and not good for other scales.

    Attributes
    ----------
    b : Float (default = 500)
        the location of the transition from linear to log-like.

    References
    ----------
    [1] Hyperlog-a flexible log-like transform for negative, zero, and positive
        valued data.
        Bagwell CB.
        Cytometry A. 2005 Mar;64(1):34-42.
        PMID: 15700280
        http://onlinelibrary.wiley.com/doi/10.1002/cyto.a.20114/abstract
    """

    id = Constant("edu.mit.synbio.cytoflow.utility.hlog")
    name = "hlog"

    experiment = Instance("xiaoyu_CF._experiment.Experiment")

    # what data do we use to compute scale parameters?  set one.
    channel = Str
    condition = Str
    statistic = Tuple(Str, Str)
    error_statistic = Tuple(Str, Str)
    data = Array

    range = Property(Float)
    b = Float(200, desc = "location of the log transition")

    def __call__(self, data):
        """
        Transforms `data` using this scale.

        Careful!  May return `NaN` if the scale domain doesn't match the data
        (ie, applying a log10 scale to negative numbers.)
        """

        f = _make_hlog_numeric(self.b, 1.0, np.log10(self.range))

        if isinstance(data, pd.Series):
            return data.apply(f)
        elif isinstance(data, np.ndarray):
            return f(data)
        elif isinstance(data, (int, float)):
            # numpy returns a 0-dim array.  wtf.
            return float(f(data))
        else:
            try:
                return list(map(f, data))
            except TypeError as e:
                raise CytoflowError("Unknown data type in HlogScale.__call__") from e


    def inverse(self, data):
        """
        Transforms 'data' using the inverse of this scale.
        """

        f_inv = lambda y, b = self.b, d = np.log10(self.range): hlog_inv(y, b, 1.0, d)

        if isinstance(data, pd.Series):
            return data.apply(f_inv)
        elif isinstance(data, np.ndarray):
            inverse = np.vectorize(f_inv)
            return inverse(data)
        elif isinstance(data, float):
            return f_inv(data)
        else:
            try:
                return list(map(f_inv, data))
            except TypeError as e:
                raise CytoflowError("Unknown data type in HlogScale.inverse") from e

    def clip(self, data):
        """
        Clips data to the range of the scale function
        """

        return data

    def norm(self):
        """
        A factory function that returns `matplotlib.colors.Normalize` instance,
        which normalizes values for a `matplotlib` color palette.
        """

        if self.channel:
            vmin = self.experiment[self.channel].min()
            vmax = self.experiment[self.channel].max()
        elif self.condition:
            vmin = self.experiment[self.condition].min()
            vmax = self.experiment[self.condition].max()
        elif self.statistic:
            stat = self.experiment.statistics[self.statistic]
            try:
                vmin = min([min(x) for x in stat])
                vmax = max([max(x) for x in stat])
            except (TypeError, IndexError):
                vmin = stat.min()
                vmax = stat.max()
        elif self.data is not None:
            vmin = self.data.min()
            vmax = self.data.max()
        else:
            raise CytoflowError("Must set one of 'channel', 'condition' "
                                "or 'statistic'.")

        class HlogNormalize(matplotlib.colors.Normalize):
            def __init__(self, vmin, vmax, scale):
                super().__init__(vmin, vmax)

                self._scale = scale
                matplotlib.colors.Normalize.__init__(self, vmin, vmax)

            def __call__(self, value, clip = None):
                # as implemented here, hlog already transforms onto a (0, 1)
                # scale
                scaled_value = self._scale(value)
                return np.ma.masked_array(scaled_value)

        return HlogNormalize(vmin, vmax, self)

    def _get_range(self):
        if self.experiment:
            if self.channel and self.channel in self.experiment.channels:
                if "range" in self.experiment.metadata[self.channel]:
                    return self.experiment.metadata[self.channel]["range"]
                else:
                    return self.experiment.data[self.channel].max()
            elif self.condition and self.condition in self.experiment.conditions:
                return self.experiment.data[self.condition].max()
            elif self.statistic and self.statistic in self.experiment.statistics:
                return self.experiment.statistics[self.statistic].max()
            elif self.data.size > 0:
                return self.data.max()
            else:
                return Undefined
        else:
            return Undefined

    def get_mpl_params(self):
        """
        Returns a dict with the traits needed to initialize an instance of
        `MatplotlibHlogScale`
        """

        return {"b" : self.b,
                "range" : self.range}


class MatplotlibHlogScale(HasTraits, matplotlib.scale.ScaleBase):
    """
    A class that inherits from `matplotlib.scale.ScaleBase`, which
    implements all the bits for `matplotlib` to use a new scale.
    """

    name = "hlog"

    b = Float
    range = Float

    def __init__(self, axis, **kwargs):
        HasTraits.__init__(self, **kwargs)  # @UndefinedVariable

    def get_transform(self):
        """
        Returns the matplotlib.transform instance that does the actual
        transformation
        """
        if self.b is Undefined:
            raise CytoflowError("You can't set a 'hlog' scale directly.")

        return MatplotlibHlogScale.HlogTransform(b = self.b, range = self.range)

    def set_default_locators_and_formatters(self, axis):
        """
        Set the locators and formatters to reasonable defaults for
        linear scaling.
        """
        axis.set_major_locator(HlogMajorLocator())
        axis.set_major_formatter(LogFormatterMathtext(10))
        axis.set_minor_locator(HlogMinorLocator())
        axis.set_minor_formatter(NullFormatter())

    class HlogTransform(HasTraits, transforms.Transform):
        """
        A class that implements the actual transformation
        """

        # There are two value members that must be defined.
        # ``input_dims`` and ``output_dims`` specify number of input
        # dimensions and output dimensions to the transformation.
        # These are used by the transformation framework to do some
        # error checking and prevent incompatible transformations from
        # being connected together.  When defining transforms for a
        # scale, which are, by definition, separable and have only one
        # dimension, these members should always be set to 1.
        input_dims = 1
        output_dims = 1
        is_separable = True
        has_inverse = True

        # the hyperlog params
        b = Float
        range = Float

        def __init__(self, **kwargs):
            transforms.Transform.__init__(self)
            HasTraits.__init__(self, **kwargs)  # @UndefinedVariable

        def transform_non_affine(self, values):

            f = _make_hlog_numeric(self.b, 1.0, np.log10(self.range))

            if isinstance(values, pd.Series):
                return values.apply(f)
            elif isinstance(values, np.ndarray):
                return f(values)
            elif isinstance(values, float):
                return f(values)
            else:
                raise CytoflowError("Unknown data type in MatplotlibHlogScale.HlogTransform.transform_non_affine")


        def inverted(self):
            return MatplotlibHlogScale.InvertedHlogTransform(b = self.b, range = self.range)

    class InvertedHlogTransform(HasTraits, transforms.Transform):
        """
        A class that implements the inverse transformation
        """

        input_dims = 1
        output_dims = 1
        is_separable = True
        has_inverse = True

        # the hyperlog params
        b = Float
        range = Float

        def __init__(self, **kwargs):
            transforms.Transform.__init__(self)
            HasTraits.__init__(self, **kwargs)  # @UndefinedVariable

        def transform_non_affine(self, values):

            f_inv = lambda y, b = self.b, d = np.log10(self.range): hlog_inv(y, b, 1.0, d)

            if isinstance(values, pd.Series):
                return values.apply(f_inv)
            elif isinstance(values, np.ndarray):
                inverse = np.vectorize(f_inv)
                return inverse(values)
            elif isinstance(values, float):
                return f_inv(values)
            else:
                raise CytoflowError("Unknown data type in MatplotlibHlogScale.InvertedHlogTransform.transform_non_affine")


        def inverted(self):
            return MatplotlibHlogScale.HlogTransform(b = self.b, range = self.range)


class HlogMajorLocator(Locator):
    """
    Determine the tick locations for hlog axes.
    Based on matplotlib.LogLocator
    """

    def set_params(self):
        """Empty"""
        pass

    def __call__(self):
        'Return the locations of the ticks'
        vmin, vmax = self.axis.get_view_interval()
        return self.tick_values(vmin, vmax)

    def tick_values(self, vmin, vmax):
        'Every decade, including 0 and negative'

        vmin, vmax = self.view_limits(vmin, vmax)
        max_decade = 10 ** np.ceil(np.log10(vmax))

        if vmin < 0:
            min_decade = -1.0 * 10 ** np.floor(np.log10(-1.0 * vmin))
            ticks = [-1.0 * 10 ** x for x in np.arange(np.log10(-1.0 * min_decade), 1, -1)]
            ticks.append(0.0)
            ticks.extend( [10 ** x for x in np.arange(2, np.log10(max_decade), 1)])
        else:
            ticks = [0.0] if vmin == 0.0 else []
            ticks.extend( [10 ** x for x in np.arange(1, np.log10(max_decade), 1)])

        return self.raise_if_exceeds(np.asarray(ticks))

    def view_limits(self, data_min, data_max):
        'Try to choose the view limits intelligently'

        if data_max < data_min:
            data_min, data_max = data_max, data_min

        # get the nearest tenth-decade that contains the data

        if data_max > 0:
            logs = np.ceil(np.log10(data_max))
            vmax = np.ceil(data_max / (10 ** (logs - 1))) * (10 ** (logs - 1))
        else:
            vmax = 100

        if data_min >= 0:
            vmin = 0
        else:
            logs = np.ceil(np.log10(-1.0 * data_min))
            vmin = np.floor(data_min / (10 ** (logs - 1))) * (10 ** (logs - 1))

        return transforms.nonsingular(vmin, vmax)

class HlogMinorLocator(Locator):
    """
    Determine the tick locations for hlog axes.
    Based on matplotlib.LogLocator
    """

    def set_params(self):
        """Empty"""
        pass

    def __call__(self):
        'Return the locations of the ticks'
        vmin, vmax = self.axis.get_view_interval()
        return self.tick_values(vmin, vmax)

    def tick_values(self, vmin, vmax):
        'Every tenth decade, including 0 and negative'

        vmin, vmax = self.view_limits(vmin, vmax)

        if vmin < 0:
            lt = [np.arange(10 ** x, 10 ** (x - 1), -1.0 * (10 ** (x-1)))
                  for x in np.arange(np.ceil(np.log10(-1.0 * vmin)), 1, -1)]

            # flatten and take the negative
            lt = [-1.0 * item for sublist in lt for item in sublist]

            # whoops! missed an endpoint
            lt.extend([-10.0])

            gt = [np.arange(10 ** x, 10 ** (x + 1), 10 ** x)
                  for x in np.arange(1, np.log10(vmax))]

            # flatten
            gt = [item for sublist in gt for item in sublist]

            ticks = lt
            ticks.extend(gt)
        else:
            vmin = max((vmin, 1))
            ticks = [np.arange(10 ** x, 10 ** (x + 1), 10 ** x)
                     for x in np.arange(np.log10(vmin), np.log10(vmax))]
            ticks = [item for sublist in ticks for item in sublist]

        return self.raise_if_exceeds(np.asarray(ticks))


def hlog_inv(y, b, r, d):
    '''
    Inverse of base 10 hyperlog transform.
    '''
    aux = 1.*d/r *y
    s = np.sign(y)
    if s.shape: # to catch case where input is a single number
        s[s==0] = 1
    elif s==0:
        s = 1
    return s*10**(s*aux) + b*aux - s

def _make_hlog_numeric(b, r, d):
    '''
    Return a function that numerically computes the hlog transformation for given parameter values.
    '''
    hlog_obj = lambda y, x, b, r, d: hlog_inv(y, b, r, d) - x
    find_inv = np.vectorize(lambda x: scipy.optimize.brentq(hlog_obj, -2*r, 2*r,
                                        args=(x, b, r, d)))
    return find_inv

def hlog(x, b, r, d):
    '''
    Base 10 hyperlog transform.

    Parameters
    ----------
    x : num | num iterable
        values to be transformed.
    b : num
        Parameter controling the location of the shift
        from linear to log transformation.
    r : num (default = 10**4)
        maximal transformed value.
    d : num (default = log10(2**18))
        log10 of maximal possible measured value.
        hlog_inv(r) = 10**d

    Returns
    -------
    Array of transformed values.
    '''
    hlog_fun = _make_hlog_numeric(b, r, d)
    if not hasattr(x, '__len__'): #if transforming a single number
        y = hlog_fun(x)
    else:
        n = len(x)
        if not n: #if transforming empty container
            return x
        else:
            y = hlog_fun(x)
    return y


# ===========================================================================
# Register default scales (linear and log are auto-registered; hlog is NOT)
# ===========================================================================

register_scale(LinearScale)
register_scale(LogScale)

# hlog_scale is NOT auto-registered because it is really slow.
# If you want it for your analysis, import register_scale and call:
#   register_scale(HlogScale)

matplotlib.scale.register_scale(MatplotlibHlogScale)


# ===========================================================================
# Density gate → polygon conversion
# ===========================================================================

def density_gate_to_polygon(gate_op):
    """Extract polygon vertices from a DensityGateOp's keep contour.

    Returns list of [x, y] vertices, or None if contour extraction fails.
    Warns if multiple contour segments are found.
    """
    import matplotlib.pyplot as _plt
    hist = gate_op._histogram.get(True)
    if hist is None or not hist.size:
        return None
    level = hist[gate_op._keep_xbins[True][-1], gate_op._keep_ybins[True][-1]]
    xbins = gate_op._xbins[:-1]
    ybins = gate_op._ybins[:-1]
    cs = _plt.contour(xbins, ybins, hist.T, [level])
    segs = cs.allsegs[0] if hasattr(cs, 'allsegs') else []
    if not segs:
        _plt.close()
        return None
    segs = [s for s in segs if len(s) >= 3]
    if not segs:
        _plt.close()
        return None
    best = max(segs, key=len)
    _plt.close()
    if len(segs) > 1:
        print(f"[WARN] Density gate produced {len(segs)} contour segments; "
              "keeping the largest. Consider manual polygon gating.")
    best = best[::4]
    return [[float(v[0]), float(v[1])] for v in best]
