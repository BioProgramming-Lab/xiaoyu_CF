"""
xiaoyu_CF - self-contained flow cytometry analysis for 96-well plates.

Import as: import xiaoyu_CF as flow
"""

__version__ = "0.3.0"

import mpl_scatter_density
import matplotlib.pyplot as plt
from matplotlib.widgets import PolygonSelector
import glob
import os.path
import numpy as np
import pandas as pd

# suppress meaningless warnings from seaborn
import warnings
warnings.filterwarnings('ignore', '.*IPython widgets are experimental.*')
warnings.filterwarnings('ignore', 'axes.color_cycle is deprecated and replaced with axes.prop_cycle')

# matplotlib logging filter
import matplotlib.text
import logging
from ._utility import MplFilter

if hasattr(matplotlib.text, "_log"):
    matplotlib.text._log.addFilter(MplFilter())

# core experiment class
from ._experiment import Experiment

# operations
from ._operations import (Tube, ImportOp, RangeOp, PolygonOp,
                          DensityGateOp, AutofluorescenceOp,
                          BleedthroughLinearOp)

# views
from ._views import (HistogramView, ScatterplotView, DensityView)

# utility functions
from ._utility import (geom_mean, geom_sd, geom_sd_range,
                       geom_sem, geom_sem_range)
from ._utility import ci, percentiles
from ._utility import set_default_scale, get_default_scale

from ._helpers import *
from ._helpers import _find_channel  # explicit: _-prefixed names skipped by *
