#!/usr/bin/env python3.8
# coding: latin-1

# (c) Massachusetts Institute of Technology 2015-2018
# (c) Brian Teague 2018-2022
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
cytoflow
--------

``cytoflow`` is a package for quantitative, reproducible analysis of flow 
cytometry data.  

Written by Brian Teague, bpteague@gmail.com

Copyright Massachusetts Institute of Technology 2015-2018

Copyright Brian Teague 2018-2021
"""

# check python version
import sys
if sys.version_info < (3, 4):
    raise Exception("Cytoflow requires Python 3.4 or later")

# suppress meaningless warnings from seaborn
import warnings
warnings.filterwarnings('ignore', '.*IPython widgets are experimental.*')
warnings.filterwarnings('ignore', 'axes.color_cycle is deprecated and replaced with axes.prop_cycle')

# and matplotlib 3.1.1 -- there's some weird interaction with seaborn here.
import matplotlib.text
import logging

        
from .utility.logging import MplFilter
if hasattr(matplotlib.text, "_log"):
    matplotlib.text._log.addFilter(MplFilter())

# keep track of whether we're running in the GUI.
# there is the occasional place where we differ in behavior
RUNNING_IN_GUI = False

# basics
from .experiment import Experiment
from .operations.import_op import ImportOp, Tube

# gates
from .operations.range import RangeOp
from .operations.polygon import PolygonOp

# TASBE
from .operations.autofluorescence import AutofluorescenceOp
from .operations.bleedthrough_linear import BleedthroughLinearOp

# data-driven
from .operations.density import DensityGateOp

# views
from .views.histogram import HistogramView
from .views.scatterplot import ScatterplotView
from .views.densityplot import DensityView

# util
from .utility.util_functions import (geom_mean, geom_sd, geom_sd_range,
                                     geom_sem, geom_sem_range)
from .utility.algorithms import (ci, percentiles)
from .utility.scale import set_default_scale, get_default_scale

from ._version import get_versions  # @UnresolvedImport
__version__ = get_versions()['version']
del get_versions
