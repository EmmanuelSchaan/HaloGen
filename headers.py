#from pylab import *


import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import special, optimize, integrate, stats
from scipy.interpolate import UnivariateSpline, RectBivariateSpline, interp1d, interp2d, BarycentricInterpolator
from time import time
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from timeit import timeit
from time import time
from copy import copy
import sys

# parallelizing "map"
# version that works when the function is a class module
from pathos.multiprocessing import ProcessingPool as Pool

import basic_functions
reload(basic_functions)
from basic_functions import *

