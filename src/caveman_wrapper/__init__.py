"""
CaVEMan wrapper to replace Perl script

This script uses multiprocessing to ensure that CPU resource
allocated by HPC schedulers are adequately used
"""
from .utils import *
from .core import *
