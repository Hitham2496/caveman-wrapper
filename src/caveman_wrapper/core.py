"""
CaVEMan wrapper to replace Perl script

Core class to run caveman with parameters provided by user

This script uses multiprocessing to ensure that CPU resource
allocated by HPC schedulers are adequately used
"""
import os
import sys
import subprocess
import time
import multiprocessing
import utils

class CaveManRunner():
    """

    """
