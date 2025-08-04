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
from .utils import CavemanFlags

class CavemanRunner():
    """
    Class to initialise, distribute resources for, and run caveman
    on HPC, based on the original wrapper keyword arguments.
    """

    HELP_MESSAGE = f"Wrapper for CaVEMan, usage: caveman.py [kwargs]\n{CavemanFlags.short_flags}"

    def __init__(self, **kwargs):
        """
        Initialises the runner from allowed key word arguments
        """
        if kwargs["help"]:
            print(self.HELP_MESSAGE)
            return

    def caveman_setup(self):
        """
        Runs CAVEMAN_SETUP from the parameters used at initialisation
        """

    def caveman_split(self):
        """
        Runs CAVEMAN_SPLIT from the parameters used at initialisation
        """

    def caveman_merge(self):
        """
        Runs CAVEMAN_MERGE from the parameters used at initialisation
        """
    
    def caveman_mstep(self):
        """
        Runs CAVEMAN_MSTEP from the parameters used at initialisation
        """
    
    def caveman_estep(self):
        """
        Runs CAVEMAN_ESTEP from the parameters used at initialisation
        """

    def caveman_merge_results(self):
        """
        Runs MERGE_CAVEMAN_RESULTS from the parameters used at initialisation
        """

    def caveman_add_vcf_ids(self):
        """
        Runs CAVEMAN_VCF_IDS from the parameters used at initialisation
        """
    
    def caveman_split_vcf(self):
        """
        Runs CAVEMAN_VCF_SPLIT from the parameters used at initialisation
        """

    def count_files(self):
        """
        Runs FILE_COUNT from the parameters used at initialisation
        """

    def caveman_flag(self):
        """
        Runs CAVEMAN_FLAG from the parameters used at initialisation
        """

    def concat_flagged(self):
        """
        Runs CAVEMAN_VCF_FLAGGED_CONCAT from the parameters used at initialisation
        """
