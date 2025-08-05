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

    def concat(self):
        """
        Runs concatenation for files from caveman running
        """

    def concat_flagged(self):
        """
        Runs CAVEMAN_VCF_FLAGGED_CONCAT from the parameters used at initialisation
        """

    def zip_flagged(self):
        """
        Zips flagged files from caveman_flag stage
        """

    def pre_cleanup_zip(self):
        """
        Zip files before the cleanup stage if cleanup option is specified
        """

    def limited_indices(self):
        """
        Checks whether the index is not greater than the limit or lower than 1
        """

    def limited_flag_indices(self):
        """
        Check limited indices for the VCF split count
        """
        
    def limited_xstep_indices(self):
        """
        Check limited indices for the split list count
        """

    def load_exclude(self):
        """
        Loads the file of excludes and returns the exclude pattern
        """

    def valid_seq_indices(self):
        """
        Checks that the sequence indices provided are valid
        """

    def extend_no_analysis(self):
        """
        Extends no analysis in caveman_merge results
        """
