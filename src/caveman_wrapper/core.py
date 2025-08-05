"""
CaVEMan wrapper to replace Perl script

Core class to run caveman with parameters provided by user

This script uses multiprocessing to ensure that CPU resource
allocated by HPC schedulers are adequately used
"""
import os
import shutil
import sys
import subprocess
import time
import multiprocessing
from .utils import *

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

        # See if the user has asked for help before doing any more checking
        # Exit after printing help message
        if "help" in kwargs:
            self.print_help_message()
            sys.exit()

        # Requesting manual should take precedence over checks too
        # Exit after opening manual
        if "man" in kwargs:
            self.open_caveman_wrapper_manual()
            sys.exit()

        # Raise an error if the caveman executabel is not in the path
        # Do this before checking version
        if not self.check_caveman_in_path():
            raise FileNotFoundError("`caveman` could not be found in $PATH")

        if "version" in kwargs:
            self.print_caveman_version()
            sys.exit()

        allowed_keys = CavemanFlags.__annotations__

        for key in kwargs:
            if key not in allowed_keys:
                raise ValueError(f"Key '{key}' is not recognised as an option for caveman")

            if not isinstance(kwargs[key], allowed_keys[key]):
                raise ValueError(f"The value provided for the parameter '{key}' should be of type: {allowed_keys[key]}")

            setattr(self, key, kwargs[key])

        self.run_caveman()

    def print_help_message(self):
        """
        Prints the help message for the wrapper usage
        """
        print(self.HELP_MESSAGE)
    
    def print_caveman_version(self):
        """
        Prints the version of caveman if requested
        """
        print(0.)

    def open_caveman_wrapper_manual(self):
        """
        Opens the caveman wrapper script manual page
        """
        print("This will become the manual")

    def check_caveman_in_path(self):
        """
        Attempts to locate the C executable for caveman in $PATH

        Returns:
        --------
        `bool` - 
            True if found, False otherwise
        """
        # TODO: Replace the 'True' with the commented out bit when done testing
        self.caveman_executable = shutil.which("caveman")
        return True # bool(self.caveman_executable)

    def run_caveman(self):
        """
        Run caveman with the parameters identified at initialisation
        """
        # TODO: This is an outline of what needs to be added here
        
        # Step 1. set up processes in the same way as caveman.pl::setup

        # Step 2. register processes i.e. if `threads` is given, use as many processes as requested

        # Step 3. create temporary reference to be cleaned if reference is in `*.gz.fai` format

        # Step 4. caveman_setup: if !process OR process == setup

        # Step 5. caveman_split: if !process OR process == split

        # Step 6. split_concat: if !process OR process == split_concat
        
        # Step 7. file_line_count: if !process OR process in [mstep, estep]

        # Step 8. caveman_mstep: if !process OR process == mstep
        
        # Step 9. caveman_estep: if !process OR process == estep
        
        # Step 10. caveman_add_vcf_ids: if !process OR process == add_ids

        # Step 11. caveman_flag: if !process OR process == flag OR !noflag

        # Step 12. cleanup: if !noclean

        # VALID_PROCESSES = ["setup", "split", "split_concat", "mstep", "merge", "estep", "merge_results", "add_ids", "flag"]
        if self.process:
            if not (self.process in CavemanConstants.VALID_PROCESSES):
                raise ValueError(f"Process '{self.process}' is not a valid caveman process")
        return

   # def caveman_setup(self):
   #     """
   #     Runs CAVEMAN_SETUP from the parameters used at initialisation
   #     """

   # def caveman_split(self):
   #     """
   #     Runs CAVEMAN_SPLIT from the parameters used at initialisation
   #     """

   # def caveman_mstep(self):
   #     """
   #     Runs CAVEMAN_MSTEP from the parameters used at initialisation
   #     """

   # def caveman_merge(self):
   #     """
   #     Runs CAVEMAN_MERGE from the parameters used at initialisation
   #     """ 
 
   # def caveman_estep(self):
   #     """
   #     Runs CAVEMAN_ESTEP from the parameters used at initialisation
   #     """

   # def caveman_merge_results(self):
   #     """
   #     Runs MERGE_CAVEMAN_RESULTS from the parameters used at initialisation
   #     """

   # def caveman_add_vcf_ids(self):
   #     """
   #     Runs CAVEMAN_VCF_IDS from the parameters used at initialisation
   #     """

   # def caveman_flag(self):
   #     """
   #     Runs CAVEMAN_FLAG from the parameters used at initialisation
   #     """
 
   # def caveman_split_vcf(self):
   #     """
   #     Runs CAVEMAN_VCF_SPLIT from the parameters used at initialisation
   #     """

   # def count_files(self):
   #     """
   #     Runs FILE_COUNT from the parameters used at initialisation
   #     """

   # def concat(self):
   #     """
   #     Runs concatenation for files from caveman running
   #     """

   # def concat_flagged(self):
   #     """
   #     Runs CAVEMAN_VCF_FLAGGED_CONCAT from the parameters used at initialisation
   #     """

   # def zip_flagged(self):
   #     """
   #     Zips flagged files from caveman_flag stage
   #     """

   # def pre_cleanup_zip(self):
   #     """
   #     Zip files before the cleanup stage if cleanup option is specified
   #     """

   # def limited_indices(self):
   #     """
   #     Checks whether the index is not greater than the limit or lower than 1
   #     """

   # def limited_flag_indices(self):
   #     """
   #     Check limited indices for the VCF split count
   #     """
     
   # def limited_xstep_indices(self):
   #     """
   #     Check limited indices for the split list count
   #     """

   # def load_exclude(self):
   #     """
   #     Loads the file of excludes and returns the exclude pattern
   #     """

   # def valid_seq_indices(self):
   #     """
   #     Checks that the sequence indices provided are valid
   #     """

   # def extend_no_analysis(self):
   #     """
   #     Extends no analysis in caveman_merge results
   #     """
