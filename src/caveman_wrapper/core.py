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

        self.setup_caveman()

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

    def setup_caveman(self):
        """
        Setup caveman with the parameters identified at initialisation
        """
        
        # i. check files (ref, tumbam, normbam, ignore)  exists
        if not self.noflag:
            self.noflag = False

        if not self.noclean:
            self.noflag = False

        # set working dir as object member
        self.working_dir = os.getcwd()

        if not os.path.isfile(self.reference):
            raise ValueError("Reference file {self.reference} could not be found")

        if not (os.path.isfile("f{self.tumour_bam}.bai") or
                os.path.isfile("f{self.tumour_bam}.csi") or
                os.path.isfile("f{self.tumour_bam}.crai")):
            raise ValueError("Tumour BAM file {self.tumour_bam} could not be found")

        if not (os.path.isfile("f{self.normal_bam}.bai") or
                os.path.isfile("f{self.normal_bam}.csi") or
                os.path.isfile("f{self.normal_bam}.crai")):
            raise ValueError("Tumour BAM file {self.normal_bam} could not be found")

        if not os.path.isfile(self.ignore_file):
            raise ValueError("Ignore file {self.ignore_file} could not be found")
        
        # ii. check (tumcn, normcn) exists if provided
        control_args = {"tumour_cn" : "tum_cn_default" , "normal_cn" : "norm_cn_default" }
        for key in control_args:
            # Check if tumour and normal control exist and are defined
            if hasattr(self, key) and (getattr(self, key) is not None):
                filename = getattr(self, key)

                # Check the file actually exists
                if os.path.isfile(filename):
                    default_defined = (not hasattr(self, control_args[key])) or (getattr(self, control_args[key]) is None)

                    # Check the default value for the control is defined if the file is empty
                    if os.path.getsize(filename) == 0 and default_defined:
                        raise ValueError(f"If file specified for {key} is empty, {control_args[key]} should be defined.")

                else:
                    raise ValueError(f"File '{filename}' for option '{key}' could not be found")
                
        # iii. delete (process, index, limit, exclude) if not provided
        for del_flag in ["process", "index", "limit", "exclude"]:
            if hasattr(self, del_flag) and (getattr(self, del_flag) is None):
                delattr(self, del_flag)

        # iv. set read-count to default unless provided
        if hasattr(self, "read_count") and (getattr(self, "read_count") is None):
           setattr(self, "read_count", CavemanConstants.SPLIT_STEP_READ_COUNT) 

        # v. check outdir, if exists throw error and quit
        # vi. check (flagconfig, flagtovcfconfig, germline-indel-bed) if provided
        # vii. check reference provided is the fasta fai file
        # viii. check protocols provided for (normprot, tumprot) otherwise set to default
        # ix. set threads to 1 unless specified
        # x. if normcont is defined, try to extract it from the file, fail if it doesn't work
        # xi. create objects for output files, starting with tmp dir in outdir
        # xii. check process - if provided - is valid, set max index if it is provided otherwise default

        # VALID_PROCESSES = ["setup", "split", "split_concat", "mstep", "merge", "estep", "merge_results", "add_ids", "flag"]
        if self.process:
            if not (self.process in CavemanConstants.VALID_PROCESSES):
                raise ValueError(f"Process '{self.process}' is not a valid caveman process")

        # xiii. make paths!


     def run_caveman(self):
         """
         Runs caveman with parameters from initialisation and setup
         """
        # Step 1. setup has been completed

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
