"""
CaVEMan wrapper to replace Perl script

Core class to run caveman with parameters provided by user

This script uses multiprocessing to ensure that CPU resource
allocated by HPC schedulers are adequately used
"""
import shutil
import subprocess
import time
import multiprocessing
import re
from .utils import *

class CavemanRunner():
    """
    Class to initialise, distribute resources for, and run caveman
    on HPC, based on the original wrapper keyword arguments.
    """

    HELP_MESSAGE = f"Wrapper for CaVEMan, usage: caveman.py [kwargs]\n{CavemanFlags.short_flags}"

    # Container for value of normal contamination if it is provided
    # so as not to overwrite the filename we get it from (!)
    normcont_value = None

    # Containers for output subdirs
    tmp_dir = None
    results_dir = None
    progress_dir = None
    log_dir = None

    # Containers for configs
    cave_cfg = None
    cave_alg = None
    cave_parr = None
    cave_carr = None

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
        
        #### i. check files (ref, tumbam, normbam, ignore) exists
        for item in ["noflag", "noclean"]: 
            if (getattr(self, item, None) is None):
                setattr(self, item, False)

        # set working dir as object member
        #setattr(self, working_dir, os.getcwd())

        if not os.path.isfile(self.reference):
            raise ValueError(f"Reference file {self.reference} could not be found")

        if not (os.path.isfile(f"{self.tumour_bam}.bai") or
                os.path.isfile(f"{self.tumour_bam}.csi") or
                os.path.isfile(f"{self.tumour_bam}.crai")):
            raise ValueError(f"Tumour BAM file {self.tumour_bam} could not be found")

        if not (os.path.isfile(f"{self.normal_bam}.bai") or
                os.path.isfile(f"{self.normal_bam}.csi") or
                os.path.isfile(f"{self.normal_bam}.crai")):
            raise ValueError(f"Tumour BAM file {self.normal_bam} could not be found")

        if not os.path.isfile(self.ignore_file):
            raise ValueError(f"Ignore file {self.ignore_file} could not be found")
        
        #### ii. check (tumcn, normcn) exists if provided
        control_args = {"tumour_cn" : "tum_cn_default" , "normal_cn" : "norm_cn_default" }
        for key in control_args:
            # Check if tumour and normal control exist and are defined
             filename = getattr(self, key, None)
             if filename is None:
                 raise ValueError(f"File '{filename}' for option '{key}' could not be found")

             # Check the file actually exists
             if os.path.isfile(filename):
                 default_undefined = not hasattr(self, control_args[key])

                 # Check the default value for the control is defined if the file is empty
                 if os.path.getsize(filename) == 0 and default_undefined:
                     raise ValueError(f"If file specified for {key} is empty, {control_args[key]} should be defined.")

        #### iii. delete (process, index, limit, exclude) if not provided
        for del_flag in ["process", "index", "limit", "exclude"]:
            if getattr(self, del_flag, None) is None:
                delattr(self, del_flag)

        #### iv. set read-count to default unless provided
        if getattr(self, "read_count", None) is None:
            setattr(self, "read_count", CavemanConstants.SPLIT_STEP_READ_COUNT) 

        #### v. check outdir, if exists and non-empty throw error and quit
        outdir = getattr(self, "outdir", None)
        check_outdir(outdir)
        
        log_dir = f"{outdir}/logs"
        if os.path.isdir(log_dir):
            raise OSError(f"Presence of {log_dir} indicates that an analysis has been completed, delete to rerun")

        #### vi. check (flagconfig, flagtovcfconfig, germline-indel-bed) if provided
        for config_arg in ["flagConfig", "flagToVcfConfig", "germindel"]:
            config_val = getattr(self, config_arg, None)
            if not os.path.isfile(config_val):
                raise ValueError(f"File '{config_val}' for option '{config_arg}' could not be found")

        #### vii. get assembly and check reference provided is the fasta fai file

        # TODO: Implement more than a stub for this function
        self.get_species_assembly_from_bam()

        for bam_key in ["species", "species_assembly"]:
            if getattr(self, bam_key, None) is None:
                raise ValueError(f"{bam_key} must be defined, see BAM header for options")
        
        if getattr(self, "seqType", None) not in CavemanConstants.VALID_SEQ_TYPES:
            raise ValueError(f"seqType must be one of: {', '.join(CavemanConstants.VALID_SEQ_TYPES)}")

        #### viii. check protocols provided for (normprot, tumprot) otherwise set to default
        for protocol in ["tumour_protocol", "normal_protocol"]:

            protocol_option = getattr(self, protocol, None)

            if protocol_option is None:
                setattr(protocol, CavemanConstants.DEFAULT_PROTOCOL)

            elif protocol_option not in CavemanConstants.VALID_PROTOCOLS:
                raise ValueError(f"{protocol} option '{protocol_option}' is not recognised, protocols "
                                 f"must be one of: {', '.join(CavemanConstants.VALID_PROTOCOLS)}")

        #### ix. set threads to 1 unless specified
        if getattr(self, "threads", None) is None:
            setattr(self, "threads", 1)

        #### x. if normcont is defined, try to extract it from the file, fail if it doesn't work
        if getattr(self, "normal_contamination", None):

            if not os.path.isfile(self.normal_contamination):
                raise ValueError(f"Reference file {self.reference} could not be found")

            elif os.path.getsize(self.normal_contamination) == 0:
                raise ValueError(f"Normal contamination file {self.normal_contamination} is empty.")

            # Search for data of the form 'NormalContamination <number>' 
            with open(self.normal_contamination, "r") as norm_cont:
                norm_cont_data = norm_cont.readlines()

            for line in norm_cont_data:
                match_line = re.match(r"^NormalContamination\s(\d+\.?\d*)$", line)
                if match_line:
                    setattr(self, "normcont_value", float(match_line.group(1)))
                    break

            if getattr(self, "normcont_value", None) is None:
                raise RuntimeError(f"Could not find NormalContamination value in {self.normal_contamination}")

            # Value needs to be between 0 and 1, negatives will not give a match
            elif self.normcont_value > 1.:
                raise ValueError(f"NormalContamination value {self.normcont_value} must be between 0. and 1.")

        else:
            setattr(self, "normcont_value", CavemanConstants.DEFAULT_NORMCONT)

        #### xi. create objects for output files, starting with tmp dir in outdir
        setattr(self, "tmp_dir", f"{self.outdir}/tmpCaveman")
        setattr(self, "results_dir", f"{self.tmp_dir}/results")
        setattr(self, "progress_dir", f"{self.tmp_dir}/progresss")

        if getattr(self, "logs", None):
            setattr(self, "log_dir", f"{self.tmp_dir}/{self.logs}")
        else:
            setattr(self, "log_dir", f"{self.tmp_dir}/logs")

        setattr(self, "cave_cfg", f"{self.tmp_dir}/{CavemanConstants.CAVEMAN_CONFIG}")
        setattr(self, "cave_alg", f"{self.tmp_dir}/{CavemanConstants.CAVEMAN_ALG_BEAN}")
        setattr(self, "cave_parr", f"{self.tmp_dir}/{CavemanConstants.CAVEMAN_PROB_ARR}")
        setattr(self, "cave_carr", f"{self.tmp_dir}/{CavemanConstants.CAVEMAN_COV_ARR}")

        #### xii. check process - if provided - is valid, set max index if it is provided otherwise default

        # VALID_PROCESSES = ["setup", "split", "split_concat", "mstep", "merge", "estep", "merge_results", "add_ids", "flag"]
        if self.process:
            if not (self.process in CavemanConstants.VALID_PROCESSES):
                raise ValueError(f"Process '{self.process}' is not a valid caveman process")

        # xiii. make paths!
        return

    def get_species_assembly_from_bam(self):
         """
         Gets species assemblies from BAM files, sets `species`, `species_assembly`
         from reference
         """
         # TODO: Figure out how to implement this one with pysam?
         # This is placeholder code for the use of the function 
         setattr(self, "species", "human")
         setattr(self, "species_assembly", "human_assembly")

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
