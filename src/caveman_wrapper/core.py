"""
CaVEMan wrapper to replace Perl script

Core class to run caveman with parameters provided by user

This script uses multiprocessing to ensure that CPU resource
allocated by HPC schedulers are adequately used
"""
import subprocess
import time
from multiprocessing import Pool
import pysam
import warnings
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
    split_list = None

    # Container for valid fai indices
    valid_fai_idx = None

    # Maximum indices for processes
    index_max = {
        "setup" : 1,
        "split" : -1,
        "split_concat" : 1,
        "mstep" : -1,
        "merge" : 1,
        "estep" : -1,
        "merge_results" : 1,
        "add_ids" : 1,
        "flag" : 1
    }

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

        self.setup_caveman_environment()
        #self.run_caveman()

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
        return True #bool(self.caveman_executable)

    def setup_caveman_environment(self):
        """
        Setup caveman with the parameters identified at initialisation
        """ 
        #### Step 1. Check files (ref, tumbam, normbam, ignore) exist
        for item in ["noflag", "noclean"]: 
            if (getattr(self, item, None) is None):
                setattr(self, item, False)

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
        
        #### Step 2. Check (tumcn, normcn) exist if provided
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

        #### Step 3. Delete (process, index, limit, exclude) if not provided
        #### This is not needed since the attributes will not exist if not provided

        #for del_flag in ["process", "index", "limit", "exclude"]:
        #    if getattr(self, del_flag, None) is None:
        #        delattr(self, del_flag)

        #### Step 4. Set read-count to default unless provided
        if getattr(self, "read_count", None) is None:
            setattr(self, "read_count", CavemanConstants.SPLIT_STEP_READ_COUNT) 

        #### Step 5. Check outdir, if exists and non-empty throw error and quit
        outdir = getattr(self, "outdir", None)
        check_outdir(outdir)
        
        log_dir = f"{outdir}/logs"
        if os.path.isdir(log_dir):
            raise OSError(f"Presence of {log_dir} indicates that an analysis has been completed, delete to rerun")

        #### Step 6. Check (flagconfig, flagtovcfconfig, germline-indel-bed) if provided
        for config_arg in ["flagConfig", "flagToVcfConfig", "germindel"]:
            config_val = getattr(self, config_arg, None)
            if not os.path.isfile(config_val):
                raise ValueError(f"File '{config_val}' for option '{config_arg}' could not be found")

        #### Step 7. Get assembly and check reference provided is the fasta fai file

        # TODO: Implement more than a stub for this function
        self.get_species_assembly_from_bam()

        for bam_key in ["species", "species_assembly"]:
            if getattr(self, bam_key, None) is None:
                raise ValueError(f"{bam_key} must be defined, see BAM header for options")
        
        if getattr(self, "seqType", None) not in CavemanConstants.VALID_SEQ_TYPES:
            raise ValueError(f"seqType must be one of: {', '.join(CavemanConstants.VALID_SEQ_TYPES)}")

        #### Step 8. Check protocols provided for (normprot, tumprot) otherwise set to default
        for protocol in ["tumour_protocol", "normal_protocol"]:

            protocol_option = getattr(self, protocol, None)

            if protocol_option is None:
                setattr(protocol, CavemanConstants.DEFAULT_PROTOCOL)

            elif protocol_option not in CavemanConstants.VALID_PROTOCOLS:
                raise ValueError(f"{protocol} option '{protocol_option}' is not recognised, protocols "
                                 f"must be one of: {', '.join(CavemanConstants.VALID_PROTOCOLS)}")

        #### Step 9. Set threads to 1 unless specified
        if getattr(self, "threads", None) is None:
            setattr(self, "threads", 1)

        #### Step 10. If normcont is defined, try to extract it from the file, fail if it doesn't work
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

        #### Step 11. create objects for output files, starting with tmp dir in outdir
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
        setattr(self, "split_list", f"{self.tmp_dir}/splitList")

        # TODO: subvcf, snpvcf, noanalysisbed dynamic output file generation will be implemented in member methods

        #### Step 12. check process - if provided - is valid, set max index if it is provided otherwise default
        if getattr(self, "process", None):
            if not (self.process in CavemanConstants.VALID_PROCESSES):
                raise ValueError(f"Process '{self.process}' is not a valid caveman process")

            if getattr(self, "index", None):
                cur_idx_max = self.index_max[self.process]
                if cur_idx_max == -1:
                    if getattr(self, "limit", None):
                        cur_idx_max = self.limit
                    else:
                        cur_idx_max = self.valid_index()

                if self.index < 1 or self.index > cur_idx_max:
                    raise ValueError(f"Based on reference and exclude option, index must be between 1 and {cur_idx_max}")

                if cur_idx_max == 0:
                    raise ValueError(f"No maximum index has been provided for the process type '{self.process}'")

                valid_index_by_factor("index", self.index, cur_idx_max, 1)
        
        elif getattr(self, "index", None):
            raise AttributeError(f"Index cannot be defined without a process")

        #### Step 13. Create paths
        for path in [self.tmp_dir, self.results_dir, self.progress_dir, self.log_dir]:
            if not os.path.isdir(path):
                os.mkdir(path)
        return

    def valid_index(self):
        """
        Checks that the index provided in options is valid, otherwise returns 0
        """
        # Process should already be provided when we call this method, but exit if not:
        if getattr(self, "process", None) is None:
            raise AttributeError(f"To calculate a valid index, the process must be provided")

        if self.process == "split":
            return file_line_count(self.reference)

        elif self.process == "mstep" or self.process == "estep":
            return file_line_count(self.split_list)

        return 0

    def get_species_assembly_from_bam(self):
         """
         Gets species assemblies from BAM files, sets `species`, `species_assembly`
         from reference
         """
         SP_ASS_MESSAGE = f"{} defined at commandline {} does not match that in the BAM file {}. Defaulting to BAM file value."
         bam = pysam.AlignmentFile(self.tumour_bam, "rb")
         header = bam.header
         sq_entries = header.get('SQ', [])

         for entry in sq_entries:
             if "AS" in entry:
                 assembly = entry["AS"]
                 if getattr(self, "species_assembly", None) != assembly:
                     warnings.warn(SP_ASS_MESSAGE.format("Assembly", opts["species-assembly"], assembly))
                 setattr(self, "species_assembly", assembly)

             if "SP" in entry:
                 species = entry["SP"]
                 if getattr(self, "species", None) != species:
                     warnings.warn(SP_ASS_MESSAGE.format("Species", opts["species"], species))
                 setattr(self, "species", species)

             # Only need to process the first sequence line.
             break

         return

    def run_caveman(self):
        """
        Runs caveman with parameters from initialisation and setup
        """
        no_process = getattr(self, "process", None) is None

        #### Step 1. setup has been completed

        #### Step 2. register processes, required for Perl, not for Python

        #### Step 3. create temporary reference to be cleaned if reference is in `*.gz.fai` format
        if re.search(".*\.gz\.fai$", self.reference):
            tmp_ref = f"{self.tmp_dir}/genome.fa"
            if not os.path.isfile(tmp_ref):
                non_indexed_reference = self.reference.replace(".fai","")
                gunzip_file(non_indexed_reference, tmp_ref)
                # I suspect this may be required for bug-compatibility, since the reference
                # is not in-place modified in the originl Perl code (caveman.pl:L104-114)
                shutil.copy(self.reference, f"{tmp_ref}.fai")

        # Step 4. caveman_setup: if !process OR process == setup
        if no_process or getattr(self, "process", None) == "setup":
            self.caveman_setup()

        # Step 5. caveman_split: if !process OR process == split
        if no_process or getattr(self, "process", None) == "split":

            valid_indices = self.valid_seq_indices()
            contig_count = len(valid_indices)
            setattr(self, "valid_fai_idx", valid_indices)

            self.caveman_split()

        # Step 6. split_concat: if !process OR process == split_concat
        if no_process or getattr(self, "process", None) == "split_concat":
            self.concat()
        
        # Step 7. file_line_count: if !process OR process in [mstep, estep]

        # Step 8. caveman_mstep: if !process OR process == mstep
        
        # Step 9. caveman_estep: if !process OR process == estep
        
        # Step 10. caveman_add_vcf_ids: if !process OR process == add_ids

        # Step 11. caveman_flag: if !process OR process == flag OR !noflag

        # Step 12. cleanup: if !noclean

    def caveman_setup(self):
        """
        Runs CAVEMAN_SETUP from the parameters used at initialisation.

        Returns:
        --------
        `bool` - True if the run is successful, False otherwise, based on
        exit code of job, and output of `touch_success`.
        """
        if success_exists(self.progress_dir, 0):
            return True
        
        # In case function is being called manually, check caveman
        # is still in the path.
        if not self.check_caveman_in_path():
            raise FileNotFoundError("`caveman` could not be found in $PATH")

        command = f"caveman setup "
                  f"-t {self.tumour_bam} "
                  f"-n {self.normal_bam} "
                  f"-r {self.reference} "
                  f"-g {self.ignore_file} "
                  f"-l {self.split_list} "
                  f"-f {self.results_dir} "
                  f"-c {self.cave_cfg} "
                  f"-a {self.alg_bean} "

        if getattr(self, "normal_cn", None):
            command += f"-j {self.normal_cn}"

        if getattr(self, "tumour_cn", None):
            command += f"-e {self.tumour_cn}"

        # Only one process is required for setup, set index to 0.
        index = 0
        final_result = worker(self.log_dir, command, index)

        if not final_result["success"]:
            return False

        return touch_success(self.progress_dir, final_result["index"])

    def caveman_split(self, index=None):
        """
        Runs CAVEMAN_SPLIT from the parameters used at initialisation,
        or optionally from a specified index (starting from 1)

        If `index` is provided, the caveman split workflow is run with
        the initialisation settings, for only the value of valid_fai_idx
        at position `index - 1`. The number of processes for the pool
        is set to 1.

        Otherwise, a pool of `self.threads` size is instantiated and
        caveman split is run asynchronously for each value of valid_fai_idx
        as calculated from the initialisation

        A bool is returned, trivially equal to True/False if the run is/is
        not successful to maintain consistency with Perl wrapper.

        Parameters:
        -----------
        `index` : `int` -
            Index of `self.valid_vai_idx`, starting from 1, to run split
            method on.

        Returns:
        --------
        `bool` - 
            True if the run is successful, False if errors arise in running
        """
        index_list = None
        num_procs = self.threads
        if index:
            index_list = [self.valid_fai_idx[index-1]]
            num_procs = 1
            if index != self.index:
                return True
            if success_exists(self.progress_dir, index):
                return True
        else:
            index_list = self.valid_fai_idx

        # In case function is being called manually, check caveman
        # is still in the path.
        if not self.check_caveman_in_path():
            raise FileNotFoundError("`caveman` could not be found in $PATH")

        errors_raised = False
        with Pool(processes=num_procs) as pool:
            async_results = []
            for value in index_list:

                command = f"caveman split -i {value} "
                          f"-f {self.cave_cfg} "
                          f"-e {self.read_count}"

                result = pool.apply_async(worker, args=(self.log_dir, command, index,))
                async_results.append(result)

            for item in async_results:
                final_result = item.get()
                if not final_result["success"]:
                    print(f"Split stage failed for {final_result['index']}", file=sys.stderr)
                    print(f"Error: {final_result['error']}", file=sys.stderr)
                    errors_raised = True
                elif not touch_success(self.progress_dir, final_result["index"]):
                    # touch_success os called, so if not successful count as an error
                    errors_raised = True
                else:
                    # Explicitly show we continue if touch_success is True
                    continue

        successful_run = not errors_raised

        return successful_run

    def caveman_mstep(self, index=None):
        """
        Runs CAVEMAN_MSTEP from the parameters used at initialisation,
        or optionally from a specified index (starting from 1).

        If `index` is provided, the caveman split workflow is run with
        the initialisation settings, for only the value of valid_fai_idx
        at position `index - 1`. The number of processes for the pool
        is set to 1.

        Otherwise, a pool of `self.threads` size is instantiated and
        caveman split is run asynchronously for each value of valid_fai_idx
        as calculated from the initialisation

        A bool is returned, trivially equal to True/False if the run is/is
        not successful to maintain consistency with Perl wrapper.

        Parameters:
        -----------
        `index` : `int` -
            Index of `self.valid_vai_idx`, starting from 1, to run split
            method on.

        Returns:
        --------
        `bool` - 
            True if the run is successful, False if errors arise in running
        """
        return False
        #index_list = None
        #num_procs = self.threads
        #if index:
        #    index_list = [self.valid_fai_idx[index-1]]
        #    num_procs = 1
        #    if index != self.index:
        #        return True
        #    if success_exists(self.progress_dir, index):
        #        return True
        #else:
        #    index_list = self.valid_fai_idx



        # In case function is being called manually, check caveman
        # is still in the path.
        if not self.check_caveman_in_path():
            raise FileNotFoundError("`caveman` could not be found in $PATH")

        errors_raised = False
        with Pool(processes=num_procs) as pool:
            async_results = []
            for value in index_list:

                command = f"caveman mstep -i {value} "
                          f"-f {self.cave_cfg} "
                          f"-e {self.read_count}"

                result = pool.apply_async(worker, args=(self.log_dir, command, index,))
                async_results.append(result)

            for item in async_results:
                final_result = item.get()
                if not final_result["success"]:
                    print(f"Split stage failed for {final_result['index']}", file=sys.stderr)
                    print(f"Error: {final_result['error']}", file=sys.stderr)
                    errors_raised = True
                elif not touch_success(self.progress_dir, final_result["index"]):
                    # touch_success os called, so if not successful count as an error
                    errors_raised = True
                else:
                    # Explicitly show we continue if touch_success is True
                    continue

        successful_run = not errors_raised

        return successful_run


    def caveman_merge(self):
        """
        Runs CAVEMAN_MERGE from the parameters used at initialisation
        """ 
        # In case function is being called manually, check caveman
        # is still in the path.
        if not self.check_caveman_in_path():
            raise FileNotFoundError("`caveman` could not be found in $PATH")

        if success_exists(self.progress_dir, 0):
            return True
        
        command = f"caveman merge "
                  f"-c {self.cave_carr} "
                  f"-p {self.cave_parr} "
                  f"-f {self.cave_cfg}"

        # Only one process is required for setup, set index to 0.
        index = 0
        final_result = worker(command, index)

        if not final_result["success"]:
            return False

        return touch_success(self.progress_dir, final_result["index"])
 
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

    def concat(self):
        """
        Runs concatenation for files from caveman running
        """
        if success_exists(self.progress_dir, 0):
            return True
        
        # Concatenate all {self.split_list}.* files to {self.split_list}
        command = f"cat {self.split_list}.* > {self.split_list}"

        # Only one process is required for setup, set index to 0.
        index = 0
        final_result = worker(self.log_dir, command, index)

        if not final_result["success"]:
            return False

        return touch_success(self.progress_dir, final_result["index"])


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

   # def limited_xstep_indices(self):
   #     """
   #     Check limited indices for the split list count
   #     """

   # def limited_indices(self):
   #     """
   #     Checks whether the index is not greater than the limit or lower than 1
   #     """

   # def limited_flag_indices(self):
   #     """
   #     Check limited indices for the VCF split count
   #     """

    def load_exclude(self):
        """
        Loads the file of excludes and returns the exclude pattern
        """
        exclude_patterns = []

        if getattr(self, "exclude", None) is None:
            return exclude_patterns

        for pattern in self.exclude.split(","):
            pattern.replace("%", ".+")
            exclude_patterns.append(pattern)

        return exclude_patterns

    def valid_seq_indices(self):
        """
        Checks that the sequence indices provided are valid
        """
        valid_indices = []

        exclude_patterns = self.load_exclude()

        with open(self.reference, "r") as ifs_reference:
            # Line numbers start from 1
            for idx, line in enumerate(ifs_reference, start=1):
                seq_name = line.split('\t')[0]
                # Check if it matches any exclude pattern
                if not any(pattern.match(seq_name) for pattern in exclude_patterns):
                    valid_indices.append(idx)

        return valid_indices

   # def extend_no_analysis(self):
   #     """
   #     Extends no analysis in caveman_merge results
   #     """
