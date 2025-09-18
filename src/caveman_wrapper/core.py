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
    alg_bean = None

    # Container for valid fai indices, vcf split counts
    for_split = None
    split_count = None
    valid_fai_idx = None
    vcf_split_counts = None

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
        if kwargs.get("help", False):
            self.print_help_message()
            sys.exit()

        # Requesting manual should take precedence over checks too
        # Exit after opening manual
        if kwargs.get("man", False):
            self.open_caveman_wrapper_manual()
            sys.exit()

        # Raise an error if the caveman executabel is not in the path
        if not self.check_exec_in_path():
            raise FileNotFoundError("`caveman` could not be found in $PATH")

        if kwargs.get("version", False):
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
        # In case function is being called manually, check caveman
        # is still in the path.
        if not self.check_exec_in_path():
            raise FileNotFoundError("`caveman` could not be found in $PATH")
        command = f"caveman --version"
        subprocess.run(command.split(" "), check=True)

    def open_caveman_wrapper_manual(self):
        """
        Opens the caveman wrapper script manual page
        """
        print("This will become the manual")

    def check_exec_in_path(self, executable="caveman"):
        """
        Attempts to locate an executable for caveman in $PATH,
        by default looks for `caveman`.

        Parameters:
        -----------
        `executable` : `str` - 
            Name of executable to look for

        Returns:
        --------
        `bool` - 
            True if found, False otherwise
        """
        executable_found = shutil.which(executable)
        #TODO: REMOVE TRUE WHEN TESTING FINISHED
        return True #bool(executable_found)

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
            if config_val and not os.path.isfile(config_val):
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

        setattr(self, "subvcf", f"{self.results_dir}/%/%.muts.vcf.gz")
        setattr(self, "snpvcf", f"{self.results_dir}/%/%.snps.vcf.gz")
        setattr(self, "noanalysisbed", f"{self.results_dir}/%/%.no_analysis.bed")

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
         bam = pysam.AlignmentFile(self.tumour_bam, "rb")
         header = bam.header
         sq_entries = header.get('SQ', [])

         for entry in sq_entries:
             if "AS" in entry:
                  assembly = entry["AS"]
                  if getattr(self, "species_assembly", None) != assembly:
                       as_warn_message = f"Assembly defined at commandline "
                       f"{self.species_assembly} does not match that in the BAM file "
                       f"{assembly}. Defaulting to BAM file value."
                       warnings.warn(as_warn_message)
                  setattr(self, "species_assembly", assembly)

             if "SP" in entry:
                  species = entry["SP"]
                  if getattr(self, "species", None) != species:
                       sp_warn_message = f"Species defined at commandline "
                       f"{self.species} does not match that in the BAM file "
                       f"{species}. Defaulting to BAM file value."
                       warnings.warn(sp_warn_message)
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
        
        split_count = None
        # Step 7. file_line_count: if !process OR process in [mstep, estep]
        if no_process or getattr(self, "process", None) in ["mstep", "estep"]:
            split_count = file_line_count(self.split_list)

        if getattr(self, "limit", None):
            split_count = self.limit

        # Step 8. caveman_mstep: if !process OR process == mstep
        if no_process or getattr(self, "process", None) == "mstep":
            self.caveman_mstep(split_count)

        # Step 9. caveman_merge: if !process OR process == merge
        if no_process or getattr(self, "process", None) == "merge":
            self.caveman_merge()

        # Step 10. caveman_estep: if !process OR process == estep
        if no_process or getattr(self, "process", None) == "estep":
            self.caveman_estep(split_count)
        
        # Step 11. caveman_merge_results: if !process OR process == merge_results
        out_file = f"{self.tumour_bam}_vs_{self.normal_bam}"

        if no_process or getattr(self, "process", None) == "merge_results":
            self.caveman_merge_results(out_file)
        
        # Step 12. caveman_add_vcf_ids: if !process OR process == add_ids
        raw_muts_file = f"{out_file}.{CavemanConstants.RAW_MUTS}"
        ids_muts_file = f"{out_file}.{CavemanConstants.IDS_MUTS}"
        raw_snps_file = f"{out_file}.{CavemanConstants.RAW_SNPS}"
        ids_snps_file = f"{out_file}.{CavemanConstants.IDS_SNPS}"

        if no_process or getattr(self, "process", None) == "add_ids":
            # First do mutations then SNPs
            self.caveman_add_vcf_ids(raw_muts_file, ids_muts_file, "muts")
            self.caveman_add_vcf_ids(raw_snps_file, ids_snps_file, "snps")

        # Step 13. caveman_flag: if !process OR process == flag OR !noflag
        flag_defined_or_main_run = (no_process or getattr(self, "process", None) == "flag")
        no_flag_undefined = getattr(self, "noflag", None)

        if flag_defined_or_main_run or not no_flag_undefined:
            # Use same additions to options from perl wrapper for simplicity
            self.for_split = ids_muts_file
            self.split_out = f"{ids_muts_file}.split"
            # No need for split_lines param as it is in CavemanConstants
            # Split VCF first
            self.caveman_split_vcf()

            self.vcf_split_count = self.count_files("f{self.split_out}.*")
            self.flagged = f"{out_file}.{CavemanConstants.FLAGGED_MUTS}"
            # Concatenate flagged files to single flagged output file
            self.concat_flagged()
            # Gzip and index output flagged file
            self.zip_flagged()

        # Step 14. cleanup: if !noclean
        if not getattr(self, "noclean", None):
            # 'or' is distributive, so combining with previous bools has no
            # effect and makes the implementaiton more readable.
            add_ids_defined = getattr(self, "add_ids", None)
            if flag_defined_or_main_run or (no_flag_defined and add_ids_defined):
                self.pre_cleanup_zip()
                self.cleanup()

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
        if not self.check_exec_in_path():
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
            print(f"Setup stage failed for {final_result['index']}", file=sys.stderr)
            print(f"Error: {final_result['error']}", file=sys.stderr)
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
        if not self.check_exec_in_path():
            raise FileNotFoundError("`caveman` could not be found in $PATH")

        errors_raised = False
        with Pool(processes=num_procs) as pool:
            async_results = []
            for idx, value in enumerate(index_list):
                # If run for current index has been done, continue
                if success_exists(self.progress_dir, idx+1):
                    continue

                command = f"caveman split -i {value} "
                f"-f {self.cave_cfg} "
                f"-e {self.read_count}"

                result = pool.apply_async(worker, args=(self.log_dir, command, idx+1,))
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
            Starting index for limited x step check.

        Returns:
        --------
        `bool` - 
            True if the run is successful, False if errors arise in running
        """
        index_list = None 
        num_procs = self.threads
        if index:
            index_list = self.limited_xstep_indices(index)
            if index != self.index:
                return True
        else:
            index_list = self.valid_fai_idx

        # If we only have one index to consider, only one thread is required.
        if len(index_list) == 1:
            num_procs = 1

        # In case function is being called manually, check caveman
        # is still in the path.
        if not self.check_exec_in_path():
            raise FileNotFoundError("`caveman` could not be found in $PATH")

        errors_raised = False
        with Pool(processes=num_procs) as pool:
            async_results = []
            for idx, value in enumerate(index_list):
                # If run for current index has been done, continue
                if success_exists(self.progress_dir, idx+1):
                    continue

                command = f"caveman mstep -i {value} "
                f"-f {self.cave_cfg}"

                result = pool.apply_async(worker, args=(self.log_dir, command, idx+1,))
                async_results.append(result)

            for item in async_results:

                final_result = item.get()

                if not final_result["success"]:
                    print(f"Split stage failed for {final_result['index']}", file=sys.stderr)
                    print(f"Error: {final_result['error']}", file=sys.stderr)
                    errors_raised = True
                elif not touch_success(self.progress_dir, final_result["index"]):
                    # touch_success is called, so if not successful count as an error
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
        if not self.check_exec_in_path():
            raise FileNotFoundError("`caveman` could not be found in $PATH")

        if success_exists(self.progress_dir, 0):
            return True
        
        command = f"caveman merge "
        f"-c {self.cave_carr} "
        f"-p {self.cave_parr} "
        f"-f {self.cave_cfg}"

        # Only one process is required for setup, set index to 0.
        index = 0
        final_result = worker(self.log_dir,command, index)

        if not final_result["success"]:
            print(f"Merge stage failed for {final_result['index']}", file=sys.stderr)
            print(f"Error: {final_result['error']}", file=sys.stderr)
            return False

        return touch_success(self.progress_dir, final_result["index"])
 
    def caveman_estep(self, index=None):
        """
        Runs CAVEMAN_ESTEP from the parameters used at initialisation,
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
            Starting index for limited x step check.

        Returns:
        --------
        `bool` - 
            True if the run is successful, False if errors arise in running
        """
        index_list = None 
        num_procs = self.threads
        if index:
            index_list = self.limited_xstep_indices(index)
            if index != self.index:
                return True
        else:
            index_list = self.valid_fai_idx

        # If we only have one index to consider, only one thread is required.
        if len(index_list) == 1:
            num_procs = 1

        # In case function is being called manually, check caveman
        # is still in the path.
        if not self.check_exec_in_path():
            raise FileNotFoundError("`caveman` could not be found in $PATH")

        errors_raised = False
        with Pool(processes=num_procs) as pool:
            async_results = []
            for idx, value in enumerate(index_list):
                # If run for current index has been done, continue
                if success_exists(self.progress_dir, idx+1):
                    continue

                command = f"caveman estep -i {value} "
                f"-k {self.normcont_value} "
                f"-g {self.cave_carr} "
                f"-o {self.cave_parr} "
                f"-v {self.species_assembly} "
                f"-w {self.species} "
                f"-f {self.cave_cfg} "
                f"-l {self.normal_protocol} "
                f"-r {self.tumour_protocol} "

                if getattr(self, "norm_cn_default", None):
                    command += f"-n {self.norm_cn_default} "
                
                if getattr(self, "tum_cn_default", None):
                    command += f"-t {self.tum_cn_default} "

                if getattr(self, "priorMut", None):
                    command += f"-c {self.priorMut} "

                if getattr(self, "priorSnp", None):
                    command += f"-d {self.priorSnp} "

                if getattr(self, "nplat", None):
                    command += f"-P {self.nplat} "

                if getattr(self, "tplat", None):
                    command += f"-T {self.tplat} "

                if getattr(self, "mpc", None):
                    command += f"-p {self.mpc} "

                if getattr(self, "spc", None):
                    command += f"-q {self.spc} "

                if getattr(self, "debug_cave", None):
                    command += f"-s "

                result = pool.apply_async(worker, args=(self.log_dir, command, idx+1,))
                async_results.append(result)

            for item in async_results:

                final_result = item.get()

                if not final_result["success"]:
                    print(f"Split stage failed for {final_result['index']}", file=sys.stderr)
                    print(f"Error: {final_result['error']}", file=sys.stderr)
                    errors_raised = True
                elif not touch_success(self.progress_dir, final_result["index"]):
                    # touch_success is called, so if not successful count as an error
                    errors_raised = True
                else:
                    # Explicitly show we continue if touch_success is True
                    continue

        successful_run = not errors_raised

        return successful_run

    def caveman_merge_results(self, out_file=None):
        """
        Runs MERGE_CAVEMAN_RESULTS from the parameters used at initialisation.

        This relies on the Perl script since threading is not important for this
        application - one job is required for each call.

        Parameters:
        -----------
        `out_file` : `str` - 
            Base name of the output files for merged results

        Returns:
        --------
        `bool` - 
            True if the run is successful, False if errors arise in running
        """
        tmp = self.tmp_dir
        # TODO: Implement PCAP::sample_name analogue to redefine Caveman::Implement::prepare
        # if possible, currently just writing to tumour_filename_vs_normal_filename.
        split_list = self.split_list

        if not out_file:
            out_file = f"{self.tumour_bam}_vs_{self.normal_bam}"

        # First do substitutions
        sub_command = f"mergeCavemanResults -s {split_list} -o {out_file}.muts.vcf -f {self.subvcf}"
        if success_exists(f"{self.progress_dir}/merge_muts", 0):
            return True

        # Only one process is required for setup, set index to 0.
        sub_index = 0
        sub_final_result = worker(self.log_dir, sub_command, sub_index)
        sub_success = touch_success(f"{self.progress_dir}/merge_muts", sub_final_result["index"])
        if not sub_success:
            print(f"Merging results failed for mutations stage failed", file=sys.stderr)
            return False

        # Next do SNPs
        snp_command = f"mergeCavemanResults -s {split_list} -o {out_file}.snps.vcf -f {self.snpvcf}"
        if success_exists(f"{self.progress_dir}/merge_snps", 0):
            return True

        # Only one process is required for setup, set index to 0.
        snp_index = 0
        snp_final_result = worker(self.log_dir, snp_command, snp_index)
        snp_success = touch_success(f"{self.progress_dir}/merge_snps", snp_final_result["index"])
        if not snp_success:
            print(f"Merging results failed for SNP stage failed", file=sys.stderr)
            return False

        # Next do no analysis region.
        if success_exists(f"{self.progress_dir}/merge_no_analysis", 0):
            return True
        else:
            no_analysis_command = f"mergeCavemanResults -s {split_list} -o {out_file}.no_analysis.bed -f {self.noanalysisbed}"

            # Only one process is required for setup, set index to 0.
            no_analysis_index = 0
            no_analysis_final_result = worker(self.log_dir, no_analysis_command, no_analysis_index)
            # Extend no analysis region
            self.extend_no_analysis(f"{out_file}.no_analysis.bed")
            no_analysis_success = touch_success(f"{self.progress_dir}/merge_no_analysis", no_analysis_final_result["index"]) 
            if not no_analysis_success:
                print(f"Merging results failed for no analysis stage failed", file=sys.stderr)
                return False

        return True

    def caveman_add_vcf_ids(self, raw_file, ids_file, snps_or_muts):
        """
        Runs CAVEMAN_VCF_IDS from the parameters used at initialisation

        Parameters:
        -----------
        `raw_file` : `str` - 
            Raw VCF file without IDS

        `ids_file` : str - 
            IDs of samples for raw file

        `snps_or_muts` : str - 
            Whether SNPs or mutation signatures are being considered

        Returns:
        --------
        `bool` - 
            True if the run is successful, False if errors arise in running
        """
        if success_exists(f"{self.progress_dir}/{snps_or_muts}"):
            return True

        # Raise an error if the perl executable is not in the path
        perl_path = "perl"
        if not self.check_exec_in_path(perl_path):
            raise FileNotFoundError(f"`{perl_path}` could not be found in $PATH")

        # Raise an error if the VCF IDS executable is not in the path
        executable = "cgpAppendIdsToVcf.pl"
        if not self.check_exec_in_path(executable):
            raise FileNotFoundError(f"`{executable}` could not be found in $PATH")

        command = f"perl {executable} "
        f"-i {raw_file} "
        f"-o {ids_file}"

        final_result = worker(self.log_dir, command, 0)
        if not touch_success(f"{self.progress_dir}/{snps_or_muts}"):
            print(f"caveman_add_vcf_ids (calling {executable}) failed", file=sys.stderr)
            print(f"Error: {final_result['error']}", file=sys.stderr)
            return False

        return True

    def caveman_flag(self, index):
        """
        Runs CAVEMAN_FLAG from the parameters used at initialisation

        Parameters:
        -----------
        `index` : `int` - 
            The required index for file flagging
        """
        index_list = None 
        num_procs = self.threads
        if index:
            if index != self.index:
                return True
        else:
            index_list = self.limited_flag_indices(index)

        # If we only have one index to consider, only one thread is required.
        if len(index_list) == 1:
            num_procs = 1

        # Raise an error if the perl executable is not in the path
        perl_path = "perl"
        if not self.check_exec_in_path(perl_path):
            raise FileNotFoundError(f"`{perl_path}` could not be found in $PATH")

        # Raise an error if the flagging executable is not in the path
        executable = "cgpFlagCaVEMan.pl"

        if not self.check_exec_in_path(executable):
            raise FileNotFoundError(f"`{executable}` could not be found in $PATH")

        errors_raised = False
        with Pool(processes=num_procs) as pool:
            async_results = []
            for idx, value in enumerate(index_list):
                # If run for current index has been done, continue
                if success_exists(self.progress_dir, idx+1):
                    continue

                command = f"perl {executable} "
                f"-i {self.split_out}.{value} "
                f"-o {self.flagged}.{value} "
                f"-s {self.species} "
                f"-m {self.tumour_bam} "
                f"-n {self.normal_bam} "
                f"-b {self.flag_bed} "
                f"-g {self.germindel} "
                f"-umv {self.unmatched_vcf} "
                f"-ref {self.reference} "
                f"-t {self.seqType} "
                f"-sa {self.species_assembly} "

                if getattr(self, "flagConfig", None):
                    command += f"-c {self.flagConfig} "

                if getattr(self, "flagToVcfConfig", None):
                    command += f"-v {self.flagToVcfConfig} "

                if getattr(self, "apid", None):
                    command += f"-p {self.apid} "

                if getattr(self, "annot_bed", None):
                    command += f"-ab {self.annot_bed} "

                result = pool.apply_async(worker, args=(self.log_dir, command, idx+1,))
                async_results.append(result)

            for item in async_results:

                final_result = item.get()

                if not final_result["success"]:
                    print(f"Flag stage failed for {final_result['index']}", file=sys.stderr)
                    print(f"Error: {final_result['error']}", file=sys.stderr)
                    errors_raised = True
                elif not touch_success(self.progress_dir, final_result["index"]):
                    # touch_success is called, so if not successful count as an error
                    errors_raised = True
                else:
                    # Explicitly show we continue if touch_success is True
                    continue

        successful_run = not errors_raised

        return successful_run
 
    def caveman_split_vcf(self):
        """
        Runs CAVEMAN_VCF_SPLIT from the parameters used at initialisation

        Returns:
        -------
        `bool` - 
            True if run was successful, False otherwise.
        """
        if success_exists(f"{self.progress_dir}", 0):
            return True

        # Raise an error if the perl executable is not in the path
        perl_path = "perl"
        if not self.check_exec_in_path(perl_path):
            raise FileNotFoundError(f"`{perl_path}` could not be found in $PATH")

        # Raise an error if the VCF split executable is not in the path
        executable = "cgpVCFSplit.pl"
        if not self.check_exec_in_path(executable):
            raise FileNotFoundError(f"`{executable}` could not be found in $PATH")

        command = f"perl {executable} "
        f"-i {self.for_split} "
        f"-o {self.split_out} "
        f"-s "
        f"-l {CavemanConstants.SPLIT_LINE_COUNT}"

        final_result = worker(self.log_dir, command, 0)
        if not touch_success(f"{self.progress_dir}", 0):
            print(f"caveman_split_vcf (calling {executable}) failed", file=sys.stderr)
            print(f"Error: {final_result['error']}", file=sys.stderr)
            return False

        return True

    def count_files(self, match_pattern):
        """
        Runs FILE_COUNT from the parameters used at initialisation

        Parameters:
        -----------
        `match_pattern` : `str` - 
            Filename pattern to match for counting

        Returns:
        --------
        `count` : `int` - 
            Number of files matching the pattern `match_pattern`
        """
        try:
            # Use glob for safety and simplicity
            import glob
            files = glob.glob(match_pattern)
            return len(files)
        except Exception as e:
            # Mimic the die clause in the perl wrapper
            raise RuntimeError(f"ERROR: ({e}) Encountered counting split files for flagging. Searching {match_pattern}")

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
            print(f"Concat failed for command: `{command}`", file=sys.stderr)
            print(f"Error: {final_result['error']}", file=sys.stderr)
            return False

        return touch_success(self.progress_dir, final_result["index"])

    def concat_flagged(self):
        """
        Runs CAVEMAN_VCF_FLAGGED_CONCAT from the parameters used at initialisation

        Returns:
        ---------
        `bool` - 
            True if the run was successful, False otherwise.
        """
        if success_exists(self.progress_dir, 0):
            return True
        
        # Concatenate all {self.split_list}.* files to {self.split_list}
        concat_command = f"vcf-concat {self.flagged}.*"
        sort_command = f"vcf-sort"

        try:
            # Safer way of running piped command: concat | sort
            concat_run = subprocess.run(concat_command.split(" "), stdout=subprocess.PIPE)
            sort_run = subprocess.run(sort_command, input=concat_run.stdout, stdout=subprocess.PIPE)
            with open(self.flagged, "w") as flagged_out:
                print(sort_run.stdout.decode(), file=flagged_out)

        except Exception as e:
            print(f"concat_flagged stage failed", file=sys.stderr)
            print(f"Error: {e}", file=sys.stderr)
            return False

        return touch_success(self.progress_dir, 0)

    def zip_flagged(self):
        """
        Zips flagged files from caveman_flag stage

        Returns:
        --------
        `bool` - 
            True if run is successful, False otherwise.
        """
        # IMPORTANT: THERE IS NO CHECK AGAINST `success_exists` IN THE PERL
        # WRAPPER, TO MAINTAIN COMPATIBILITY I HAVE NOT IMPLEMENTED ONE HERE.
        vcf_gz = f"{self.flagged}.gz"

        # Raise an error if bgzip executable is not in the path
        bgzip = "bgzip"
        if not self.check_exec_in_path(bgzip):
            raise FileNotFoundError(f"`bgzip` could not be found in $PATH")

        # Raise an error if the tabix executable is not in the path
        tabix = "tabix"
        if not self.check_exec_in_path(tabix):
            raise FileNotFoundError(f"`tabix` could not be found in $PATH")

        bgzip_command = f"{bgzip} -c {self.flagged} > {vcf_gz}"
        tabix_command = f"{tabix} -p vcf {vcf_gz}"
        commands = [bgzip_command, tabix_command]

        for command in commands:
            index = 0
            result = worker(self.log_dir, command, index)

            if not result["success"]:
                print(f"Zip flagging stage failed for command: `{command}`", file=sys.stderr)
                print(f"Error: {result['error']}", file=sys.stderr)
                return False

        return touch_success(self.progress_dir, 0)

    def pre_cleanup_zip(self):
        """
        Zip files before the cleanup stage if cleanup option is specified
        """
        if success_exists(self.progress_dir, 0):
            return True

        # Raise an error if bgzip executable is not in the path
        bgzip = "bgzip"
        if not self.check_exec_in_path(perl_path):
            raise FileNotFoundError(f"`bgzip` could not be found in $PATH")

        # Raise an error if the tabix executable is not in the path
        tabix = "tabix"
        if not self.check_exec_in_path(tabix):
            raise FileNotFoundError(f"`tabix` could not be found in $PATH")

        vcf_muts_gz = f"{self.ids_muts_file}.gz"
        bgzip_muts_command = f"{bgzip} -c {self.ids_muts_file} > {vcf_muts_gz}"
        tabix_muts_command = f"{tabix} -p vcf {vcf_muts_gz}"

        vcf_snps_gz = f"{self.ids_snps_file}.gz"
        bgzip_snps_command = f"{bgzip} -c {self.ids_snps_file} > {vcf_snps_gz}"
        tabix_snps_command = f"{tabix} -p vcf {vcf_snps_gz}"

        commands = [bgzip_muts_command, tabix_muts_command, bgzip_snps_command, tabix_snps_command]

        for command in commands:
            index = 0
            result = worker(self.log_dir, command, index)

            if not result["success"]:
                print(f"Zip flagging stage failed for command: `{command}`", file=sys.stderr)
                print(f"Error: {result['error']}", file=sys.stderr)
                return False

        return touch_success(self.progress_dir, 0)

    def cleanup(self):
        """
        Runs cleanup method after all processes have concluded.
        """
        return

    def limited_indices(self, index, count):
        """
        Returns a list of indices starting from `index`,
        incrementing by `self.limit` if it exists, up to and
        including `count`, but not exceeding it.

        Parameters:
        -----------
        `index` : `int` - 
            Starting index to be checked from.

        `count` : `int` - 
            The number of indices to check.

        Returns:
        --------
        `indices` : `list` - 
            List of indices with limit imposed.
        """
        indices = []
        if getattr(self, "limit", None) is None:
            indices.append(index)
            return indices

        if index < 1 or (getattr(self, "limit", None) and self.limit <= 0):
            return indices

        base = index
        while base <= count:
            indices.append(base)
            base += self.limit

        return indices

    def limited_xstep_indices(self, index):
        """
        Check limited indices for the split list count.

        Parameters:
        -----------
        `index` : `int` - 
            Starting index to be checked from.

        Returns:
        --------
        `indices` : `list` - 
            List of limited indices by xstep.
        """
        split_count = file_line_count(self.split_list)
        indices = self.limited_indices(index, split_count)
        return indices

    def limited_flag_indices(self, index):
        """
        Check limited indices for the VCF split count

        Parameters:
        -----------
        `index` : `int` - 
            Starting index to be checked from.

        Returns:
        --------
        `indices` : `list` - 
            List of limited indices by flag.
        """
        indices = self.limited_indices(index, self.vcf_split_count)

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

    def extend_no_analysis(self, no_analysis_file):
        """
        Extends no analysis in caveman_merge_results

        Parameters:
        -----------
        `no_analysis_file` : `str` - 
            Filename to use to extend no analysis region.
        """
        exclude_patterns = self.load_exclude()

        if not exclude_patterns:
            return

        with open(no_analysis_file, "w") as output_file:
            with open(self.reference, "r") as input_reference:
                for line in input_reference:
                    seq, length = line.strip().split('\t')[0:2]
                    print(f"{seq}\t0\t{length}", file=output_file)

