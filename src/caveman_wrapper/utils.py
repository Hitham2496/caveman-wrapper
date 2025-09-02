"""
CaVEMan wrapper to replace Perl script

utilities classes for caveman params and consts
"""
import argparse
import os
import sys
from dataclasses import dataclass
import gzip
import shutil
import re


def gunzip_file(filename_in, filename_out):
    """
    Void method to unzip a '.gz.*' file

    Parameters:
    ----------
    `filename_in` : `str` - 
        Name of gzipped file to be unzipped

    `filename_out` : `str` - 
        Name of file to receive output of `gunzip`
    """
    if not re.search(".*\.gz[^a-z,A-Z,0-9].*$", filename_in):
        raise ValueError(f"File {filename_in} not compressed in gzip format")

    with gzip.open(filename_in, "rb") as f_in:
        with open(filename_out, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

def file_line_count(filename):
    """
    Returns a contig count based on the number of lines in a seuqence file

    Parameters:
    ----------
    `filename` : `str` - 
        File to count the lines of, no default

    Returns:
    -------
    `contig_count` : `int` -
        Lines in the file provided, equal to the number of contigs
    """
    contig_count = 0

    with open(filename) as filestream:
        for line in filestream:
            contig_count += 1

    return contig_count

def valid_index_by_factor(opt_name, opt_val, base, proc_factor=1):
    """
    Simplified void implementation of PCAP::Cli::valid_index_by_factor
    which raises an error if the option value is not in the required
    region.

    Parameters
    ----------
    `opt_name` - `str` - 
        Option name for the error message

    `opt_val` - `str` - 
        Value of the option provided

    `base` : `int` - 
        Base for maximum value

    `proc_factor` : `int` - 
        Multiplicative factor for the maximum calculation
    """
    max_val = base * (proc_factor if proc_factor else 1)
    if not (1 <= opt_val <= max_val):
        raise ValueError(f"Option '{opt_name}' needs to be between 1 and {max_val}")

def check_outdir(directory_name):
    """
    Simplified void implementation of PCAP::Cli::out_dir_check,
    which does not ask to overwrite an existing directory.

    Parameters
    ----------
    `directory_name` : `str` - 
        Name of the directory to use as the output dir
    """
    if directory_name is None:
        raise ValueError("Output directory name must be provided")

    if os.path.isdir(directory_name):
        if len(os.listdir(directory_name)) != 0:
            raise OSError("Directory {directory_name} exists and is non-empty")

        try:
            filestream = open(f"{directory_name}/test_file.txt", 'w')
        except IOError:
            print(f"Directory {directory_name} is write-protected", file=sys.stderr)
            sys.exit(1)

        print(f"Directory {directory_name} exists, output files will be written here", file=sys.stderr)

    elif os.path.isfile(directory_name) or os.path.islink(directory_name):
        raise OSError(f"Non-directory object already exists at location f{directory_name}")

    else:
        print(f"Directory {directory_name} does not exist, creating now", file=sys.stderr)
        try:
            os.mkdir(directory_name)
        except Exception as e:
            print(e, file=sys.stderr)
            print(f"Error creating output directory {directory_name}", file=sys.stderr)

def get_marker_filename(tmp, *indices):
    """
    Get the filename for PCAP::Threaded progress tracking utility
    functions `success_exists` and `touch_success`
    """
    frame = inspect.stack()[1]
    caller = f"{frame.frame.f_globals['__name__']}_{frame.function}".replace('.', '_')
    suffix = '.'.join(str(i) for i in indices)
    return Path(tmp) / f"{caller}.{suffix}"

def success_exists(tmp, *indices):
    """
    Utility function to replicate progress checking of
    PCAP::Threaded::success_exists.
    """
    marker = get_marker_filename(tmp, *indices)
    if marker.exists():
        print(f"Skipping {marker.name} as previously successful")
        return True
    return False

def touch_success(tmp, *indices):
    """
    Utility function to replicate progress logging of
    PCAP::Threaded::touch_success.
    """
    marker = get_marker_filename(tmp, *indices)
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.touch()
    return True

# Class for constants, the valid params for caveman running
class CavemanConstants():
    """
    Bookkeeping for constants to provide a slightly more robust guard against mutability
    """

    # Checks for valid parameters
    VALID_PROCESSES = ["setup", "split", "split_concat", "mstep", "merge", "estep", "merge_results", "add_ids", "flag"]
    VALID_PROTOCOLS = "WGS", "WXS", "RNA", "AMPLICON", "TARGETED", "RNA-Seq"
    VALID_SEQ_TYPES = ["pulldown", "exome", "genome", "genomic", "followup", "targeted", "rna_seq"]

    # Configuration settings
    CAVEMAN_CONFIG = "caveman.cfg.ini"
    CAVEMAN_ALG_BEAN = "alg_bean"
    CAVEMAN_PROB_ARR = "prob_arr"
    CAVEMAN_COV_ARR = "cov_arr"
    DEFAULT_PROTOCOL = "WGS"
    DEFAULT_NORMCONT = 0.1
    SPLIT_LINE_COUNT = 25000
    SPLIT_STEP_READ_COUNT = 500000

    # Suffixes
    RAW_MUTS = ".muts.vcf"
    IDS_MUTS = ".muts.ids.vcf"
    FLAGGED_MUTS = ".flagged.muts.vcf"
    FLAGGED_MUTS_GZ = ".flagged.muts.vcf.gz"
    FLAGGED_MUTS_TBI = ".flagged.muts.vcf.gz.tbi"
    RAW_SNPS = ".snps.vcf"
    IDS_SNPS = ".snps.ids.vcf"
    IDS_SNPS_GZ = ".snps.ids.vcf.gz"
    IDS_SNPS_TBI = ".snps.ids.vcf.gz.tbi"
    IDS_MUTS_GZ = ".muts.ids.vcf.gz"
    IDS_MUTS_TBI = ".muts.ids.vcf.gz.tbi"
    NO_ANALYSIS = ".no_analysis.bed"
    # TODO: Figure out where this should live later
    #SP_ASS_MESSAGE = qq{%s defined at commandline (%s) does not match that in the BAM file (%s). Defaulting to BAM file value.\n}


@dataclass
class CavemanFlags():
    """
    Caveman runner flags to generate from arguments
    provided at runtime, naming maximally reflective
    of Perl script naming convention
    """
    man: bool = False
    version: bool = False
    reference: str = None
    outdir: str = None
    tumour_bam: str = None
    normal_bam: str = None
    ignore_file: str = None
    tumour_cn: str = None
    normal_cn: str = None
    threads: int = None
    normal_contamination: str = None
    species: str = None
    species_assembly: str = None
    process: str = None
    logs: str = None
    index: int = None
    limit: int = None
    flag_bed_files: str = None
    annot_bed: str = None
    germindel: str = None
    unmatched_vcf: str = None
    normal_protocol: str = None
    tumour_protocol: str = None
    tum_cn_default: int = None
    norm_cn_default: int = None
    flagConfig: str = None
    flagToVcfConfig: str = None
    priorMut: float = None
    priorSnp: float = None
    apid: int = None
    tplat: str = None
    nplat: str = None
    seqType: str = None
    noflag: bool = False
    noclean: bool = False
    mut_probability_cutoff: float = None
    snp_probability_cutoff: float = None
    read_count: int = None
    exclude: str = None
    debug_cave: bool = False

    # Dictionary of shortened flags for the argument parser
    short_flags = {
       "help" : "h",
       "man" : "m",
       "version" : "v",
       "reference" : "r",
       "outdir" : "o",
       "tumour_bam" : "tb",
       "normal_bam" : "nb",
       "ignore_file" : "ig",
       "tumour_cn" : "tc",
       "normal_cn" : "nc",
       "threads" : "t",
       "normal_contamination" : "k",
       "species" : "s",
       "species_assembly" : "sa",
       "process" : "p",
       "logs" : "g",
       "index" : "i",
       "limit" : "l",
       "flag_bed_files" : "b",
       "annot_bed" : "ab",
       "germindel" : "in",
       "unmatched_vcf" : "u",
       "normal_protocol" : "np",
       "tumour_protocol" : "tp",
       "tum_cn_default" : "td",
       "norm_cn_default" : "nd",
       "flagConfig" : "c",
       "flagToVcfConfig" : "f",
       "priorMut" : "pm",
       "priorSnp" : "ps",
       "apid" : "a",
       "tplat" : "TP",
       "nplat" : "NP",
       "seqType" : "st",
       "noflag" : "noflag",
       "noclean" : "noclean",
       "mut_probability_cutoff" : "mpc",
       "snp_probability_cutoff" : "spc",
       "read_count" : "e",
       "exclude" : "x",
       "debug_cave" : "dbg"
    }
    
    @classmethod
    def parser(cls):
        """
        An argument parser to prepare flags for the CavemanFlags class

        Returns:
        --------
        `parser` : argpaerse.ArgumentParser - 
            Argument parsing object based on the valid arguments of caveman
        """
        parser = argparse.ArgumentParser(description="Usage: caveman.py [kwargs]")
        for key in cls.__annotations__.keys():
            value = cls.__annotations__[key]
            long_flag = f"--{key.replace('_', '-')}"
            short_flag = f"-{cls.short_flags[key]}"

            # Help is added automatically
            if key == "help":
                continue

            if value == bool:
                parser.add_argument(long_flag, action=argparse.BooleanOptionalAction)
                continue

            parser.add_argument(long_flag, short_flag, type=value)
        return parser

