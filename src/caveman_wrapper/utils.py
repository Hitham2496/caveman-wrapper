"""
CaVEMan wrapper to replace Perl script

utilities classes for caveman params and consts
"""
import argparse
from dataclasses import dataclass


def file_line_count(filename, contig_count = 0):
    """
    Increments a contig count based on the number of lines in a seuqence file

    Parameters:
    ----------
    `filename` : `str` - 
        File to count the lines of, no default

    `contig_count` : `int` - 
        Count of contigous regions, defaults to zero

    Returns:
    -------
    `contig_count` : `int` -
        Updated `contig_count` variable
    """

    if not isinstance(contig_count, int):
        raise ValueError("contig_count needs to be provided as an integer")

    with open(filename) as filestream:
        for line in filestream:
            contig_count++

    return contig_count


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
       "help" : "h"
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

