"""
CaVEMan wrapper to replace Perl script

This script uses multiprocessing to ensure that CPU resource
allocated by HPC schedulers are adequately used
"""
import os
import subprocess
import argparse
import time
import multiprocessing

# Class for constants, the valid params for caveman
class CavemanParameters():

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
