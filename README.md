#  CaVEMan Wrapper

A Python wrapper for the [CaVEMan (Cancer Variants Through Expectation Maximization)](https://github.com/cancerit/CaVEMan) variant caller, designed for efficient execution on HPC clusters using multiprocessing where Perl threading was employed by the previous wrapper to avoid underutilising platform resources (i.e. cores provided by LSF/other schedulers).

## Overview

This wrapper script replaces the original Perl implementation with multiprocessed Python version to distribute work across HPC resources without threading. Full compatibility with short flag names is ensured with the original.
- Setting up and managing CaVEMan configuration
- Splitting jobs for parallel processing
- Managing variant calling workflow stages
- Handling results merging and flagging
- Providing cleanup and resource management

## Prerequisites
- Python 3.5+
- CaVEMan and its dependencies installed and available in PATH (e.g. in a singularity container)
- Required non-standard Python packages:
  - pysam

The following executables must be in your PATH:
- `caveman`
- `bgzip`
- `tabix`
- `perl` (for some legacy components)
- `mergeCavemanResults.pl`
- `cgpAppendIdsToVcf.pl`
- `cgpFlagCaVEMan.pl`
- `cgpVCFSplit.pl`

## Installation

```bash
git clone https://github.com/Hitham2496/caveman-wrapper.git
cd caveman-wrapper
pip install -e .
```

## Usage

Create a new entry point script `caveman.py`:

```python
#!/usr/bin/env python3

from caveman_wrapper.core import *
from caveman_wrapper.utils import *

def main():
    # Parse command line arguments
    parser = CavemanFlags.parser()
    args = parser.parse_args()
    
    # Initialise and run CaVEMan
    runner = CavemanRunner(**vars(args))
    runner.run_caveman()

if __name__ == "__main__":
    main()
```

### Basic Usage

```bash
python caveman.py \
    --reference /path/to/reference.fa \
    --outdir /path/to/output \
    --tumour-bam /path/to/tumor.bam \
    --normal-bam /path/to/normal.bam \
    --ignore-file /path/to/ignore.bed
```

### Advanced Usage

```bash
python caveman.py \
    --reference /path/to/reference.fa \
    --outdir /path/to/output \
    --tumour-bam /path/to/tumor.bam \
    --normal-bam /path/to/normal.bam \
    --ignore-file /path/to/ignore.bed \
    --threads 16 \
    --normal-contamination 0.1 \
    --species "human" \
    --species-assembly "GRCh38" \
    --seqType "genome"
```

## Key Parameters

| Parameter | Short Flag | Description |
|-----------|------------|-------------|
| --reference | -r | Reference genome FASTA file |
| --outdir | -o | Output directory |
| --tumour-bam | -tb | Tumor sample BAM file |
| --normal-bam | -nb | Normal sample BAM file |
| --ignore-file | -ig | BED file of regions to ignore |
| --threads | -t | Number of threads to use |
| --normal-contamination | -k | Normal contamination value (0-1) |
| --species | -s | Species (e.g., "human") |
| --species-assembly | -sa | Assembly version (e.g., "GRCh38") |
| --seqType | -st | Sequencing type (genome/exome/targeted) |

For a complete list of parameters:
```bash
python caveman.py --help
```
which will call the native argument parsing method to produce the complete list of commands


## Workflow Stages

The wrapper manages the following CaVEMan stages:
1. Setup - Configure CaVEMan parameters
2. Split - Divide samples into chunks to be processed subsequently
3. Split concatenation - Concatenate split lists into one file - **uses multiple process as given by `threads`**
4. M-Step - Maximise parameters - **uses multiple process as given by `threads`**
5. Merge - Combine M-step results
6. E-Step - Expectation step for variant calling - **uses multiple process as given by `threads`**
7. Results Merge - Combine all results
8. Add IDs - Add variant IDs
9. Flagging - Apply filters and flags
10. Cleanup - Clean the files in the temporary directory

## Output Files

The wrapper generates the following key output files in the specified output directory:

```
outdir/
├── logs/                    # Log files for each process
├── tmpCaveman/             # Temporary processing files
│   ├── results/            # Intermediate results
│   └── progress/           # Progress tracking files
├── {tumor}_vs_{normal}.muts.vcf.gz      # Final mutations VCF
├── {tumor}_vs_{normal}.snps.vcf.gz      # Final SNPs VCF
└── {tumor}_vs_{normal}.no_analysis.bed  # Regions without analysis
```
