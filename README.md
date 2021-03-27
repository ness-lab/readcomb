# `readcomb` - fast detection of recombinant reads in BAMs

[![PyPI version](https://badge.fury.io/py/readcomb.svg)](https://badge.fury.io/py/readcomb)

`readcomb` is a collection of command line and Python tools for fast detection
of recombination events in pooled high-throughput sequencing data. `readcomb`
searches for changes in parental haplotype phase across individual reads and classifies
recombination events based on various properties of the observed recombinant haplotypes.

`readcomb` was designed for use with the model alga _Chlamydomonas reinhardtii_ and
currently only supports haploids. Although the means of specifically detecting gene
conversion are more specific to _C. reinhardtii_, everything else in `readcomb` is
generalizable to the detection of recombination events in any haploid species. 

## Installation

```bash
pip install readcomb
```

## Dependencies
- [cyvcf2](https://brentp.github.io/cyvcf2/) - Fast retrieval and filtering of VCF files and VCF objects written in C
- [pysam](https://pysam.readthedocs.io/en/latest/api.html) - Interface for SAM and BAM files and provides SAM and BAM objects
- [pandas](https://pandas.pydata.org/) - Support for data tables
- [tqdm](https://tqdm.github.io/) - Provides updating progress bars for command line programs
- [samtools](http://www.htslib.org/) - Used for preprocessing of VCF and BAM files

# Usage:

## bamprep 

Command line preprocessing script for BAM files. `bamprep` will prepare an
index file, filter out unusuable reads, and output a BAM _sorted by read name_.
`readcomb` requires BAMs sorted by read name for fast parsing and filtering.


```bash
readcomb-bamprep --bam [bam_filepath] --out [outdir]
```

Optional parameters:

- `--samtools` - Path to samtools binary
- `--threads [int]` - Number of threads samtools should use (default 1)
- `--index_csi` - Create CSI index instead of BAI
- `--no_progress` - Disable index creation - this will speed up `bamprep` but
  will mean no progress bars when filtering

## vcfprep

Command line preprocessing script for VCF files

```bash
readcomb-vcfprep --vcf [vcf_filepath] --out [output_filepath]
```

Optional arguments
- `--snps_only` - Keep only SNPs
- `--indels_only` - Keep only indels
- `--no_hets` - Remove heterozygote calls
- `--min_GQ [int]` - Minimum genotype quality at both sites (default 30)

## filter

Command line multiprocessing script for identification of bam sequences with phase changes

```bash
readcomb-filter --bam [bam_filepath] --vcf [vcf_filepath]
```

Optional arguments:
- `-p, --processes [processes]`, Number of processes available for filter (default 4)
- `-m, --mode [phase_change|no_match]`, Filtering mode (default `phase_change`)
- `-l, --log [log_filepath]`, Filename for log metric output
- `-o, --out [output_filepath]`, File to write filtered output to (default `recomb_diagnosis`)

## classification

Python module for detailed classification of sequences containing phase changes

```python
>>> import readcomb.classification as rc
>>> from cyvcf2 import VCF

>>> bam_filepath = 'data/example_sequences.bam'
>>> vcf_filepath = 'data/example_variants.vcf.gz'
>>> pairs = rc.pairs_creation(bam_filepath, vcf_filepath)     # generate list of Pair objects
>>> cyvcf_object = VCF(vcf_filepath)                          # cyvcf2 file object

>>> print(pairs[0])
Record name: chromosome_1-199370 
Read1: chromosome_1:499417-499667 
Read2: chromosome_1:499766-500016 
VCF: data/example_variants.vcf.gz

>>> pairs[0].classify(cyvcf_object)                           # run classification algorithm
>>> print(pairs[0])
Record name: chromosome_1-199370 
Read1: chromosome_1:499417-499667 
Read2: chromosome_1:499766-500016 
VCF: data/example_variants.vcf.gz
Unmatched Variant(s): False 
Condensed: [['CC2936', 499417, 499626], ['CC2935', 499626, 499736], ['CC2936', 499736, 500016]] 
Call: gene_conversion 
Condensed Masked: [['CC2936', 499487, 499626], ['CC2935', 499626, 499736], ['CC2936', 499736, 499946]] 
Call Masked: gene_conversion 
```
## License

GNU General Public License v3 (GPLv3+)

## Development

Currently in alpha

[Source code](https://github.com/ness-lab/readcomb)

[Development repo](https://github.com/ness-lab/recombinant-reads)
