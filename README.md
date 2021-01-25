# `Readcomb` - fast detection of recombinant reads in BAMs

## What is it?

Readcomb is a hybrid command line script and Python module for the filtering and classification of bam sequences based on their phase change properties. 

## Installation

``` Python
pip install readcomb
```

## Dependencies
- [cyvcf2](https://brentp.github.io/cyvcf2/) - Fast retrieval and filtering of vcf files and vcf objects written in C
- [pysam](https://pysam.readthedocs.io/en/latest/api.html) - Interface for SAM and BAM files and provides SAM and BAM objects
- [pandas](https://pandas.pydata.org/) - Support for data tables
- [tqdm](https://tqdm.github.io/) - Provides updating progress bars for command line programs
- [samtools](http://www.htslib.org/)

# Usage:

## bamprep 
Command line preprocessing script for bam files
```bash

```

Optional parameters:

- 

## vcfprep
Command line preprocessing script for vcf files
```bash
readcomb-vcfprep --vcf [vcf_filepath] --out [output_filepath]
```

Optional arguments
- `--snps_only`, Keep only SNPs
- `--indels_only`, Keep only indels
- `--no_hets`, Remove heterozygote calls
- `--min_GQ [quality]`, Minimum genotype quality at both sites (default is 30)

## filter
Command line multiprocessing script for identification of bam sequences with phase changes
```bash
readcomb-filter --bam [bam_filepath] --vcf [vcf_filepath]
```

Optional arguments:
- `-p, --processes [processes]`, Number of processes available for filter (default is 4)
- `-m, --mode [phase_change|no_match]`, Filtering mode (default phase_change)
- `-l, --log [log_filepath]`, Filename for log metric output
- `-o, --out [output_filepath]`, File to write filtered output to (default recomb_diagnosis)

## classification
Python module for detailed classification of sequences containing phase changes
```python
from readcomb.classification import rc

# generate list of bam read pairs
pairs = rc.pairs_creation(bam_filepath, vcf_filepath)

# call each of the pairs to analyse and classify them
# map and lambda function
map(lambda x:x.call(), pairs)

# or use a for loop
for pair in pairs:
    pair.call() 

# get classification of first read pair
pairs[0].classify
# > gene_conversion
```
## License

GNU General Public License v3 (GPLv3+)

## Development

Currently in alpha

[Source code](https://github.com/ness-lab/readcomb)

[Development repo](https://github.com/ness-lab/recombinant-reads)