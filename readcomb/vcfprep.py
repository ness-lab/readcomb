#!/usr/bin/env python3
"""
vcfprep.py - prepare parental VCF for phase change detection and apply further
filters if needed
"""

import subprocess
import argparse
import sys
from cyvcf2 import VCF
from cyvcf2 import Writer
from tqdm import tqdm

class ReadcombParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write(f'error: {message}\n')
        self.print_help()
        sys.exit(1)

def arg_parser():
    parser = ReadcombParser(
        description='prep parental VCF for phase change detection',
        usage='readcomb-vcfprep [options]')

    parser.add_argument('-v', '--vcf', required=True,
        type=str, help='File to filter (.vcf.gz)')
    parser.add_argument('--snps_only', required=False,
        action='store_true', help='Keep only SNPs [optional]')
    parser.add_argument('--indels_only', required=False,
        action='store_true', help='Keep only indels [optional]')
    parser.add_argument('--no_hets', required=False,
        action='store_true', help='Remove heterozygote calls [optional]')
    parser.add_argument('--min_GQ', required=False, default=30,
        type=int, help='Min GQ at both sites (default 30)')
    parser.add_argument('-o', '--out', required=True,
        type=str, help='File to write to. If .gz, script will bgzip and tabix file.')
    parser.add_argument('--version', action='version', version='readcomb 0.1.4')

    return parser

def vcfprep(args):
    """
    Iterate through parental VCF and remove records by given filters.

    All records where parent1 allele is the same as parent2
    ('uninformative sites') are automatically removed.

    Raises ValueError if >2 samples in VCF (should just be 2 parents)

    Parameters
    ----------
    args : Namespace
        Namespace containing all user given arguments compiled by arg_parse()

    Returns
    -------
    total_count : int
        total records considered
    kept_count : int
        number of records that passed filters
    """
    vcf_in = VCF(args.vcf)

    if len(vcf_in.samples) > 2:
        raise ValueError('more than 2 parental samples in input VCF')
    if args.snps_only and args.indels_only:
        raise ValueError('both --snps_only and --indels_only provided. pick one or neither!')

    if args.out.endswith('.gz'):
        outfile = args.out.replace('.gz', '')
    vcf_out = Writer(outfile, vcf_in)
    vcf_out.write_header()

    total_count, kept_count = 0, 0

    for record in tqdm(vcf_in):
        total_count += 1

        if not len(record.ALT) > 0:
            continue

        # SNP filter
        if args.snps_only and not record.is_snp:
            continue

        # indel filter
        if args.indels_only and not record.is_indel:
            continue

        # heterozygote call filter
        if args.no_hets and record.num_het != 0:
            continue

        # genotype quality filter
        if not all(record.gt_quals >= args.min_GQ):
            continue

        # ensure parental alleles differ
        if record.gt_bases[0] == record.gt_bases[1]:
            continue

        # remove any calls with deleted alleles
        if '*' in ' '.join(record.gt_bases):
            continue

        # remove uninformative indels
        # e.g. A/ATC - either both match, or we have a no match
        if not args.snps_only:
            if record.is_indel:
                allele_1 = record.gt_bases[0].split('/')[0].split('|')[0]
                allele_2 = record.gt_bases[1].split('/')[0].split('|')[0]
                if allele_1 in allele_2 or allele_2 in allele_1:
                    continue

        # only passes if record not caught in above filters
        vcf_out.write_record(record)
        kept_count += 1

    return total_count, kept_count

def main():
    parser = arg_parser()
    args = parser.parse_args()

    print(f'[readcomb] Filtering {args.vcf}')
    if args.min_GQ < 30:
        print('[readcomb] WARNING: min GQ below 30 selected')
    total_count, kept_count = vcfprep(args)
    print('[readcomb] Complete.')
    print(f'[readcomb] {kept_count} of {total_count} records retained.')
    if args.out.endswith('.gz'):
        print('[readcomb] compressing outfile...')
        proc_bgzip = subprocess.Popen(['bgzip', args.out.replace('.gz', '')])
        stdout, stderr = proc_bgzip.communicate()
        print('[readcomb] creating tabix index...')
        proc_tabix = subprocess.Popen(['tabix', args.out])
        stdout, stderr = proc_tabix.communicate()
        print('[readcomb] done.')

if __name__ == '__main__':
    main()
