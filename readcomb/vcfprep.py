#!/usr/bin/env python3
"""
vcfprep.py - prepare parental VCF for phase change detection and apply further
filters if needed
"""

import subprocess
import argparse
from cyvcf2 import VCF
from cyvcf2 import Writer
from tqdm import tqdm

def arg_parser():
    parser = argparse.ArgumentParser(
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

    args = parser.parse_args()

    return args.vcf, args.snps_only, args.indels_only, \
            args.no_hets, args.min_GQ, args.out

def vcfprep(vcf, snps_only, indels_only, no_hets, min_GQ, outfile):
    """
    Iterate through parental VCF and remove records by given filters. 
    
    All records where parent1 allele is the same as parent2
    ('uninformative sites') are automatically removed.

    Raises ValueError if >2 samples in VCF (should just be 2 parents)

    Parameters
    ----------
    vcf : str
        path to input VCF
    snps_only : bool
        whether or not to only keep SNPs
    indels_only : bool
        whether or not to only keep indels
    no_hets : bool
        whether or not to remove het calls
    min_GQ : int
        minimum GQ value required to consider sites
    out : str
        path to file to write to

    Returns
    -------
    total_count : int
        total records considered
    kept_count : int
        number of records that passed filters
    """
    vcf_in = VCF(vcf)

    if len(vcf_in.samples) > 2:
        raise ValueError('more than 2 parental samples in input VCF')
    if snps_only and indels_only:
        raise ValueError('both --snps_only and --indels_only provided. pick one or neither!')

    if outfile.endswith('.gz'):
        outfile = outfile.replace('.gz', '')
    vcf_out = Writer(outfile, vcf_in)
    vcf_out.write_header()

    total_count, kept_count = 0, 0

    for record in tqdm(vcf_in):
        total_count += 1

        if not len(record.ALT) > 0:
            continue

        # SNP filter
        if snps_only and not record.is_snp:
            continue

        # indel filter
        if indels_only and not record.is_indel:
            continue

        # heterozygote call filter
        if no_hets and record.num_het != 0:
            continue

        # genotype quality filter
        if not all(record.gt_quals >= min_GQ):
            continue

        # ensure parental alleles differ
        if record.gt_bases[0] == record.gt_bases[1]:
            continue

        # only passes if record not caught in above filters
        vcf_out.write_record(record)
        kept_count += 1

    return total_count, kept_count


def main():
    vcf, snps_only, indels_only, no_hets, min_GQ, out = arg_parser()
    print('[readcomb] filtering {}'.format(vcf))
    if min_GQ < 30:
        print('[readcomb] WARNING: min GQ below 30 selected')
    total_count, kept_count = vcfprep(vcf, snps_only, indels_only, no_hets, min_GQ, out)
    print('[readcomb] Complete.')
    print('[readcomb] {} of {} records retained.'.format(kept_count, total_count))
    if out.endswith('.gz'):
        print('[readcomb] compressing outfile...')
        proc_bgzip = subprocess.Popen(['bgzip', out.replace('.gz', '')])
        stdout, stderr = proc_bgzip.communicate()
        print('[readcomb] creating tabix index...')
        proc_tabix = subprocess.Popen(['tabix', out])
        stdout, stderr = proc_tabix.communicate()
        print('[readcomb] done.')

if __name__ == '__main__':
    main()

        

