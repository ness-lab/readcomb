'''
phase_change_filter.py - filter BAM for reads containing phase changes

usage:
    python3.5 phase_change_filter.py \
            --filename input.bam \
            --vcf parental.vcf.gz \
            --out output.bam
'''

import os
import argparse
import pysam
from cyvcf2 import VCF
from tqdm import tqdm
from Bio import SeqIO

def args():
    parser = argparse.ArgumentParser(
        description='filter BAM for reads containing phase changes', 
        usage='python3.5 phase_change_filter.py [options]')

    parser.add_argument('-f', '--filename', required=True,
                        type=str, help='BAM to filter')
    parser.add_argument('-v', '--vcf', required=True,
                        type=str, help='VCF containing parents')
    parser.add_argument('-l', '--log', required=False,
                        type=str, help='Log metrics to provided filename')
    parser.add_argument('-o', '--out', required=True,
                        type=str, help='File to write to')

    args = parser.parse_args()

    return args.filename, args.vcf, args.log, args.out

def check_snps(f_name, chromosome, left_bound, right_bound):
    vcf_in = VCF(f_name)
    # 1 is added to record.reference_start and the following parameter because vcf is 1 indexed
    # in order to keep code consistent
    region = '{c}:{l}-{r}'.format(c=chromosome, l=left_bound+1, r=right_bound+1)
    records = [rec for rec in vcf_in(region)]
    return records


def parse_cigar(cigar_tuples, query_sequence):
    # TODO: improve cigar parsing
    segment = ''

    # if read is clean 150 matches:
    if cigar_tuples == [(0, 150)]:
        segment = query_sequence
        return segment
    
    # otherwise -
    index = 0
    
    for cigar_tuple in cigar_tuples:
        # 4 = soft clipping, the record.query_sequence has the portion that is soft clipping
        # so we need to skip it with index
        if cigar_tuple[0] == 4:
            index += cigar_tuple[1]
            
        # 5 = hard clipping, record.query_sequence does not have the portion that is
        # hard clipping so we don't skip it and we don't add anything
        elif cigar_tuple[0] == 5:
            continue
        
        # if it is a match (0) or an insertion, then just add it onto the segment 
        elif cigar_tuple[0] in [0, 1]:
            segment += query_sequence[index:index+cigar_tuple[1]]
            index += cigar_tuple[1]
        
        else:
            print('oops forgot to consider this: ' + str(cigar_tuple))
            print(cigar_tuples)

    return segment


def parse_phase_changes(snps, cigar_tuples, segment, record):

    no_match_flag = False
    phase_change_flag = False
    previous_strand = None
    snp_lst = []
    
    for snp in snps:
        start = snp.start - record.reference_start # snp.start is 0-based
        
        current_tuple = 0
        current_base = 0
        
        while current_base < start and current_tuple < len(cigar_tuples):
            if cigar_tuples[current_tuple][0] == 1:
                # shift the start over by the amount of insertion to compensate for it
                start += cigar_tuples[current_tuple][1]
            
            current_base += cigar_tuples[current_tuple][1]
            current_tuple += 1
    
        # indexing for VCF seems to be a bit weird and will sometimes be -1
        if start < 0:
            raise Exception('VCF indexing is off. Check SNP at {}'.format(snp))
        
        strand1 = snp.gt_bases[0][0]
        strand2 = snp.gt_bases[1][0]
        
        if start >= len(segment):
            break
        
        if segment[start] == strand1:
            snp_lst[start] = '1'
        elif segment[start] == strand2:
            snp_lst[start] = '2'
        else:
            snp_lst[start] = 'N'
            
    snp_str = ''.join(snp_lst)
    if '1' in snp_str and '2' in snp_str:
        phase_change_flag = True
    if 'N' in snp_str:
        no_match_flag = True
            
    return phase_change_flag, no_match_flag


def main():
    bam, vcf, log, out = args()

    with open(out, 'w') as f_out:
        bam_file_obj = pysam.AlignmentFile(bam, 'r')
        f_obj = pysam.AlignmentFile(out + '.sam', 'wh', template=bam_file_obj) 
        
        # instantiate counters
        total_reads, phase_change_reads = 0, 0
        no_match_reads, no_snp_reads, no_match_phase_change_reads = 0, 0, 0
        for record in tqdm(bam_file_obj):
            total_reads += 1

            ### check for variation in region spanning read
            chrom = record.reference_name
            start, end = record.reference_start, \
                record.reference_start + record.query_alignment_length
            snps = check_snps(vcf, chrom, start, end)

            if len(snps) <= 1:
                no_snp_reads += 1
                continue # ignore reads without enough SNPs

            ### parse cigar string
            cigar_tuples = record.cigartuples
            segment = parse_cigar(cigar_tuples, record.query_sequence)

            ### scan for phase changes
            phase_change_flag, no_match_flag = parse_phase_changes(snps,
                    cigar_tuples, segment, record)
            if phase_change_flag:
                phase_change_reads += 1
                f_obj.write(record)
            if no_match_flag:
                no_match_reads += 1
            if phase_change_flag and no_match_flag:
                no_match_phase_change_reads += 1

    print('''
    Done.
    {} phase change reads extracted from {} total ({}%)
    {} reads had no-match variants, and {} of these were phase-change reads.
    {} reads did not have enough SNPs (>= 2) to call ({}%)
    '''.format(phase_change_reads, total_reads, phase_change_reads /
        total_reads, no_match_reads, no_match_phase_change_reads, no_snp_reads,
        no_snp_reads / total_reads))
    if log:
        needs_header = True
        if os.path.isfile(log):
            needs_header = False
        with open(log, 'a') as f:
            if needs_header:
                fieldnames = ['phase_change_reads', 'total_reads',
                'no_match_reads', 'no_match_phase_change_reads',
                'no_snp_reads']
                out_values = [phase_change_reads, total_reads, no_match_reads,
                        no_match_phase_change_reads, no_snp_reads]
                f.write(','.join(fieldnames) + '\n')
                f.write(','.join([str(n) for n in out_values]) + '\n')


if __name__ == '__main__':
    main()

        

