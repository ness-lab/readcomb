import argparse
import collections
import os
import time
import datetime

import numpy as np
import pysam
from Bio import SeqIO
from cyvcf2 import VCF

from phase_change_filter import cache_pairs, check_snps


def arg_parser():
    parser = argparse.ArgumentParser(
        description='filter BAM for reads containing phase changes and outputs data in a human readable text file', 
        usage='python3.5 phase_change_filter.py [options]')

    parser.add_argument('-b', '--bam', required=True,
                        type=str, help='BAM to filter')
    
    parser.add_argument('-v', '--vcf', required=True,
                        type=str, help='VCF containing parents')

    parser.add_argument('-c', '--chrom', required=False,
                        type=str, default='all', help='Specify which chromsome sequences to run, default is all')
    
    parser.add_argument('-r', '--reference', required=True,
                        type=str, help='FASTA reference file to align bam reads to')
    
    parser.add_argument('-m', '--mode', required=False,
                        type=str, default='phase_change no_match', help='Mode to execute the program')
    
    parser.add_argument('-l', '--log', required=False,
                        type=str, help='Log metrics to provide filename')
    
    parser.add_argument('-o', '--out', default='recomb_diagnosis.txt', required=False,
                        type=str, help='File to write to')

    parser.add_argument('-i', '--improper', action='store_true',
                        help='Include unpaired/improper reads in filter results')

    args = parser.parse_args()

    return {
        "bam": args.bam, 
        "vcf": args.vcf, 
        "chrom": args.chrom, 
        "ref": args.reference, 
        "mode": args.mode, 
        "log": args.log, 
        "out": args.out,
        "improper": args.improper,
        }

def human_cigar(record, ref):
    '''
    Parameters:
    record: bam record from pysam alignment file
    ref: reference sequence that is aligned with the query sequence obtained from a fasta file
    
    Returns:
    ref: modified reference sequence that was given with gaps if there are insertions in the query sequence
    segment: dna sequence that is built using the bam cigar tuples and query segment
    
    Function for human readable version of mate pair sequence analysis. Takes in a bam record and a reference
    sequence and builds the query segment using the index and cigar tuple given by the bam record
    '''
    cigar_tuples = record.cigartuples
    
    # initialize segment for building
    segment = []
    
    # index to keep track of where we are in the query_segment
    query_segment = record.query_sequence
    index = 0
    ref_id = 0 # reference segment does not include soft clips
    
    # no alignment, cigar_tuple is equal to None
    if cigar_tuples == None:
        return ref, segment

    for cigar_tuple in cigar_tuples:
        # 4 = soft clipping, the record.query_sequence has the portion that is soft clipping
        # so we need to skip it with index
        if cigar_tuple[0] == 4:
            index += cigar_tuple[1]

        # 5 = hard clipping, record.query_sequence does not have the portion that is
        # hard clipping so we don't skip it and we don't add anything
        elif cigar_tuple[0] == 5:
            continue

        # 0 is a match, just add it onto the segment 
        elif cigar_tuple[0] == 0:
            segment.append(query_segment[index:index+cigar_tuple[1]])
            index += cigar_tuple[1]
            ref_id += cigar_tuple[1]

        # 1 is an insertion, we will add gaps to the reference and snip off extras on the end
        elif cigar_tuple[0] == 1:
            segment.append(query_segment[index:index+cigar_tuple[1]])
            
            ref = ref[:ref_id] + '-' * cigar_tuple[1] + ref[ref_id:]
            
            index += cigar_tuple[1]
            ref_id += cigar_tuple[1]

        # 2 is an deletion, add a gap to the query to realign it to the reference
        elif cigar_tuple[0] == 2:
            segment.append('-' * cigar_tuple[1])
            ref_id += cigar_tuple[1]

        else:
            raise Exception('No condition for tuple ' + str(cigar_tuples.index(cigar_tuple)) + ' of ' + str(cigar_tuples))

    return ref, ''.join(segment)

def human_phase_detection(snps, segment, record):
    
    '''
    snps: list of snps in the area of the sequence and record given by the function check_snps(),
        each snp is a cyvcf2 variant object
    segment: string sequence that is the segment built by the function 
    '''
    
    snp_lst = [' '] * len(segment)
    
    cigar_tuples = record.cigartuples
    
    # check for no alignment
    if cigar_tuples == None:
        return ''.join(snp_lst)
    
    for snp in snps:
        # Using SNP.start and record.reference_start since they are both 0 based
        # SNP.start grabs vcf positions in 0 index while vcfs are 1 indexed

        start = snp.start - record.reference_start

        # extra calculations to realign start if there is an insertion or deletion
        current_tuple = 0
        current_base = 0

        while current_base < start and current_tuple < len(cigar_tuples):
            if cigar_tuples[current_tuple][0] == 1:
                # shift the start to the right by the amount of insertion to compensate for it
                start += cigar_tuples[current_tuple][1]
            
            # the segment doesn't include soft clipped and hard clipped parts anymore
            if cigar_tuples[current_tuple][0] != 4 and cigar_tuples[current_tuple][0] != 5:
                current_base += cigar_tuples[current_tuple][1]
            
            current_tuple += 1
        
        # indexing for VCF seems to be a bit weird and will sometimes be -1
        # For indels, this can be valid
        #if start < 0:
        #    raise Exception('VCF indexing is off. Check SNP at {}'.format(snp))

        strand1 = snp.gt_bases[0].split('/')[0]
        strand2 = snp.gt_bases[1].split('/')[0]

        # weird anomaly
        if strand1 == '.' or strand2 == '.':
            continue
  
        # skip uninformative sites
        if strand1 == strand2:
            continue

        if snp.is_indel:
            # snip off the first position if it's a pure indel
            if strand1[0] == strand2[0]:
                strand1 = strand1[1:]
                strand2 = strand2[1:]
                start += 1
            
            # check if indel is outside of segment
            if start + max(len(strand1), len(strand2)) > len(segment):
                continue

            # check longer strand first
            if len(strand1) > len(strand2):
                if segment[start:start+len(strand1)] == strand1:
                    snp_lst[start] = '1'
                elif segment[start:start+len(strand2)] == strand2:
                    snp_lst[start] = '2'
                else:
                    snp_lst[start] = 'N'
            else:
                if segment[start:start+len(strand2)] == strand2:
                    snp_lst[start] = '2'
                elif segment[start:start+len(strand1)] == strand1:
                    snp_lst[start] = '1'
                else:
                    snp_lst[start] = 'N'
        else:
            if start >= len(segment):
                break

            if segment[start] == strand1:
                snp_lst[start] = '1'

            elif segment[start] == strand2:
                snp_lst[start] = '2'

            else:
                snp_lst[start] = 'N'

    return ''.join(snp_lst)

def human_write_recomb(pair_info, f_obj, snp_str, record, ref, segment):
                    
    f_obj.write('{info}: {name} \n'.format(info=pair_info, name=record.query_name))
    
    if record.reference_end:
        f_obj.write('Ref Alignment: {start} - {end} \n'.format(start=record.reference_start,
                                                              end=record.reference_end))
    else:
        f_obj.write('Ref Alignment: {start} - {end} \n'.format(start=record.reference_start, 
                                                   end=record.reference_start+record.query_alignment_length))

    f_obj.write('Cigar Tuples: ' + str(record.cigartuples) + '\n')

    f_obj.write(str(ref.seq) + '\n' + segment + '\n' + snp_str +  '\n \n')

    return

def human_read_matepairs_recomb():
    '''
    Given a bam file object, return a human-readable file with aligned reference, quality string, 
        bam sequence, and SNPs    

    If mode='no_match', create a file with only no_match sequences.
    If mode='all', create a file with all sequences
    If mode='phase_change', create a file with only sequences with phase changes
    
    Makes assumptions on the following file locations:
        VCF: vcf_files_path + 'parental_filtered.vcf.gz'
        Reference: reference_files_path + 'chlamy.5.3.w_organelles_mtMinus.fasta'
    '''
    # start timer
    start = time.time()

    args = arg_parser()

    counters = {"no_match": 0,
                "phase_change": 0,
                "phase_change_mate_pair": 0,
                "seq_with_snps": 0,
                "seq": 0}

    # create bam pysam alignment file object
    bam_file_obj = pysam.AlignmentFile(args['bam'], 'r')

    vcf_file_obj = VCF(args['vcf'])

    f_obj = open(args['out'], 'w')

    f_obj.write('Key: \nSequence: Start - End \nReference sequence \nPhred Scale Quality \nQuery alignment sequence \n')
    f_obj.write('1: CC2935, 2: CC2936, N: Does not match SNP \n\n')
    
    # get reference segment
    seq_obj = SeqIO.parse(args['ref'], 'fasta')
    
    # grab chromosome 1 from generator seq_obj
    chrom_1 = next(seq_obj)
        
    #counters
    no_match_counter = 0
    phase_change_counter = 0
    phase_change_mate_pair_counter = 0
    seq_with_snps_counter = 0
    all_seq_counter = 0
    
    # get mate pairs
    pairs, counters['paired'], counters['unpaired'] = cache_pairs(bam_file_obj, args)
    
    for query_name in pairs:
        
        all_seq_counter += 1
        
        # initialize first record
        record = pairs[query_name][0]
        
        # qualities string
        #qualities_str = ''.join([chr(quality + 33) for quality in record.query_alignment_qualities])
        
        # reference sequence
        ref = chrom_1[record.reference_start:record.reference_start + record.query_alignment_length] 
        
        # analyze cigar string
        ref, segment = human_cigar(record, ref)
                    
        snps = check_snps(vcf_file_obj, record.reference_name, 
                          record.reference_start,
                          record.reference_start + record.query_alignment_length)
        
        snp_str = human_phase_detection(snps, segment, record)
        
        if len(snps) > 0:
            seq_with_snps_counter += 1
        
        # initialize second record if there is a matepair
        if pairs[query_name][1]:
            
            counters['seq'] += 1
            
            record2 = pairs[query_name][1]
            
            #qualities_str2 = ''.join([chr(quality + 33) for quality in record2.query_alignment_qualities])
            
            ref2 = chrom_1[record2.reference_start:record2.reference_start + record2.query_alignment_length] 
            
            ref2, segment2 = human_cigar(record2, ref2)
            
            snps2 = check_snps(vcf_file_obj, record2.reference_name, 
                              record2.reference_start,
                              record2.reference_start + record2.query_alignment_length)            
            
            snp_str2 = human_phase_detection(snps2, segment2, record2)
            
            if len(snps2) > 0:
                counters['seq_with_snps'] += 1
            
        # no_match_counter update - moved this out here because it's needed for both below modes
        if 'N' in snp_str:
            counters['no_match'] += 1

            if 'no_match' in args['mode']:                        
                human_write_recomb('Unpaired', f_obj, snp_str, record, ref, segment)
                
        # pair 2 doesn't exist and pair 1 has more than 2 snps
        if pairs[query_name][1] == None:

            #phase_change_counter update
            if '1' in snp_str and '2' in snp_str:
                counters['phase_change'] += 1

                if 'phase_change' in args['mode']:                        
                    human_write_recomb('Unpaired', f_obj, snp_str, record, ref, segment)


            if args['mode'] == 'all':
                human_write_recomb('Unpaired', f_obj, snp_str, record, ref, segment)

        # both pairs exist        
        else:
            
                        
            if args['mode'] == 'all':
                human_write_recomb('Pair 1', f_obj, snp_str, record, ref, segment)

                human_write_recomb('Pair 2', f_obj, snp_str2, record2, ref2, segment2)
            
            # no match in second pair
            if 'N' in snp_str2:
                counters['no_match'] += 1

                if 'no_match' in args['mode']:                        
                    human_write_recomb('Unpaired', f_obj, snp_str2, record2, ref2, segment2)
                
                
            # phase change across mate pairs
            if ('1' in (snp_str + snp_str2) and '2' in (snp_str + snp_str2)):
                counters['phase_change_mate_pair'] += 1
                
                if 'phase_change' in args['mode']:
                    human_write_recomb('Pair 1', f_obj, snp_str, record, ref, segment)

                    human_write_recomb('Pair 2', f_obj, snp_str2, record2, ref2, segment2)

    # end timer
    end = time.time()
    runtime = str(datetime.timedelta(seconds=round(end - start)))

    print('[readcomb] Done.')
    print('[readcomb] Run stats:')
    print('{} phase change reads from {} total unpaired'.format(counters['phase_change'], 
        counters['unpaired']))
    print('{} phase change reads across mate pairs from {} total read pairs'.format(
        counters['phase_change_mate_pair'], counters['paired']))
    print('{} reads had no-match variants'.format(counters['no_match']))
    print('{} reads did not have enough SNPs (> 0) to call'.format(
        counters['seq'] - counters['seq_with_snps']))
    print('time taken: {}'.format(runtime))

    if args["log"]:
        needs_header = True
        if os.path.isfile(args["log"]):
            needs_header = False
        with open(args["log"], 'a') as f:
            if needs_header:
                fieldnames = ['run', 'phase_change_reads', 'unpaired_reads',
                              'phase_change_across_mate_pair', 'read_pairs',
                              'no_match_reads', 'no_snp_reads', 'total_reads',
                              'time_taken']
                fieldnames.extend(sorted(args.keys()))
                f.write(','.join(fieldnames) + '\n')
            out_values = [datetime.datetime.now(), counters["phase_change"], counters["unpaired"],
                          counters["phase_change_mate_pair"], counters["paired"],
                          counters["no_match"], counters["seq"] - counters["seq_with_snps"],
                          counters["seq"], runtime]
            out_values.extend([args[arg] for arg in sorted(args.keys())])
            f.write(','.join([str(n) for n in out_values]) + '\n')

if __name__ == '__main__':
    human_read_matepairs_recomb()
