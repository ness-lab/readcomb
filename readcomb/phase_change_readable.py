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
    reference: modified reference sequence that was given with gaps if there are insertions in the query sequence
    segment: dna sequence that is built using the bam cigar tuples and query segment
    
    Function for human readable version of mate pair sequence analysis. Takes in a bam record and a reference
    sequence and builds the query segment using the index and cigar tuple given by the bam record
    '''
    cigar_tuples = record.cigartuples
    
    # initialize segment and reference for building
    segment = []
    reference = []
    
    # index to keep track of where we are in the query_segment
    query_segment = record.query_sequence
    query_id = 0
    ref_id = 0 # reference segment does not include soft clips
    
    # no alignment, cigar_tuple is equal to None
    if cigar_tuples == None:
        return ref, segment

    for cigar_tuple in cigar_tuples:
        # 4 = soft clipping, the record.query_sequence has the portion that is soft clipping
        # so we need to skip it with index
        if cigar_tuple[0] == 4:
            query_id += cigar_tuple[1]

        # 5 = hard clipping, record.query_sequence does not have the portion that is
        # hard clipping so we don't skip it and we don't add anything
        elif cigar_tuple[0] == 5:
            continue

        # 0 is a match, just add it onto the segment and reference
        elif cigar_tuple[0] == 0:
            segment.append(query_segment[query_id:query_id+cigar_tuple[1]])
            reference.append(ref[ref_id:ref_id+cigar_tuple[1]])

            query_id += cigar_tuple[1]
            ref_id += cigar_tuple[1]

        # 1 is an insertion, we will add gaps to the reference and snip off extras on the end
        elif cigar_tuple[0] == 1:
            segment.append(query_segment[query_id:query_id+cigar_tuple[1]])
            reference.append('-' * cigar_tuple[1])
            
            query_id += cigar_tuple[1]

        # 2 is an deletion, add a gap to the query to realign it to the reference
        elif cigar_tuple[0] == 2:
            segment.append('-' * cigar_tuple[1])
            reference.append(ref[ref_id:ref_id+cigar_tuple[1]])
            
            ref_id += cigar_tuple[1]

        else:
            raise Exception('No condition for tuple ' + str(cigar_tuples.index(cigar_tuple)) + ' of ' + str(cigar_tuples))

    return ''.join(reference), ''.join(segment)

def human_phase_detection(snps, segment, record):
    
    '''
    Parameters:
    snps: list of snps/indels in the area of the sequence and record given by the function check_snps(),
        each snp/indel is a cyvcf2 variant object
    segment: string sequence that is the segment built by the function human_cigar()

    Returns:
    snp_str: string of spaces, 1, 2, and N that align with the segment and reference that
        signify which parent each SNP/indel belongs to

    Takes in a list of snps that are in the range of the segment/record. The function
    then takes each snp/indel, realigns them if there are insertions, then compares 
    the snp/indel of each parent to the segment, outputting the results to the snp_str
    '''
    
    snp_lst = [' '] * len(segment)
    
    cigar_tuples = record.cigartuples
    
    # check for no alignment
    if cigar_tuples == None:
        return ''.join(snp_lst)
    
    for snp in snps:
        # Using SNP.start and record.reference_start since they are both 0 based
        # SNP.start grabs vcf positions in 0 index while vcfs are 1 indexed

        idx = snp.start - record.reference_start

        if idx < 0:
            continue

        # extra calculations to realign idx if there is an insertion
        current_tuple = 0
        current_base = 0

        while current_base < idx and current_tuple < len(cigar_tuples):
            if cigar_tuples[current_tuple][0] == 1:
                # shift the idx to the right by the amount of insertion to compensate for it
                idx += cigar_tuples[current_tuple][1]
            
            # the segment doesn't include soft clipped and hard clipped parts anymore
            if cigar_tuples[current_tuple][0] != 4 and cigar_tuples[current_tuple][0] != 5:
                current_base += cigar_tuples[current_tuple][1]
            
            current_tuple += 1

        parent1 = snp.gt_bases[0].split('/')[0]
        parent2 = snp.gt_bases[1].split('/')[0]
  
        # skip uninformative sites
        if parent1 == parent2:
            continue

        if snp.is_indel:
            # check if indel is outside of segment
            if idx + max(len(parent1), len(parent2)) > len(segment):
                continue

            # check longer strand first
            if len(parent1) > len(parent2):
                if segment[idx:idx+len(parent1)] == parent1:
                    snp_lst[idx] = '1'
                elif segment[idx:idx+len(parent2)] == parent2:
                    snp_lst[idx] = '2'
                else:
                    snp_lst[idx] = 'N'
            else:
                if segment[idx:idx+len(parent2)] == parent2:
                    snp_lst[idx] = '2'
                elif segment[idx:idx+len(parent1)] == parent1:
                    snp_lst[idx] = '1'
                else:
                    snp_lst[idx] = 'N'
        else:
            if idx >= len(segment):
                break

            if segment[idx] == parent1:
                snp_lst[idx] = '1'

            elif segment[idx] == parent2:
                snp_lst[idx] = '2'

            else:
                snp_lst[idx] = 'N'

    return ''.join(snp_lst)

def human_parent_squence(snps, segment, reference, snp_str, record):
    '''
    Parameters:
    snps: list of snps/indels in the area of the sequence and record given by the function check_snps(),
        each snp/indel is a cyvcf2 variant object
    segment: string sequence that is the segment built by the function human_cigar()
    reference: modified reference sequence that is realigned with segment by the function human_cigar()
    snp_str: string of spaces, 1, 2, and N that align with the segment and reference that
        signify which parent each SNP/indel belongs to from the function human_phase_detection()
    record: bam record from pysam alignment file

    Returns:
    reference: reference segment with extra gaps to align with parent1_seq and parent2_seq
    segment: segment with extra gaps to align with parent1_seq and parent2_seq
    snp_str: snp_str with extra spaces to align with parent1_seq and parent2_seq
    parent1_seq: modified reference sequence with all snps and indels of parent1
    parent2_seq: modified reference sequence with all snps and indels of parent2

    Takes in a list of snps/indels in the range of the segment/reference/record and
    uses the reference to generate sequences for both parents. Also returns a modified
    segment, reference, and snp_str with extra gaps if parent1 or parent2 contains
    a snp/indel that needs more space 

    '''
    cigar_tuples = record.cigartuples
    parent1_seq = [i for i in reference]
    parent2_seq = [i for i in reference]
    
    # index to avoid SNPs in the middle of deletion dashes
    deletion_override_idx = 0

    for snp in snps:

        idx = snp.start - record.reference_start

        if idx < 0:
            continue

        parent1 = snp.gt_bases[0].split('/')[0]
        parent2 = snp.gt_bases[1].split('/')[0]

        # extra calculations to realign idx if there is an insertion
        current_tuple = 0
        current_base = 0

        while current_base < idx and current_tuple < len(cigar_tuples):
            if cigar_tuples[current_tuple][0] == 1:
                # shift the idx to the right by the amount of insertion to compensate for it
                idx += cigar_tuples[current_tuple][1]
            
            # the segment doesn't include soft clipped and hard clipped parts anymore
            if cigar_tuples[current_tuple][0] != 4 and cigar_tuples[current_tuple][0] != 5:
                current_base += cigar_tuples[current_tuple][1]
            
            current_tuple += 1
        
        if snp.is_indel:
            # update deletion override if it is a deletion
            if snp.is_deletion:
                deletion_override_idx = idx + max(len(parent1), len(parent2))
            
            # check if indel is outside of segment
            if idx + max(len(parent1), len(parent2)) > len(segment):
                continue

            if len(parent1) > len(parent2):
                # make reference and segment longer if needed
                if reference[idx:idx + len(parent1)] != parent1 and segment[idx:idx + len(parent1)] != parent1:
                    reference = reference[:idx] + '-' * len(parent1) + reference[idx:]
                    segment = segment[:idx] + '-' * len(parent1) + segment[idx:]
                    snp_str = snp_str[:idx + 1] + ' ' * (len(parent1) - 1) + snp_str[idx + 1:]
                    parent1_seq.insert(idx, parent1)
                    parent2_seq.insert(idx, parent2 + '-' * (len(parent1) - len(parent2)))
                else:
                    parent1_seq[idx:idx + len(parent1)] = parent1
                    parent2_seq[idx:idx + len(parent1)] = parent2 + '-' * (len(parent1) - len(parent2))
            else:
                # make reference and segment longer if needed
                if reference[idx:idx + len(parent2)] != parent2 and segment[idx:idx + len(parent2)] != parent2:
                    reference = reference[:idx] + '-' * len(parent2) + reference[idx:]
                    segment = segment[:idx] + '-' * len(parent2) + segment[idx:]
                    snp_str = snp_str[:idx + 1] + ' ' * (len(parent2) - 1) + snp_str[idx + 1:]
                    parent1_seq.insert(idx, parent1 + '-' * (len(parent2) - len(parent1)))
                    parent2_seq.insert(idx, parent2)
                else:
                    parent1_seq[idx:idx + len(parent2)] = parent1 + '-' * (len(parent2) - len(parent1))
                    parent2_seq[idx:idx + len(parent2)] = parent2 
            

        # simple snp base switching for the parents
        else:
            # skip SNP if it's in the middle of a deletion
            if deletion_override_idx > idx:
                continue

            if idx >= len(parent1_seq) or idx >= len(parent2_seq):
                break

            parent1_seq[idx] = parent1
            parent2_seq[idx] = parent2

    return reference, segment, snp_str, ''.join(parent1_seq), ''.join(parent2_seq)
            

def human_write_recomb(pair_info, f_obj, snp_str, record, ref, segment, parent1, parent2):
                    
    f_obj.write('{info}: {name} \n'.format(info=pair_info, name=record.query_name))
    
    if record.reference_end:
        f_obj.write('Ref Alignment: {start} - {end} \n'.format(start=record.reference_start,
                                                              end=record.reference_end))
    else:
        f_obj.write('Ref Alignment: {start} - {end} \n'.format(start=record.reference_start, 
                                                   end=record.reference_start+record.query_alignment_length))

    f_obj.write('Cigar Tuples: ' + str(record.cigartuples) + '\n')

    f_obj.write(parent1 + '\n' + parent2 + '\n' + ref + '\n' + segment + '\n' + snp_str +  '\n \n')

    return

def human_read_matepairs_recomb():
    '''
    Given a bam file object, return a human-readable file with aligned reference, quality string, 
        bam sequence, and SNPs    

    If mode='no_match', create a file with only no_match sequences.
    If mode='all', create a file with all sequences
    If mode='phase_change', create a file with only sequences with phase changes
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

    f_obj.write('Key: \nSequence: Start - End \nParent1 \nParent2 \nReference sequence \nQuery alignment sequence \n')
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
        ref = str(chrom_1[record.reference_start:record.reference_start + record.query_alignment_length].seq) 
        
        # analyze cigar string
        ref, segment = human_cigar(record, ref)
                    
        snps = check_snps(vcf_file_obj, record.reference_name, 
                          record.reference_start,
                          record.reference_start + record.query_alignment_length)
        
        snp_str = human_phase_detection(snps, segment, record)

        ref, segment, snp_str, first_parent1, first_parent2 = human_parent_squence(snps, segment, ref, snp_str, record)
        
        if len(snps) > 0:
            seq_with_snps_counter += 1
        
        # initialize second record if there is a matepair
        if pairs[query_name][1]:
            
            counters['seq'] += 1
            
            record2 = pairs[query_name][1]
            
            #qualities_str2 = ''.join([chr(quality + 33) for quality in record2.query_alignment_qualities])
            
            ref2 = str(chrom_1[record2.reference_start:record2.reference_start + record2.query_alignment_length].seq) 
            
            ref2, segment2 = human_cigar(record2, ref2)
            
            snps2 = check_snps(vcf_file_obj, record2.reference_name, 
                              record2.reference_start,
                              record2.reference_start + record2.query_alignment_length)            
            
            snp_str2 = human_phase_detection(snps2, segment2, record2)

            ref2, segment2, snp_str2, second_parent1, second_parent2 = human_parent_squence(snps2, segment2, ref2, snp_str2, record2)
            
            if len(snps2) > 0:
                counters['seq_with_snps'] += 1
            
        # no_match_counter update - moved this out here because it's needed for both below modes
        if 'N' in snp_str:
            counters['no_match'] += 1

            if 'no_match' in args['mode']:                        
                human_write_recomb('Unpaired', f_obj, snp_str, record, ref, segment, first_parent1, first_parent2)
                
        # pair 2 doesn't exist and pair 1 has more than 2 snps
        if pairs[query_name][1] == None:

            #phase_change_counter update
            if '1' in snp_str and '2' in snp_str:
                counters['phase_change'] += 1

                if 'phase_change' in args['mode']:                        
                    human_write_recomb('Unpaired', f_obj, snp_str, record, ref, segment, first_parent1, first_parent2)


            if args['mode'] == 'all':
                human_write_recomb('Unpaired', f_obj, snp_str, record, ref, segment, first_parent1, first_parent2)

        # both pairs exist        
        else:
            
                        
            if args['mode'] == 'all':
                human_write_recomb('Pair 1', f_obj, snp_str, record, ref, segment, first_parent1, first_parent2)

                human_write_recomb('Pair 2', f_obj, snp_str2, record2, ref2, segment2, second_parent1, second_parent2)
            
            # no match in second pair
            if 'N' in snp_str2:
                counters['no_match'] += 1

                if 'no_match' in args['mode']:                        
                    human_write_recomb('Unpaired', f_obj, snp_str2, record2, ref2, segment2, second_parent1, second_parent2)
                
                
            # phase change across mate pairs
            if ('1' in (snp_str + snp_str2) and '2' in (snp_str + snp_str2)):
                counters['phase_change_mate_pair'] += 1
                
                if 'phase_change' in args['mode']:
                    human_write_recomb('Pair 1', f_obj, snp_str, record, ref, segment, first_parent1, first_parent2)

                    human_write_recomb('Pair 2', f_obj, snp_str2, record2, ref2, segment2, second_parent1, second_parent2)

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
