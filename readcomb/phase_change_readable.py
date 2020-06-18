import pandas as pd
import numpy as np
import pysam   
import collections
import os
from Bio import SeqIO
from cyvcf2 import VCF

reference_files_path = os.path.dirname(__file__) + '../tests/'

def check_snps(f_name, chromosome, left_bound, right_bound):
    '''
    Generate all SNPs on given chromsome and VCF file within the left_bound and right_bound using cyvcf2
    
    f_name - string: filepath of a VCF file
    chromsome - string: chromosome name
    left_bound/right_bound - integer: 0-based index of reference sequence

    Return list of cyvcf Variant objects
    '''

    vcf_in = VCF(f_name)
    # 1 is added to record.reference_start and the following parameter because vcf is 1 indexed
    # in order to keep code consistent
    region = '{c}:{l}-{r}'.format(c=chromosome, l=left_bound+1, r=right_bound+1)
    records = [rec for rec in vcf_in(region)]
    return records


def cache_pairs(bam_file_obj):
    '''
    Iterates through a bam file to find mate pairs and cache them together in a dictionary

    bam_file_obj: pysam alignment file object

    Returns a dictionary with unique sequence read id as the key and a tuple pair of bam 
    records as the value. If there is no mate pair, the second object in the tuple is None.    
    '''
    
    cache = {}
    
    paired = 0
    unpaired = 0
    
    for record in bam_file_obj:
        name = record.query_name
        
        if name not in cache:
            cache[name] = [record,None]
            unpaired += 1
        else:
            cache[name][1] = record
            paired += 1
            unpaired -= 1
            
    print('Paired: {paired}, Unpaired: {unpaired}'.format(paired=paired, unpaired=unpaired))    
    return cache

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

        # 1 is an insertion, we will add gaps to the reference and snip off extras on the end
        elif cigar_tuple[0] == 1:
            segment.append(query_segment[index:index+cigar_tuple[1]])
            
            ref = ref[:index] + '-' * cigar_tuple[1] + ref[index:]
            ref = ref[:-cigar_tuple[1]]
            
            index += cigar_tuple[1]

        # 2 is an deletion, add a gap to the query to realign it to the reference
        elif cigar_tuple[0] == 2:
            segment.append('-' * cigar_tuple[1])

        else:
            raise Exception('No condition for tuple ' + str(cigar_tuples.index(cigar_tuple)) + ' of ' + str(cigar_tuples))

    return ref, ''.join(segment)

def phase_detection(snps, segment, record):
    
    '''
    snps: list of snps in the area of the sequence and record given by the function check_snps(),
        each snp is a cyvcf2 variant object
    segment: string sequence that is the segment built by the function 
    '''
    
    snp_lst = [' '] * record.query_alignment_length
    
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

    return ''.join(snp_lst)

def write_recomb(pair_info, f_obj, snp_str, record, ref, segment):
                    
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

def human_read_matepairs_recomb(bam_file_obj, vcf, mode='no_match', output_filename='recomb_diagnosis'):
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
    
    f_obj = open(output_filename + '.txt', 'w')

    f_obj.write('Key: \nSequence: Start - End \nReference sequence \nPhred Scale Quality \nQuery alignment sequence \n')
    f_obj.write('1: CC2935, 2: CC2936, N: Does not match SNP \n\n')
    
    # get reference segment
    seq_obj = SeqIO.parse(reference_files_path + 'chlamy.5.3.w_organelles_mtMinus.fasta', 'fasta')
    
    # grab chromosome 1 from generator seq_obj
    chrom_1 = next(seq_obj)
        
    #counters
    no_match_counter = 0
    phase_change_counter = 0
    phase_change_mate_pair_counter = 0
    seq_with_snps_counter = 0
    all_seq_counter = 0
    
    # get mate pairs
    pairs = cache_pairs(bam_file_obj)
    
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
                    
        snps = check_snps(vcf, record.reference_name, 
                          record.reference_start,
                          record.reference_start + record.query_alignment_length)
        
        snp_str = phase_detection(snps, segment, record)
        
        if len(snps) > 0:
            seq_with_snps_counter += 1
        
        # initialize second record if there is a matepair
        if pairs[query_name][1]:
            
            all_seq_counter += 1
            
            record2 = pairs[query_name][1]
            
            #qualities_str2 = ''.join([chr(quality + 33) for quality in record2.query_alignment_qualities])
            
            ref2 = chrom_1[record2.reference_start:record2.reference_start + record2.query_alignment_length] 
            
            ref2, segment2 = human_cigar(record2, ref2)
            
            snps2 = check_snps(vcf, record2.reference_name, 
                              record2.reference_start,
                              record2.reference_start + record2.query_alignment_length)            
            
            snp_str2 = phase_detection(snps2, segment2, record2)
            
            if len(snps2) > 0:
                seq_with_snps_counter += 1
            
        #no_match_counter update - moved this out here because it's needed for both below modes
        if 'N' in snp_str:
            no_match_counter += 1

            if 'no_match' in mode:                        
                write_recomb('Unpaired', f_obj, snp_str, record, ref, segment)
                
        # pair 2 doesn't exist and pair 1 has more than 2 snps
        if pairs[query_name][1] == None:

            #phase_change_counter update
            if '1' in snp_str and '2' in snp_str:
                phase_change_counter += 1

                if 'phase_change' in mode:                        
                    write_recomb('Unpaired', f_obj, snp_str, record, ref, segment)


            if mode == 'all':
                write_recomb('Unpaired', f_obj, snp_str, record, ref, segment)

        # both pairs exist        
        else:
            
                        
            if mode == 'all':
                write_recomb('Pair 1', f_obj, snp_str, record, ref, segment)

                write_recomb('Pair 2', f_obj, snp_str2, record2, ref2, segment2)
            
            # no match in second pair
            elif 'N' in snp_str2:
                no_match_counter += 1

                if 'no_match' in mode:                        
                    write_recomb('Unpaired', f_obj, snp_str2, record2, ref2, segment2)
                
                
            # phase change across mate pairs
            elif ('1' in (snp_str + snp_str2) and '2' in (snp_str + snp_str2)):
                phase_change_mate_pair_counter += 1
                
                if mode == 'phase_change':
                    write_recomb('Pair 1', f_obj, snp_str, record, ref, segment)

                    write_recomb('Pair 2', f_obj, snp_str2, record2, ref2, segment2)

    
    # print counters
    print('Sequences with phase changes: ' + str(phase_change_counter))
    print('Sequences with a no match: ' + str(no_match_counter))
    print('Sequences with phase change across mate pairs: ' + str(phase_change_mate_pair_counter))
    print('Sequences with SNP(s): ' + str(seq_with_snps_counter))
    print('No. of sequences: ' + str(all_seq_counter))

if __name__ == '__main__':
    bam_file_obj = pysam.AlignmentFile('../tests/mock.sam', 'r')
    human_read_matepairs_recomb(bam_file_obj, '../tests/mock.vcf.gz', mode='all')
        
    
