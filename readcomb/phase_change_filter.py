import os
import argparse
import pysam
from cyvcf2 import VCF

def args():
    parser = argparse.ArgumentParser(
        description='filter BAM for reads containing phase changes', 
        usage='python3.5 phase_change_filter.py [options]')

    parser.add_argument('-f', '--filename', required=True,
                        type=str, help='BAM to filter')
    parser.add_argument('-v', '--vcf', required=True,
                        type=str, help='VCF containing parents')
    parser.add_argument('-m', '--mode', required=False,
                        type=str, default='phase_change', help='Mode to execute the program')
    parser.add_argument('-o', '--out', required=True,
                        type=str, help='File to write to')

    args = parser.parse_args()

    return args.filename, args.vcf, args.mode, args.out

def check_snps(f_name, chromosome, left_bound, right_bound):
    '''
    Using cyvcf, go through each of the breakpoints and check the number of SNPs per segment.
    Return the SNP locations in a list.
    '''
    vcf_in = VCF(f_name)
    # 1 is added to record.reference_start and the following parameter because vcf is 1 indexed
    # in order to keep code consistent
    region = '{c}:{l}-{r}'.format(c=chromosome, l=left_bound+1, r=right_bound+1)
    records = [rec for rec in vcf_in(region)]
    return records


def cache_pairs(bam_file_obj):

    '''
    Function to cache mate pairs of a bam file

    Returns a dictionary of the cached reads in the format:

    {read_name: (bam_record_1, bam_record_2), ...}
    '''
    
    cache = {}
    
    paired = 0
    total = 0
    
    for record in bam_file_obj:
        name = record.query_name
        total += 1
        
        if name not in cache:
            cache[name] = [record,None]
        else:
            cache[name][1] = record
            paired += 1
            
    print('Paired: {paired}, Unpaired: {unpaired}'.format(paired=paired, unpaired=total-paired))    
    return cache
        
    

def cigar(record):
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
    segment = ''
    
    # index to keep track of where we are in the query_segment
    query_segment = record.query_sequence
    index = 0
    
    # no alignment, cigar_tuple is equal to None
    if cigar_tuples == None:
        return segment

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
            segment += query_segment[index:index+cigar_tuple[1]]
            index += cigar_tuple[1]

        # 1 is an insertion, treated same as a match
        elif cigar_tuple[0] == 1:
            segment += query_segment[index:index+cigar_tuple[1]]
            index += cigar_tuple[1]\

        # 2 is an insertion, add a gap to the query
        elif cigar_tuple[0] == 2:
            segment += '-' * cigar_tuple[1]
            index += cigar_tuple[1]

        else:
            print('forgot to consider this: ' + str(cigar_tuple))
            print(cigar_tuples)
    
    return segment


def phase_detection(snps, segment, record):
    
    '''
    snps: list of snps in the area of the sequence and record given by the function check_snps(),
        each snp is a cyvcf2 variant object
    segment: string sequence that is the segment built by the function 
    '''
    
    snp_lst = []
    
    cigar_tuples = record.cigartuples
    
    # check for no alignment
    if cigar_tuples == None:
        return snp_lst
    
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
            
            elif cigar_tuples[current_tuple][0] == 2:
                # shift the start to the left by the amount of deletion to compensate for it
                start -= cigar_tuples[current_tuple][1]

            current_base += cigar_tuples[current_tuple][1]
            current_tuple += 1


        # indexing for VCF seems to be a bit weird and will sometimes be -1
        # this has been fixed but left check in anyways
        if start < 0:
            raise Exception('VCF indexing is off. Check SNP at {}'.format(snp))

        strand1 = snp.gt_bases[0][0]
        strand2 = snp.gt_bases[1][0]

        if start >= len(segment):
            break

        if segment[start] == strand1:
            snp_lst.append('1')

        elif segment[start] == strand2:
            snp_lst.append('2')

        else:
            snp_lst.append('N')

    return snp_lst

def matepairs_recomb():
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

    bam, vcf, mode, output_filename = args()

    bam_file_obj = pysam.AlignmentFile(bam, 'r')
    
    f_obj = pysam.AlignmentFile(output_filename + '.sam', 'wh', template=bam_file_obj)
    
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
        
        # analyze cigar string
        segment = cigar(record)
                    
        snps = check_snps(vcf, record.reference_name, 
                        record.reference_start,
                        record.reference_start + record.query_alignment_length)
        
        snp_lst = phase_detection(snps, segment, record)
        
        # initialize second record if there is a matepair
        if pairs[query_name][1]:
            
            seq_with_snps_counter += 1
            
            record2 = pairs[query_name][1]
            
            segment2 = cigar(record2)
            
            snps2 = check_snps(vcf, record2.reference_name, 
                            record2.reference_start,
                            record2.reference_start + record2.query_alignment_length)            
            
            snp_lst2 = phase_detection(snps2, segment2, record2)
            
        #no_match_counter update - moved this out here because it's needed for both unpaired and paired mates
        if 'N' in snp_lst:
            no_match_counter += 1

            if 'no_match' in mode:                        
                f_obj.write(record)
                
        # pair 2 doesn't exist
        if pairs[query_name][1] == None:

            #phase_change_counter update
            if '1' in snp_lst and '2' in snp_lst:
                phase_change_counter += 1

                if 'phase_change' in mode:                        
                    f_obj.write(record)


            if mode == 'all':
                f_obj.write(record)

        # both pairs exist        
        else:
            
            # no match in second pair
            if 'N' in snp_lst2:
                no_match_counter += 1

                if 'no_match' in mode:                        
                    f_obj.write(record2)
                
                
            # phase change across mate pairs
            elif ('1' in (snp_lst + snp_lst2) and '2' in (snp_lst + snp_lst2)):
                phase_change_mate_pair_counter += 1
                
                if mode == 'phase_change':
                    f_obj.write(record)

                    f_obj.write(record2)
    
    print('''
    Done.
    {} phase change reads extracted from {} total ({}%)
    {} phase change reads across mate pairs
    {} reads had no-match variants.
    {} reads did not have enough SNPs (>= 2) to call ({}%)
    '''.format(phase_change_counter, all_seq_counter, phase_change_counter / all_seq_counter, 
        phase_change_mate_pair_counter,
        no_match_counter, all_seq_counter - seq_with_snps_counter,
        (all_seq_counter - seq_with_snps_counter) / all_seq_counter)
    )


if __name__ == '__main__':
    matepairs_recomb()
