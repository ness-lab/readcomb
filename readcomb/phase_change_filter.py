import os
import argparse
import pysam
import time
import datetime
from cyvcf2 import VCF
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(
        description='filter BAM for reads containing phase changes', 
        usage='python3.5 phase_change_filter.py [options]')

    parser.add_argument('-b', '--bam', required=True,
                        type=str, help='BAM to filter, required')
                        
    parser.add_argument('-v', '--vcf', required=True,
                        type=str, help='VCF containing parents, required')

    parser.add_argument('-c', '--chrom', required=False,
                        type=str, default='all', help='Specify which chromsome sequences to run, default is all')

    parser.add_argument('-m', '--mode', required=False,
                        type=str, default='phase_change', help='Mode to execute the program, default is phase_change')

    parser.add_argument('-l', '--log', required=False,
                        type=str, help='Log metrics to provide filename, default is False')

    parser.add_argument('-o', '--out', required=False,
                        type=str, default='recomb_diagnosis', help='File to write to, default is recomb_diagnosis')

    args = parser.parse_args()

    return args.bam, args.vcf, args.chrom, args.mode, args.log, args.out

def check_snps(vcf_file_obj, chromosome, left_bound, right_bound):
    '''
    Generate all SNPs on given chromsome and VCF file within the left_bound and right_bound using cyvcf2
    
    vcf_file_obj - cyvcf2 VCF file object
    chromsome - string: chromosome name
    left_bound/right_bound - integer: 0-based index of reference sequence

    Return list of cyvcf Variant objects
    '''

    # 1 is added to record.reference_start and the following parameter because vcf is 1 indexed
    # in order to keep code consistent
    region = '{c}:{l}-{r}'.format(c=chromosome, l=left_bound+1, r=right_bound+1)\

    # list comp with if statement to only include SNPs
    records = [rec for rec in vcf_file_obj(region) if rec.is_snp and len(rec.ALT) > 0]
    return records


def cache_pairs(bam_file_obj, chromosome):
    '''
    Iterates through a bam file to find mate pairs and cache them together in a dictionary

    bam_file_obj: pysam alignment file object
    chromosome: string

    Returns a dictionary with unique sequence read id as the key and a tuple pair of bam 
    records as the value. If there is no mate pair, the second object in the tuple is None.    

    Also returns the number of unpaired reads (unpaired) and the number of mate pairs (paired)
    '''
    
    print('Caching reads for ' + chromosome + ' sequences')

    cache = {}
    
    paired = 0
    unpaired = 0
    
    for record in bam_file_obj:
        
        # error checking
        if record.is_secondary or record.is_supplementary:
            continue
        
        # chromosome argument checking
        if chromosome != record.reference_name and chromosome != 'all':
            continue 

        # check if query_name and reference_name exist
        if record.query_name == None or record.reference_name ==  None:
            continue

        name = record.query_name + record.reference_name
        
        if name not in cache:
            cache[name] = [record,None]
            unpaired += 1
        else:
            cache[name][1] = record
            paired += 2
            unpaired -= 1
            
    print('Number of unpaired sequences: {}, read pairs: {}'.format(unpaired, paired))
    return cache, paired, unpaired
        
    

def cigar(record):
    ''' 
    Build the query segment using the cigar tuple given by the bam record so that it aligns with
    indexes of SNPs in the VCF that are aligned to the reference sequence

    record: bam record from pysam alignment file
    
    returns dna sequence string that is built using the bam cigar tuples and query segment
    '''
    cigar_tuples = record.cigartuples
    
    # initialize segment for building
    segment = []
    
    # index to keep track of where we are in the query_segment
    query_segment = record.query_sequence
    index = 0
    
    # no alignment, cigar_tuple is equal to None
    if cigar_tuples == None:
        return segment

    for cigar_tuple in cigar_tuples:
      
        # 1 is an insertion to query segment, skip it because SNPs are aligned to reference and do not exist in this region
        # 4 is soft clipping, query sequence includes this but is not aligned to reference
        if cigar_tuple[0] == 4:
            index += cigar_tuple[1]

        elif cigar_tuple[0] == 1:
            segment.append(query_segment[index:index+cigar_tuple[1]])
            index += cigar_tuple[1]

        # 5 = hard clipping, record.query_sequence does not have the portion that is
        # hard clipping so skip it and don't add anything
        elif cigar_tuple[0] == 5:
            continue

        # 0 is a match, add it onto the segment 
        elif cigar_tuple[0] == 0:
            segment.append(query_segment[index:index+cigar_tuple[1]])
            index += cigar_tuple[1]
            
        # 2 is an deletion to query segment, add a gap to realign it to reference
        # don't add index to not moving further through the query segment
        elif cigar_tuple[0] == 2:
            segment.append('-' * cigar_tuple[1])

        else:
            raise Exception('No condition for tuple ' + str(cigar_tuples.index(cigar_tuple)) + ' of ' + str(cigar_tuples))
        
    return ''.join(segment)


def phase_detection(snps, segment, record):
    '''
    snps - list of variants: list of snps in the area of the sequence and record given by the function check_snps(),
        can include snps that are outside of the area
    segment - string: sequence that is the segment built by the function cigar()
    '''
    
    snp_lst = []
    
    cigar_tuples = record.cigartuples
    
    # check for no alignment
    if not cigar_tuples:
        return snp_lst
    
    for snp in snps:
        # Using SNP.start and record.reference_start since they are both 0 based
        # SNP.start grabs vcf positions in 0 index while vcfs are 1 indexed

        start = snp.start - record.reference_start

        # extra calculations to realign start if there is an insertion
        current_tuple = 0
        current_base = 0

        while current_base < start and current_tuple < len(cigar_tuples):
            if cigar_tuples[current_tuple][0] == 1:
                # shift the start to the right by the amount of insertion to compensate for it
                start += cigar_tuples[current_tuple][1]

            current_base += cigar_tuples[current_tuple][1]
            current_tuple += 1

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
    Parses arguements and filters bam records using SNPS from a vcf

    bam - string: bam filepath
    vcf - string: vcf filepath
    mode - string:
        phase_change - only write to output bam records that have phase changes
        no_match - only write to output bam records that have a base that does not match either variation of a SNP
    log - boolean: when true, logs sequence counts to a log file
    output_filename - string: file to write the filtered bams to, by default the function writes to recomb_diagnosis.sam
    '''    
    # start timer
    start = time.time()

    bam, vcf, chromosome, mode, log, output_filename = args()

    bam_file_obj = pysam.AlignmentFile(bam, 'r')

    vcf_file_obj = VCF(vcf)
    
    f_obj = pysam.AlignmentFile(output_filename + '.sam', 'wh', template=bam_file_obj)
    
    #counters
    no_match_counter = 0
    phase_change_counter = 0
    phase_change_mate_pair_counter = 0
    seq_with_snps_counter = 0
    all_seq_counter = 0
    
    pairs, paired, unpaired = cache_pairs(bam_file_obj, chromosome)

    print('Beginning phase change analysis')
    
    for query_name in tqdm(pairs):

        snp_lst = []
        for record in pairs[query_name]:
            if record:
                all_seq_counter += 1
                
                # analyze cigar string
                segment = cigar(record)
                            
                snps = check_snps(vcf_file_obj, record.reference_name, 
                                record.reference_start,
                                record.reference_start + record.query_alignment_length)
                
                snp_lst += phase_detection(snps, segment, record)

                if len(snps) > 0:
                    seq_with_snps_counter += 1
            
                

        if '1' in snp_lst and '2' in snp_lst:
            
            if pairs[query_name][1]:
                phase_change_mate_pair_counter += 1
            else:
                phase_change_counter += 1

            if 'phase_change' in mode:                   
                f_obj.write(pairs[query_name][0])

                if pairs[query_name][1]:                                        
                    f_obj.write(pairs[query_name][1])
        
        if 'N' in snp_lst:
            no_match_counter += 1

            if 'no_match' in mode:                        
                f_obj.write(record)

    # end timer
    end = time.time()
    runtime = str(datetime.timedelta(seconds=round(end - start)))
    
    print('''
    Done.
    {} phase changes reads from {} total unpaired ({}%)
    {} phase changes reads across mate pairs from {} total read pairs ({}%)
    {} reads had no-match variants.
    {} reads did not have enough SNPs (> 0) to call ({}%)
    time taken: {}
    '''.format(phase_change_counter, unpaired, round(phase_change_counter / unpaired * 100, 2), 
        phase_change_mate_pair_counter, paired, round(phase_change_mate_pair_counter / paired * 100, 2),
        no_match_counter, all_seq_counter - seq_with_snps_counter,
        round((all_seq_counter - seq_with_snps_counter) / all_seq_counter * 100, 2),
        runtime)
    )

    if log:
        needs_header = True
        if os.path.isfile(log):	
            needs_header = False	
        with open(log, 'a') as f:	
            if needs_header:	
                fieldnames = ['phase_change_reads', 'unpaired_reads',	
                'phase_change_across_mate_pairs', 'read_pairs',
                'no_match_reads'	
                'no_snp_reads', 'total_reads',
                'time_taken']	
                out_values = [phase_change_counter, unpaired, 
                        phase_change_mate_pair_counter, paired,	
                        no_match_counter, 
                        all_seq_counter - seq_with_snps_counter, all_seq_counter,
                        runtime]	
                f.write(','.join(fieldnames) + '\n')	
                f.write(','.join([str(n) for n in out_values]) + '\n')


if __name__ == '__main__':
    matepairs_recomb()
