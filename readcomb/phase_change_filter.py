import os
import argparse
import pysam
import time
import datetime
import threading
import itertools
from cyvcf2 import VCF
from tqdm import tqdm

def arg_parser():
    parser = argparse.ArgumentParser(
        description='filter BAM for reads containing phase changes', 
        usage='python3.5 phase_change_filter.py [options]')

    parser.add_argument('-b', '--bam', required=True,
                        type=str, help='BAM to filter, required')
                        
    parser.add_argument('-v', '--vcf', required=True,
                        type=str, help='VCF containing parents, required')

    parser.add_argument('-t', '--threads', required=False,
                        type=int, default=1, help='Number of threads to run readcomb filter on, default is 1')

    parser.add_argument('-c', '--chrom', required=False,
                        type=str, default='all', help='Specify which chromsome sequences to run, default is all')

    parser.add_argument('-m', '--mode', required=False,
                        type=str, default='phase_change', help='Mode to execute the program, default is phase_change')

    parser.add_argument('-l', '--log', required=False,
                        type=str, help='Log metrics to provide filename, default is False')

    parser.add_argument('-o', '--out', required=False,
                        type=str, default='recomb_diagnosis', help='File to write to, default is recomb_diagnosis')

    args = parser.parse_args()

    return {"bam": args.bam,
            "vcf": args.vcf,
            "threads": args.threads,
            "chrom": args.chrom,
            "mode": args.mode,
            "log": args.log,
            "out": args.out,}

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
      
        # 0 is a match, add it onto the segment 
        if cigar_tuple[0] == 0:
            segment.append(query_segment[index:index+cigar_tuple[1]])
            index += cigar_tuple[1]

        # 1 is an insertion to query segment, skip it because SNPs are aligned to reference and do not exist in this region
        elif cigar_tuple[0] == 1:
            segment.append(query_segment[index:index+cigar_tuple[1]])
            index += cigar_tuple[1]               
            
        # 2 is an deletion to query segment, add a gap to realign it to reference
        # don't add index to not moving further through the query segment
        elif cigar_tuple[0] == 2:
            segment.append('-' * cigar_tuple[1])

        # 4 is soft clipping, query sequence includes this but is not aligned to reference
        elif cigar_tuple[0] == 4:
            index += cigar_tuple[1]

        # 5 = hard clipping, record.query_sequence does not have the portion that is
        # hard clipping so skip it and don't add anything
        elif cigar_tuple[0] == 5:
            continue

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

class matepairs_thread(threading.Thread):
    def __init__(self, counters, args, pairs, index):
        threading.Thread.__init__(self)
        self.counters = counters
        self.args = args
        # cyvcf2 VCF file reader
        self.vcf_file_obj = VCF(args["vcf"])
        self.pairs = pairs
        self.index = index

        self.filtered_pairs = {}

    def run(self):
        
        # self.pairs is a list of dictionaries and we use self.index to get the 
        # subset of cached pairs assigned to this thread
        for query_name in tqdm(self.pairs[self.index]):
            snp_lst = []
            pair = self.pairs[self.index][query_name]
            for record in pair:
                if record:
                    self.counters["seq"] += 1
                    
                    segment = cigar(record)
                                
                    snps = check_snps(self.vcf_file_obj, record.reference_name, 
                                    record.reference_start,
                                    record.reference_start + record.query_alignment_length)
                    
                    snp_lst += phase_detection(snps, segment, record)

                    if len(snps) > 0:
                        self.counters["seq_with_snps"] += 1
                
                    

            if '1' in snp_lst and '2' in snp_lst:
                
                if pair[1]:
                    self.counters["phase_change_mate_pair"] += 1
                else:
                    self.counters["phase_change"] += 1

                if 'phase_change' in self.args["mode"]:                   
                    self.filtered_pairs[query_name] = pair
            
            if 'N' in snp_lst:
                self.counters["no_match"] += 1

                if 'no_match' in self.args["mode"]:                        
                    self.filtered_pairs[query_name] = pair

        self.pairs[self.index] = self.filtered_pairs

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

    # dictionary of arguements
    args = arg_parser()

    # pysam alignment file object
    bam_file_obj = pysam.AlignmentFile(args["bam"], 'r')
    
    # pysam alignment file with input bam as filter
    f_obj = pysam.AlignmentFile(args["out"] + '.sam', 'wh', template=bam_file_obj)
    
    #counters
    counters = {"no_match": 0,
                "phase_change": 0,
                "phase_change_mate_pair": 0,
                "seq_with_snps": 0,
                "seq": 0}
    
    pairs, counters["paired"], counters["unpaired"] = cache_pairs(bam_file_obj, args["chrom"])

    print('Beginning phase change analysis')

    split_pairs = []
    split_index = len(pairs) // args["threads"]
    pairs_iter = iter(pairs.items())

    threads = []

    
    for thread_index in range(args["threads"]):
        # split pairs dictionary into multiple dictionaries and store them in the list split_pairs to be used by threads
        if thread_index < args["threads"] - 1:
            split_pairs.append(dict(itertools.islice(pairs_iter, split_index)))
        else:
            split_pairs.append(dict(pairs_iter))

        # initiate and start threads
        threads.append(matepairs_thread(counters, args, split_pairs, thread_index))
        threads[thread_index].start()

    # end threads
    for thread in threads:
        thread.join()

    # write filtered data onto output sam file
    for pairs in split_pairs:
        for query_name in pairs:
            f_obj.write(pairs[query_name][0])

            if pairs[query_name][1]:
                f_obj.write(pairs[query_name][0])


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
    '''.format(counters["phase_change"], counters["unpaired"], round(counters["phase_change"] / counters["unpaired"] * 100, 2), 
        counters["phase_change_mate_pair"], counters["paired"], round(counters["phase_change_mate_pair"] / counters["paired"] * 100, 2),
        counters["no_match"], counters["seq"] - counters["seq_with_snps"],
        round((counters["seq"] - counters["seq_with_snps"]) / counters["seq"] * 100, 2),
        runtime)
    )

    if args["log"]:
        needs_header = True
        if os.path.isfile(log):	
            needs_header = False	
        with open(log, 'a') as f:	
            if needs_header:	
                fieldnames = ['phase_change_reads', 'unpaired_reads',	
                'phase_change_across_mate_pair', 'read_pairs',
                'no_match_reads'	
                'no_snp_reads', 'total_reads',
                'time_taken']
                f.write(','.join(fieldnames) + '\n')		
            out_values = [counters["phase_change"], counters["unpaired"], 
                    counters["phase_change_mate_pair"], counters["paired"],	
                    counters["no_match"], 
                    counters["seq"] - counters["seq_with_snps"], counters["seq"],
                    runtime]	
            f.write(','.join([str(n) for n in out_values]) + '\n')


if __name__ == '__main__':
    matepairs_recomb()
