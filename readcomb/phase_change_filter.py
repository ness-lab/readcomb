#!/usr/bin/env python3
"""
Main readcomb filtering script
"""
import os
import argparse
import time
import datetime
import threading
import itertools
import pysam
from cyvcf2 import VCF
from tqdm import tqdm

def arg_parser():
    """
    Parse command line args.
    """
    parser = argparse.ArgumentParser(
        description='filter BAM for reads containing phase changes',
        usage='python3.5 phase_change_filter.py [options]')

    parser.add_argument('-b', '--bam', required=True,
                        type=str, help='BAM to filter, required')

    parser.add_argument('-v', '--vcf', required=True,
                        type=str, help='VCF containing parents, required')

    parser.add_argument('-t', '--threads', required=False, type=int, default=1,
                        help='Number of threads to run readcomb filter on, default is 1')

    parser.add_argument('-c', '--chrom', required=False, type=str, default='all',
                        help='Chromosome sequences in the bam file to run (default all)')

    parser.add_argument('-m', '--mode', required=False, type=str, default='phase_change',
                        help='Read filtering mode (default phase_change) \
                             [phase_change|no_match|all]')

    parser.add_argument('-l', '--log', required=False,
                        type=str, help='Filename for log metric output [optional]')

    parser.add_argument('-o', '--out', required=False, type=str, default='recomb_diagnosis.sam',
                        help='File to write to (default recomb_diagnosis)')

    parser.add_argument('-i', '--improper', action='store_true',
                        help='Include unpaired/improper reads in filter results')

    args = parser.parse_args()

    return {
        "bam": args.bam,
        "vcf": args.vcf,
        "threads": args.threads,
        "chrom": args.chrom,
        "mode": args.mode,
        "log": args.log,
        "out": args.out,
        "improper": args.improper,
        }

def check_snps(vcf_file_obj, chromosome, left_bound, right_bound):
    """
    Generate all SNPs on given chromosome and VCF file within the
    left_bound and right_bound using cyvcf2

    SNPs must have all calls with GQ >= 30 and no heterozygous calls.
    (The heterozygous calls filter might need to be revisited for compatibility
    with diploid species)

    vcf_file_obj: cyvcf2 VCF file object
    chromosome: string: chromosome name
    left_bound/right_bound: integer: 0-based index of reference sequence

    Return list of cyvcf Variant objects
    """

    # 1 is added to record.reference_start and the following parameter because vcf is 1 indexed
    # in order to keep code consistent
    region = '{c}:{l}-{r}'.format(c=chromosome, l=left_bound+1, r=right_bound+1)

    # list comp with if statement to only include SNPs
    records = [rec for rec in vcf_file_obj(region) if rec.is_snp and len(rec.ALT) > 0 
               and rec.num_het == 0 and all(rec.gt_quals >= 30)
               and rec.gt_bases[0] != rec.gt_bases[1]]
    return records


def cache_pairs(bam_file_obj, args):
    """
    Iterates through a bam file to find mate pairs and cache them together in a dictionary

    bam_file_obj: pysam alignment file object
    chromosome: string

    Returns a dictionary of tuples (cache), the number of unpaired reads (unpaired),
    and the number of mate pairs (paired).

    Cache is a dictionary with a unique sequence read id as the key and a tuple pair of bam
    records as the value. If there is no mate pair, the second object in the tuple is None
    """

    print('Caching reads for ' + args['chrom'] + ' sequences')

    cache = {}

    paired = 0
    unpaired = 0

    for record in bam_file_obj:

        # filter out the read if it is in a improper read 
        # unless specified not to in arguement
        if not record.is_proper_pair and not args['improper']:
            continue

        # filter on the read if it is a secondary or supplementary mate pair
        if record.is_secondary or record.is_supplementary:
            continue

        # filter out read if it isn't the chromosome specified in arguement
        if args['chrom'] != record.reference_name and args['chrom'] != 'all':
            continue

        # check if query_name and reference_name exist
        if record.query_name is None or record.reference_name is None:
            continue

        # key is unique combination of query name and the chromosome
        name = record.query_name + record.reference_name

        if name not in cache:
            cache[name] = [record, None]
            unpaired += 1
        elif cache[name][1] is None:
            cache[name][1] = record
            paired += 2
            unpaired -= 1
        else:
            raise ValueError('More than 2 sequences for mate pairs ' + record.query_name \
                              + ' in chromosome ' + record.reference_name)

    print('Number of unpaired sequences: {}, read pairs: {}'.format(unpaired, paired))
    return cache, paired, unpaired

def cigar(record):
    '''
    Build the query segment using the cigar tuple given by the bam record so that it aligns with
    indexes of SNPs in the VCF and the reference sequence

    record: bam record from pysam alignment file

    Returns a dna sequence string
    '''
    cigar_tuples = record.cigartuples

    # initialize segment for building
    segment = []

    # index to keep track of where we are in the query_segment
    query_segment = record.query_sequence
    index = 0

    # no alignment, cigar_tuple is equal to None
    if cigar_tuples is None:
        return segment

    for cigar_tuple in cigar_tuples:

        # 0 is a match, add it onto the segment
        if cigar_tuple[0] == 0:
            segment.append(query_segment[index:index+cigar_tuple[1]])
            index += cigar_tuple[1]

        # 1 is an insertion to query segment
        # skip it because SNPs are aligned to reference and do not exist in this region
        elif cigar_tuple[0] == 1:
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
            raise ValueError(
                'No condition for tuple {} of {}'.format(
                    cigar_tuples.index(cigar_tuple), cigar_tuples))

    return ''.join(segment)

def phase_detection(snps, segment, record):
    '''
    Takes a segment and a list in the region of the segment and generates a list of
    strings to represent in order the parent of each variant on the segment

    snps: list of cyvcf2 variants in the region of the segment given by the function check_snps(),
        can include variants that are outside of the region when there is soft clipping
    segment: string sequence that is the segment built by the function cigar()
    record: bam record from pysam alignment file

    Returns a list of strings ('1', '2', or 'N')
    '''

    snp_lst = []

    cigar_tuples = record.cigartuples

    # check for no alignment
    if not cigar_tuples:
        return snp_lst

    for snp in snps:
        # Using SNP.start and record.reference_start since they are both 0 based
        # SNP.start grabs vcf positions in 0 index while vcfs are 1 indexed
        idx = snp.start - record.reference_start

        current_tuple = 0
        current_base = 0

        if idx < 0:
            raise ValueError('VCF indexing is off. Check SNP at {}'.format(snp))

        parent1 = snp.gt_bases[0][0]
        parent2 = snp.gt_bases[1][0]

        if idx >= len(segment):
            break

        if segment[idx] == parent1:
            snp_lst.append('1')

        elif segment[idx] == parent2:
            snp_lst.append('2')

        else:
            snp_lst.append('N')

    return snp_lst

class MatepairsThread(threading.Thread):
    """
    Inherited class from python's threading.Thread module. This class can be called
    to create a parallel thread that will filter cached bam sequences

    counters: dictionary of counters shared across all threads
    args: dictionary of arguments parsed from the command line call of this program
    pairs: list of dictionaries
    index: int, the index of a dictionary in pairs assigned to this thread

    The thread does not return anything but rather overwrites the pairs[index] dictionary
    assigned to this thread with the filtered results
    """
    def __init__(self, counters, args, pairs, index):
        # initiate threading.Thread super class
        threading.Thread.__init__(self)
        self.counters = counters
        self.args = args

        # cyvcf2 VCF file object creation
        self.vcf_file_obj = VCF(args["vcf"])

        # pairs is a list of dictionaries
        # self.pairs[self.index] is the dictionary assigned to this thread
        self.pairs = pairs
        self.index = index

        # output dictionary
        self.filtered_pairs = {}

    def run(self):
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
    """
    Main function of phase_change_filter.py

    Parses arguments, generates file objects, and creates threads for sequence filtering
    then waits for the thread to finish before writing the results to a file and outputting
    counters to the console and a log file if chosen
    """
    # start timer
    start = time.time()

    # dictionary of arguments
    args = arg_parser()

    # pysam alignment file object
    bam_file_obj = pysam.AlignmentFile(args["bam"], 'r')

    # pysam alignment file with input bam as filter
    f_obj = pysam.AlignmentFile(args["out"], 'wh', template=bam_file_obj)

    #counters
    counters = {"no_match": 0,
                "phase_change": 0,
                "phase_change_mate_pair": 0,
                "seq_with_snps": 0,
                "seq": 0}

    pairs, counters["paired"], counters["unpaired"] = cache_pairs(bam_file_obj, args)

    print('Beginning phase change analysis')

    split_pairs = []
    split_index = max(len(pairs) // args["threads"], 1)
    pairs_iter = iter(pairs.items())

    threads = []

    for thread_index in range(args["threads"]):
        # split pairs dictionary into multiple dictionaries
        # and store them in the list split_pairs to be used by threads
        if thread_index < args["threads"] - 1:
            split_pairs.append(dict(itertools.islice(pairs_iter, split_index)))
        else:
            split_pairs.append(dict(pairs_iter))

        # initiate and start threads
        threads.append(MatepairsThread(counters, args, split_pairs, thread_index))
        threads[thread_index].start()

    # end threads
    for thread in threads:
        thread.join()

    # write filtered data onto output sam file
    for pairs in split_pairs:
        for query_name in pairs:
            f_obj.write(pairs[query_name][0])

            if pairs[query_name][1]:
                f_obj.write(pairs[query_name][1])


    # end timer
    end = time.time()
    runtime = str(datetime.timedelta(seconds=round(end - start)))

    print('''
    Done.
    {} phase changes reads from {} total unpaired
    {} phase changes reads across mate pairs from {} total read pairs
    {} reads had no-match variants.
    {} reads did not have enough SNPs (> 0) to call
    time taken: {}
    '''.format(counters["phase_change"], counters["unpaired"],
               counters["phase_change_mate_pair"], counters["paired"],
               counters["no_match"], counters["seq"] - counters["seq_with_snps"],
               runtime)
         )

    if args["log"]:
        needs_header = True
        if os.path.isfile(args["log"]):
            needs_header = False
        with open(args["log"], 'a') as f:
            if needs_header:
                fieldnames = ['phase_change_reads', 'unpaired_reads',
                              'phase_change_across_mate_pair', 'read_pairs',
                              'no_match_reads' 'no_snp_reads', 'total_reads',
                              'time_taken']
                f.write(','.join(fieldnames) + '\n')
            out_values = [counters["phase_change"], counters["unpaired"],
                          counters["phase_change_mate_pair"], counters["paired"],
                          counters["no_match"], counters["seq"] - counters["seq_with_snps"],
                          counters["seq"], runtime]
            f.write(','.join([str(n) for n in out_values]) + '\n')


if __name__ == '__main__':
    matepairs_recomb()
