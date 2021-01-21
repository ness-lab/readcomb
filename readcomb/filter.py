#!/usr/bin/env python3
"""
Readcomb filtering script utilizing multiprocessing
"""
import os
import sys
import argparse
import time
import datetime
import pysam
import subprocess
import pandas as pd
from io import BytesIO
from multiprocessing import Queue, Process
from cyvcf2 import VCF
from tqdm import tqdm

def arg_parser():
    """
    Parse command line args.
    """
    parser = argparse.ArgumentParser(
        description='filter BAM for reads containing phase changes',
        usage='readcomb-filter [options]')

    parser.add_argument('-b', '--bam', required=True,
                        type=str, help='BAM to filter, required')

    parser.add_argument('-v', '--vcf', required=True,
                        type=str, help='VCF containing parents, required')

    parser.add_argument('-p', '--processes', required=False, type=int, default=1,
                        help='Number of processes to run readcomb filter on, default is 1')

    parser.add_argument('-m', '--mode', required=False, type=str, default='phase_change',
                        help='Read filtering mode (default phase_change) \
                             [phase_change|no_match]')

    parser.add_argument('-l', '--log', required=False,
                        type=str, help='Filename for log metric output [optional]')

    parser.add_argument('-o', '--out', required=False, type=str, default='recomb_diagnosis.sam',
                        help='File to write to (default recomb_diagnosis)')

    args = parser.parse_args()

    return {
        'bam': args.bam,
        'vcf': args.vcf,
        'processes': args.processes,
        'mode': args.mode,
        'log': args.log,
        'out': args.out,
        }

class SilentVCF:
    def __enter__(self):
        self._original_stderr = sys.stderr
        sys.stderr = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stderr.close()
        sys.stderr = self._original_stderr

def check_variants(vcf_file_obj, chromosome, left_bound, right_bound):
    """
    Generate all variants on given chromosome and VCF file within the
    left_bound and right_bound using cyvcf2

    Variants must have all calls with GQ >= 30 and no heterozygous calls.
    Please use the filter script to preprocess parental VCF prior to
    phase change detection.

    vcf_file_obj: cyvcf2 VCF file object
    chromosome: string: chromosome name
    left_bound/right_bound: integer: 0-based index of reference sequence

    Return list of cyvcf Variant objects
    """

    # 1 is added to record.reference_start and the following parameter because vcf is 1 indexed
    # in order to keep code consistent
    region = '{c}:{l}-{r}'.format(c=chromosome, l=left_bound+1, r=right_bound+1)

    with SilentVCF():
        records = [rec for rec in vcf_file_obj(region)]
    return records

def cigar(record):
    """
    Build the query segment using the cigar tuple given by the bam record so that it aligns with
    indexes of variants in the VCF and the reference sequence

    record: bam record from pysam alignment file

    Returns a dna sequence string
    """
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
        # skip it because variants are aligned to reference and do not exist in this region
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

def phase_detection(variants, segment, record):
    """
    Takes a segment and a list in the region of the segment and generates a list of
    strings to represent in order the parent of each variant on the segment

    variants: list of cyvcf2 variants in the region of the segment given by the function check_variants(),
        can include variants that are outside of the region when there is soft clipping
    segment: string sequence that is the segment built by the function cigar()
    record: bam record from pysam alignment file

    Returns a list of strings ('1', '2', or 'N')
    """

    variant_lst = []

    cigar_tuples = record.cigartuples

    # check for no alignment
    if not cigar_tuples:
        return variant_lst

    for variant in variants:
        # Using variant.start and record.reference_start since they are both 0 based
        # variant.start grabs vcf positions in 0 index while vcfs are 1 indexed
        idx = variant.start - record.reference_start

        parent1 = variant.gt_bases[0].split('/')[0]
        parent2 = variant.gt_bases[1].split('/')[0]

        if parent1 == parent2: # ignore uninformative variants
            continue
        
        # ignore if variant is before sequence
        if idx < 0:
            continue
        
        if variant.is_indel:
            # check if indel is outside of segment
            if idx + max(len(parent1), len(parent2)) > len(segment):
                continue

            parent1_match = segment[idx:idx + len(parent1)] == parent1
            parent2_match = segment[idx:idx + len(parent2)] == parent2

            if len(parent1) > len(parent2):
                if parent1_match:
                    variant_lst.append('1')
                elif parent2_match:
                    variant_lst.append('2')
                else:
                    variant_lst.append('N')
            else:
                if parent2_match:
                    variant_lst.append('2')
                elif parent1_match:
                    variant_lst.append('1')
                else:
                    variant_lst.append('N')
        else:
            if idx >= len(segment):
                break

            if segment[idx] == parent1:
                variant_lst.append('1')

            elif segment[idx] == parent2:
                variant_lst.append('2')

            else:
                variant_lst.append('N')

    return variant_lst

class Processor(Process):
    """
    Inherited class of multiprocessing process that takes in bam sequences from main
    scheduler and does phase change analysis, outputting bams that fit the arguement
    criteria to the writer process

    input_queue -- multiprocessing.Queue(), get bam sequences in string form from main
                    scheduler
    counter_queue -- multiprocessing.Queue(), put sequence information to be tallied
    writer_queue -- multiprocessing.Queue(), put sequences that fit the argument
                    criteria to the writer process
    args -- dictionary of args from arg_parse()
    **kwargs -- extra arguements passed onto the super class process
    """

    def __init__(self, input_queue, counter_queue, 
                writer_queue, args, **kwargs):
        Process.__init__(self, **kwargs)
        self.input_queue = input_queue
        self.counter_queue = counter_queue
        self.writer_queue = writer_queue
        self.vcf_file_obj = VCF(args['vcf'])
        # header used for creating records
        self.header = pysam.AlignmentFile(args['bam'], 'r').header
        self.args = args

    def run(self):
        pair = self.input_queue.get(block=True)

        while pair:
            record_1 = pysam.AlignedSegment.fromstring(pair[0], self.header)
            record_2 = pysam.AlignedSegment.fromstring(pair[1], self.header)

            variants_1 = check_variants(self.vcf_file_obj, record_1.reference_name,
                                    record_1.reference_start,
                                    record_1.reference_start + record_1.query_alignment_length)
            variants_2 = check_variants(self.vcf_file_obj, record_2.reference_name,
                                    record_2.reference_start,
                                    record_2.reference_start + record_2.query_alignment_length)

            self.counter_queue.put('seq')
            if len(variants_1) + len(variants_2) > 1: 
                self.counter_queue.put('seq_with_variants')
            # skip if not enough variants
            else:
                pair = self.input_queue.get(block=True)
                continue

            segment_1 = cigar(record_1)
            segment_2 = cigar(record_2)

            variant_lst = phase_detection(variants_1, segment_1, record_1) + phase_detection(variants_2, segment_2, record_2)

            if '1' in variant_lst and '2' in variant_lst:
                self.counter_queue.put('phase_change')

                if 'phase_change' in self.args['mode']:
                    self.writer_queue.put(pair)
                
            if 'N' in variant_lst:
                self.counter_queue.put('no_match')

                if 'no_match' in self.args['mode']:
                    self.writer_queue.put(pair)

            pair = self.input_queue.get(block=True)


    
class Counter(Process):
    """
    Inherited class of multiprocessing process that receives sequence information
    from Processor processes and tallies them to print to STDOUT or a log file

    input_queue -- multiprocessing.Queue() get sequence information from all Processor processes
    args -- dictionary of arguements from arg_parse()
    **kwargs -- extra args that is passed onto super class process
    """

    def __init__(self, input_queue, args, **kwargs):
        Process.__init__(self, **kwargs)
        self.input_queue = input_queue
        self.counters = {
            'no_match': 0,
            'phase_change': 0,
            'seq_with_variants': 0,
            'seq': 0
            }
        self.args = args

    def run(self):

        # start timer
        start = time.time()

        # create tqdm iteration counter if no bars
        if 'pair_count' in self.args:
            self.progress = tqdm(total=self.args['pair_count'])

        count = self.input_queue.get(block=True)

        while count:
            self.counters[count] += 1
            # update iteration counter if no bars
            if count == 'seq' and 'pair_count' in self.args:
                self.progress.update(n=1)
            
            count = self.input_queue.get(block=True)
        
        self.progress.close()

        # end timer
        end = time.time()
        runtime = str(datetime.timedelta(seconds=round(end - start)))

        # print counters to STDOUT
        print('{} phase change reads pairs from total {} read pairs'.format(
            self.counters['phase_change'], self.counters['seq']))
        print('{} reads had no-match variants'.format(self.counters['no_match']))
        print('{} reads did not have enough variants (> 0) to call'.format(
            self.counters['seq'] - self.counters['seq_with_variants']
        ))
        print('time taken: {}'.format(runtime))
        
        if self.args["log"]:
            # determine if we are adding to file or 
            needs_header = False if os.path.isfile(self.args['log']) else True

            with open(self.args["log"], 'a') as f:
                if needs_header:
                    fieldnames = ['time',
                                'phase_change_reads',
                                'no_match_reads', 
                                'seq_with_variants', 
                                'no_variant_reads', 
                                'seqs',
                                'time_taken']
                    fieldnames.extend(sorted(self.args.keys()))
                    f.write(','.join(fieldnames) + '\n')

                out_values = [datetime.datetime.now(), 
                            self.counters['phase_change'], 
                            self.counters['no_match'],
                            self.counters['seq_with_variants'], 
                            self.counters['seq'] - self.counters['seq_with_variants'],
                            self.counters["seq"],
                            runtime]
                out_values.extend([self.args[key] for key in sorted(self.args.keys())])
                f.write(','.join([str(n) for n in out_values]) + '\n')


class Writer(Process):
    """
    Inherited class of multiprocessing process that receives bam sequences that fit the
    critera of the arguements and writes them to a new pysam alignment file

    input_queue -- multiprocessing.Queue() gets filtered bam sequences from Processor processes
    args -- dictionary of arguements from arg_parse
    **kwargs -- extra args that is passed onto super class process
    """
    def __init__(self, input_queue, args, **kwargs):
        Process.__init__(self, **kwargs)
        self.input_queue = input_queue

        pysam_obj = pysam.AlignmentFile(args['bam'], 'r')
        self.out = pysam.AlignmentFile(args['out'], 'wh', 
                                        template=pysam_obj)
        self.header = pysam.AlignmentFile(args['bam'], 'r').header

    def run(self):
        pair = self.input_queue.get(block=True)
        while pair:
            for str_record in pair:
                record = pysam.AlignedSegment.fromstring(str_record, self.header)
                self.out.write(record)

            pair = self.input_queue.get(block=True)
        self.out.close()




def matepair_process():
    """
    Main process that manages Processors, Counter, and Writer processes and then divides
    up the bam file sequences equally into the processors to analyze.

    Gets arguements from command line parsing
    """
    # dictionary of arguements
    args = arg_parser()

    # idxstats on bam/bai file
    if not os.path.isfile(args['bam'] + '.bai'):
        print('Bai file not found, continuing without progress bars')
    else:
        stats = subprocess.check_output(['samtools', 'idxstats', args['bam']])
        stats_table = pd.read_csv(BytesIO(stats), sep='\t', 
                                    names=['chrom', 'length', 'map_reads', 'unmap_reads'])
        pairs_sum = stats_table['map_reads'].sum()
        args['pair_count'] = int(pairs_sum / 2) # divide 2 to get pairs

        if pairs_sum % 2 != 0:
            raise ValueError('Preprocessing of bam file went wrong')


    # pysam argument object
    bam = pysam.AlignmentFile(args['bam'], 'r')

    print('Creating processes')
    # set up counter and writer
    count_input = Queue()
    counter = Counter(count_input, args, daemon=True)
    counter.start()

    write_input = Queue()
    writer = Writer(write_input, args, daemon=True)
    writer.start()

    processes = []
    input_queues = []
    # set up processes
    for i in range(args['processes']):
        input_queue = Queue()
        input_queues.append(input_queue)
        processes.append(Processor(input_queue, count_input, 
                            write_input, args, 
                            daemon=True, name='Process ' + str(i)))
        processes[i].start()

    prev_record = None
    process_idx = 0

    print('Processes created')

    for record in bam:
        # check if record is in a proper pair
        if not record.is_proper_pair or record.is_secondary or record.is_supplementary:
            raise ValueError('Preprocessing of bam file went wrong')

        if not prev_record:
            prev_record = record
        else:
            if prev_record.query_name != record.query_name:
                raise ValueError('Preprocessing of bam file went wrong')

            pair = [prev_record.to_string(), record.to_string()]

            input_queues[process_idx].put(pair)

            # give record to next process
            if process_idx < args['processes'] - 1:
                process_idx += 1
            else:
                process_idx = 0 

            # clear the pair
            prev_record = None

    # shut down processes
    for i in range(len(processes)):
        input_queues[i].put(None)
        input_queues[i].close()
        processes[i].join()

    # shut down counter and writer
    count_input.put(None)
    count_input.close()
    counter.join()

    write_input.put(None)
    write_input.close()
    writer.join()

if __name__ == '__main__':
    matepair_process()
    
        

            


    

    


