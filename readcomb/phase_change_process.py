#!/usr/bin/env python3
"""
Readcomb filtering script utilizing multiprocessing
"""
import os
import sys
import argparse
import time
import datetime
import itertools
import pysam
from multiprocessing import Queue, Process
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

def check_snps(vcf_file_obj, chromosome, left_bound, right_bound):
    """
    Generate all SNPs on given chromosome and VCF file within the
    left_bound and right_bound using cyvcf2

    SNPs must have all calls with GQ >= 30 and no heterozygous calls.
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

    # list comp with if statement to only include SNPs
    with SilentVCF():
        records = [rec for rec in vcf_file_obj(region)]
    return records

def cigar(record):
    """
    Build the query segment using the cigar tuple given by the bam record so that it aligns with
    indexes of SNPs in the VCF and the reference sequence

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
    """
    Takes a segment and a list in the region of the segment and generates a list of
    strings to represent in order the parent of each variant on the segment

    snps: list of cyvcf2 variants in the region of the segment given by the function check_snps(),
        can include variants that are outside of the region when there is soft clipping
    segment: string sequence that is the segment built by the function cigar()
    record: bam record from pysam alignment file

    Returns a list of strings ('1', '2', or 'N')
    """

    snp_lst = []

    cigar_tuples = record.cigartuples

    # check for no alignment
    if not cigar_tuples:
        return snp_lst

    for snp in snps:
        # Using SNP.start and record.reference_start since they are both 0 based
        # SNP.start grabs vcf positions in 0 index while vcfs are 1 indexed
        idx = snp.start - record.reference_start

        parent1 = snp.gt_bases[0].split('/')[0]
        parent2 = snp.gt_bases[1].split('/')[0]

        if parent1 == parent2: # ignore uninformative SNPs
            continue
        
        # ignore if snp is before sequence
        if idx < 0:
            continue
        
        if snp.is_indel:
            # check if indel is outside of segment
            if idx + max(len(parent1), len(parent2)) > len(segment):
                continue

            parent1_match = segment[idx:idx + len(parent1)] == parent1
            parent2_match = segment[idx:idx + len(parent2)] == parent2

            if len(parent1) > len(parent2):
                if parent1_match:
                    snp_lst.append('1')
                elif parent2_match:
                    snp_lst.append('2')
                else:
                    snp_lst.append('N')
            else:
                if parent1_match:
                    snp_lst.append('1')
                elif parent2_match:
                    snp_lst.append('2')
                else:
                    snp_lst.append('N')
        else:
            if idx >= len(segment):
                break

            if segment[idx] == parent1:
                snp_lst.append('1')

            elif segment[idx] == parent2:
                snp_lst.append('2')

            else:
                snp_lst.append('N')

    return snp_lst

class Processor(Process):

    def __init__(self, input_queue, counter_queue, 
                writer_queue, vcf_file_name, args, **kwargs):
        Process.__init__(self, **kwargs)
        self.input_queue = input_queue
        self.counter_queue = counter_queue
        self.writer_queue = writer_queue
        self.vcf_file_obj = VCF(vcf_file_name)
        # header used for creating records
        self.header = pysam.AlignmentFile(args['bam'], 'r').header
        self.args = args
        self.counter = 0

    def run(self):     
        pair = self.input_queue.get(block=True)
        while pair:
            record_1 = pysam.AlignedSegment.fromstring(pair[0], self.header)
            record_2 = pysam.AlignedSegment.fromstring(pair[1], self.header)

            snps_1 = check_snps(self.vcf_file_obj, record_1.reference_name,
                                    record_1.reference_start,
                                    record_1.reference_start + record_1.query_alignment_length)
            snps_2 = check_snps(self.vcf_file_obj, record_2.reference_name,
                                    record_2.reference_start,
                                    record_2.reference_start + record_2.query_alignment_length)
            
            # updating counters
            self.counter += 1
            '''
            if self.counter % 1000 == 0:
                print(self.name + ' at ' 
                + str(self.counter) + ' iterations')
            '''

            self.counter_queue.put('seq')
            if len(snps_1) + len(snps_2) > 1: 
                self.counter_queue.put('seq_with_snps')
            # skip if not enough snps
            else:
                pair = self.input_queue.get(block=True)
                continue

            segment_1 = cigar(record_1)
            segment_2 = cigar(record_2)

            snp_lst = phase_detection(snps_1, segment_1, record_1) + phase_detection(snps_2, segment_2, record_2)

            if '1' in snp_lst and '2' in snp_lst:
                self.counter_queue.put('phase_change')

                if 'phase_change' in self.args['mode']:
                    self.writer_queue.put(pair)
                
            if 'N' in snp_lst:
                self.counter_queue.put('no_match')

                if 'no_match' in self.args['mode']:
                    self.writer_queue.put(pair)

            pair = self.input_queue.get(block=True)


    
class Counter(Process):
    def __init__(self, input_queue, **kwargs):
        Process.__init__(self, **kwargs)
        self.input_queue = input_queue
        self.counters = {
            'no_match': 0,
            'phase_change': 0,
            'seq_with_snps': 0,
            'seq': 0
            }

    def run(self):

        # start timer
        start = time.time()

        count = self.input_queue.get(block=True)
        while count:
            self.counters[count] += 1
            count = self.input_queue.get(block=True)
        
        # end timer
        end = time.time()
        runtime = str(datetime.timedelta(seconds=round(end - start)))

        # print counters to STDOUT
        print('{} phase change reads pairs from total {} read pairs'.format(
            self.counters['phase_change'], self.counters['seq']))
        print('{} reads had no-match variants'.format(self.counters['no_match']))
        print('{} reads did not have enough SNPs (> 0) to call'.format(
            self.counters['seq'] - self.counters['seq_with_snps']
        ))
        print('time taken: {}'.format(runtime))
        
        # TODO: add logging support


class Writer(Process):
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
                sp = str_record.split('\t')
                record = pysam.AlignedSegment.fromstring(str_record, 
                                                        self.header)
                try:
                    assert len(record.query_sequence) == len(record.qual)
                except:
                    raise ValueError(str_record + 'bonked')
                self.out.write(record)
            pair = self.input_queue.get(block=True)




def matepair_process():
    """
    """
    # dictionary of arguements
    args = arg_parser()

    # pysam argument object
    bam = pysam.AlignmentFile(args['bam'], 'r')

    print('Creating processes')
    # set up counter and writer
    count_input = Queue()
    counter = Counter(count_input, daemon=True)
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
                            write_input, args['vcf'], args, 
                            daemon=True, name='Process ' + str(i)))
        processes[i].start()

    prev_record = None
    process_idx = 0

    print('Processes created')

    for record in tqdm(bam):
        # check if record is in a proper pair
        if not record.is_proper_pair or record.is_secondary or record.is_supplementary:
            continue

        if not prev_record:
            prev_record = record
        else:
            if prev_record.reference_name != record.reference_name:
                raise ValueError('Pre-sorting of bam file went wrong')

            pair = [prev_record.to_string(), record.to_string()]

            input_queues[process_idx].put(pair)

            # give record to next process
            if process_idx < args['processes'] - 1:
                process_idx += 1
            else:
                process_idx = 0 

            # clear the pair
            prev_record = None

    print('Waiting for processes to complete')

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
    time.sleep(1)
    write_input.close()
    writer.join()

if __name__ == '__main__':
    matepair_process()
    
        

            


    

    



