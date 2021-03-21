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
from multiprocessing import Queue
from multiprocessing import Process
from cyvcf2 import VCF
from tqdm import tqdm

def arg_parser():
    parser = argparse.ArgumentParser(
        description='filter BAM for reads containing phase changes',
        usage='readcomb-filter [options]')

    parser.add_argument('-b', '--bam', required=True,
                        type=str, help='BAM to filter, required')

    parser.add_argument('-v', '--vcf', required=True,
                        type=str, help='VCF containing parents, required')

    parser.add_argument('-p', '--processes', required=False, type=int, default=4,
                        help='Number of processes to run readcomb filter on, default is 1')

    parser.add_argument('-m', '--mode', required=False, type=str, default='phase_change',
                        help='Read filtering mode (default phase_change) \
                             [phase_change|no_match]')

    parser.add_argument('-l', '--log', required=False,
                        type=str, help='Filename for log metric output [optional]')

    parser.add_argument('-o', '--out', required=False, type=str, default='recomb_diagnosis',
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
    # helper class for check_variants that silences cyvcf2 warnings
    def __enter__(self):
        self._original_stderr = sys.stderr
        sys.stderr = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stderr.close()
        sys.stderr = self._original_stderr

def check_variants(vcf_file_obj, chromosome, left_bound, right_bound):
    """
    Fetch variants from vcf file in given bounds.

    Parameters
    ----------
    vcf_file_obj : cyvcf2.VCF
        ``cyvcf2`` VCF file object
    chromosome : str
        reference sequence name of bound
    left_bound : int
        0-based index of sequence left bound
    right_bound : int
        0-based index of sequence right bound

    Returns
    -------
    list of cyvcf.Variant
    """

    # 1 is added to record.reference_start and the following parameter because vcf is 1 indexed
    # in order to keep code consistent
    region = '{c}:{l}-{r}'.format(c=chromosome, l=left_bound+1, r=right_bound+1)
    
    with SilentVCF():
        records = [rec for rec in vcf_file_obj(region)]
    return records

def cigar(record):
    """
    Realign the given query segment to its reference.

    Rebuild the sequence using the query sequence and cigar tuple from the bam record to realign
    it with the reference sequence prior to variant matching

    Parameters
    ----------
    record : pysam.AlignedSegment
        bam record from pysam alignment file

    Returns
    -------
    str
        a sequence of bases representing the realigned DNA sequence
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

    # terminology on consume query/reference from samtools SAM specifications
    for cigar_tuple in cigar_tuples:

        # 0 is a match, add it onto the segment
        # consume query and reference
        if cigar_tuple[0] == 0:
            segment.append(query_segment[index:index+cigar_tuple[1]])
            index += cigar_tuple[1]

        # 1 is an insertion in query segment
        # consume query by increasing index
        elif cigar_tuple[0] == 1:
            index += cigar_tuple[1]

        # 2 is an deletion to query segment
        # consume reference by adding gaps to segment
        elif cigar_tuple[0] == 2:
            segment.append('-' * cigar_tuple[1])

        # 4 is soft clipping
        # consume query by increasing index
        elif cigar_tuple[0] == 4:
            index += cigar_tuple[1]

        # 5 = hard clipping
        # does NOT consume query or reference
        elif cigar_tuple[0] == 5:
            continue

        else:
            # currently does not support padding (6),
            # sequence match (7), and sequence mismatch (8)
            raise ValueError(
                'No condition for tuple {} of {}'.format(
                    cigar_tuples.index(cigar_tuple), cigar_tuples))

    return ''.join(segment)

def phase_detection(variants, segment, record):
    """
    Detect haplotype of variants in the given DNA sequence.

    Takes a ``segment`` and a list of ``cyvcf2`` variants in the region of the segment 
    to detect the haplotype of the sequence and generate a list of {'1','2','N'} to represent it
    in the order of the variants on the segment

    Parameters
    ----------
    variants : list of cyvcf2.Variant
        ``cyvcf2`` variants in the bounds of ``segment`` given by ``check_variants()``
    segment : str
        realigned sequence built by ``cigar()``
    record : pysam.AlignedSegment
        bam record from pysam alignment file

    Returns
    -------
    list of {'1','2','N'}
        list of haplotypes of variants in the bounds of the segment
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

        # ignore uninformative variants
        if parent1 == parent2:
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

            # always check parent haplotype that is longer first as
            # most include the shorter haplotype thus the longer haplotype
            # is harder to match
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

        else: # variant is a SNP
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
    def __init__(self, input_queue, counter_queue, 
                writer_queue, args, **kwargs):
        """
        Initialize the read pair processor.

        Processor is inherited from ``multiprocessing.Process``; when started it receives pairs of
        bam sequences through the ``input_queue`` from the scheduler and does phase change analysis, outputting 
        through ``writer_queue`` bams that fit the user given arguement criteria to ``Writer``.
        
        When the scheduler has parsed and distributed all bam sequences, the ``Processor`` receives
        a ``None`` through the ``input_queue`` and ends its loop and prepares to join.

        Parameters
        ----------
        input_queue : multiprocessing.Queue
            input queue of bam read pairs in string form from scheduler
        counter_queue : multiprocessing.Queue
            output queue for result of phase change analysis
        writer_queue : multiprocessing.Queue
            output queue for bam read pairs in string form that fulfill the user given criteria
        args : dictionary of strings
            dictionary containing all user given arguements compiled by arg_parse()
        **kwargs 
            extra parameters passed onto super class ``multiprocessing.Queue()``
        """
        Process.__init__(self, **kwargs)
        self.input_queue = input_queue
        self.counter_queue = counter_queue
        self.writer_queue = writer_queue
        self.vcf_file_obj = VCF(args['vcf'])
        # header used for converting bam string to bam object
        self.header = pysam.AlignmentFile(args['bam'], 'r').header
        self.args = args

    def run(self):
        """
        Start the ``Processor``

        ``run()`` is automatically called when the scheduler runs ``Processor.start()``,
        this begins the operations of both the ``multiprocessing.Process()`` super class and
        ``Processor``. 
        
        See ``check_variants()``, ``cigar()``, and ``phase_detection`` for
        more information on the analysis process.
        """
        pair = self.input_queue.get(block=True)

        while pair:
            # use bam header to convert bam string back to pysam bam object
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
    def __init__(self, input_queue, args, **kwargs):
        """
        Intializes the ``Counter`` process.

        Inherited from ``multiprocessing.Process``, ``Counter`` receives sequence information
        from Processor processes through ``input_queue`` and tallies them to generate a ``tqdm``
        progress bar in ``STDOUT``. 
        
        After being passed a ``None`` from the ``Scheduler, ``Counter``
        either writes tallied results to ``STDOUT`` or a log depending on the user's arguements.
        
        Parameters
        ----------
        input_queue : multiprocessing.Queue
            get sequence information from all ``Processor`` processes
        args : dictionary of str
            dictionary containing all user given arguements compiled by arg_parse()
        **kwargs
            extra args that are passed onto super class ``Process``
        """
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
        """
        Start the ``Counter``

        ``run()`` is automatically called when the scheduler runs ``Counter.start()``,
        this begins the operations of both the ``multiprocessing.Process()`` super class and
        ``Counter``
        """

        # start timer
        start = time.time()

        # create tqdm iteration counter if no bars
        if 'pair_count' in self.args:
            self.progress = tqdm(total=self.args['pair_count'])
        else:
            self.progress = tqdm()

        count = self.input_queue.get(block=True)

        while count:
            self.counters[count] += 1
            # update iteration counter if sequence
            if count == 'seq':
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
        
        # write to log if argument defined
        if self.args["log"]:
            # determine if we are adding to file or creating new log file
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
    def __init__(self, input_queue, args, **kwargs):
        """
        Intializes the ``Writer`` process.

        Inherited from ``multiprocessing.Process``, ``Writer`` receives bam records
        in string form from ``Processor`` processes that fulfill the filters given by the user.
        ``Writer`` converts bam records in string form to ``pysam.AlignedSegment`` objects and
        writes them to a SAM file. When ``Writer`` receives a ``None`` from the Scheduler, the
        writing loop is ended and the output file is closed.
        

        Parameters
        ----------
        input_queue : multiprocessing.Queue
            queue to receive bam record strings from all ``Processor`` processes
        args : dictionary of str
            dictionary containing all user given arguements compiled by arg_parse()
        **kwargs
            extra args that are passed onto super class ``Process``
        """
        Process.__init__(self, **kwargs)
        self.input_queue = input_queue

        pysam_obj = pysam.AlignmentFile(args['bam'], 'r')
        self.out = pysam.AlignmentFile(args['out'] + '.sam', 'wh', 
                                        template=pysam_obj)
        self.header = pysam.AlignmentFile(args['bam'], 'r').header

    def run(self):
        """
        Start the ``Writer``

        ``run()`` is automatically called when the scheduler runs ``Writer.start()``,
        this begins the operations of both the ``multiprocessing.Process()`` super class and
        ``Writer``
        """
        pair = self.input_queue.get(block=True)
        while pair:
            for str_record in pair:
                record = pysam.AlignedSegment.fromstring(str_record, self.header)
                self.out.write(record)

            pair = self.input_queue.get(block=True)
        # very important to close or last piece of information
        # in buffer will not write to file
        self.out.close()




def matepair_process():
    """
    Start processes, queues, and distribute bam records to ``Processors``

    ``matepair_process`` calls ``arg_parse()`` to parse user given arguments.

    It then initializes and starts ``Processors``, ``Counter``, and ``Writer`` processes/daemons 
    and then equally distribute the bam records in string form to ``Processors``. 
    
    After all bam records have been distributed, ``None`` is passed to ``Processor`` 
    processes, ``Counter process, and ``Writer`` process to signal their stop and
    ``close()`` the processes.
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
    # daemon causes other sub-processes to auto shutdown if main process is interrupted
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
            
            # convert bam sequence to string so it can be pickled and sent in queue
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
        # pass None to let processes know all bam sequences have been scheduled
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
    
        

            


    

    


