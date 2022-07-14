"""
Script to filter out false positives caused by misalignment and/or
sketchy regions of the genome

Three options:
1. Midpoint - creates/uses bed file of false positive phase changes to compare
against
2. Overlap - compares all variants shared by false positive read and putative
recombinant read
3. Nuclear - removes any read that overlaps a false positive read at all

"""

import os
import sys
import csv
import argparse
import ast
import subprocess
import datetime

import pysam
from tqdm import tqdm

try:
    from readcomb.classification import pairs_creation
    from readcomb.filter import ReadcombParser
except ImportError as e:
    print('WARNING: readcomb is not installed')
    print('Use command: pip install readcomb')
    from classification import pairs_creation
    from filter import ReadcombParser 

def arg_parser():
    parser = ReadcombParser(
        description='filter putative rcmb events for alignment-derived false positives', 
        usage='readcomb-fp [options]')

    parser.add_argument('-f', '--fname', required=True, type=str, 
        help='Cross SAM output by readcomb-filter')
    parser.add_argument('-fp', '--false_plus', required=False, type=str, 
        help='False positives from plus parent')
    parser.add_argument('-fm', '--false_minus', required=False, type=str, 
        help='False positives from minus parent')
    parser.add_argument('-v', '--vcf', required=True, type=str, 
        help='VCF with calls for cross of interest')
    parser.add_argument('-m', '--method', required=True, type=str, 
        help='Filtering method to use (midpoint/overlap/nuclear) - see docs')
    parser.add_argument('--false_bed_out', required=False, type=str, 
        help='Where to save false phase change bed file [optional]')
    parser.add_argument('--false_bed_in', required=False, type=str, 
        help='Bed file of false phase changes to use, if already generated [optional]')
    parser.add_argument('--log', required=True, type=str, default='fp.log',
        help='Path to log file')
    parser.add_argument('-o', '--out', required=True, type=str,
        help='File to write filtered reads to')
    parser.add_argument('--version', action='version', version='readcomb 0.3.15')

    return parser

class BedGenerator():
    def __init__(self, false_plus, false_minus, vcf_fname, method, bed_fname):
        """Generate BED file of false positive midpoints

        BED format will depend on what filtering method is requested -

        - ``midpoint`` filtering will generate rows where positions correspond to the two
          variants that denote a phase change
        - ``overlap`` filtering will generate rows that also contain the output of
          `classification.Pair.detection()` - lists of tuples denoting sites and haps
        - ``nuclear`` filtering will only include positions of reads - overlapping
          reads will have their reads reduced

        Parameters
        ----------
        false_plus : str
            path to readcomb-filter output SAM file when run on parent 1
        false_minus : str
            path to readcomb-filter output SAM file when run on parent 2
        vcf_fname : str
            path to VCF used to call recombination events
        method : str
            one of 'midpoint', 'overlap', or 'nuclear'
        bed_fname : str
            path to bed outfile
        """
        # paths to the plus and minus false positive files
        self.false_plus = false_plus
        self.false_minus = false_minus
        self.vcf_fname = vcf_fname

        self.method = method
        self.bed_fname = bed_fname
        if not self.bed_fname: # if fname not provided - generate temp bed
            self.temp = True
            self.bed_fname = 'readcomb.parental_fp.bed.temp'
        
    def generate_bed(self):
        """main method to invoke specific BED creation methods
        """
        if self.method == 'midpoint':
            self._generate_midpoint_bed()
        elif self.method == 'overlap':
            self._generate_overlap_bed()
        elif self.method == 'nuclear':
            self._generate_nuclear_bed()
        else:
            raise ValueError(
                f'method must be one of midpoint, overlap, or nuclear, not {self.method}')

    def _generate_midpoint_bed(self):
        """generate BED file for midpoint method

        the format is as follows:
        #chrom	start	end	hap1	hap2	parental_file

        where the start and end columns correspond to the two closest variants
        surrounding a phase change, and parental_file denotes whether this
        false positive was in the plus or minus sequence
        """
        print('[readcomb] generating bed for midpoint method filtering...')
        plus_reader = pairs_creation(self.false_plus, self.vcf_fname)
        minus_reader = pairs_creation(self.false_minus, self.vcf_fname)

        f = open(self.bed_fname, 'w')
        fieldnames = ['#chrom', 'start', 'end', 'hap1', 'hap2', 'parental_file']
        bed_writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
        bed_writer.writeheader()

        # iterate through false positives
        for i, reader in enumerate([plus_reader, minus_reader]):
            for pair in tqdm(reader, desc=f'parental false positive bed generation {i+1}'):
                pair.classify(masking=0)
                if pair.call == 'no_phase_change': # ignore if no phase change
                    continue
                else:
                    detection_filt = [
                        variant for variant in pair.detection 
                        if variant[0] != 'N']
                    prev_hap, prev_pos = detection_filt[0][:2] # set to first hap
                    for variant in detection_filt:
                        hap, pos, base = variant
                        if hap == prev_hap:
                            prev_pos = pos # update variable
                            continue
                        elif hap != prev_hap:
                            bed_writer.writerow({
                                '#chrom': pair.rec_1.reference_name, 
                                'start': prev_pos, 'end': pos, 
                                'hap1': prev_hap, 'hap2': hap, 
                                'parental_file': 'plus' if i == 0 else 'minus'}
                                )
                            prev_hap = hap
                            prev_pos = pos
        f.close()
    

    def _generate_overlap_bed(self):
        """generate BED file for overlap method

        the format is as follows:
        #chrom	start	end	parental_file	call	detection

        where the detection column contains the literal output of
        ``Pair.detection``, for comparison with putative recombinant pairs of
        interest in full, while call contains the recombination event the false
        positive was called as
        """
        print('[readcomb] generating bed for overlap method filtering...')
        plus_reader = pairs_creation(self.false_plus, self.vcf_fname)
        minus_reader = pairs_creation(self.false_minus, self.vcf_fname)

        f = open(self.bed_fname, 'w')
        fieldnames = ['#chrom', 'start', 'end', 'parental_file', 'call', 'detection']
        bed_writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
        bed_writer.writeheader()

        # iterate through false positives
        for i, reader in enumerate([plus_reader, minus_reader]):
            for pair in tqdm(reader, desc=f'parental false positive bed generation {i+1}/2'):
                pair.classify(masking=0)
                if pair.call == 'no_phase_change': # ignore if no phase change
                    continue
                else:
                    bed_writer.writerow({
                        '#chrom': pair.rec_1.reference_name,
                        'start': pair.condensed[0][1],
                        'end': pair.condensed[-1][2],
                        'parental_file': 'plus' if i == 0 else 'minus',
                        'call': pair.call, 'detection': pair.detection}
                        )
        f.close()


    def _generate_nuclear_bed(self):
        """generate bed file for nuclear method

        the format is as follows:
        #chrom	start	end

        this method just stores intervals overlapped by false positive reads in
        either parent - which means the bed first needs to be generated, and
        then 'reduced'

        the reduction is done by calling on self.sort_tabix_bed without
        actually running tabix, and using the temp sorted file to reduce
        overlapping intervals. the final output file is then named the original
        filename provided to self.bed_fname
        """
        print('[readcomb] generating bed for nuclear method filtering...')
        plus_reader = pairs_creation(self.false_plus, self.vcf_fname)
        minus_reader = pairs_creation(self.false_minus, self.vcf_fname)

        f = open(self.bed_fname, 'w')
        fieldnames = ['#chrom', 'start', 'end']
        bed_writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
        bed_writer.writeheader()
        
        # iterate through false positives
        for i, reader in enumerate([plus_reader, minus_reader]):
            for pair in tqdm(reader, desc=f'parental false positive bed generation {i+1}/2'):
                pair.classify(masking=0)
                if pair.call == 'no_phase_change': # ignore if no phase change
                    continue
                else:
                    bed_writer.writerow({
                        '#chrom': pair.rec_1.reference_name,
                        'start': pair.condensed[0][1],
                        'end': pair.condensed[-1][2]}
                        )
        f.close()

        # reduce
        print('[readcomb] initial bed generated - reducing...')
        self.sort_tabix_bed(tabix_bed=False) # creates sorted bed without tabixing
        with open(self.bed_fname + 'reduced', 'w') as f_out:
            writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=fieldnames)
            writer.writeheader()
            with open(self.bed_fname, 'r') as f_in:
                reader = csv.reader(f_in, delimiter='\t')
                _ = next(reader) # skip header
                prev_chrom, left_bound, prev_end = next(reader)
                left_bound, prev_end = int(left_bound), int(prev_end)

                # iterate through sorted lines
                for line in tqdm(reader, desc='reducing overlapping intervals'):
                    chrom = line[0]
                    start, end = [int(n) for n in line[1:]]
                    if not prev_chrom:
                        prev_chrom = chrom

                    # new chr started on line
                    if chrom != prev_chrom:
                        writer.writerow({
                            '#chrom': prev_chrom,
                            'start': left_bound,
                            'end': prev_end})
                        left_bound = start
                        prev_end = end
                        prev_chrom = None
                        continue
                    # current read overlaps with previous - extend
                    elif start <= prev_end <= end: 
                        prev_end = end
                    # current read does not overlap - write prev interval
                    elif prev_end < start: 
                        writer.writerow({
                            '#chrom': chrom,
                            'start': left_bound,
                            'end': prev_end})
                        left_bound = start
                        prev_end = end

        # fix up names
        os.remove(self.bed_fname)
        os.rename(self.bed_fname + 'reduced', self.bed_fname)


    def sort_tabix_bed(self, tabix_bed=True):
        """sort, uniq, and tabix bed file

        invokes subprocess for all three steps. the same method can be used on
        all three bed types

        tabix arg can be set to false if only sorting (used in nuclear bed
        generation)
        """
        # 'cat {in} | (sed -u 1q; sort -k1,1 -k2n,3n) | uniq > {out}'
        # write header to separate file
        header_file = open('header.temp', 'w')
        header_proc = subprocess.run(
            ['head', '-n', '1', self.bed_fname], stdout=header_file)
        header_file.close()

        # sort + uniq rest of file
        sorted_file = open('sorted.temp', 'w')
        header_remove_proc = subprocess.Popen(['tail', '-n', '+2', self.bed_fname], 
            stdout=subprocess.PIPE)
        sort_proc = subprocess.Popen(['sort', '-k1,1', '-k2n,3n'], 
            stdin=header_remove_proc.stdout,
            stdout=subprocess.PIPE)
        uniq_proc = subprocess.run(['uniq'], 
            stdin=sort_proc.stdout, stdout=sorted_file)
        sorted_file.close()

        # combine header and sorted file
        final_bed = open(self.bed_fname, 'w')
        cat_proc = subprocess.run(['cat', 'header.temp', 'sorted.temp'],
            stdout=final_bed)
        final_bed.close()
        print('[readcomb] created bed file')
        print('[readcomb] clearing temp files...')

        # clear temp files
        os.remove('header.temp')
        os.remove('sorted.temp')

        # bgzip and tabix file
        if tabix_bed:
            print('[readcomb] bgzip and tabix on bed file...')
            bgzip_proc = subprocess.run(['bgzip', self.bed_fname])
            tabix_proc = subprocess.run(['tabix', self.bed_fname + '.gz'])
            print('[readcomb] complete.')


class FalsePositiveFilterer():
    def __init__(self, fname, method, vcf_fname, tabix_reader, out, log):
        """Use generated BED file to filter out false positive reads

        Three available methods:

        - ``midpoint`` filtering will filter out putative recombinant reads that
          share at least one exact phase change with a false positive parental read
        - ``overlap`` filtering will compare variants common to a putative
          recombinant read and all overlapping false positive parental reads.
          even one instance of shared haplotypes will cause read filtering
        - ``nuclear`` filtering will remove all putative recombinant reads that
          have any overlap at all with a false positive read, regardless of the
          false positive's haplotype

        Parameters
        ----------
        fname : str
            path to putative recombinant reads from readcomb-filter
        method : str
            one of 'midpoint', 'overlap', or 'nuclear'
        vcf_fname : str
            path to VCF used to call recombination events
        tabix_reader : pysam.libctabix.TabixFile
            path to bed outfile
        out : str
            path to output SAM to write filtered reads to
        log : str
            path to log file
        """
        self.fname = fname
        self.method = method
        self.vcf_fname = vcf_fname
        self.tabix_reader = tabix_reader
        self.bed_fname = tabix_reader.filename.decode('utf-8')
        self.log = log

        # instantiate writer
        bam_template = pysam.AlignmentFile(self.fname, 'r')
        self.writer = pysam.AlignmentFile(out, 'wh', template=bam_template)

    def filter_false_positives(self):
        print(f'[readcomb] filtering {self.fname} with {self.method} method')
        print(f'[readcomb] using {self.bed_fname} as false positive lookup')

        # set up log
        needs_header = bool(not os.path.isfile(self.log))
        log_f = open(self.log, 'a', newline='')
        log_writer = csv.DictWriter(
            log_f, delimiter='\t', 
            fieldnames=[
                'time', 'total_pairs', 'pairs_kept', 
                'pairs_removed', 'no_phase_change', 
                'perc_kept', 'method', 'fname'])
        if needs_header:
            log_writer.writeheader()
        log_row = {
            'time': datetime.datetime.now(), 'total_pairs': 0, 'pairs_kept': 0, 
            'pairs_removed': 0, 'no_phase_change': 0, 'perc_kept': 0.0, 
            'method': self.method, 'fname': self.fname}

        # filter
        reader = pairs_creation(self.fname, self.vcf_fname)

        for pair in tqdm(reader, desc=f'filtering {self.fname}'):
            log_row['total_pairs'] += 1
            pair.classify(masking=0)
            if pair.call == 'no_phase_change':
                log_row['no_phase_change'] += 1
                continue
            chrom = pair.rec_1.reference_name
            start, end = pair.detection[0][1], pair.detection[-1][1]

            # handle tabix issues if organelle/contig region iterator can't be made
            try:
                false_lookup = list(set([
                    line for line in self.tabix_reader.fetch(chrom, start, end)]))
            except ValueError as e:
                if any(region in e.args[0] for region in ['cpDNA', 'mtDNA', 'contig']):
                    false_lookup = []
                else:
                    raise e

            if not false_lookup: # no false phase changes at all
                log_row['pairs_kept'] += 1 
                self.writer.write(pair.rec_1)
                self.writer.write(pair.rec_2)
                continue # no need to apply a filter at all

            # if false reads in vicinity - filter
            if self.method == 'midpoint':
                passed = self._filter_midpoint(pair, chrom, start, end, false_lookup)
            elif self.method == 'overlap':
                passed = self._filter_overlap(pair, chrom, start, end, false_lookup)
            elif self.method == 'nuclear':
                log_row['pairs_removed'] += 1
                continue # nuclear option means indiscriminately ignoring all overlapping reads
            else:
                raise ValueError(
                    f'method must be one of midpoint, overlap, or nuclear, not {self.method}')

            # if filter was applied and read did not match false read - write
            if passed:
                log_row['pairs_kept'] += 1
                self.writer.write(pair.rec_1)
                self.writer.write(pair.rec_2)
                continue
            elif not passed:
                log_row['pairs_removed'] += 1

        log_row['perc_kept'] = round(log_row['pairs_kept'] / log_row['total_pairs'], 3)
        log_writer.writerow(log_row)
                

    def _filter_midpoint(self, pair, chrom, start, end, false_lookup):
        # reduce pair.detection output for comparison
        detection_reduced = [(hap, pos) for hap, pos, base in pair.detection]
        # iterate over false reads
        for false_read in false_lookup:
            chrom, pos1, pos2, hap1, hap2, parent = false_read.split('\t')
            false_detection = [(hap1, int(pos1)), (hap2, int(pos2))]
            if all(var in detection_reduced for var in false_detection):
                return False # false phase change match found - exit loop
            else:
                continue
        return True # only returned if all false phase changes aren't matched


    def _filter_overlap(self, pair, chrom, start, end, false_lookup):
        for false_read in false_lookup:
            chrom, pos1, pos2, hap1, hap2, fp_detection = false_read.split('\t')
            false_detection = ast.literal_eval(fp_detection)
            overlap_sites = set(
                [v[1] for v in pair.detection]).intersection(
                    [v[1] for v in false_detection])
            detection_filtered = [
                var for var in pair.detection if var[1] in overlap_sites]
            if all(var in false_detection for var in pair.detection):
                return False # all overlapping sites are the same in the false read
            else:
                continue
        return True # only returned if all false reads in region aren't matched

            
def handle_bed(args):
    """Function that generates/handles BED file

    This function should be used to handle BED files instead of invoking
    BedGenerator directly.

    Parameters
    ----------
    args : argparse.ArgumentParser
        command line arguments

    Returns
    -------
    tabix_reader : pysam.libctabix.TabixFile
        TabixFile object for querying
    """
    if args.false_bed_in: # bed file already provided
        tabix_reader = pysam.TabixFile(args.false_bed_in)
        return tabix_reader

    # if not - generate lookup bed file
    if args.false_bed_out.endswith('.gz'):
        args.false_bed_out = args.false_bed_out.rstrip('.gz')
    bed_generator = BedGenerator(
        false_plus=args.false_plus, 
        false_minus=args.false_minus,
        vcf_fname=args.vcf,
        method=args.method,
        bed_fname=args.false_bed_out)

    bed_generator.generate_bed() # will use method specified above
    bed_generator.sort_tabix_bed()
    tabix_reader = pysam.TabixFile(args.false_bed_out + '.gz')
    return tabix_reader


def main():
    parser = arg_parser()
    args = parser.parse_args()
    tabix_reader = handle_bed(args)
    filterer = FalsePositiveFilterer(
        fname=args.fname,
        method=args.method,
        vcf_fname=args.vcf,
        tabix_reader=tabix_reader,
        out=args.out,
        log=args.log)
    filterer.filter_false_positives()
    print(f'[readcomb] written to {args.out}')


if __name__ == '__main__':
    main()

        

