#!/usr/bin/env python3
"""
Tools for recombination event classification.

Unlike the other modules in readcomb, which should be run at the command line,
classification is intended to have its contents imported and used in a Python
environment.
"""

import itertools
import pysam
from cyvcf2 import VCF

try:
    from readcomb.filter import check_variants
    from readcomb.filter import cigar
except ImportError as e:
    print('WARNING: readcomb is not installed')
    print('Use command: pip install readcomb')
    from filter import check_variants
    from filter import cigar

def downstream_phase_detection(variants, segment, record, quality):
    """
    Detect haplotype of variants in the given DNA sequence.

    This function is exactly the same as ``phase_detection`` utilized in ``readcomb-filter``
    except it also returns the index of variants.

    ``downstream_phase_detection()`` takes a ``segment`` and a list of ``cyvcf2``
    variants in the region of the segment to detect the haplotype of the sequence and
    generate a list of tuples in the form ``[({'1', '2', 'N'}, variant_index), ...]`` to represent
    haplotypes in the order of the variants on the segment.

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
    list of tuples
        list of haplotypes of variants and their index in the bounds of the segment
    """

    detection_results = []

    cigar_tuples = record.cigartuples

    # check for no alignment
    if not cigar_tuples:
        return detection_results

    for variant in variants:
        # Using SNP.start and record.reference_start since they are both 0 based
        # SNP.start grabs vcf positions in 0 index while vcfs are 1 indexed
        idx = variant.start - record.reference_start

        parent1 = variant.gt_bases[0].split('/')[0]
        parent2 = variant.gt_bases[1].split('/')[0]

        if parent1 == parent2: # ignore uninformative SNPs
            continue

        # ignore if snp is before sequence
        if idx < 0:
            continue

        # ignore variant if quality of sequencing at that base is below threshold
        query_qualities = record.query_qualities.tolist()
        # realign query_qualities if there are insertions
        count = 0
        insertions = 0
        for tupl in cigar_tuples:
            if count > idx:
                break
            if tupl[0] == 2:
                insertions += tupl[1]
            count += tupl[1]

        if query_qualities[idx - insertions] < quality:
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
                    detection_results.append(('1', variant.start))
                elif parent2_match:
                    detection_results.append(('2', variant.start))
                else:
                    detection_results.append(('N', variant.start))
            else:
                if parent2_match:
                    detection_results.append(('2', variant.start))
                elif parent1_match:
                    detection_results.append(('1', variant.start))
                else:
                    detection_results.append(('N', variant.start))

        else: # variant is a SNP
            if idx >= len(segment):
                break

            if segment[idx] == parent1:
                detection_results.append(('1', variant.start))

            elif segment[idx] == parent2:
                detection_results.append(('2', variant.start))

            else:
                detection_results.append(('N', variant.start))

    return detection_results

class Pair():
    def __init__(self, record1, record2, vcf_filepath):
        """
        Initialize a ``Pair`` object

        It is highly recommended to use ``classification.pairs_creation`` to create pairs by
        parsing through a SAM/BAM file.

        ``record1`` and ``record2`` are reordered so that they are in increasing sequential order

        Parameters
        ----------
        record1: pysam.AlignedSegment
            one of the bam sequences in a read pair
        record2: pysam.AlignedSegment
            bam sequence in a read pair corresponding to the one in ``record1``
        vcf_filepath: str
            filepath to VCF file that contains variants for ``record1`` and ``record2``
        """
        # put pairs in order by reference_start of the bam sequences
        if record1.reference_start < record2.reference_start:
            self.rec_1 = record1
            self.rec_2 = record2
        else:
            self.rec_1 = record2
            self.rec_2 = record1

        self.vcf_filepath = vcf_filepath

        # descriptive attributes for pairs
        self.midpoint = None
        self.segment_1 = cigar(self.rec_1)
        self.segment_2 = cigar(self.rec_2)
        self.variants_1 = None
        self.variants_2 = None

        # recombination event calling
        self.detection_1 = None
        self.detection_2 = None
        self.no_match = None
        self.condensed = None
        self.masked_condensed = None
        self.call = None
        self.masked_call = None
        self.gene_conversion_len = None

        # quality
        self.variants_per_haplotype = None
        self.min_variants_in_haplotype = None

    def __str__(self):
        """
        Generate string representation of ``Pair`` class.

        This function overwites the built-in ``__str__`` function of ``Pair``.
        The ``__str__`` function is called whenever a user prints a ``Pair`` object or
        converts the object to a ``str``. ``Pair`` has two different states, before and
        after calling ``classify()``, so two different string representations are
        needed.

        Returns
        -------
        string : str
            string representation of the ``Pair`` object
        """
        if hasattr(self, 'call'):
            string = f'Record name: {self.rec_1.query_name} \n' + \
                    f'Read1: {self.rec_1.reference_name}:{self.rec_1.reference_start}' + \
                    f'-{self.rec_1.reference_start + self.rec_1.query_alignment_length} \n' + \
                    f'Read2: {self.rec_2.reference_name}:{self.rec_2.reference_start}' + \
                    f'-{self.rec_2.reference_start + self.rec_2.query_alignment_length} \n' + \
                    f'VCF: {self.vcf_filepath} \n' + \
                    f'Unmatched Variant(s): {self.no_match} \n' + \
                    f'Condensed: {self.condensed} \n' + \
                    f'Call: {self.call} \n' + \
                    f'Condensed Masked: {self.masked_condensed} \n' + \
                    f'Masked Call: {self.masked_call} \n' + \
                    f'Midpoint: {self.get_midpoint()} \n' + \
                    f'Variants Per Haplotype: {self.variants_per_haplotype} \n' + \
                    f'Gene Conversion Length: {self.gene_conversion_len}'
        else:
            string = f'Record name: {self.rec_1.query_name} \n' + \
                    f'Read1: {self.rec_1.reference_name}:{self.rec_1.reference_start}' + \
                    f'-{self.rec_1.reference_start + self.rec_1.query_alignment_length} \n' + \
                    f'Read2: {self.rec_2.reference_name}:{self.rec_2.reference_start}' + \
                    f'-{self.rec_2.reference_start + self.rec_2.query_alignment_length} \n' + \
                    f'VCF: {self.vcf_filepath}'

        return string

    def package(self):
        """
        Convert ``self.rec_1`` and ``self.rec_2`` from ``pysam.AlignedSegment`` to ``str``

        The user may want to implement multiprocessing to decrease the amount
        of time to classify all reads in a SAM/BAM file. ``self.rec_1`` and
        self.``rec_2`` and all ``pysam.AlignedSegment`` objects are not pickleable and
        cannot be passed through a ``multiprocessing.Queue``. Instead of directly
        handling BAM/SAM strings, users can choose to create a ``Pair``, call
        ``package`` to convert the records to strings using the ``to_string()``
        function from ``pysam`` and pass the ``Pair`` object through a ``Queue``.
        """
        if (
            type(self.rec_1) == type(pysam.AlignedSegment()) and
            type(self.rec_2) == type(pysam.AlignedSegment())
        ):
            self.rec_1 = self.rec_1.to_string()
            self.rec_2 = self.rec_2.to_string()

    def unpackage(self, header):
        """
        Convert ``self.rec_1`` and ``self.rec_2`` from ``str`` to ``pysam.AlignedSegment``

        If the user is implementing multiprocessing to decrease the amount of
        time to classify all reads in a SAM/BAM file, the ``package()`` function should
        be called on ``Pair`` to convert ``self.rec_1`` and ``self.rec_2`` to ``str``
        and passed to a subprocess. After the subprocess receives the packaged
        ``Pair``, it should be unpackaged using this function.

        ``pysam.AlignedSegment.fromstring()`` also requires the ``pysam.AlignmentHeader``
        object which can be obtained from calling ``.header`` on the
        ``pysam.AlignmentFile`` object that this ``Pair`` was parsed from using
        ``pair_creation()``.

        Parameters
        ----------
        header : pysam.AlignmentHeader
            ``pysam.AlignmentHeader`` of the ``pysam.AlignmentFile`` associated with this ``Pair``
        """
        if isinstance(self.rec_1, str) and isinstance(self.rec_1, str):
            self.rec_1 = pysam.AlignedSegment.fromstring(self.rec_1, header)
            self.rec_2 = pysam.AlignedSegment.fromstring(self.rec_2, header)

    def classify(self, masking=70, quality=30, vcf=None):
        """
        Determine the type of recombination event that occursed in this read pair.

        ``classify()`` does phase change analysis similar to ``readcomb-filter`` with addtional
        steps to classify the type of recombination event. See ``readcomb.check_variants()``,
        ``readcomb.cigar()``, and ``downstream_phase_detection`` for more information on the
        analysis process.

        Classification is done on the haplotypes present in the full read pair and also when the
        beginning and end are shortened by ``masking`` bases as phase changes close to the start
        and end of a read pair are difficult to call as crossovers or gene conversions. The two
        different classifications allows a more nuanced call on the phase change type present in
        the pair.

        ``classify()`` also takes in optional ``vcf`` for a ``cyvcf2.VCF`` object. It is highly
        recommended to pass in a pre-initialized ``VCF`` object when large amounts of ``Pair``
        objects are being classified as the creation of ``VCF`` parsers greatly affects
        performance and tends to be less reliable when large amounts of parsers are created.

        Do not use the same ``cyvcf2.VCF`` parser across multiple processes/threads as it leads to
        errors involving file access permissions. Use a dedicated ``cyvcf2.VCF`` parser for
        each process/thread instead.

        Parameters
        ----------
        masking : int, default 70
            number of bases to be ignored when determing the ``masked_classification``
        quality : int, default 30
            filter quality for individual bases in a sequence
        vcf: cyvcf2.VCF, optional
            pre-initialized ``cyvcf2.VCF`` for this ``Pair`` to use
            when classifying, highly recommended
        """
        # create new VCF object if none is provided
        if not vcf:
            vcf = VCF(self.vcf_filepath)

        self.variants_1 = check_variants(
            vcf, self.rec_1.reference_name, self.rec_1.reference_start, 
            self.rec_1.reference_start + self.rec_1.query_alignment_length)
        self.variants_2 = check_variants(
            vcf, self.rec_2.reference_name, self.rec_2.reference_start, 
            self.rec_2.reference_start + self.rec_2.query_alignment_length)

        self.detection_1 = downstream_phase_detection(
            self.variants_1, self.segment_1, self.rec_1, quality)
        self.detection_2 = downstream_phase_detection(
            self.variants_2, self.segment_2, self.rec_2, quality)

        # set no_match variable if there are unmatched variants
        self.no_match = bool('N' in self.detection_1 or 'N' in self.detection_2)

        # simplification of results
        # [(haplotype, beginning, end), ...]
        self.condensed = []

        for variant in self.detection_1 + self.detection_2:
            haplotype = variant[0]
            location = variant[1]

            # first variant
            if len(self.condensed) == 0:
                self.condensed.append(
                    [haplotype, self.rec_1.reference_start, self.rec_1.reference_start])

            # different haplotype
            elif self.condensed[-1][0] != haplotype:
                # middle of previous variant location and current variant
                # gonna round down so it's not a decimal
                midpoint = int((self.condensed[-1][2] + location) // 2)
                self.condensed[-1][2] = midpoint
                self.condensed.append([haplotype, midpoint, midpoint])

            # last variant
            if len(self.detection_2) == 0:
                if variant == self.detection_1[-1]:
                    self.condensed[-1][2] = self.rec_1.reference_start + \
                        self.rec_1.query_alignment_length
            else:
                if variant == self.detection_2[-1]:
                    self.condensed[-1][2] = self.rec_2.reference_start + \
                        self.rec_2.query_alignment_length

        # create list of just haplotype information no range from condensed
        haplotypes = [tupl[0] for tupl in self.condensed if tupl[0] != 'N']

        # classify condensed
        if len(haplotypes) == 2:
            self.call = 'cross_over'
        elif len(haplotypes) == 3:
            self.call = 'gene_conversion'
        elif len(haplotypes) > 3:
            self.call = 'complex'
        else:
            self.call = 'no_phase_change'

        if self.rec_2.reference_start + self.rec_2.query_alignment_length \
            - self.rec_1.reference_start < masking * 2:

            print('Masking size too large for pair: ' + self.rec_1.query_name)
            self.masked_condensed = None
            self.masked_call = None
            return

        mask_start = self.rec_1.reference_start + masking
        mask_end = self.rec_2.reference_start + self.rec_2.query_alignment_length - masking

        # create masked_condensed from condensed
        # [(haplotype, beginning, end), ...]
        self.masked_condensed = []

        for phase in self.condensed:
            # one phase contains both mask start and end
            if phase[1] < mask_start and mask_end < phase[2]:
                self.masked_condensed.append([phase[0], mask_start, mask_end])
            # phase contains only mask start
            elif phase[1] < mask_start < phase[2]:
                self.masked_condensed.append([phase[0], mask_start, phase[2]])
            # phase contains only mask end
            elif phase[1] < mask_end < phase[2]:
                self.masked_condensed.append([phase[0], phase[1], mask_end])
            # phase is in the middle
            elif mask_start < phase[1] and phase[2] < mask_end:
                self.masked_condensed.append(phase)

        # create list of just haplotype information no range from condensed
        haplotypes = [tupl[0] for tupl in self.masked_condensed if tupl[0] != 'N']

        # classify masked_condensed
        if len(haplotypes) == 2:
            self.masked_call = 'cross_over'
        elif len(haplotypes) == 3:
            self.masked_call = 'gene_conversion'
        elif len(haplotypes) > 3:
            self.masked_call = 'complex'
        else:
            self.masked_call = 'no_phase_change'

        # length of gene_conversion
        if self.call == 'gene_conversion':
            self.gene_conversion_len = self.condensed[-1][1] - self.condensed[0][2]
        else:
            self.gene_conversion_len = 'N/A'

        # calculate average number of variants per haplotype
        self.variants_per_haplotype = len(self.variants_1 + self.variants_2) / \
            max(len(haplotypes), 1)

        # get the lowest number of variants a haplotype has -
        # splits variant list (e.g. ['1', '1', '2', '1']) and gets min variant count across haps
        if 'no_phase_change' not in [self.call, self.masked_call]:
            self.min_variants_in_haplotype = min(
                len(list(grouper)) for value, grouper
                in itertools.groupby([hap for hap, position in self.detection_1 + self.detection_2])
                )

        # convert haplotypes 1/2/N to VCF sample names
        samples = vcf.samples

        for tupl in self.condensed:
            if tupl[0] == '1':
                tupl[0] = samples[0]
            elif tupl[0] == '2':
                tupl[0] = samples[1]

        for tupl in self.masked_condensed:
            if tupl[0] == '1':
                tupl[0] = samples[0]
            elif tupl[0] == '2':
                tupl[0] = samples[1]

    def get_midpoint(self):
        """
        Determine and return the midpoint of ``Pair``

        The midpoint of a read pair with no phase changes is halfway between the start
        of the first read and end of second read.

        The midpoint of a read pair with a phase change is halfway between the two closest
        variants that signify a phase change event. For gene conversions, it is halfway between
        the two outer variants of the group of 3. This logic extends to read pairs
        with complex haplotypes.

        Returns
        -------
        midpoint : int
            the middle point of the phase change event of the ``Pair``
        """
        # give error if Pair object is packaged
        if isinstance(self.rec_1, str) or isinstance(self.rec_2, str):
            raise TypeError('cannot classify if Pair object is not unpackaged yet')

        # return midpoint if it's already been called
        if self.midpoint:
            return self.midpoint

        # classify if read has not been already
        if not hasattr(self, 'call'):
            self.classify()

        # simplification of results
        # [(haplotype, beginning, end), ...]

        if self.call == 'no_phase_change':
            # midpoint is middle of two paired reads if no phase change
            start = self.rec_1.reference_start
            end = self.rec_2.reference_start + self.rec_2.query_alignment_length
            self.midpoint = int((start + end) / 2)
        elif self.call == 'phase_change':
            # (X, begin, end), (Y, begin, end): end of X = beginning of Y = midpoint
            self.midpoint = self.condensed[0][2]
        else:
            # kind of a shortcut using midpoints to find the middle
            start = self.condensed[0][2]
            end = self.condensed[-1][1]
            self.midpoint = (start + end) / 2

        return self.midpoint


def pairs_creation(bam_filepath, vcf_filepath):
    """
    Parses through a BAM/SAM file and generates ``Pair`` objects for read pairs.

    Given a read name sorted BAM/SAM file, ``pairs_creation`` creates ``Pair``
    objects with read pairs.

    It is recommended to use this function to create ``Pair`` objects for an entire
    BAM/SAM file instead of manually creating the ``Pair``. It is also recommended
    to use ``pairs_creation`` on SAM files that have already been filtered by
    ``readcomb-filter`` instead of full BAMs/SAMs.

    Parameters
    ----------
    bam_filepath : str
        filepath of BAM/SAM file processed with ``readcomb-bamprep``
    vcf_filepath : str
        filepath of VCF files that contain variants for reads in ``bam_filepath``

    Returns
    -------
    pairs : generator
        generator that yields reads from ``bam_filepath`` in the form of ``Pair`` objects
    """
    bam = pysam.AlignmentFile(bam_filepath, 'r')

    prev_rec = None
    for rec in bam:
        # first record
        if not prev_rec:
            prev_rec = rec
        # check if query_name pairs exist
        elif rec.query_name == prev_rec.query_name:
            yield Pair(rec, prev_rec, vcf_filepath)
        elif rec.query_name != prev_rec.query_name:
            prev_rec = rec
            continue
