#!/usr/bin/env python3
"""
Tools for recombination event classification.

Unlike the other modules in readcomb, which should be run at the command line,
classification is intended to have its contents imported and used in a Python
environment.
"""

import itertools
import re
import numpy as np
import pysam
from cyvcf2 import VCF

try:
    from readcomb.filter import check_variants
    from readcomb.filter import cigar
    from readcomb.filter import qualities_cigar
except ImportError as e:
    print('WARNING: readcomb is not installed')
    print('Use command: pip install readcomb')
    from filter import check_variants
    from filter import cigar
    from filter import qualities_cigar

__version__ = '0.4.3'

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

    # get realigned base qualities
    query_qualities = qualities_cigar(record)

    for variant in variants:
        # Using SNP.start and record.reference_start since they are both 0 based
        # SNP.start grabs vcf positions in 0 index while vcfs are 1 indexed
        idx = variant.start - record.reference_start

        parent1 = variant.gt_bases[0].split('/')[0].split('|')[0]
        parent2 = variant.gt_bases[1].split('/')[0].split('|')[0]

        if parent1 == parent2: # ignore uninformative SNPs
            continue

        # ignore if snp is before sequence
        if idx < 0:
            continue

        # ignore variant if quality of sequencing at that base is below threshold
        if not query_qualities[idx]:
            continue # variant site is deleted in read

        if query_qualities[idx] < quality:
            continue

        if variant.is_indel:
            # check if indel is outside of segment
            if idx + max(len(parent1), len(parent2)) > len(segment):
                continue

            parent1_match = segment[idx:idx + len(parent1)] == parent1
            parent2_match = segment[idx:idx + len(parent2)] == parent2

            if parent1_match and not parent2_match:
                detection_results.append(('1', variant.start, parent1))
            elif parent2_match and not parent1_match:
                detection_results.append(('2', variant.start, parent2))
            elif parent1_match and parent2_match:
                continue
            else:
                detection_results.append(('N', variant.start, None))

        else: # variant is a SNP
            if idx >= len(segment):
                break

            if segment[idx] == parent1:
                detection_results.append(('1', variant.start, parent1, query_qualities[idx]))

            elif segment[idx] == parent2:
                detection_results.append(('2', variant.start, parent2, query_qualities[idx]))

            else:
                detection_results.append(('N', variant.start, None, query_qualities[idx]))

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
        elif (
            record1.reference_start == record2.reference_start and
            record1.reference_end > record2.reference_end # rare, but possible
        ):
            self.rec_1 = record2
            self.rec_2 = record1
        else:
            self.rec_1 = record2
            self.rec_2 = record1

        self.vcf_filepath = vcf_filepath

        # descriptive attributes for pairs
        self.segment_1 = cigar(self.rec_1)
        self.segment_2 = cigar(self.rec_2)
        self.variants_1 = None
        self.variants_2 = None
        self.variants_all = None
        self.variants_filt = None
        self.midpoint = None
        self.relative_midpoint = -1
        self.location = ''.join(f'{self.rec_1.reference_name}:\
        {self.rec_1.reference_start}-\
        {self.rec_2.reference_start + len(self.segment_2)}'.split())

        # recombination event calling
        self.detection_1 = None
        self.detection_2 = None
        self.detection = None
        self.masked_detection = None
        self.no_match = None
        self.condensed = None
        self.masked_condensed = None
        self.call = None
        self.masked_call = None
        self.phase_change_variants = None
        self.gene_conversion_len = None

        # quality metrics
        self.overlap = None
        self.overlap_disagree = None
        self.variants_per_haplotype = -1
        self.min_variants_in_haplotype = -1
        self.outer_bound = -1
        self.min_end_proximity = -1
        self.variant_skew = -1
        self.mismatch_variant_ratio = -1
        self.indels = None
        self.indel_proximity = -1
        self.proximate_indel_length = -1
        self.min_base_qual = -1
        self.mean_base_qual = -1
        self.variant_counts = None
        self.converted_haplotype = None
        self.converted_variants = None
        self.conversion_tract = None

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
        if self.call:
            string = f'Record name: {self.rec_1.query_name} \n' + \
                    f'Read1: {self.rec_1.reference_name}:{self.rec_1.reference_start}' + \
                    f'-{self.rec_1.reference_start + len(self.segment_1)} \n' + \
                    f'Read2: {self.rec_2.reference_name}:{self.rec_2.reference_start}' + \
                    f'-{self.rec_2.reference_start + len(self.segment_2)} \n' + \
                    f'VCF: {self.vcf_filepath} \n' + \
                    f'Unmatched Variant(s): {self.no_match} \n' + \
                    f'Condensed: {self.condensed} \n' + \
                    f'Call: {self.call} \n' + \
                    f'Condensed Masked: {self.masked_condensed} \n' + \
                    f'Masked Call: {self.masked_call} \n' + \
                    f'Midpoint: {self.get_midpoint()} \n' + \
                    f'Variants Per Haplotype: {self.variants_per_haplotype} \n' + \
                    f'Variant Skew: {self.variant_skew} \n' + \
                    f'Gene Conversion Length: {self.gene_conversion_len} \n' + \
                    f'Min Variants In Haplotype: {self.min_variants_in_haplotype} \n' + \
                    f'Mismatch/Variant Ratio: {self.mismatch_variant_ratio} \n'
        else:
            string = f'Record name: {self.rec_1.query_name} \n' + \
                    f'Read1: {self.rec_1.reference_name}:{self.rec_1.reference_start}' + \
                    f'-{self.rec_1.reference_start + len(self.segment_1)} \n' + \
                    f'Read2: {self.rec_2.reference_name}:{self.rec_2.reference_start}' + \
                    f'-{self.rec_2.reference_start + len(self.segment_2)} \n' + \
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

    def classify(self, masking=0, quality=0, vcf=None):
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
        masking : int, default 0
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
            self.rec_1.reference_start + len(self.segment_1))
        self.variants_2 = check_variants(
            vcf, self.rec_2.reference_name, self.rec_2.reference_start, 
            self.rec_2.reference_start + len(self.segment_2))
        vcf.close()

        self.detection_1 = downstream_phase_detection(
            self.variants_1, self.segment_1, self.rec_1, quality)
        self.detection_2 = downstream_phase_detection(
            self.variants_2, self.segment_2, self.rec_2, quality)

        # deal with 'het pairs' and 'duplicate variants' - see _resolve_overlap() documentation
        # updates the read-specific detection methods
        self._resolve_overlap()
        self.detection = sorted(self.detection_1 + self.detection_2,
            key=lambda variant: variant[1]) # sort by position

        # create condensed variant list for easy look up
        # this hideous listcomp removes duplicates and sorts variants
        self.variants_all = [
            [variant for variant in self.variants_1 + self.variants_2 if variant.POS == pos][0]
            for pos in sorted(set(variant.POS for variant in self.variants_1 + self.variants_2))]

        # create filtered variant list only containing variants used for detection
        self.variants_filt = [
            variant for variant in self.variants_all
            if variant.POS in [v[1] + 1 for v in self.detection]]
            
        # set no_match variable if there are unmatched variants
        self.no_match = any([v[0] == 'N' for v in self.detection])

        # simplification of results - creating condensed hap view
        # [(haplotype, beginning, end), ...]
        self.condensed = self._create_condensed(self.detection)

        # create list of just haplotype information without range from condensed
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

        if self.rec_2.reference_start + len(self.segment_2) \
            - self.rec_1.reference_start < masking * 2:
            # masking size is too large for pair
            self.masked_condensed = None
            self.masked_call = None
            return

        mask_start = self.rec_1.reference_start + masking
        mask_end = self.rec_2.reference_start + len(self.segment_2) - masking
        self.masked_detection = [
            v for v in self.detection if mask_start < v[1] < mask_end]

        # create masked_condensed from masked detection
        # [(haplotype, beginning, end), ...]
        self.masked_condensed = self._create_condensed(self.masked_detection)

        # create list of just haplotype information no range from condensed
        masked_haplotypes = [tupl[0] for tupl in self.masked_condensed if tupl[0] != 'N']

        # classify masked_condensed
        if len(masked_haplotypes) == 2:
            self.masked_call = 'cross_over'
        elif len(masked_haplotypes) == 3:
            self.masked_call = 'gene_conversion'
        elif len(masked_haplotypes) > 3:
            self.masked_call = 'complex'
        else:
            self.masked_call = 'no_phase_change'

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

        # set descriptive + quality attributes post-classification
        self._describe(haplotypes, vcf)

        # if gene conversion - set descriptive info
        if 'gene_conversion' in [self.call, self.masked_call]:
            self._describe_gene_conversion(vcf)
        else:
            self.gene_conversion_len = 'N/A'
            self.converted_haplotype = 'N/A'
            self.converted_variants = 'N/A'
            self.conversion_tract = 'N/A'

        self.midpoint, self.relative_midpoint = self.get_midpoint()

    def _create_condensed(self, detection):
        """
        Internal method to generate haplotype representation in
        ``Pair.condensed``.

        Parameters
        ----------
        detection : list
            list containing variants in downstream_phase_detection format
        """
        condensed = []

        for variant in detection:
            haplotype, location, base, qual = variant
            if haplotype == 'N':
                continue

            # first variant
            if len(condensed) == 0:
                condensed.append(
                    [haplotype, self.rec_1.reference_start, self.rec_1.reference_start])

            # if same haplotype - extend segment for correct eventual midpoint calc
            elif condensed[-1][0] == haplotype:
                condensed[-1][2] = location

            # different haplotype
            elif condensed[-1][0] != haplotype:
                # middle of previous variant location and current variant
                # rounding down so it's not a decimal
                midpoint = int((condensed[-1][2] + location) // 2)
                condensed[-1][2] = midpoint
                condensed.append([haplotype, midpoint, midpoint])

            # last variant
            if len(self.detection_2) == 0:
                if variant == self.detection_1[-1]:
                    condensed[-1][2] = self.rec_1.reference_start + \
                        len(self.segment_1)
            else:
                if variant == self.detection_2[-1]:
                    condensed[-1][2] = self.rec_2.reference_start + \
                        len(self.segment_2)

        # this occurs during false positive phase changes when reads overlap
        if any(start == end for hap, start, end in condensed):
            for false_hap in [h for h in condensed if h[1] == h[2]]:
                condensed.remove(false_hap)
        
        return condensed

    def _describe(self, haplotypes, vcf):
        """
        Internal method to calculate various attributes of ``Pair``. Sets
        the following attributes:

        - mismatch_variant_ratio - number of mismatches to reference over the
          number of variants - see ``Pair._get_mismatch_variant_ratio()``
          documentation for details
        - variants_per_haplotype - average number of variants supporting each
          haplotype
        - min_variants_in_haplotype - number of variants in least
          well-supported haplotype
        - outer_bound - relative location of 'most outer' variant involved
          in a phase change
        - min_end_proximity - lowest distance in bp between a variant involved in a
          phase change and the end of the read it's on
        - phase_change_variants - list of phase change variant positions
        - indels - list of indels
        - indel_proximity - lowest distance in bp between a variant involved in
          a phase change and an indel
        - proximate_indel_length - if indel_proximity != -1, size in bp of
          nearest indel
        - min_base_qual - minimum base qual for a phase change variant
        - mean_base_qual - mean base qual across phase change variants

        Parameters
        ----------
        haplotypes : list
            haplotype listing generated from self.masked_condensed
        vcf : cyvcf2.VCF
            VCF file object generated by ``classify()``
        """
        if not self.call:
            self.classify()

        # use helper function to calculate unexpected mismatch count
        self.mismatch_variant_ratio = self._get_mismatch_variant_ratio()

        # calculate average number of variants per haplotype and create dict of counts
        self.variants_per_haplotype = len(self.variants_filt) / max(len(haplotypes), 1)
        self.variant_counts = {
            hap: [t[0] for t in self.detection].count(hap) for hap in set(haplotypes)}
        if len(self.variant_counts.values()) == 2: # only calc if both haps present
            self.variant_skew = max(self.variant_counts.values()) / min(self.variant_counts.values())

        # subsequent metrics can't be calculated for non-phase-change pairs
        if self.call == 'no_phase_change':
            return

        # get the lowest number of variants a haplotype has -
        # splits variant list (e.g. ['1', '1', '2', '1']) and gets min variant count across haps
        self.min_variants_in_haplotype = min(
            len(list(grouper)) for value, grouper
            in itertools.groupby([hap for hap, position, allele, qual in self.detection])
            )

        # get variants involved in phase change(s)
        # used for outer_bound, min_end_proximity, indel_proximity, proximate_indel_length
        phase_change_variants = [
            [self.detection[i][1], self.detection[i+1][1]]
            for i in range(len(self.detection) - 1)
            if self.detection[i][0] != self.detection[i+1][0]]
        self.phase_change_variants = sorted(list(set(
            [pos for var_pair in phase_change_variants for pos in var_pair]
        )))

        # only check leftmost and rightmost variants to assign outer bound
        if (
            self.phase_change_variants[0] - self.condensed[0][1] < \
            abs(self.phase_change_variants[-1] - self.condensed[-1][2])
        ):
            self.outer_bound = round(
                (self.phase_change_variants[0] - self.rec_1.reference_start) / \
                (self.condensed[-1][2] - self.condensed[0][1]), 3)
        else:
            self.outer_bound = round(
                (self.phase_change_variants[-1] - self.rec_1.reference_start) / \
                (self.condensed[-1][2] - self.condensed[0][1]), 3)

        # get min_end_proximity
        read_bounds = [
            self.rec_1.reference_start, self.rec_1.reference_start + len(self.segment_1),
            self.rec_2.reference_start, self.rec_2.reference_start + len(self.segment_2)]
        self.min_end_proximity = min(
            [abs(pos - bound) for bound in read_bounds for pos in self.phase_change_variants])

        # get indels for calculation of indel attributes
        indels = []
        for rec in [self.rec_1, self.rec_2]:
            idx = rec.reference_start
            for op, basecount in rec.cigartuples:
                if op == 0:
                    idx += basecount
                elif op == 1:
                    # type, ref_start_idx, ref_end_idx, indel_length
                    indels.append(['ins', idx, idx, basecount])
                elif op == 2:
                    indels.append(['del', idx, idx + basecount, basecount])
        # overlapping reads may have the same indel register twice
        self.indels = list(indel for indel, _ in itertools.groupby(indels)) 

        if indels:
            self.indel_proximity = min(abs(np.concatenate(
                [np.array(self.phase_change_variants) - pos
                 for indel in self.indels
                 for pos in indel[1:3]])))

            # length of most proximate indel
            for indel in self.indels:
                diffs = np.concatenate(
                    [abs(np.array(self.phase_change_variants) - pos)
                    for pos in indel[1:3]])
                if self.indel_proximity in diffs:
                    self.proximate_indel_length = indel[-1]
                    break
        else:
            self.indel_proximity = -1
            self.proximate_indel_length = -1

        # base qual metrics
        self.min_base_qual = min(site[-1] for site in self.detection)
        self.mean_base_qual = round(
            sum(site[-1] for site in self.detection) / len(self.detection), 2)


    def _describe_gene_conversion(self, vcf):
        """
        Internal method to calculate various attributes of ``Pair`` if the
        pair is called as a gene conversion. Sets the following attributes:

        The following attributes are only set if the read pair is called as a
        gene conversion:
        - gene_conversion_len - length of gene conversion segment, if applicable,
          as determined by midpoint method
        - converted_haplotype - which haplotype was 'converted to'
        - converted_variants - which variants were converted in the tract
        - conversion_tract - 2-length list with bounds of conversion tract
        """
        self.gene_conversion_len = self.condensed[-1][1] - self.condensed[0][2]

        # converted haplotype and variants
        self.converted_haplotype = [
            hap for hap, count in self.variant_counts.items() 
            if count == min(self.variant_counts.values())][0]
        self.converted_variants = [
            variant for variant in self.detection if variant[0] == self.converted_haplotype]

        # set converted haplotype to sample name
        if self.converted_haplotype == '1':
            self.converted_haplotype = vcf.samples[0]
        elif self.converted_haplotype == '2':
            self.converted_haplotype = vcf.samples[1]

        # get conversion tract
        self.conversion_tract = np.array([
            tract for tract in self.condensed 
            if tract[0] == self.converted_haplotype][0][1:])

    def _get_mismatch_variant_ratio(self):
        """
        Helper function to calculate mismatch_variant_ratio quality metric. Will
        use MD tags in reads to get the number of mismatches to the reference
        per read, and then compare with the number of variants. Ratios above 1
        indicate variation that was either deemed too low quality to be
        included in the VCF. It is possible some sites could have diverged from
        the reference and resulted in the same allele in both parents, but this
        is less likely than the region simply being prone to paralogous
        alignment if the ratio is particularly high.

        The MD tag also includes deletions, designated as '^ATG' if 'ATG' is a
        3 bp deletion from the reference in the read. These deletions are
        removed from consideration.
        """
        if len(self.variants_1) + len(self.variants_2) == 0:
            return -1
        if len(self.variants_filt) == 0:
            return -1

        if not self.overlap:
            rec_1_mismatches = len(
                re.findall(
                    r'[ACGT]', re.sub(r'\^[ACGT]+', '', self.rec_1.get_tag('MD'))
                    )
                )
            rec_2_mismatches = len(
                re.findall(
                    r'[ACGT]', re.sub(r'\^[ACGT]+', '', self.rec_2.get_tag('MD'))
                    )
                )
            return (rec_1_mismatches + rec_2_mismatches) / len(self.variants_filt)

        elif self.overlap:
            # if rec2 is entirely 'contained' by rec1 - just do rec1
            if (
                self.rec_1.reference_start <= self.rec_2.reference_start and
                self.rec_2.reference_end < self.rec_1.reference_end
            ):
                m1 = self.rec_1.get_tag('MD')
                mismatches = len(
                    re.findall(
                        r'[ACGT]', re.sub(r'\^[ACGT]+', '', self.rec_1.get_tag('MD'))
                        )
                    )
                return mismatches / len(self.variants_filt)

            # otherwise - handle overlap as normal
            m1 = self.rec_1.get_tag('MD')
            m2 = self.rec_2.get_tag('MD')
            left_region = None
            right_region = None

            # get overlap
            for i, _ in enumerate(m1):
                # remove up to 3 digits from end
                # since these may differ
                subseq = re.sub(r'[0-9]{0,3}$', '', m1[i:])
                # test subseqs
                if subseq in m2 and len(subseq) != 1: # not last possible base
                    overlap_seq = subseq
                    break
            if not overlap_seq: # happens if the segment removed in re.sub above was the overlap
                if m2 in m1[-3:]:
                    overlap_seq = m2
                else:
                    # overlap still not detected - just keep both
                    left_region = m1
                    right_region = m2

            # get regions on either side of overlap if not assigned just above
            if not left_region and not right_region:
                left_region = m1[:m1.index(overlap_seq)]
                right_region = re.sub(
                    r'^[0-9]{0,3}', '',
                    m2[m2.index(overlap_seq) + len(overlap_seq):]
                )
            mismatches = 0
            for region in [left_region, overlap_seq, right_region]:
                mismatches += len(
                    re.findall(
                        r'[ACGT]', re.sub(r'\^[ACGT]+', '', region)
                        )
                    )
            return mismatches / len(self.variants_filt)

    def _resolve_overlap(self):
        """
        Helper function for ``Pair.classify()``. 

        Resolves cases where the same site is sequenced in both reads (i.e. the
        reads overlap) but the calls at the two sites differ. This will
        otherwise be incorrectly flagged as a phase change when it's more
        likely to be a sequencing error or paralogous alignment. This method
        filters the variants_* and detection_* attributes to remove any
        variants where this is the case.

        If the variants do have the same call, the variant in the second (3')
        read will be removed so that the call doesn't show multiple times in
        ``Pair.detection``.

        Will set the self.overlap attribute to True if variants overlap.
        Finally, will check whether there are discrepancies in base calls
        between overlapping reads at all (and will set overlap_disagree
        accordingly).
        """
        # check variant site overlap first
        overlapping_sites = list(set(var[1] for var in self.detection_1).intersection(
            set(var[1] for var in self.detection_2)))

        # set self.overlap to True if variant sites overlap
        if overlapping_sites and not self.overlap:
            self.overlap = True

        for site in overlapping_sites:
            read_1_call = [(hap, pos, base, qual) for hap, pos, base, qual in self.detection_1 if pos == site][0]
            read_2_call = [(hap, pos, base, qual) for hap, pos, base, qual in self.detection_2 if pos == site][0]

            # if parental hap calls at same site disagree, remove variant from consideration
            if read_1_call[0] != read_2_call[0]:
                self.detection_1.remove(read_1_call)
                self.detection_2.remove(read_2_call)
            elif read_1_call[0] == read_2_call[0]:
                self.detection_2.remove(read_2_call) # only remove from right read

        # check if there is any overlap at all
        if not self.overlap:
            if self.rec_2.reference_start < self.rec_1.reference_start + len(self.segment_1):
                self.overlap = True
            else:
                self.overlap = False

        # if there's overlap, check whether overlapping sites anywhere in read disagree
        if self.overlap:
            # identify overlapping region
            overlap_size = self.rec_1.reference_start + len(self.segment_1) - self.rec_2.reference_start
            idx = self.rec_2.reference_start - self.rec_1.reference_start
            region = [idx, idx + overlap_size] # idx for segment 1
            segment_1 = self.segment_1[region[0]:region[1]]
            segment_2 = self.segment_2[0:overlap_size]
            if segment_1 != segment_2: # compare full strings
                self.overlap_disagree = False
            else:
                self.overlap_disagree = True

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
        relative_midpoint : float
            the location of the midpoint relative to the length of the ``Pair``
        """
        # give error if Pair object is packaged
        if isinstance(self.rec_1, str) or isinstance(self.rec_2, str):
            raise TypeError('cannot classify if Pair object is not unpackaged yet')

        # return midpoint if it's already been called
        if self.midpoint:
            return self.midpoint, self.relative_midpoint

        # classify if read has not been already
        if not hasattr(self, 'call'):
            self.classify()

        # simplification of results
        # [(haplotype, beginning, end), ...]

        if self.call == 'no_phase_change':
            # midpoint is middle of two paired reads if no phase change
            start = self.rec_1.reference_start
            end = self.rec_2.reference_start + len(self.segment_2)
            self.midpoint = int((start + end) / 2)
            self.relative_midpoint = -1
        elif self.call == 'cross_over':
            # (X, begin, end), (Y, begin, end): end of X = beginning of Y = midpoint
            self.midpoint = self.condensed[0][2]
            # get relative midpoint
            start = self.condensed[0][1]
            end = self.condensed[-1][2]
            self.relative_midpoint = round(
                (self.midpoint - self.rec_1.reference_start) / (end - start), 3)
        elif self.call == 'gene_conversion':
            # (X, begin, end), (Y, begin, end), (X', begin, end)
            # midpoint of Y tract
            self.midpoint = sum([
                [start, end] for hap, start, end in self.condensed
                if hap == self.converted_haplotype][0]) / 2
            start = self.condensed[0][1]
            end = self.condensed[-1][2]
            self.relative_midpoint = round(
                (self.midpoint - self.rec_1.reference_start) / (end - start), 3)
        else:
            # kind of a shortcut using midpoints to find the middle
            start = self.condensed[0][2]
            end = self.condensed[-1][1]
            self.midpoint = (start + end) / 2
            self.relative_midpoint = 0.5

        # return in samtools format
        self.midpoint = f'{self.rec_1.reference_name}:{round(self.midpoint)}'

        return self.midpoint, self.relative_midpoint


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
            yield Pair(prev_rec, rec, vcf_filepath)
        elif rec.query_name != prev_rec.query_name:
            prev_rec = rec
            continue
