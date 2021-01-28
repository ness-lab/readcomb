#!/usr/bin/env python3

import pysam
from cyvcf2 import VCF
from readcomb.filter import check_variants
from readcomb.filter import cigar
    
def downstream_phase_detection(variants, segment, record):
    """
    Takes a segment and a list in the region of the segment and generates a list of
    strings to represent in order the parent of each variant on the segment

    snps: list of cyvcf2 variants in the region of the segment given by the function check_variants(),
        can include variants that are outside of the region when there is soft clipping
    segment: string sequence that is the segment built by the function cigar()
    record: bam record from pysam alignment file

    Returns a list of strings ('1', '2', or 'N')
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
    '''
    bam read pair object
    '''
    def __init__(self, record1, record2, vcf_filepath):

        # put pairs in order by reference_start of the bam sequences
        if record1.reference_start < record2.reference_start:
            self.rec_1 = record1
            self.rec_2 = record2
        else:
            self.rec_1 = record2
            self.rec_2 = record1
        
        self.vcf_filepath = vcf_filepath
    
    def package(self):
        '''
        convert bam sequences to strings so that they can be passed through queues/pipes
        '''
        self.rec_1 = self.rec_1.to_string()
        self.rec_2 = self.rec_2.to_string()

    def unpackage(self, header):
        '''
        convert bam strings to pysam bam objects, requires the header of the bam file that
        the string/sequence came from
        '''
        self.rec_1 = pysam.AlignedSegment.fromstring(self.rec_1, header)
        self.rec_2 = pysam.AlignedSegment.fromstring(self.rec_2, header)

    def classify(self, masking=70, vcf=None):
        '''
        Does phase change analysis on bam pair associated with this pair object and classifies
        them depending on the results, providing both a classification on the whole sequence
        and masked classification when bases of argument masking length is ignored at both ends

        Takes in optional vcf arguement to take in a cyvcf2.VCF object for use in multiprocessing/threading
        Note: Do not use the same VCF object across threads/processes
        '''
        self.segment_1 = cigar(self.rec_1)
        self.segment_2 = cigar(self.rec_2)

        # create new VCF object if none is provided
        if not vcf:
            vcf = VCF(self.vcf_filepath)

        self.variants_1 = check_variants(vcf, self.rec_1.reference_name, 
                                self.rec_1.reference_start,
                                self.rec_1.reference_start + self.rec_1.query_alignment_length)
        self.variants_2 = check_variants(vcf, self.rec_2.reference_name, 
                                self.rec_2.reference_start,
                                self.rec_2.reference_start + self.rec_2.query_alignment_length)

        self.detection_1 = downstream_phase_detection(self.variants_1, self.segment_1, self.rec_1)
        self.detection_2 = downstream_phase_detection(self.variants_2, self.segment_2, self.rec_2)

        # simplification of results
        # [(haplotype, beginning, end), ...]
        self.condensed = []

        for variant in self.detection_1 + self.detection_2:
            haplotype = variant[0]
            location = variant[1]

            # first variant
            if len(self.condensed) == 0:
                self.condensed.append([haplotype, self.rec_1.reference_start, None])
            
            # different haplotype
            if self.condensed[-1][0] != haplotype:
                self.condensed[-1][2] = location
                self.condensed.append([haplotype, location, None]) 
            
            # last variant
            if len(self.detection_2) == 0:
                if variant == self.detection_1[-1]:
                    self.condensed[-1][2] = self.rec_1.reference_start + self.rec_1.query_alignment_length
            else:
                if variant == self.detection_2[-1]:
                    self.condensed[-1][2] = self.rec_2.reference_start + self.rec_2.query_alignment_length

        # create condensed
        haplotypes = [tupl[0] for tupl in self.condensed]

        # classify condensed
        if 'N' in haplotypes:
            self.call = 'no_match'
        elif len(haplotypes) == 2:
            self.call = 'ambiguous_cross_over'
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

        # created masked_condensed
        # [(haplotype, beginning, end), ...]
        self.masked_condensed = []

        for phase in self.condensed:
            # one phase contains both mask start and end
            if phase[1] < mask_start and mask_end < phase[2]:
                self.masked_condensed.append((phase[0], mask_start, mask_end))
            # phase contains only mask start
            elif phase[1] < mask_start and mask_start < phase[2]:
                self.masked_condensed.append((phase[0], mask_start, phase[2]))
            # phase contains only mask end
            elif phase[1] < mask_end and mask_end < phase[2]:
                self.masked_condensed.append((phase[0], phase[1], mask_end))
            # phase is in the middle
            elif mask_start < phase[1] and phase[2] < mask_end:
                self.masked_condensed.append(phase)
                     

        # classify masked_condensed
        haplotypes = [tupl[0] for tupl in self.masked_condensed]

        if 'N' in haplotypes:
            self.masked_call = 'no_match'
        elif len(haplotypes) == 2:
            self.masked_call = 'unambiguous_cross_over'
        elif len(haplotypes) == 3:
            self.masked_call = 'gene_conversion'
        elif len(haplotypes) > 3:
            self.masked_call = 'complex'
        else:
            self.masked_call = 'no_phase_change'
                




def pairs_creation(bam_filepath, vcf_filepath):

    bam = pysam.AlignmentFile(bam_filepath, 'r')

    pairs = []
    prev_rec = None
    for rec in bam:
        # first record
        if not prev_rec:
            prev_rec = rec
        # check if query_name pairs exist
        elif rec.query_name == prev_rec.query_name:
            pairs.append(Pair(rec, prev_rec, vcf_filepath))
            prev_rec = None

    return pairs

    
