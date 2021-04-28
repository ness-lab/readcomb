import pytest
import pysam
import os
import cyvcf2
from cyvcf2 import VCF
from readcomb import filter

@pytest.fixture
def vcf_obj():
    root_dir = os.path.dirname(os.path.abspath(__file__))
    vcf_filepath = os.path.join(root_dir, 'vcfs/filtered_full.vcf.gz')
    return VCF(vcf_filepath)

class TestCheckVariants:
    def test_no_variants(self, vcf_obj):
        variants = filter.check_variants(vcf_obj, 'chromosome_1', 1, 10)
        assert variants == []

    def test_wrong_reference(self, vcf_obj):
        variants = filter.check_variants(vcf_obj, 'nonexistent_chrom', 1, 1000)
        assert variants == []

    def test_one_variant(self, vcf_obj):
        locations = [14525]
        variants = filter.check_variants(vcf_obj, 'chromosome_1', 14500, 14600)

        for variant in variants:
            assert variant.start in locations
            assert type(variant) == cyvcf2.Variant

    def test_multiple_variants(self, vcf_obj):
        locations = [14525, 14746, 14753, 14759, 14767, 14769]
        variants = filter.check_variants(vcf_obj, 'chromosome_1', 14500, 15000)

        for variant in variants:
            assert variant.start in locations
            assert type(variant) == cyvcf2.Variant

class TestCigar:
    def test_match(self):
        input = {'query_name': 'match',
                'reference_start': 1,
                'mapping_quality': 60,
                'cigarstring': '20M',
                'query_sequence': 'ATCGATCGATCGATCGATCG'
                }

        record = pysam.AlignedSegment()
        for key in input:
            record.__setattr__(key, input[key])
        
        segment = filter.cigar(record)
        assert segment == 'ATCGATCGATCGATCGATCG'

    def test_insert(self):
        input = {'query_name': 'insert',
                'reference_start': 1,
                'mapping_quality': 60,
                'cigarstring': '10M2I8M',
                'query_sequence': 'ATCGATCGATCGATCGATCG'
                }

        record = pysam.AlignedSegment()
        for key in input:
            record.__setattr__(key, input[key])
        
        segment = filter.cigar(record)
        assert segment == 'ATCGATCGATATCGATCG'

    def test_deletion(self):
        input = {'query_name': 'deletion',
                'reference_start': 1,
                'mapping_quality': 60,
                'cigarstring': '10M2D10M',
                'query_sequence': 'ATCGATCGATCGATCGATCG'
                }

        record = pysam.AlignedSegment()
        for key in input:
            record.__setattr__(key, input[key])
        
        segment = filter.cigar(record)
        assert segment == 'ATCGATCGAT--CGATCGATCG'

    def test_soft_clipping(self):
        input = {'query_name': 'soft_clipping',
                'reference_start': 1,
                'mapping_quality': 60,
                'cigarstring': '10S10M',
                'query_sequence': 'ATCGATCGATCGATCGATCG'
                }

        record = pysam.AlignedSegment()
        for key in input:
            record.__setattr__(key, input[key])
        
        segment = filter.cigar(record)
        assert segment == 'CGATCGATCG'

    def test_hard_clipping(self):
        input = {'query_name': 'hard_clipping',
                'reference_start': 1,
                'mapping_quality': 60,
                'cigarstring': '10H10M',
                'query_sequence': 'CGATCGATCG'
                }

        record = pysam.AlignedSegment()
        for key in input:
            record.__setattr__(key, input[key])
        
        segment = filter.cigar(record)
        assert segment == 'CGATCGATCG'

class TestPhaseDetection:
    def test_no_variants(self):
        input = {'query_name': 'match',
                'reference_start': 1,
                'mapping_quality': 60,
                'cigarstring': '20M',
                'query_sequence': 'ATCGATCGATCGATCGATCG'
                }

        record = pysam.AlignedSegment()
        for key in input:
            record.__setattr__(key, input[key])
        
        segment = 'ATCGATCGATCGATCGATCG'
        variants = []

        haplotypes = filter.phase_detection(variants, segment, record)
        assert haplotypes == []


