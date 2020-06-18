import pysam
import os
from cyvcf2 import VCF
from cyvcf2 import Writer

# os.path.dirname(__file__) is the directory of mock_create.py
current_dir = os.path.dirname(__file__)
bam_files_path = current_dir + '/mock_sams/'
vcf_files_path = current_dir + '/'


def mock_sequence(snp_tuples, seq_len=150):
    '''
    creates mock_sequences that are all cytosine bases except bases defined in snp_tuples
    
    snp_tuples is a list of tuples that contains tuples in the form [(index, base), ...]
    where index is an integer and base is a string
    
    returns a string of seq_len (default 150)
    '''
    seq = ['C'] * seq_len
    
    for snp in snp_tuples:
        seq[snp[0]] = snp[1] 
        
    return ''.join(seq)

def mock_vcf():
    # generate gz and tabix file
    # bgzip -c file.vcf > file.vcf.gz
    # tabix -p vcf file.vcf.gz

    # VCFs: 1-based 
    data = [
        # CHROM	POS	ID REF ALT QUAL FILTER INFO	FORMAT CC2935 CC2936
        ['chromosome_1', '50', '.', 'A', 'T', '100', '.', 'AC=2;AF=0.5;AN=4', 'GT:AD:DP:GQ', '0/0:100,0:100:99', '1/1:100,0:100:99'],
        ['chromosome_1', '100', '.', 'A', 'T', '100', '.', 'AC=2;AF=0.5;AN=4', 'GT:AD:DP:GQ', '0/0:100,0:100:99', '1/1:100,0:100:99'],
        ['chromosome_1', '250', '.', 'A', 'T', '100', '.', 'AC=2;AF=0.5;AN=4', 'GT:AD:DP:GQ', '0/0:100,0:100:99', '1/1:100,0:100:99'],
        ['chromosome_1', '300', '.', 'A', 'T', '100', '.', 'AC=2;AF=0.5;AN=4', 'GT:AD:DP:GQ', '0/0:100,0:100:99', '1/1:100,0:100:99'],
    ]
    
    template = VCF(vcf_files_path + 'parental_filtered.vcf.gz')
    mock = Writer(current_dir + '/mock.vcf', template)
    mock.write_header()

    for snp in data:
        mock.write_record(mock.variant_from_string('\t'.join(snp)))

    mock.close()

    os.system('bgzip -f ' + current_dir + '/mock.vcf')
    os.system('tabix -p vcf ' + current_dir + '/mock.vcf.gz')


def mock_sam():
    template = pysam.AlignmentFile(bam_files_path + 'mock_1.sam')
    mock = pysam.AlignmentFile(current_dir + '/mock.sam', 'wh', template=template)

    '''
    QNAME/query_name - str
    FLAG
    RNAME/reference_name - str
    POS/reference_start - int
    MAPQ/mapping_quality - int 
    CIGAR/cigarstring - str 
    RNEXT/next_reference_id - str
    PNEXT/next_reference_start - int 
    TLEN/template_length - int
    SEQ/query_sequence - str
    '''

    unpaired_data = [
        # phase change 1 to 2
        {'query_name': 'seq1', 
        'reference_name': 'chromosome_1', 
        'reference_start': 0,
        'mapping_quality': 60, 
        'cigarstring': '150M', 
        'query_sequence': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
        },
        # phase change 2 to 1
        {'query_name': 'seq2', 
        'reference_name': 'chromosome_1', 
        'reference_start': 0,
        'mapping_quality': 60, 
        'cigarstring': '150M', 
        'query_sequence': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
        },
        # insertion before both SNPs, note that the reference_start is 20 instead of 0
        {'query_name': 'seq3',
         'reference_name': 'chromosome_1',
         'reference_start': 20,
         'mapping_quality': 60,
         'cigarstring': '20I130M',
         'query_sequence': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
        },
        # deletion before both SNPs
        {'query_name': 'seq4',
         'reference_name': 'chromosome_1',
         'reference_start': 0,
         'mapping_quality': 60,
         'cigarstring': '10M10D140M',
         'query_sequence': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'             
        },
        # insertion between the SNPs
        {'query_name': 'seq5',
         'reference_name': 'chromosome_1',
         'reference_start': 0,
         'mapping_quality': 60,
         'cigarstring': '60M10I80M',
         'query_sequence': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
        },
        # deletion between the SNPs
        {'query_name': 'seq6',
         'reference_name': 'chromosome_1',
         'reference_start': 0,
         'mapping_quality': 60,
         'cigarstring': '60M10D90M',
         'query_sequence': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
        },
        # soft clipping at beginning and end
        {'query_name': 'seq7',
         'reference_name': 'chromosome_1',
         'reference_start': 20,
         'mapping_quality': 60,
         'cigarstring': '20S110M20S',
         'query_sequence': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
        },
        # hard clipping at begging and end
        {'query_name': 'seq8',
         'reference_name': 'chromosome_1',
         'reference_start': 20,
         'mapping_quality': 60,
         'cigarstring': '20H110M20H',
         'query_sequence': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
        },
    ]

    for seq in unpaired_data:
        segment = pysam.AlignedSegment(template.header)
        segment.query_name = seq['query_name']
        segment.reference_name = seq['reference_name']
        segment.reference_start = seq['reference_start']
        segment.mapping_quality = seq['mapping_quality']
        segment.cigarstring = seq['cigarstring']
        segment.query_sequence = seq['query_sequence']

        mock.write(segment)
        
    paired_data = [
        # 2 mates that are both 150M and lineage 1 to lineage 2
        {'query_name': 'mate1',
         'reference_name': 'chromosome_1',
         'reference_start': 0,
         'mapping_quality': 60,
         'cigarstring': '150M',
         'query_sequence': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
         'next_reference_id': 0,
         'next_reference_start': 200,
         'template_length': 351            
        },
        {'query_name': 'mate1',
         'reference_name': 'chromosome_1',
         'reference_start': 200,
         'mapping_quality': 60,
         'cigarstring': '150M',
         'query_sequence': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
         'next_reference_id': 0,
         'next_reference_start': 0,
         'template_length': -349            
        },
        # 2 mates that are both 150M and lineage 2 to lineage 1
        {'query_name': 'mate2',
         'reference_name': 'chromosome_1',
         'reference_start': 0,
         'mapping_quality': 60,
         'cigarstring': '150M',
         'query_sequence': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
         'next_reference_id': 0,
         'next_reference_start': 200,
         'template_length': 351            
        },
        {'query_name': 'mate2',
         'reference_name': 'chromosome_1',
         'reference_start': 200,
         'mapping_quality': 60,
         'cigarstring': '150M',
         'query_sequence': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
         'next_reference_id': 0,
         'next_reference_start': 0,
         'template_length': -349            
        },
    ]
    
    for seq in paired_data:
        segment = pysam.AlignedSegment(template.header)
        segment.query_name = seq['query_name']
        segment.reference_name = seq['reference_name']
        segment.reference_start = seq['reference_start']
        segment.mapping_quality = seq['mapping_quality']
        segment.cigarstring = seq['cigarstring']
        segment.query_sequence = seq['query_sequence']
        segment.next_reference_id = seq['next_reference_id']
        segment.next_reference_start = seq['next_reference_start']
        segment.template_length = seq['template_length']

        mock.write(segment)

if __name__ == '__main__':
    mock_vcf()
    mock_sam()
