#!/usr/bin/env python3
"""
bamprep.py - prepare SAM/BAM file for phase change detection

TODO: add explicit check for samtools in userpath
-> raise informative error if not present

this script will require a path to your samtools binary - 
this can be found by running `which samtools` in the command line
(or, alternatively, you could symlink samtools into a bin/ folder
in your project dir to keep paths clean, and then provide `./bin/samtools`
for the path - recommended)
"""

import os
import re
import argparse
import subprocess
import pysam
from tqdm import tqdm

def arg_parser():
    parser = argparse.ArgumentParser(
        description='prep SAM/BAM for phase change detection', 
        usage='python3.5 bamprep.py [options]')

    parser.add_argument('-f', '--fname', required=True,
        type=str, help='File to prepare')
    parser.add_argument('-s', '--samtools', required=True,
        type=str, help='Path to samtools binary')
    parser.add_argument('-t', '--threads', required=False,
        type=int, help='Number of threads to use [optional]')
    parser.add_argument('--index_csi', required=False,
        action='store_true', help='Use CSI index instead of BAI [optional]')
    parser.add_argument('-o', '--outdir', required=False,
        type=str, help='Directory to write to [optional] (default: current dir)')

    args = parser.parse_args()

    return args.fname, args.samtools, args.threads, args.index_csi, args.outdir

def bamprep(fname, samtools, threads, index_csi, outdir):
    """
    Prepare BAM for phase change detection with readcomb.
    Script will:
    1. Apply filters (proper pair, not secondary, not supplementary)
    2. Sort by read name instead of position (samtools sort)
    3. Create .bai file (samtools index)

    The output file will be named similarly to the input - only
    the extension will have been changed. 
    
    Parameters
    -------
    fname : str
        path to input SAM/BAM
    outdir : str
        path to directory to write to (optional)

    Returns
    -------
    None
    """
    # check extension
    match = re.search('\.[sb]am$', fname)
    if match:
        extension = match.group()
        basename = fname[:fname.find(extension)]
    else:
        raise ValueError('input file is not SAM or BAM - make sure extension is correct')

    # format args
    if outdir:
        if not outdir.endswith('/'):
            outdir += '/'
    elif not outdir:
        outdir = ''
        print('[readcomb] no outdir specified, saving to current wd')

    # filter out non proper pair, supplementary, and secondary reads
    print('[readcomb] Filtering...') 
    filter_pos_cmd = '{s} view -O bam '.format(s=samtools)
    filter_pos_cmd += '-f 0x2 -F 0x100 -F 0x800 '
    filter_pos_cmd += '-o {base}.temp.filtered.bam {fname}'.format(
            s=samtools, out=outdir, base=basename, fname=fname)
    proc = subprocess.run(filter_pos_cmd.split(' '))
    proc.check_returncode() # raises CalledProcessError if proc failed

    # TODO: make positional sorting + bai creation optional - ie if no progress bars wanted
    # sort positionally for bai file creation
    print('[readcomb] Positionally sorting for bai creation...')
    sort_pos_cmd = '{s} sort -O bam '.format(s=samtools)
    if threads:
        sort_pos_cmd += '-@{n} '.format(n=threads)
    sort_pos_cmd += '-o {base}.sorted.bam {base}.temp.filtered.bam'.format(base=basename)
    proc = subprocess.run(sort_pos_cmd.split(' '))
    proc.check_returncode()

    # create bai file
    if index_csi:
        idx_type = '.csi'
    else:
        idx_type = '.bai'
    print('[readcomb] Creating {} index file...'.format(idx_type))
    index_cmd = '{s} index {base}.sorted.bam {out}{base}.sorted{idx}'.format(
            s=samtools, out=outdir, base=basename, idx=idx_type)
    proc = subprocess.run(index_cmd.split(' '))
    proc.check_returncode()

    # remove positionally sorted file
    try:
        os.remove('{base}.sorted.bam'.format(base=basename))
    except OSError as e:
        raise FileNotFoundError('Unable to delete {} - {}'.format(e.filename, e.strerror))

    # sort by name and ensure file is in BAM format
    print('[readcomb] Creating read name sorted file...')
    sort_cmd = '{s} sort -n -O bam '.format(s=samtools)
    if threads:
        sort_cmd += '-@{n} '.format(n=threads)
    sort_cmd += '-o {out}{base}.sorted.bam {base}.temp.filtered.bam'.format(
            out=outdir, base=basename)
    proc = subprocess.run(sort_cmd.split(' '))
    proc.check_returncode() 

    # remove filtered unsorted file
    print('[readcomb] Removing temp files...')
    try:
        os.remove('{base}.temp.filtered.bam'.format(base=basename))
    except OSError as e:
        raise FileNotFoundError('Unable to delete {} - {}'.format(e.filename, e.strerror))


def main():
    fname, samtools, threads, index_csi, outdir = arg_parser()
    if not os.path.isfile(samtools):
        raise ValueError('No samtools binary found at {}'.format(samtools))
    bamprep(fname, samtools, threads, index_csi, outdir)
    print('[readcomb] Complete.')

if __name__ == '__main__':
    main()

        

