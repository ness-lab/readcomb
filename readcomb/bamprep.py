#!/usr/bin/env python3
"""
bamprep.py - prepare SAM/BAM file for phase change detection

TODO: add explicit check for samtools in userpath
-> raise informative error if not present

this script will automatically try to detect your samtools binary by default - 
the one it uses can be found by running `which samtools` in the command line
(alternatively, you could symlink your preferred version of samtools into a bin/ folder
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
        usage='readcomb-bamprep [options]')

    parser.add_argument('-b', '--bam', required=True,
        type=str, help='BAM file to prepare')
    parser.add_argument('-s', '--samtools', required=False,
        type=str, help='Path to samtools binary [optional]')
    parser.add_argument('-t', '--threads', required=False,
        type=int, help='Number of threads to use [optional]')
    parser.add_argument('--index_csi', required=False,
        action='store_true', help='Use CSI index instead of BAI [optional]')
    parser.add_argument('--no_progress', required=False,
        action = 'store_true', help='Skip positional sorting (faster, but disables progress \
        bars in filtering step).')
    parser.add_argument('-o', '--outdir', required=False,
        type=str, help='Directory to write to [optional] (default: current dir)')

    args = parser.parse_args()

    return args.bam, args.samtools, args.threads, args.index_csi, \
        args.no_progress, args.outdir

def bamprep(fname, samtools, threads, index_csi, no_progress, outdir):
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

    # get samtools path if not specified
    if not samtools:
        samtools_check_cmd = 'which samtools'
        proc = subprocess.Popen(samtools_check_cmd.split(' '),
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate()
        if not out:
            raise RuntimeError("""ERROR: samtools not found. Please
                    either ensure samtools is in your PATH or specify a
                    samtools binary using the --samtools argument""")
        else:
            samtools = out.decode('utf-8')

    # format args
    if outdir:
        if not outdir.endswith('/'):
            outdir += '/'
    elif not outdir:
        outdir = ''
        print('[readcomb] no outdir specified, saving to current working directory')

    # filter out non proper pair, supplementary, and secondary reads
    print('[readcomb] Filtering...') 
    filter_pos_cmd = '{s} view -O bam '.format(s=samtools)
    filter_pos_cmd += '-f 0x2 -F 0x100 -F 0x800 '
    filter_pos_cmd += '-o {base}.temp.filtered.bam {fname}'.format(
            s=samtools, out=outdir, base=basename, fname=fname)
    proc = subprocess.run(filter_pos_cmd.split(' '))
    proc.check_returncode() # raises CalledProcessError if proc failed

    if not no_progress:
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
    bam, samtools, threads, index_csi, no_progress, outdir = arg_parser()
    if samtools:
        if not os.path.isfile(samtools):
            raise ValueError('ERROR: No samtools binary found at {}'.format(samtools))
    bamprep(bam, samtools, threads, index_csi, no_progress, outdir)
    print('[readcomb] Complete.')

if __name__ == '__main__':
    main()

        

