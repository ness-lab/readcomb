#!/usr/bin/env python3
"""
bamprep.py - prepare SAM/BAM file for phase change detection

this script will automatically try to detect your samtools binary by default -
the one it uses can be found by running `which samtools` in the command line
(alternatively, you could symlink your preferred version of samtools into a bin/ folder
in your project dir to keep paths clean, and then provide `./bin/samtools`
for the path - recommended)
"""

import os
import re
import sys
import argparse
import subprocess

class ReadcombParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write(f'error: {message}\n')
        self.print_help()
        sys.exit(1)

def arg_parser():
    parser = ReadcombParser(
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
    parser.add_argument('--version', action='version', version='readcomb 0.1.3')

    return parser

def bamprep(args):
    """
    Prepare BAM for phase change detection with readcomb.

    Script will:
    1. Apply filters (proper pair, not secondary, not supplementary)
    2. Sort by read name instead of position (samtools sort)
    3. Create .bai file (samtools index)

    The output file will be named similarly to the input - only
    the extension will have been changed.

    Parameters
    ----------
    args : Namespace
        Namespace containing all user given arguements compiled by arg_parse()
    """
    # check extension
    match = re.search(r'\.[sb]am$', args.bam)
    if match:
        extension = match.group()
        basename = os.path.basename(args.bam[:args.bam.find(extension)])
    else:
        raise ValueError('input file is not SAM or BAM - make sure extension is correct')

    # get samtools path if not specified
    if not args.samtools:
        samtools_check_cmd = 'which samtools'
        proc = subprocess.Popen(samtools_check_cmd.split(' '),
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, _ = proc.communicate()
        if not out:
            raise RuntimeError("""ERROR: samtools not found. Please
                    either ensure samtools is in your PATH or specify a
                    samtools binary using the --samtools argument""")
        else:
            args.samtools = out.decode('utf-8').rstrip('\n')

    # format args
    if args.outdir:
        if not args.outdir.endswith('/'):
            args.outdir += '/'
    elif not args.outdir:
        args.outdir = './'
        print('[readcomb] no outdir specified, saving to current working directory')

    # filter out non proper pair, supplementary, and secondary reads
    print('[readcomb] Filtering...')
    filter_pos_cmd = f'{args.samtools} view -O bam '
    filter_pos_cmd += '-f 0x2 -F 0x100 -F 0x800 '
    filter_pos_cmd += f'-o {args.outdir}{basename}.temp.filtered.bam {args.bam}'
    print(f'[readcomb] cmd: {filter_pos_cmd}')
    proc = subprocess.run(filter_pos_cmd.split(' '), check=True)

    if not args.no_progress:
    # sort positionally for bai file creation
        print('[readcomb] Positionally sorting for bai creation...')
        sort_pos_cmd = f'{args.samtools} sort -O bam '
        if args.threads:
            sort_pos_cmd += f'-@{args.threads} '
        sort_pos_cmd += f'-o {args.outdir}{basename}.sorted.bam '
        sort_pos_cmd += f'{args.outdir}{basename}.temp.filtered.bam'
        print('[readcomb] Sorting positionally for .bai creation...')
        print(f'[readcomb] cmd: {sort_pos_cmd}')
        proc = subprocess.run(sort_pos_cmd.split(' '), check=True)

        # create bai file
        if args.index_csi:
            idx_type = '.csi'
        else:
            idx_type = '.bai'
        print(f'[readcomb] Creating {idx_type} index file...')
        index_cmd = f'{args.samtools} index {args.outdir}{basename}.sorted.bam '
        index_cmd += f'{args.outdir}{basename}.sorted{idx_type}'
        print(f'[readcomb] cmd: {index_cmd}')
        proc = subprocess.run(index_cmd.split(' '), check=True)

        # remove positionally sorted file
        try:
            os.remove(f'{args.outdir}{basename}.sorted.bam')
        except OSError as e:
            raise FileNotFoundError(f'Unable to delete {e.filename} - {e.strerror}') from e

    # sort by name and ensure file is in BAM format
    print('[readcomb] Creating read name sorted file...')
    sort_cmd = f'{args.samtools} sort -n -O bam '
    if args.threads:
        sort_cmd += f'-@{args.threads} '
    sort_cmd += f'-o {args.outdir}{basename}.sorted.bam {args.outdir}{basename}.temp.filtered.bam'
    print(f'[readcomb] cmd: {sort_cmd}')
    proc = subprocess.run(sort_cmd.split(' '), check=True)

    # remove filtered unsorted file
    print('[readcomb] Removing temp files...')
    try:
        os.remove(f'{args.outdir}{basename}.temp.filtered.bam')
    except OSError as e:
        raise FileNotFoundError(f'Unable to delete {e.filename} - {e.strerror}')

def main():
    parser = arg_parser()
    # Namespace of arguments
    args = parser.parse_args()

    if args.samtools:
        if not os.path.isfile(args.samtools):
            raise ValueError(f'ERROR: No samtools binary found at {args.samtools}')
    bamprep(args)
    print('[readcomb] Complete.')


if __name__ == '__main__':
    main()
