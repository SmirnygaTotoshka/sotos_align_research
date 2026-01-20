#!/usr/bin/env python3

import os
import sys
import pysam
import argparse


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Clean single reads from paired name-sorted BAM file'
    )
    
    parser.add_argument('input', help='Input name-sorted BAM file with paired reads with singletons after filtering')
    parser.add_argument('output', help='Output BAM file')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {file_path}")
        sys.exit(1)
            
    with pysam.AlignmentFile(args.input, 'rb') as input_bamfile:
        with pysam.AlignmentFile(args.output, 'wb', template = input_bamfile) as output_bamfile:
            
            assert input_bamfile.header.as_dict()["HD"]["SO"] == "queryname", "Input file doesn`t name-sorted"        
            prev_read = None
            i = 0
            for read in input_bamfile.fetch(until_eof=True):
                assert read.is_paired, "This script only for paired reads"
                i+=1
                if i % 2 == 0:  # Process two lines at a time
                    if prev_read.query_name == read.query_name:
                        output_bamfile.write(prev_read)
                        output_bamfile.write(read)
                    else:
                        print(f'Error: Read IDs do not match at lines {i-1} and {i}')
                        print(f'Line {i-1}: {prev_read}')
                        print(f'Line {i}: {read}')
                        i-=1
                prev_read = read
            print(f"Total written {i} reads")

