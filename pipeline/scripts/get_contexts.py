#!/usr/bin/env python3
import os
import sys
import pysam
import argparse
import pandas as pd

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Create BED files with cytosine coords in CpG and non-CpG context'
    )
    
    parser.add_argument('reference', help='Input fasta reference file')
    parser.add_argument('panel', help='Input BED for panel coordinates.')
    parser.add_argument('cpg', help='Output CpG BED')
    parser.add_argument('conversion', help='Output non-CpG cytosines BED')
    
    args = parser.parse_args()

    if not os.path.exists(args.reference) or not os.path.exists(args.panel):
        not_existing_path = args.reference if not os.path.exists(args.reference) else args.panel
        print(f"Error: Input file not found: {not_existing_path}")
        sys.exit(1)
            
    with pysam.FastaFile(args.reference) as reference:
        panel = pd.read_csv(args.panel, sep = "\t", header = None)
        j = 0
        k = 0
        cpg_context = pd.DataFrame(columns = ["chr","start","end"])
        conversion_context = pd.DataFrame(columns = ["chr","start","end"])
        for i in panel.index:
            chr = panel.iloc[i,0]
            start = panel.iloc[i,1]
            end = panel.iloc[i,2]
            seq = reference.fetch(chr, start, end) #See BED specification. Region coords are first 3 columns
            for pos in range(0,len(seq)-1):
                if seq[pos] == 'C':
                    if seq[pos+1] == 'G' or seq[max(pos-1,0)] == 'G':
                        cpg_context.loc[j] = [chr, start+pos, start+pos+1]
                        j += 1
                    else:
                        conversion_context.loc[k] = [chr, start+pos, start+pos+1]
                        k += 1
    cpg_context.to_csv(args.cpg,sep = "\t", header = False, index = False)
    conversion_context.to_csv(args.conversion,sep = "\t", header = False, index = False)
