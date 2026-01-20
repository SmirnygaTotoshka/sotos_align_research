#!/usr/bin/env python3

import os
import sys
import pysam
import argparse
import numpy as np
import pandas as pd
from collections import namedtuple

def get_seq_metrics(bam_path, contig, start, end, depth_threshold = 30):
    with pysam.AlignmentFile(bam_path,"rb") as bam_file:
        alignments = bam_file.fetch(contig, start, end)
        total_positions = end - start
        depths = np.zeros(total_positions, dtype=int)
        total_reads = 0
        
        for a in alignments:
            if a.is_unmapped or a.is_qcfail or a.is_supplementary or a.is_secondary:
                continue
            total_reads += 1
            read_positions = a.get_reference_positions(full_length=False)
            for ref_pos in read_positions:
                if start <= ref_pos < end:
                    depths[ref_pos - start] += 1

        Metrics = namedtuple('SeqMetrics', ['mean_depth', "coverage_width", "uniformity", "total_reads"])
        mean_depth = round(depths.mean(),2)
        return Metrics(
            mean_depth,
            round(np.sum(depths > depth_threshold) / (total_positions) * 100,2),
            round(np.sum(depths > 0.2 * mean_depth) / (total_positions) * 100,2),
            total_reads
        )
    
def get_methylation_metrics(bam_path, cpg_context_path, conversion_context_path, contig, start, end):
    cpg_context = pd.read_csv(cpg_context_path,sep = "\t",header = None)
    cpg_context = cpg_context.loc[(cpg_context[0] == contig)&(cpg_context[1]>=start)&(cpg_context[1] < end)]
    conversion_context = pd.read_csv(conversion_context_path,sep = "\t",header = None)
    conversion_context = conversion_context.loc[(conversion_context[0] == contig)&(conversion_context[1]>=start)&(conversion_context[1] < end)]
    cpg_context["T_count"] = 0 
    cpg_context["C_count"] = 0
    conversion_context["T_count"] = 0 
    conversion_context["C_count"] = 0
    with pysam.AlignmentFile(bam_path,"rb") as bam_file:
        pileupcolumns = bam_file.pileup(contig, start, end, truncate=True, stepper='all')
        for pileupcolumn in pileupcolumns:
            sequences = pileupcolumn.get_query_sequences()
            if pileupcolumn.reference_pos in cpg_context[1].values:
                cpg_context.loc[cpg_context[1] == pileupcolumn.reference_pos,"C_count"] = sequences.count('c') + sequences.count('C')
                cpg_context.loc[cpg_context[1] == pileupcolumn.reference_pos,"T_count"] = sequences.count('t') + sequences.count('T')
            elif pileupcolumn.reference_pos in conversion_context[1].values:
                conversion_context.loc[conversion_context[1] == pileupcolumn.reference_pos,"C_count"] = sequences.count('c') + sequences.count('C')
                conversion_context.loc[conversion_context[1] == pileupcolumn.reference_pos,"T_count"] = sequences.count('t') + sequences.count('T')
    mean_methylation = round((cpg_context["C_count"] / (cpg_context["C_count"] + cpg_context["T_count"])).mean(skipna=True) * 100,2)
    mean_conversion = round((conversion_context["T_count"] / (conversion_context["C_count"] + conversion_context["T_count"])).mean() * 100,2)
    Metrics = namedtuple('MethylMetrics', ['mean_methylation', "mean_conversion"])
    return Metrics(mean_methylation, mean_conversion)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Calculate statistics by locus after bsbolt align and methylation calling.'
    )
    
    parser.add_argument('--input', required=True, help='Input fasta reference file')
    parser.add_argument('--panel', required=True, help='Input BED for panel coordinates.')
    parser.add_argument('--cpg', required=True, help='Output CpG BED')
    parser.add_argument('--conversion', required=True, help='Output non-CpG cytosines BED')
    parser.add_argument('--cgmap', required=True, help='Output non-CpG cytosines BED')
    parser.add_argument('--depth_threshold', default=30, type=int, help='Output non-CpG cytosines BED')
    parser.add_argument('--output', default="statistics.csv", help='Output non-CpG cytosines BED')
    
    args = parser.parse_args()
    
    for path in [args.input,args.panel,args.cpg,args.conversion,args.cgmap]:
        if not os.path.exists(path):
            print(f"Error: Input file not found: {path}")
            sys.exit(1)

    panel = pd.read_csv(args.panel,sep = "\t",header = None)
    panel["region"] = panel[0] + ":" + panel[1].astype(str) + "-" + panel[2].astype(str)

    for i in panel.index:
        chr = panel.iloc[i,0]
        start = panel.iloc[i,1]
        end = panel.iloc[i,2]
        metrics = get_seq_metrics(args.input, chr, start, end,depth_threshold=args.depth_threshold)
        methylation_metrics = get_methylation_metrics(args.input, args.cpg, args.conversion, chr, start, end)
        panel.loc[i,"mean_depth"] = metrics.mean_depth
        panel.loc[i,"coverage_width"] = metrics.coverage_width
        panel.loc[i,"uniformity"] = metrics.uniformity
        panel.loc[i,"total_reads"] = metrics.total_reads
        panel.loc[i,"mean_methylation"] = methylation_metrics.mean_methylation
        panel.loc[i,"mean_conversion"] = methylation_metrics.mean_conversion
    
    cgmap = pd.read_csv(args.cgmap, sep = "\t", header = None)
    cgmap.columns = ["chr",'nucleotide','pos','context','subcontext','value','n_meth_bases','depth']
    for i in cgmap.index:
        for j in panel.index:
            if cgmap.loc[i,"chr"] == panel.iloc[j,0] and cgmap.loc[i,"pos"] >= panel.iloc[j,1] and cgmap.loc[i,"pos"] < panel.iloc[j,2]:
                cgmap.loc[i,"region"] = panel.loc[j,"region"]
                
    meth = cgmap.query("subcontext == 'CG'")
    conversion = cgmap.query("subcontext != 'CG'")
    m1 = meth.groupby('region').agg(meth_by_caller = pd.NamedAgg('value', aggfunc='mean'))
    m2 = conversion.groupby('region').agg(conversion_by_caller = pd.NamedAgg('value', aggfunc='mean'))
    p1 = pd.merge(panel, m1, on='region', how='inner')
    p2 = pd.merge(p1, m2, on='region', how='inner')
    p2.to_csv(args.output, sep = ";", index = False)
