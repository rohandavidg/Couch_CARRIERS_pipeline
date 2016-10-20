#!/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bcbio/anaconda/bin/python

"""
CARRIERS script for primer performance and
ratio of high to low coverage
"""


import argparse
import csv
import pandas as pd
import numpy as np 
from Bio.Seq import Seq
import pylab as pl
import matplotlib.pyplot as plt
from Bio.SeqUtils import GC
from Bio import SeqIO
import os
import subprocess
from collections import defaultdict
import pprint
from pyfaidx import Fasta
import pysam
import pybedtools
from pybedtools import BedTool
from csv import DictWriter

samtools = "/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bin/samtools"
genes = Fasta('/data2/bsi/reference/sequence/human/ncbi/hg19/allchr.fa')
sortbed = "/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bin/sortBed"

def main():
    args = parse_args()
    run(args.mayo_panel_primer, args.list_bam, args.annotated_bed)


def run(mayo_panel_primer, list_bam, annotated_bed):
    samtools_loc5_position, samtools_faidx_dict =  parse_primer(mayo_panel_primer)    
    primer_performance_dict, run_sample_names = calculate_loc5_count(list_bam,
                                                                     samtools_loc5_position)
    loc5_region_dict = find_upstream_loc5(samtools_faidx_dict)
    region_pos_annotation_dict = parse_annotated_bed_file(samtools_loc5_position, annotated_bed)
    intermidiate_file =  generate_output(run_sample_names, primer_performance_dict,
                                         loc5_region_dict, region_pos_annotation_dict)
    primer_count_file = write_out_primer_count()

def parse_args():
    """
    parsing arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-p', dest='mayo_panel_primer',
                       help='mayo panel primer excel file',
                       required=True)
    parser.add_argument('-l', dest='list_bam', type=argparse.FileType('r'),
                        help="list of bams in the run",
                        required=True)
    parser.add_argument('-b', dest='annotated_bed',
                        help='CARRIERS annoated bed file', required=True)
    args = parser.parse_args()
    return args


def parse_primer(primer_file):
    df = pd.read_excel(primer_file, 0)
    df_forward = df[(df['strand'] == 0)]
    df_reverse = df[(df['strand'] == 1)]
    df_forward['loc5'] = df_forward.loc[:,('loc5')] + 1
    df_reverse['loc5'] = df_reverse.loc[:,('loc5')] + 1
    df_forward['loc3'] = df_forward['loc5'] + 150
    df_reverse['loc3'] = df_reverse['loc5'] - 150
    df_forward['loc5'] = df_forward['loc5'].astype(int)
    df_reverse['loc5'] = df_reverse['loc5'].astype(int)
    df_forward['loc3'] = df_forward['loc3'].astype(int)
    df_reverse['loc3'] = df_reverse['loc3'].astype(int)
    df_out = pd.concat([df_forward, df_reverse], axis=0)
    gc_count_dict = {i:GC(i) for i in df_out.loc[:,('primer')]}
    df_out['primer_gc_count'] = df_out['primer'].map(gc_count_dict)
    result = df_out.loc[:,('chrom', 'loc5', 'loc3', 'strand','primer', 'primer_gc_count')]
    mayo_loc5 = result.loc[:,['chrom', 'loc5', 'strand']]
    loc5_querry = mayo_loc5['chrom'].map(str) + ":" + mayo_loc5['loc5'].map(str)
    loc5_list = loc5_querry.tolist()
    samtools_query = {str(k):samp_querry(a,k,c).encode('ascii', 'ignore') for a,k,c in zip(result['chrom'],result['loc5'],result['loc3'])}
    result.to_csv('qiagen_primer_24_samples_gc.tsv', index=False, sep='\t', encoding='utf-8')
    return loc5_list, samtools_query
#    mayo_loc5.to_csv('mayo_primer_loc5.bed', index=False, sep='\t', encoding='utf-8')


def samp_querry(a,b,c):
    if b < c:
        value = a + ":" + str(b) + "-" + str(c)
        return value
    else:
        value = a + ':' + str(c) + "-" + str(b)
        return value


def calculate_loc5_count(bam_list, loc5_list):
    sample_loc5_count = defaultdict(list)
    sample_names = []
    for bam in bam_list.readlines():
        new_bam = bam.strip()
        sample_name = new_bam.split("/")[-1].split(".")[0].encode('ascii', 'ignore')
        sample_names.append(sample_name)
        for position in loc5_list:
            key = position.split(":")[1]
            count_line = pysam.view("-c", new_bam, position)
            for line in count_line.splitlines():
                if line.strip().encode('ascii', 'ignore'):
                    sample_loc5_count[key].append({sample_name:line})
    return sample_loc5_count, sample_names


def find_upstream_loc5(samtools_query):
    region_seq_dict = {}
    for key, value in samtools_query.items():
        chrom = value.split(":")[0]
        start = int(value.split(":")[1].split("-")[0])
        stop = int(value.split(":")[1].split("-")[1])
        region = genes[chrom][start:stop]
        region_gc_content = GC(str(region))
        region_seq_dict[key] =  [region.seq.encode('ascii', 'ignore'), region_gc_content]
    return region_seq_dict


def create_temp_bed_file(loc5_list, annotated_bed):
    with open("temp.bed", 'wa') as fout:
        for pos in loc5_list:
            chrom = pos.split(":")[0]
            start = pos.split(":")[1]
            stop =  pos.split(":")[1]
            fout.write(chrom + "\t" + start + "\t" + stop + "\n")


def parse_annotated_bed_file(loc5_list, annotated_bed):
    pos_annotate_dict = {}
    if os.path.isfile("temp.bed"):
        os.remove("temp.bed")
        create_temp_bed = create_temp_bed_file(loc5_list, annotated_bed)
        cmd = sortbed + " -i " + 'temp.bed' + " > temp.sort.bed"
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        create_sort_bed,err2 = p.communicate()
    else:
        create_temp_bed = create_temp_bed_file(loc5_list, annotated_bed)
        cmd = sortbed + " -i " + 'temp.bed' + " > temp.sort.bed"
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        create_sort_bed,err2 = p.communicate()


    temp_sort_bed = os.path.abspath("temp.sort.bed")

    a = BedTool(temp_sort_bed)
    b = a.closest(annotated_bed, stream=True)
    for region in b:
        pos = region[1].encode('ascii', 'ignore')
        regions = region[6].encode('ascii', 'ignore')
        pos_annotate_dict[pos] = regions
    return pos_annotate_dict

 
def generate_output(sample_names, sample_loc5_count, region_seq_dict, pos_annotate_dict):
    with open ('outfile.txt', 'wb') as fout:
        fout.write("loc5\ttarget_seq\ttarget_gc\ttarget_gene\t" + "\t".join(str(i) for i in sample_names) + "\n") 
        for key, value in pos_annotate_dict.items():
#        if region_seq_dict[key]:# and sample_loc5_count[key]:
            if region_seq_dict[key]:
                primer_count =  sample_loc5_count[key]
                fout.write(str(key) + "\t" + str(region_seq_dict[key][0]) + "\t" +  str(region_seq_dict[key][1]) + "\t" +  str(value) + "\t" +  "\t".join([str(v) for k in primer_count for key, v in k.items() if key in sample_names]) + "\n")



def write_out_primer_count():
    primer_count_result =  os.path.abspath('outfile.txt')
    primer_metrics = os.path.abspath('qiagen_primer_24_samples_gc.tsv')
    df_a = pd.DataFrame.from_csv(primer_count_result, sep='\t', index_col=None)
#    print df_a
    df_b =  pd.DataFrame.from_csv(primer_metrics, sep='\t',index_col=None)
    headers_a = df_a.dtypes.index
    headers_b = df_b.dtypes.index

    new_df = pd.merge(df_b, df_a, on='loc5')
    new_df.to_csv('results.tsv', index=False, sep='\t', encoding='utf-8')   


if __name__ == "__main__":
    main()

