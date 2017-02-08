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
import re
from random import randint

samtools = "/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bin/samtools"
genes = Fasta('/data2/bsi/reference/sequence/human/ncbi/hg19/allchr.fa')
sortbed = "/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bin/sortBed"

def main():
    args = parse_args()
    run(args.mayo_panel_primer, args.bam_path, args.annotated_bed, args.sample_info_list, args.bam_list)

def run(mayo_panel_primer, bam_path, annotated_bed, sample_info_list, bam_list):
    samtools_loc5_position, temp_outfile =  parse_primer(mayo_panel_primer)
    sample_name_lanes_dict = get_lanes_from_sample_info(sample_info_list)
    region_pos_annotation_dict = parse_annotated_bed_file(samtools_loc5_position, annotated_bed)
    if bam_path:
        primer_performance_dict, run_sample_name = calculate_loc5_count(bam_path,
                                                                     samtools_loc5_position)    
        print primer_performance_dict
        intermidiate_file =  generate_output(run_sample_name, primer_performance_dict,
                                             region_pos_annotation_dict, sample_name_lanes_dict)
        primer_count_file = write_out_primer_count(run_sample_name, sample_name_lanes_dict, temp_outfile)
    else:
        primer_performance_dict, run_sample_name = calculate_loc5_count_all(bam_list,
                                                                            samtools_loc5_position)
        intermidiate_file =  generate_output_all(run_sample_name, primer_performance_dict,
                                                 region_pos_annotation_dict, sample_name_lanes_dict)
        primer_count_file = write_out_primer_count_all(temp_outfile)



def parse_args():
    """
    parsing arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-p', dest='mayo_panel_primer',
                       help='mayo panel primer excel file',
                       default="/data5/bsi/epibreast/m087494.couch/Couch/Huge_Breast_VCF/Mayo_panel_primer.xlsx")
    parser.add_argument('-l', dest='bam_path',
                        help="bam path")
    parser.add_argument('-b', dest='annotated_bed',
                        help='CARRIERS annoated bed file', default="/data5/bsi/epibreast/m087494.couch/Couch/Huge_Breast_VCF/CARRIERS_PANC.targets.annotated.bed")
    parser.add_argument('-s', dest='sample_info_list', nargs='+',
                        help='sample info from the ggps run')
    parser.add_argument('-B', dest='bam_list', type=argparse.FileType('r'),
                       help="bam file with list of bams")
    args = parser.parse_args()
    return args


def parse_primer(primer_file):
    df = pd.read_excel(primer_file, 0)
    df_forward = df[(df['strand'] == 0)]
    df_reverse = df[(df['strand'] == 1)]
    df_forward.loc[:,'loc5'] = df_forward.loc[:,('loc5')] + 1
    df_reverse.loc[:,'loc5'] = df_reverse.loc[:,('loc5')] + 1
    df_forward.loc[:,'loc3'] = df_forward['loc5'] + 150
    df_reverse.loc[:,'loc3'] = df_reverse['loc5'] - 150
    df_forward.loc[:,'loc5'] = df_forward['loc5'].astype(int)
    df_reverse['loc5'] = df_reverse['loc5'].astype(int)
    df_forward['loc3'] = df_forward['loc3'].astype(int)
    df_reverse['loc3'] = df_reverse['loc3'].astype(int)
    df_out = pd.concat([df_forward, df_reverse], axis=0)
#    gc_count_dict = {i:GC(i) for i in df_out.loc[:,('primer')]}
#    df_out['primer_gc_count'] = df_out['primer'].map(gc_count_dict)
    result = df_out.loc[:,('chrom', 'loc5', 'loc3', 'strand','primer')]
    mayo_loc5 = result.loc[:,['chrom', 'loc5', 'strand']]
    tsv_outfile = 'qiagen_primer_24_samples_gc' + str(randint(0,100000)) + ".tsv"
    loc5_querry = mayo_loc5['chrom'].map(str) + ":" + mayo_loc5['loc5'].map(str)
    result.to_csv(tsv_outfile, index=False, sep='\t', encoding='utf-8')
    loc5_list = loc5_querry.tolist()
    return loc5_list, tsv_outfile


def calculate_loc5_count(bam_path, loc5_list):
    sample_loc5_count = defaultdict(list)
    new_bam = bam_path.strip()
    sample_name = new_bam.split("/")[-1].split(".")[0].encode('ascii', 'ignore')
    for position in loc5_list:
        key = position.split(":")[1]
        count_line = pysam.view("-c", new_bam, position)
        for line in count_line.splitlines():
            if line.strip().encode('ascii', 'ignore'):
                sample_loc5_count[key].append({sample_name:line})
    return sample_loc5_count, sample_name


def calculate_loc5_count_all(bam_list, loc5_list):
    sample_loc5_count = defaultdict(list)
    sample_names = set()
    for bam in bam_list.readlines():
        new_bam = bam.strip()
        sample_name = new_bam.split("/")[-1].split(".")[0].encode('ascii', 'ignore')
        sample_names.add(sample_name)
        for position in loc5_list:
            key = position.split(":")[1]
            count_line = pysam.view("-c", new_bam, position)
            for line in count_line.splitlines():
                if line.strip().encode('ascii', 'ignore'):
                    sample_loc5_count[key].append({sample_name:line})
    return sample_loc5_count, sample_names



def create_temp_bed_file(loc5_list, annotated_bed):
    a = randint(0,100000)
    temp_file = "temp_" + str(a) + ".bed"
    with open(temp_file, 'wa') as fout:
        for pos in loc5_list:
            chrom = pos.split(":")[0]
            start = pos.split(":")[1]
            stop =  pos.split(":")[1]
            fout.write(chrom + "\t" + start + "\t" + stop + "\n")
    return temp_file


def parse_annotated_bed_file(loc5_list, annotated_bed):
    if os.path.isfile("temp.bed"):
        os.remove("temp.bed")
        create_temp_bed = create_temp_bed_file(loc5_list, annotated_bed)
        output_file = create_temp_bed + ".sort.bed"
        cmd = sortbed + " -i " + temp_file + " > " + output_file 
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        create_sort_bed,err2 = p.communicate()
        chr_pos_annotate_dict = run_closedBet(output_file, annotated_bed)
        return chr_pos_annotate_dict
    else:
        create_temp_bed = create_temp_bed_file(loc5_list, annotated_bed)
        output_file = create_temp_bed + ".sort.bed"
        cmd = sortbed + " -i " + create_temp_bed + " > " + output_file 
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        create_sort_bed,err2 = p.communicate()
        chr_pos_annotate_dict = run_closedBet(output_file, annotated_bed)
        return chr_pos_annotate_dict        


def run_closedBet(output_file, annotated_bed):
    pos_annotate_dict = {}
    temp_sort_bed = os.path.abspath(output_file)
    a = BedTool(output_file)
    b = a.closest(annotated_bed, stream=True)
    for region in b:
        pos = region[1].encode('ascii', 'ignore')
        regions = region[6].encode('ascii', 'ignore')
        pos_annotate_dict[pos] = regions
    return pos_annotate_dict


def get_lanes_from_sample_info(sample_info_list):
    sample_name_lane_dict = {}
    for sample_info in sample_info_list:
        with open(sample_info) as fin:
            for raw_line in fin:
                line = raw_line.strip().split("\t")
                sample_fastq = line[0].split(":")[1]
                sample_name = sample_fastq.split("=")[0]
                fastq = sample_fastq.split("=")[1]
                flowcell_lane_index =  fastq.split(".")[1]
                m = re.search('_L(.+?)_', flowcell_lane_index)
                if m:
                    lane = m.group(1)
                    sample_name_lane_dict[sample_name] = lane
    return sample_name_lane_dict


def generate_output(sample_names, sample_loc5_count, pos_annotate_dict, sample_name_lane_dict):
    with open ('outfile_' + sample_names + "_L" + sample_name_lane_dict[sample_names] +'.txt', 'wb') as fout:
        fout.write("loc5\ttarget_gene\tlane\t" + sample_names + "\n")
        for key, value in pos_annotate_dict.items():
            if sample_loc5_count[key]:
                primer_count =  sample_loc5_count[key]
                fout.write(str(key) + "\t" + str(value) + "\t" +  sample_name_lane_dict[sample_names] + "\t" + "\t".join([str(v) for k in primer_count for key, v in k.items() if key in sample_names]) + "\n")


def write_out_primer_count(sample_name, sample_name_lane_dict, tsv_outfile):
    primer_count_result =  os.path.abspath('outfile_' + sample_name + "_L" + sample_name_lane_dict[sample_name]+ '.txt')
    print primer_count_result
    primer_metrics = os.path.abspath(tsv_outfile)
    df_a = pd.DataFrame.from_csv(primer_count_result, sep='\t', index_col=None)
    print df_a
    df_b =  pd.DataFrame.from_csv(primer_metrics, sep='\t',index_col=None)
#    print df_b
    headers_a = df_a.dtypes.index
    headers_b = df_b.dtypes.index
    new_df = pd.merge(df_b, df_a, on='loc5')
    new_df.to_csv('CARRIERS_'+ sample_name + '_L' + sample_name_lane_dict[sample_name] + '_results.tsv', index=False, sep='\t', encoding='utf-8')
#    os.remove('qiagen_primer_24_samples_gc.tsv')
#    os.remove('outfile.txt')
#    os.remove('temp.bed')


def generate_output_all(sample_names, sample_loc5_count,  pos_annotate_dict, sample_name_lane_dict):
    with open ('outfile.txt', 'wb') as fout:
        fout.write("loc5\ttarget_seq\ttarget_gene\t" + "\t".join(str(i) for i in sample_names) + "\n")
        for key, value in pos_annotate_dict.items():
            if sample_loc5_count[key]:
                primer_count = sample_loc5_count[key]
#                print str(key) + "\t" + str(value) + "\t" + "\t".join([str(v) for k in primer_count for key, v in k.items() if key in sample_names]) + "\n"
                fout.write(str(key) + "\t" + str(value) + "\t" + "\t".join([str(v) for k in primer_count for key, v in k.items() if key in sample_names]) + "\n")
                

def write_out_primer_count_all(tsv_outfile):
    primer_count_result =  os.path.abspath('outfile.txt')
    primer_metrics = os.path.abspath(tsv_outfile)
    df_a = pd.DataFrame.from_csv(primer_count_result, sep='\t', index_col=None)
#    print df_a
    df_b =  pd.DataFrame.from_csv(primer_metrics, sep='\t',index_col=None)
    headers_a = df_a.dtypes.index
    headers_b = df_b.dtypes.index

    new_df = pd.merge(df_b, df_a, on='loc5')
    new_df.to_csv('results.tsv', index=False, sep='\t', encoding='utf-8')


if __name__ == "__main__":
    main()
