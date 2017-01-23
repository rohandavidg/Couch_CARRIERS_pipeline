#!/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bcbio/anaconda/bin/python

"""
script to generate a table from the align out files in the 
number directory
"""

import os
import argparse
import csv
import pprint
import subprocess
import collections
from collections import defaultdict
import pandas as pd
import matplotlib 
from CARRIERS_calculate_median_absolute_deviation import mad
from CARRIERS_calculate_median_absolute_deviation import get_score
from decimal import Decimal
from CARRIERS_calculate_median_absolute_deviation import parse_file
import re

samtools = "/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bin/samtools"

def main():
    args = parse_args()
    run(args.align_out_list, args.CNV_txt_list, args.index_html_table, args.sample_info)


def run(align_out_list, CNV_txt_list, index_html_table, sample_info):
    if align_out_list:
        coverage_table_dict, carrier_run_lanes = align_out_parse(align_out_list)
        sample_mad_scrore_dict = compute_mad(CNV_txt_list)
        create_dataframe_output = create_table(coverage_table_dict, sample_mad_scrore_dict)
    else:
        reads_mapped_dict, sample_reads_mapped_dict, reads_mapped_list = parse_file(index_html_table)
        sample_lane_index_dict = get_lanes_from_sample_info(sample_info)
        read_metrics_dict = merge_lane_reads(sample_lane_index_dict, sample_reads_mapped_dict)
        sample_mad_scrore_dict = compute_mad(CNV_txt_list)
        create_dataframe_output = create_table(read_metrics_dict, sample_mad_scrore_dict)
#    compute_mad = create_dataframe()


def parse_args():
    """
    parsing arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-l', dest='align_out_list', type=argparse.FileType('r'),
                        help="list of align out files")
    parser.add_argument('-c', dest='CNV_txt_list', type=argparse.FileType('r'),
                        help='list of CNV txt files present in variants directory',
                        required=True)
    parser.add_argument('-i', dest='index_html_table',
                       help='coverage table copied from the index.html')
    parser.add_argument('-s', dest='sample_info', type=argparse.FileType('r'),
                       help='sample info for ggps run')
    args = parser.parse_args()
    return args


def align_out_parse(align_out_list):
    metrics_table_dict = defaultdict(list)
    run_lanes = set()
    for raw_path in align_out_list:
        raw_line = raw_path.strip()
        carrier_sample = get_sample_name(raw_line)
        with open(raw_line) as fin:
            line = fin.readlines()
            lanes = line[1].strip()
            run_lanes.add(lanes)
            Indexes = line[3].strip()
            TotalReads = line[5].strip()
            TotalMappedReads = line[7].strip()
            DuplicateReads = line[9].strip()
            RealignedReads = line[11].strip()
            ReadsInCaptureRegion = line[13].strip()
            metrics_table_dict[carrier_sample] = [lanes, Indexes, TotalReads, TotalMappedReads,
                                                  DuplicateReads, RealignedReads, ReadsInCaptureRegion]
    return metrics_table_dict, run_lanes


def get_sample_name(some_path):
    path_list = some_path.split("/")[-1]
    s_sample = path_list.split(".")[0]
    sample_name = s_sample.split("_")[1]
    return sample_name


def get_lanes_from_sample_info(sample_info):
    sample_name_lane_dict = {}
    for raw_line in sample_info:
        line = raw_line.strip().split("\t")
        sample_fastq = line[0].split(":")[1]
        sample_name = sample_fastq.split("=")[0]
        fastq = sample_fastq.split("=")[1]
        flowcell_lane_index =  fastq.split(".")[1]
        index_string = flowcell_lane_index.split("_")[3]
        index = index_string[1:]
        m = re.search('_L(.+?)_', flowcell_lane_index)
        if m:
            lane = m.group(1)
            sample_name_lane_dict[sample_name] = [lane, index]
    return sample_name_lane_dict


def merge_lane_reads(sample_name_lane_dict, sample_reads_mapped_dict):
    merge_sample_lane_index_dict = {}
    for key, value in sample_name_lane_dict.items():
        if sample_reads_mapped_dict[key]:
            sample = key.split("_")[1]
            merge_sample_lane_index_dict[sample] = [value[0], value[1],  int(sample_reads_mapped_dict[key][0].replace(',','')), 
                                                    int(sample_reads_mapped_dict[key][1].replace(',','')),int(sample_reads_mapped_dict[key][3].replace(',','')), 
                                                    int(sample_reads_mapped_dict[key][5].replace(',',''))]
    return merge_sample_lane_index_dict
            
# [total_reads, TotalMappedReads, total_reads_per,RealignedReads,RealignedReads_per, ReadsInCaptureRegion, ReadsInCaptureRegion_per]

def compute_mad(CNV_txt_list):
    CNV_mad_score = {}
    for cnv_file in CNV_txt_list:
        filename = cnv_file.strip()
        base = os.path.basename(filename).split("(")[0].strip()
        sample_name = base.split("_")[1]
        df = pd.read_csv(filename, sep='\t')
        df['absolute_deviation'] = abs(df['CNV.log2ratio'] - df['CNV.log2ratio'].median())
        median_absolute_deviation = df['absolute_deviation'].median()
        CNV_mad_score[sample_name] = round(median_absolute_deviation,4)
    return CNV_mad_score


def create_table(metrics_table_dict, CNV_mad_score):
    with open("CARRIERS_coverage_metric.tsv", "wa+") as fout:
        fout.write("sample\tlanes\tIndexes\tTotalReads\tTotalMappedReads\tDuplicateReads\tRealignedReads\tReadsInCaptureRegion\tAbsolute_median_deviation\n")
        for key, value in metrics_table_dict.items():
            if CNV_mad_score[key]:
                fout.write(key + "\t" + "\t".join(str(i) for i in value) + '\t' + str(CNV_mad_score[key]) + "\n")

            
if __name__ == "__main__":
    main()
