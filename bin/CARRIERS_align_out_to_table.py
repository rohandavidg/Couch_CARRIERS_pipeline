#!/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bcbio/anaconda/bin/python

"""
script to generate a table from the align out files in the 
number directory
"""

import argparse
import csv
import pprint
import subprocess
import collections
from collections import defaultdict
import pandas as pd
import matplotlib 


samtools = "/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bin/samtools"

def main():
    args = parse_args()
    run(args.align_out_list)


def run(align_out_list):
    coverage_table_dict, carrier_run_lanes = align_out_parse(align_out_list)
    write_out_bam_list_by_lane = generate_lane_bam_list(coverage_table_dict, carrier_run_lanes)
    create_dataframe_output = create_table(coverage_table_dict)

def parse_args():
    """
    parsing arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-l', dest='align_out_list', type=argparse.FileType('r'),
                        help="list of align out files",
                        required=True)
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


def generate_lane_bam_list(metrics_table_dict, run_lanes):
    for key, value in metrics_table_dict.items():
        if value[0] == run_lanes[0]:
            

def create_table(metrics_table_dict):
    with open("CARRIERS_coverage_metric.tsv", "wa+") as fout:
        fout.write("sample\tlanes\tIndexes\tTotalReads\tTotalMappedReads\tDuplicateReads\tRealignedReads\tReadsInCaptureRegion\n")
        for key, value in metrics_table_dict.items():
            fout.write(key + "\t" + "\t".join(i for i in value) + "\n")
            

if __name__ == "__main__":
    main()
