#!/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bcbio/anaconda/bin/python


"""
This script creates a dataframe with index and mad scores
"""

import sys
import numpy as np
from numpy import mean, absolute
import pandas as pd
import os
import argparse
import openpyxl
import csv

def main():
    args = parse_args()
    run(args.input_file, args.submission_excel_list)


def run(input_file, submission_excel_list):
    reads_mapped_dict, sample_reads_mapped_dict, reads_mapped_list = parse_file(input_file)
    create_submssion_csv_list = create_csv_from_excel(submission_excel_list)
    sample_index_dict = map_sample_index(create_submssion_csv_list)
#    print sample_index_dict
    compute_mad = get_score(reads_mapped_list)
    out_file = join_files(sample_reads_mapped_dict, sample_index_dict)
    compute_mad = create_dataframe()


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', dest='input_file',
                        help="coverage table in the index.html in tsv format",
                        required=True)
    parser.add_argument('-l', dest='submission_excel_list', nargs="+",
                        help="submission excel manifest submitted to BAP",
                        required=True)
    args = parser.parse_args()
    return args


def parse_file(input_file):
    sample_reads_mapped_dict = {}
    reads_mapped_dict  = {}
    reads_mapped_list = []
    with open(input_file) as fin:
        next(fin)
        for raw_line in fin:
            line = raw_line.strip().split("\t")
            reads_maped_to_target = line[-1].split(" ")[0]
            reads_maped = reads_maped_to_target.replace(",","")
            reads_mapped_list.append(reads_maped)
            sample = line[0]
            total_reads = line[1]
            TotalMappedReads = line[2].split(" ")[0]
            total_reads_per =  line[2].split(" ")[1]            
            RealignedReads = line[4].split(" ")[0]
            RealignedReads_per = line[4].split(" ")[1]
            ReadsInCaptureRegion = line[5].split(" ")[0]
            ReadsInCaptureRegion_per = line[5].split(" ")[1]
            sample_reads_mapped_dict[sample] = [total_reads, TotalMappedReads, total_reads_per,RealignedReads,
                                                RealignedReads_per, ReadsInCaptureRegion, ReadsInCaptureRegion_per]
            reads_mapped_dict[sample] = reads_maped
    return reads_mapped_dict, sample_reads_mapped_dict, reads_mapped_list


def parse_info(info_file):
    sample_info_dict = {}
    with open(info_file) as fin:
        next(fin)
        for raw_line in fin:
            line = raw_line.strip().split("\t")
            sample = line[1]
            new_sample =  "s_" + sample.replace("/", "_")
            sample_info_dict[new_sample] = line[2:]
    return sample_info_dict


def create_csv_from_excel(submission_excel_list):
    index_tmp_file = []
    for i, submission_excel in enumerate(submission_excel_list):
        csv_filename = "carriers_pancreas" + str(i + 1) +"_tmp.csv"
        to_csv = convert_excel(submission_excel, csv_filename)
        index_tmp_file.append(csv_filename)
    return index_tmp_file


def convert_excel(input_xl, csv_file):
    df = pd.read_excel(input_xl, 0, header=5, index_col=None)
    df = df.dropna(axis=1, how='all')
    header = list(df.columns.values)
    df_sample = df[header[0]] + "," + df[header[1]] + "," + df[header[5]] + "-" +  df[header[6]]
    df_sample.to_csv(csv_file, sep='\t', header=None, index=False)
    

def map_sample_index(index_tmp_files):
    index_dict = {}
    for tmp_file in index_tmp_files:
        with open(tmp_file) as fin:
            raw_lines = fin.readlines()
            for raw_line in raw_lines:
                line = raw_line.strip()
                value = line.split(',')
                sample = "s_" + value[0]                                        
                sample_type = value[1]
                index = value[2]
                index_dict[sample] = [sample_type, index]
    return index_dict


def mad(data):
    reads_mapped_int = map(int, data)
    reads_mapped_array = np.asarray(reads_mapped_int)
    arr = np.ma.array(reads_mapped_array).compressed()
    med = np.median(reads_mapped_array)
    per_med = np.abs((reads_mapped_array - med))
    return np.median(np.abs(reads_mapped_array - med))


def get_score(value):
    show_score = mad(value)
    return show_score


def join_files(sample_reads_mapped_dict, index_dict):
    with open( 'coverage_outfile.tsv', 'wb') as fout:
        fout.write("sample\tTotal_Reads\tTotal_Reads_Mapped\tPercent_Total_ReadsMapped\tRealigned_reads\tpercent_realigned_reads\treads_mapped_to_target\tpercent_reads_mapped_to_target\tsample_type\tindex\n")
        for key, value in sample_reads_mapped_dict.items():
            if index_dict[key]:
#                print value, index_dict[key]
                out =  value + index_dict[key]
#                print out
#                out = "\t".join(str(i) for i in value) + "\t" + str(index_dict[key])
                fout.write(key + "\t" + "\t".join(out) + "\n")
                

def create_dataframe():
    outfile =  os.path.abspath('coverage_outfile.tsv')
    df = pd.DataFrame.from_csv(outfile, sep='\t',index_col=None)
    df['reads_mapped_to_target'] =  df['reads_mapped_to_target'].str.replace(",","")
    df['reads_mapped_to_target'] = df['reads_mapped_to_target'].astype(int)
    df['Percent_Total_ReadsMapped'] = df['Percent_Total_ReadsMapped'].str.strip("()")
    df['percent_realigned_reads'] = df['percent_realigned_reads'].str.strip("()")
    df['percent_reads_mapped_to_target'] = df['percent_reads_mapped_to_target'].str.strip("()")
    df['absolute_median'] = abs( df['reads_mapped_to_target'] -  df['reads_mapped_to_target'].median())
    df['standard_deviation'] = abs
    df.to_csv('carriers_coverage_absolute_median_results.tsv', index=False, sep='\t', encoding='utf-8')
    print df['reads_mapped_to_target'].median()


if __name__ == "__main__":
    main()
