#!/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bcbio/anaconda/bin/python

"""
This script filters and removes unnecessary columns and also
annotates the targets with gene
"""


import csv
import sys
import os
import pandas as pd
import pprint
import argparse

CARRIERS_annotated_bed="/data5/bsi/epibreast/m087494.couch/Couch/Huge_Breast_VCF/CARRIERS_PANC.targets.annotated.bed"

def main():
    args = parse_args()
    run(args.CARRIERS_annotated_bed, args.CARRIERS_interval_stats_file, args.output_dir, args.output_name)


def parse_args():
    """
    parsing arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-c', dest='CARRIERS_interval_stats_file',
                        help="interval stat file from depth of coverage", required=True)
    parser.add_argument('-o', dest='output_dir',
                        help='output dir to write output to',
                        required=True)
    parser.add_argument('-n', dest='output_name',
                       help='output file name')
    args = parser.parse_args()
    return args


def run(CARRIERS_annotated_bed, CARRIERS_interval_stats_file, output_dir, output_name):
    outpath = output_dir + '/' + output_name
    annot_bed_dict = bed_dict(CARRIERS_annotated_bed)
    stop_interval_dict = parse_interval_stats(CARRIERS_interval_stats_file, annot_bed_dict)
    remove_granular = remove_granular_columns(outpath)


def bed_dict(CARRIERS_annotated_bed):
    some_bed_dict ={}
    with open(CARRIERS_annotated_bed) as f:
        for raw_line in f:
            line = raw_line.strip().split()
            stop = line[2]
            start = line[1]
            chrm = line[0]
            annot = line[3]
            some_bed_dict[stop] = annot
    return some_bed_dict


def parse_interval_stats(CARRIERS_interval_stats_file, some_bed_dict):
    with open(CARRIERS_interval_stats_file) as csvfile, open('temp.tsv','wb') as out:
        reader = csv.DictReader(csvfile, delimiter='\t',skipinitialspace=True)
        headers = reader.fieldnames
        new_header = [headers[0], 'primers'] + headers[1:]
        writer = csv.DictWriter(out, new_header, extrasaction='ignore', delimiter='\t')
        writer.writeheader()
        for row in reader:
            stop = row['Target'].split("-")[1]
            primer = some_bed_dict[stop]
            some_out = [row[headers[0]], primer] + [row[i] for i in headers[1:]] 
            out.write("\t".join(i for i in some_out) + "\n")


def remove_granular_columns(OUTFILE):
    df = pd.read_csv('temp.tsv', sep='\t', low_memory=False)
    headers = list(df.columns.values)
    keep_header = [ i for i in headers if 'granular' not in i]
#    print keep_header
    df2 = df[keep_header]
    df2.to_csv(OUTFILE, sep='\t', index=False)
    os.remove('temp.tsv')

if __name__ == "__main__":
    main()
