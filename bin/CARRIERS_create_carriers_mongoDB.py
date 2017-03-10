#!/Users/m149947/anaconda/bin//python

"""
this scripts set up a carries 
mongo db
"""

import sys
import pymongo
#from pymongo import Connection
from pymongo import MongoClient
from bson.json_util import dumps
from bson.json_util import loads
import csv
import json
import os
import datetime
import pprint
import bson
import argparse
import pandas as pd
from collections import defaultdict
#from collections import OrderedDict
import collections
import pprint

def main():
    args = parse_args()
    run(args.study, args.primer_performance, args.mad_file, args.coverage_file)


def run(study, primer_performance, mad_file, coverage_file):
    sample_mad_dict = parse_mad_file(mad_file)
    reformat_primer = eval_primer_performance(primer_performance)
    reformat_coverage = eval_coverage_performace(coverage_file)
#    pprint.pprint(reformat_coverage)
    csv_to_json = parse_csv_json(sample_mad_dict, study)


def parse_args():
    """
    parsing arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-s', dest='study',
                        help='CARRIERS study name')
    parser.add_argument('-p', dest='primer_performance',
                        help='primer performance output file',
                        required=True)
    parser.add_argument('-c', dest='coverage_file',
                        help='coverage file for each target',
                        required=True)
    parser.add_argument('-m', dest='mad_file', type=argparse.FileType('r'),
                        help='coverage MAD output file')
    args = parser.parse_args()
    return args


def create_db():
    client = MongoClient('rcf-mongo-shr1-dev', 27017 )
    db = client['hart']
    db.authenticate("m087494", "mayo72841")    
    #client.the_database.autheticate('m087494', 'mayo72841')
    #ppcollection = db["couch_primer_performance"]
    return db["couch_carriers_bionformatics"]

def parse_mad_file(mad_file):
    mad_dict = {}
    reader = csv.DictReader(mad_file, delimiter='\t')
    for row in reader:
        mad_dict[row['sample']] = [row['lanes'], row['flowcell'],
                                   row['Indexes'], int(row['TotalReads']),  
                                   int(row['RealignedReads']), int(row['ReadsInCaptureRegion']), 
                                   float(row['Absolute_median_deviation']),
                                   row['sample']]

    return mad_dict


def eval_primer_performance(primer_file):
    CARRIERS_primer_range_dict = {}
    df = pd.read_csv(primer_file, sep='\t')
    df['gene'] = df['target_gene'].str.split('_').str[0]
    df['target'] = df['chrom'] +':' + df['loc5'].astype(str) + '-' + df['loc3'].astype(str) + '-' + df['gene']
    headers = df.columns.tolist()
    new_df = df.drop(df.columns[[0,1,2,3,4,5,6,7]], axis=1)
    new_headers = new_df.columns.tolist()
    rearrange_headers = new_headers[-1:] + new_headers[:-2]
    filtered_new_df = new_df[rearrange_headers]
    transporse_df = filtered_new_df.transpose()
    transporse_df.columns = transporse_df.iloc[0]
    transporse_df = transporse_df[1:]
    transporse_df = transporse_df.reset_index()
    transporse_df = transporse_df.rename(columns ={'index' : 'sample'})
    transporse_df['sample'] = transporse_df['sample'].str.split('_').str[-1]
    transporse_df.to_csv('temp_primer_file.csv', index=False)


def eval_coverage_performace(coverage_file):
    df = pd.read_csv(coverage_file, sep='\t',low_memory=False)
    headers = df.columns.tolist()
    new_header = [i for i in headers if 'granular' not in i]
    new_df = df[new_header]
    new_df['gene'] = df['primers'].str.split('_').str[0]
    new_df['Target'] = new_df['Target'] + '-' + new_df['gene']
    new_df = new_df.drop(new_df.columns[[1,2,3,-1]], axis=1)
    transpose_df = new_df.transpose()
    transpose_df.columns = transpose_df.iloc[0]
    transpose_df = transpose_df[1:]
    reset_transpose_df = transpose_df.reset_index()
    reset_transpose_df = reset_transpose_df.rename(columns ={'index' : 'sample'})
    reset_transpose_df['new_sample'] = reset_transpose_df['sample'].str.split('_').str[1]
    reset_transpose_df['percentage_above_x'] = reset_transpose_df['sample'].str.split('_').str[-2] + '_' + reset_transpose_df['sample'].str.split('_').str[-1]
    some_header = reset_transpose_df.columns.tolist()
    out_header = some_header[-2:] + some_header[1:-2]
    coverage_out_df = reset_transpose_df[out_header]
    coverage_out_df.to_csv('temp_coverage_file.csv', index=False) 


def generate_coverage_dict():
    new_data_dict = {}
    with open('temp_coverage_file.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        headers = reader.fieldnames
        chr_headers = headers[2:]
        for row in reader:
            item = new_data_dict.get(row['new_sample'], dict())
            item[row['percentage_above_x']] = [{i:float(row[i])} for i in chr_headers]
            new_data_dict[row['new_sample']] = item
    return new_data_dict


def parse_csv_json(mad_dict, study):
    coverage_performance_dict = generate_coverage_dict()
#    pprint.pprint(coverage_performance_dict)
#    connect_mongo = check_connection()
#    primer_performance = []
    ppcollection = create_db()
    pprint.pprint(ppcollection)
    with open('temp_primer_file.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        header = reader.fieldnames
        for row in reader:
            new_row = [{i:int(row[i])} for i in header if 'sample' not in i]
            try:
                if mad_dict[row['sample']] and coverage_performance_dict[row['sample']]:
                    value = {'primer_performance':new_row,
                             'coverage_performance':coverage_performance_dict[row['sample']],
                             'study': study, 'lane':mad_dict[row['sample']][0],
                             'flowcell':mad_dict[row['sample']][1],
                             'Indexes':mad_dict[row['sample']][2],
                             'TotalReads':float(mad_dict[row['sample']][3]),
                             'RealignedReads': float(mad_dict[row['sample']][4]),
                             'ReadsInCaptureRegion': float(mad_dict[row['sample']][5]),
                             'Absolute_median_deviation': float(mad_dict[row['sample']][6]),
                             'CARRIERS_ID':mad_dict[row['sample']][7]}                             
                    ppcollection.insert(loads(dumps(value)))
            except KeyError:
                pass




if __name__ == '__main__':
    main()
