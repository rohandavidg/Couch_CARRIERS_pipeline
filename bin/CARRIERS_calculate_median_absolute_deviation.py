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
import gzip
from collections import defaultdict
import sys
import logging
import time
import datetime
import pprint

def main():
    args = parse_args()
    run(args.submission_excel_list, args.stats_file, args.ggps_coverage_file)


def run(submission_excel_list, stats_file, ggps_coverage_file):
    logger = configure_logger()
    reads_mapped_dict, sample_reads_mapped_dict, reads_mapped_list = parse_file(ggps_coverage_file)
    if stats_file:
        coverage_metrics_list = parse_stats_file(stats_file, logger)
        stats_dict = index_coverage_stats(coverage_metrics_list, logger)
        sample_index_dict = map_sample_index(create_submssion_csv_list, logger)
        sample_stats_index_dict = find_samples_from_index(stats_dict, sample_index_dict, logger)
    else:
        create_submssion_csv_list = create_csv_from_excel(submission_excel_list, logger)
        sample_index_dict = map_sample_index(create_submssion_csv_list, logger)
        compute_mad = get_score(reads_mapped_list)
        out_file = join_files(sample_reads_mapped_dict, sample_index_dict)
        compute_mad = create_dataframe()


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', dest='ggps_coverage_file',
                        help="coverage table in the index.html in tsv format",
                        required=True)
    parser.add_argument('-l', dest='submission_excel_list', nargs="+",
                        help="submission excel manifest submitted to BAP",
                        required=True)
    parser.add_argument('-s', dest='stats_file',
                        help='path to ggps stats file stats.tar.gz')
    args = parser.parse_args()
    return args


def configure_logger():
    """
    setting up logging
    """
    logger = logging.getLogger('CARRIERS_coverage_MAD')
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(time.strftime("CARRIERS_coverage_MAD-%Y%m%d.log"))
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s'\t'%(name)s'\t'%(levelname)s'\t'%(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger


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
            try:
                ReadsInCaptureRegion = line[5].split(" ")[0]
                ReadsInCaptureRegion_per = line[5].split(" ")[1]
                sample_reads_mapped_dict[sample] = [total_reads, TotalMappedReads, total_reads_per,RealignedReads,
                                                    RealignedReads_per, ReadsInCaptureRegion, ReadsInCaptureRegion_per]
                reads_mapped_dict[sample] = reads_maped            
            except IndexError:
                ReadsInCaptureRegion = "-"
                ReadsInCaptureRegion_per = "-"
                sample_reads_mapped_dict[sample] = [total_reads, TotalMappedReads, total_reads_per,RealignedReads,
                                                    RealignedReads_per, ReadsInCaptureRegion, ReadsInCaptureRegion_per]
                reads_mapped_dict[sample] = reads_maped
    return reads_mapped_dict, sample_reads_mapped_dict, reads_mapped_list


def parse_stats_file(stats_file, logger):
    coverage_values = []
    with gzip.open(stats_file, 'r') as fin:
        logger.info('parsing out only read metrics for stats.gz and breaking before GATKTable')
        for raw_line in fin:
            line = raw_line.strip()
            if line.startswith('#:GATKTable:'):
                break
            else:
                coverage_values.append(line)
    return coverage_values
    

def index_coverage_stats(coverage_values, logger):
    sample_index = []
    sample_lanes = []
    sample_TotalReads = []
    sample_TotalMappedReads = []
    sample_DuplicateReads = []
    sample_RealignedReads = []
    sample_ReadsInCaptureRegion = []
    sample_stat_dict = defaultdict(list)
    for i, x in enumerate(coverage_values[:-1]):
        if i % 14 == 3:
            index = set([a for a in x.split(",")])
            sample_index.append(index)
            continue
        elif i % 14 == 1:
            lanes = x
            sample_lanes.append(lanes)
            continue
        elif i % 14 == 5:
            TotalReads = int(x)
            if isinstance(TotalReads, int):
                sample_TotalReads.append(TotalReads)
            else:
                logger.debug('Total Reads:{0} is not a valid value, on line {1}'.format(TotalReads,i))
                sys.exit(1)
        elif i % 14 == 7:
            TotalMappedReads = int(x)
            if isinstance(TotalMappedReads, int):
                sample_TotalMappedReads.append(TotalMappedReads)
            else:
                logger.debug('Total Mapped Reads:{0} is not a valid value, on line {1}'.format(TotalMappedReads, i))
                sys.exit(1)
        elif i % 14 == 9:
            DuplicateReads = x
            sample_DuplicateReads.append(DuplicateReads)
            continue
        elif i % 14 == 11:
            RealignedReads = int(x)
            if isinstance(RealignedReads, int):
                sample_RealignedReads.append(RealignedReads)
            else:
                logger.debug('RealignedReads:{0} is not a valid value, on line {1}'.format(str(RealignedReads), i))
                sys.exit(1)
        elif i % 14 == 13:
            ReadsInCaptureRegion = x
            sample_ReadsInCaptureRegion.append(ReadsInCaptureRegion)
            continue
        else:
            ReadsInCaptureRegion = x
            logger.debug('ReadsInCaptureRegion: {0} is not a valid value, on line {1}'.format(ReadsInCaptureRegion, i))
                ##TODO:Not sure how to handle this error- re-run count reads or pass it
    pprint.pprint(sample_index)
#    sample_stat_dict = dict(("".join(z[0]), list(z[1:])) for z in zip(sample_index, sample_lanes, sample_TotalReads,
#                                                                      sample_TotalMappedReads, sample_DuplicateReads,
#                                                                      sample_RealignedReads, sample_ReadsInCaptureRegion))
#    sample_stat_dict[sample_index].append([sample_lanes, sample_TotalReads, sample_TotalMappedReads, sample_DuplicateReads,
#                                           sample_RealignedReads, sample_ReadsInCaptureRegion])

#    print sample_index
#    return sample_stat_dict


#/usr/local/biotools/java/jdk1.7.0_67/bin/java -Xmx24g -jar /data5/bsi/bictools/alignment/gatk/3.5/GenomeAnalysisTK.jar -T CountReads -L /data5/bsi/epibreast/m087494.couch/Couch/Huge_Breast_VCF/CARRIERS_PANC.targets.bed -I /data2/delivery/Couch_Fergus_coucf/161024_K00316_0040_AHFV37BBXX/secondary/bam/s_CA00005790.bam -R /data2/bsi/reference/sequence/human/ncbi/37.1/indexed/allchr.fa  --filter_bases_not_stored 2>&1 | grep CountReads | tail -1 | cut -d '-' -f2 | tr ' ' '\n' | awk '/^[0-9]/' > s_CA00005790.captureReads.txt


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


def create_csv_from_excel(submission_excel_list, logger):
    index_tmp_file = []
    for i, submission_excel in enumerate(submission_excel_list):
        logger.info('creating a csv from {0} excel file'.format(submission_excel))
        csv_filename = "carriers_pancreas" + str(i + 1) +"_tmp.csv"
        to_csv = convert_excel(submission_excel, csv_filename)
        logger.info('succesfully created csv from chunlings excel file containing sample and index information submitted for sequencing')
        index_tmp_file.append(csv_filename)
    return index_tmp_file


def convert_excel(input_xl, csv_file):
    df = pd.read_excel(input_xl, 0, header=5, index_col=None)
    df = df.dropna(axis=1, how='all')
    header = list(df.columns.values)
    df_sample = df[header[0]] + "," + df[header[1]] + "," + df[header[5]] + "-" +  df[header[6]]
    df_sample.to_csv(csv_file, sep='\t', header=None, index=False)
    

def map_sample_index(index_tmp_files, logger):
    index_dict = {}
    for tmp_file in index_tmp_files:
        logger.info('parsing the csv file {0}'.format(tmp_file))
        with open(tmp_file) as fin:
            raw_lines = fin.readlines()
            for raw_line in raw_lines:
                line = raw_line.strip()
                value = line.split(',')
                sample = "s_" + value[0]
                sample_type = value[1]
                index = value[2]
                new_index = index.split('-')
                if len(new_index) < 2:
                    logger.debug('index {0} in missing either a i5 or i7 in the excel sheet'.format(index))
                    print "check log file"
                    sys.exit(1)
                else:
                    pass
                index_dict[sample] = [sample_type, index]
    return index_dict


def find_samples_from_index(sample_stat_dict, index_dict, logger):
    logger.info('writing out temp file after mapping index to sample name')
    print "writing out temp file this may contain '-' value in reads in capture region"
    with open( 'coverage_outfile.tsv', 'wb') as fout:
        fout.write("sample\tTotal_Reads\tTotal_Reads_Mapped\tPercent_Total_ReadsMapped\tRealigned_reads\tpercent_realigned_reads\treads_mapped_to_target\tpercent_reads_mapped_to_target\tsample_type\tindex\n")
        for key, value in index_dict.items():
            if sample_stat_dict[value[1]]:
                sample = key 
                index = value[1]
                sample_type = value[0]
                total_reads = sample_stat_dict[value[1]][1]
                total_TotalMappedReads = sample_stat_dict[value[1]][2]
                percent_total_mapped = (float(total_TotalMappedReads)/float(total_reads)) * 100
                total_DuplicateReads = sample_stat_dict[value[1]][3]
                total_RealignedReads = sample_stat_dict[value[1]][4]
                percent_total_RealignedReads = (float(total_RealignedReads)/float(total_reads)) * 100
                total_ReadsInCaptureRegion =  sample_stat_dict[value[1]][5]
                try:
                    percent_total_ReadsInCaptureRegion = (float(total_ReadsInCaptureRegion)/float(total_reads)) * 100
                    out = [sample, total_reads, total_TotalMappedReads, percent_total_mapped, 
                           total_RealignedReads, percent_total_RealignedReads, total_ReadsInCaptureRegion,
                           percent_total_ReadsInCaptureRegion, sample_type, index]
                    fout.write("\t".join(str(i) for i in out) + "\n")
                except ValueError:
                    logger.debug('total_ReadsInCaptureRegion:{0} has a non-integer value'.format(total_ReadsInCaptureRegion))
                    out = [sample, total_reads, total_TotalMappedReads, percent_total_mapped, 
                           total_RealignedReads, percent_total_RealignedReads, total_ReadsInCaptureRegion,
                           "-", sample_type, index]
                    fout.write("\t".join(str(i) for i in out) + "\n")


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
#    print df
    df['reads_mapped_to_target'] =  df['reads_mapped_to_target'].str.replace(",","")
#    print df['reads_mapped_to_target'].applymap(np.isreal)
    df['reads_mapped_to_target'] = df['reads_mapped_to_target'].astype(int)
    df['Percent_Total_ReadsMapped'] = df['Percent_Total_ReadsMapped'].str.strip("()")
    df['percent_realigned_reads'] = df['percent_realigned_reads'].str.strip("()")
    df['percent_reads_mapped_to_target'] = df['percent_reads_mapped_to_target'].str.strip("()")
    df['absolute_median'] = abs( df['reads_mapped_to_target'] -  df['reads_mapped_to_target'].median())
#    df['standard_deviation'] = abs
    df.to_csv('carriers_coverage_absolute_median_results.tsv', index=False, sep='\t', encoding='utf-8')
#    print df['reads_mapped_to_target'].median()


if __name__ == "__main__":
    main()
