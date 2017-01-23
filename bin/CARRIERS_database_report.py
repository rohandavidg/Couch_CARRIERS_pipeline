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
import shlex


#Databse_header_dict = ['NGS portal name (bio_ngs_portal_name)': ', 'Sample BAM name (bio_sample_bam_name)', 'NGS Lane Number (bio_ngs_lane)', 'NGS Flowcell (bio_ngs_flowcell)', 'Workflow Version (bio_workflow_vs)', 'BAM Local (bio_bam_local)', 'BAM Remote (FTP) (bio_bam_ftp)', 'VCF Local (bio_vcf_local)', 'VCF Remote (bio_vcf_remote)', 'gVCF Local (bio_gvcf_local)', 'gVCF Remote (bio_gvcf_remote)', 'Quality control pass/fail (bio_qc)', 'Quality control failed reason (bio_qc_reason)', 'Bioinformatics comments (bio_comments)', 'Bioinformatics complete (bio_complete)']


def main():
    args = parse_args()
    run(args.run_info, args.sample_info, args.database_csv_template, args.mad_table)


def parse_args():
    """
    parsing arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-r', dest='run_info', type=argparse.FileType('r'),
                        help="GGPS run info file",
                        required=True)
    parser.add_argument('-s', dest='sample_info', type=argparse.FileType('r'),
                        help="GGPS run info file",
                        required=True)
    parser.add_argument('-d',dest='database_csv_template',
                        help='database csv template',
                        required=True)
    parser.add_argument('-m', dest='mad_table', 
                        help='MAD table from align out files in GGPS', 
                        required=True)
    args = parser.parse_args()
    return args


def run(run_info, sample_info, database_csv_template, mad_table):
    delivery_folder, NGS_portal = get_run_info_values(run_info)
    get_bam_gvcf_dict = parse_sample_info(sample_info, delivery_folder)
    CARRIERS_sample_mad_dict = parse_mad_table(mad_table)
#    print CARRIERS_sample_mad_dict
    database_info = parse_database_template(database_csv_template, get_bam_gvcf_dict, CARRIERS_sample_mad_dict)


def get_run_info_values(run_info):
    NGS_portal = []
    delivery_folder = []
    for raw_line in run_info:
        try:
            (key, val) = raw_line.split("=")
            if key == 'DELIVERY_FOLDER':
                NGS_portal.append(shlex.split(val))
            if key == 'PROJECTNAME':
                delivery_folder.append(shlex.split(val))
        except ValueError:
            pass
    return ''.join(NGS_portal[0]), ''.join(delivery_folder[0])


def parse_sample_info(sample_info, run_delivery_folder):
    sample_bam_gvcf_dict = {}
    for raw_line in sample_info:
        some_line = raw_line.strip().split('\t')[0].split(':')[1].split("=")[0]
        sample = some_line.split("_")[1]
        sample_bam_gvcf_dict[sample] = [run_delivery_folder + '/' + 'bam' + "/" + some_line +".bam", 
                                        run_delivery_folder + "/" + 'variants/gVCF' + some_line + '.g.vcf.gz'] 
    return sample_bam_gvcf_dict


def parse_mad_table(mad_table):
    mad_dict = {}
    with open (mad_table) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            mad_dict[row['sample']] = [row['Indexes'], row['TotalReads'], row['lanes'], row['ReadsInCaptureRegion'], row['absolute_median']]
    return mad_dict


def parse_database_template(database_csv_template, sample_bam_gvcf_dict, mad_dict):
    with open (database_csv_template) as csvfile:
        headers = [header.strip() for header in next(csvfile).split(",")]
#        print headers
        reader = csv.DictReader(csvfile, fieldnames=headers)
        for row in reader:
            try:
                mad_dict[row['CARRIERS ID (carriers_id)']]
                print row
            except KeyError:
                pass

if __name__== "__main__":
    main()

