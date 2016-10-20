#!/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bcbio/anaconda/bin/python

# cython: profile=True

#%load_ext cythonmagic

"""
This script will create a csv file with the results
of the carriers run to upload into
redcap database
"""

import cython
import csv
import os
import pprint
import argparse
import logging
import datetime
import time


reguired_fields = ['NGS portal name (bio_ngs_portal_name)', 'Sample BAM name (bio_sample_bam_name)', 'NGS Lane Number (bio_ngs_lane)', 'NGS Flowcell (bio_ngs_flowcell)',
'Workflow Version (bio_workflow_vs)', 'BAM Local (bio_bam_local)', 'BAM Remote (FTP) (bio_bam_ftp)', 'VCF Local (bio_vcf_local)', 'VCF Remote (bio_vcf_remote)', 'gVCF Local (bio_gvcf_local)', 'gVCF Remote (bio_gvcf_remote)', 'Quality control pass/fail (bio_qc)', 'Quality control failed reason (bio_qc_reason)', 'Bioinformatics comments (bio_comments)', 'Bioinformatics complete (bio_complete)']


def main():
    args = parse_args()
    run(args.red_cap_cnv, args.run_info)


def run(red_cap_cnv, run_info):
    carriers_database_id_dict = parse_cnv(red_cap_cnv)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-c', dest='red_cap_cnv',
                        help="red cap cnv for bioinformatics that has information on the database id",
                        required=True)
    parser.add_argument('-r', dest='run_info',
                        help="run info txt file used in the GGPS run",
                        required=True)
    args = parser.parse_args()
    return args


def parse_cnv(red_cap_cnv):
    with open(red_cap_cnv, "rb") as csvfile:
        header = csv.reader(csvfile).next()
        new_header = [i.replace(" ","_") for i in header]
        reader = csv.DictReader(csvfile, fieldnames=new_header)
        for row in reader:
            print row


if __name__ == "__main__":
    main()
