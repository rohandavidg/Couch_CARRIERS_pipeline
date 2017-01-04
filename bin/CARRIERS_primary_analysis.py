#!/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bcbio/anaconda/bin/python

"""
script to generate read metrics from primary
bams
"""

import argparse
import csv
import pprint
import subprocess

samtools = "/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bin/samtools"

def main():
    args = parse_args()
    run(args.bam_list)


def run(bam_list):
    show_flagstat = flag_stat(bam_list) 


def parse_args():
    """
    parsing arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-l', dest='list_bam', type=argparse.FileType('r'),
                        help="list of bams in the run",
                        required=True)
    args = parser.parse_args()
    return args


def flag_stat(bam_list):
    for bam_file in bam_list:
        cmd = samtools + " flagstat " + bam_file + " head -8 " + '|' + " tail -3 "
        print cmd
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        flagstat, err2 = p.communicate()
        print flagstat


if __name__ == '__main__':
    main()
