#!/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bcbio/anaconda/bin/python

"""
This program is used to cleave the gene
specific primer, Common sequencing element
and N12 from CARRIERS qiagen fastq
"""

from Bio import SeqIO
import argparse
import os
from collections import Counter
from Bio.SeqRecord import SeqRecord
import regex as re
from collections import defaultdict
import multiprocessing
import logging
import sys
from Bio.Seq import Seq
import pandas as pd
import csv
import time
import gzip
from Bio import pairwise2
from Bio.SeqIO.QualityIO import FastqGeneralIterator


CSE='ATTGGAGTCCT'
REV_CSE='AGGACTCCAAT'


def main():
    logger = configure_logger()
    args = parse_args()
    run(args.read1_fastq, args.read2_fastq, args.primer_excel)


def run(read1_fastq, read2_fastq, primer_excel):
    forward_primer_list = convert_to_csv(primer_excel)
    reverse_primer_list = rev_compliment(forward_primer_list)
    read1_fastq = search_fastq(read1_fastq, read2_fastq)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-r1', dest='read1_fastq',
                        help="path to read1 fastq", required=True)
    parser.add_argument('-r2', dest='read2_fastq',
                        help="path to read2 fastq", required=True)
    parser.add_argument('-p', dest='primer_excel', default='/data5/bsi/epibreast/m087494.couch/Couch/Huge_Breast_VCF/Mayo_panel_primer.xlsx',
                       help='path to the mayo panel primer excel file')
    args = parser.parse_args()
    return args


def configure_logger():
    """
    setting up logging
    """
    logger = logging.getLogger('CARRIERS_primer_cleave')
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(time.strftime("CARRIERS_primer_cleave-%Y%m%d.log"))
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s'\t'%(name)s'\t'%(levelname)s'\t'%(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger


class fastq(object):
    
    def __init__(self, read1_fastq, read2_fastq):
        self.read1_fastq = read1_fastq
        self.read2_fastq = read2_fastq
        

    def __iter__(self):
        return self

    
    def fastq_record(self, fastq):
        self.fastq = fastq
        with gzip.open(self.fastq, 'rU') as fa:
            for record in SeqIO.parse(fa, "fastq"):
                yield record


    def fastq_record_id(self, records):
        for record in records:
            yield record.id            


    def fastq_record_seq(self, records):
        for record in records:
            yield record.seq


        
def convert_to_csv(primer_excel):
    """
    using the original source primer excel file
    """
    df = pd.read_excel(primer_excel)
    forward_primer = list(df.primer)
    return forward_primer


def rev_compliment(forward_sequence_list):
    rev_compliment_list = []
    for seq in forward_sequence_list:
        rev_compliment_list.append(Seq(seq).reverse_complement())
    return rev_compliment_list


def find_sequence_forward(sequence):
    if re.findall(CSE+'{e<=1}',sequence):
        index = re.finditer(CSE+'{e<=1}', sequence)
        for m in index:
            return m.start(0), m.end(0)


def find_sequence_reverse(sequence):
    if re.findall(REV_CSE+'{e<=1}',sequence):
        index = re.finditer(REV_CSE+'{e<=1}', sequence)
        for m in index:
            return m.start(0), m.end(0)


def search_fastq(read1_fastq, read2_fastq):
    check_fastq = fastq(read1_fastq, read2_fastq)
    R1_fastq_record = check_fastq.fastq_record(read1_fastq)
    for record in R1_fastq_record:
        read1_fastq = str(record.seq)
        search_seq = find_sequence_forward(read1_fastq)
        if search_seq:
            print CSE, record.seq, search_seq
        else:
            rev_search = find_sequence_reverse(read1_fastq)
            if rev_search:
                print REV_CSE, record.seq, rev_search


if __name__ == '__main__':
    main()
