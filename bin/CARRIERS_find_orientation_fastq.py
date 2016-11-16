#!/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bcbio/anaconda/bin/python

"""
script to split a fastq into specific strand orientation
"""

import os
from Bio import SeqIO
import argparse
import gzip
import multiprocessing
import Levenshtein
from Levenshtein import distance
from Levenshtein import hamming
from Levenshtein import jaro
#import re
from collections import Counter
from Bio.SeqRecord import SeqRecord
import regex as re 
from collections import defaultdict

CSE = 'ATTGGAGTCCT'
REV_CSE = 'AGGACTCCAAT'

def main():
    args = parse_args()
    run(args.read1_fastq, args.read2_fastq)


def run(read1_fastq, read2_fastq):
    fastq_check = handle_fastq(read1_fastq, read2_fastq, CSE)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-r1', dest='read1_fastq',
                        help="path to read1 fastq", required=True)
    parser.add_argument('-r2', dest='read2_fastq',
                        help="path to read2 fastq", required=True)
    args = parser.parse_args()
    return args


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
    

def handle_fastq(read1_fastq, read2_fastq, CSE):
    check_fastq = fastq(read1_fastq, read2_fastq)
    read_2_records = check_fastq.fastq_record(read2_fastq)
    read_1_records = check_fastq.fastq_record(read1_fastq)
    read_2_forward_index_count_dict, read_2_reverse_index_count_dict,\
total_record, unaccounted_read = check_fastq_record(read_2_records, CSE, REV_CSE)
    read_1_forward_index_count_dict, read_1_reverse_index_count_dict,\
total_record, unaccounted_read = check_fastq_record(read_1_records, CSE, REV_CSE)
#    read2_stats = show_stats(read_2_forward_index_count_dict, read_2_reverse_index_count_dict, total_record, unaccounted_read, "read2")
#    read1_stats = show_stats(read_1_forward_index_count_dict, read_1_reverse_index_count_dict, total_record, unaccounted_read, "read1")


def check_fastq_record(some_record, CSE, REV_CSE):
    record_id = set()
    total_record = 0
    unaccounted_read = 0
    forward_index_count_dict = defaultdict(int)
    reverse_index_count_dict = defaultdict(int)
    for record in some_record:            
        total_record += 1 
        rev_seq = record.seq.reverse_complement()
        check_forward = find_sequence_forward(str(record.seq))
        check_reverse = find_sequence_reverse(str(record.seq))
        if check_forward:
            window = [str(record.seq)[i:i+11] for i in range(len(str(record.seq))-5)]
            forward_score = [compute_edit_distance_forward(CSE, i) for i in window]
            check_score_for = [True if i == 0 or i == 1 or i == 2 else False for i in forward_score]
            forward_index = check_score_for.index(True)
            forward_index_count_dict[forward_index] += 1
        elif check_reverse:
            window = [str(record.seq)[i:i+11] for i in range(len(str(record.seq))-5)]
            reverse_score = [compute_edit_distance_reverse(REV_CSE, i) for i in window]
            check_score_rev = [True if i == 0 or i == 1 or i == 2 else False for i in reverse_score]
            try:
                reverse_index = check_score_rev.index(True)
                reverse_index_count_dict[reverse_index] += 1
            except ValueError:
                print reverse_score
        else:
            unaccounted_read += 1
    print reverse_index_count_dict
    return forward_index_count_dict, reverse_index_count_dict, total_record, unaccounted_read
    

def show_stats(forward_index_count_dict, reverse_index_count_dict,
               total_record, forward_count, reverse_count, unaccounted_read, side):
    print ("Total number of reads in {0} fastq: {1}").format(side, total_record)
    print ("Total number of reads on {0} forward strand: {1}").format(side, forward_count)
    print ("Total number of reads on {0} reverse strand: {1}").format(side, reverse_count)
    print ("Total number of reads on {0} unaccounted for: {1}").format(side, unaccounted_read)
    if forward_index_count_dict:
        for k, v in forward_index_count_dict.items():
            print ("Total count of CSE starting at {0} forward strand index {1} is {2}").format(side, k,v)
    else:
        pass
    if reverse_index_count_dict:
        for k, v in reverse_index_count_dict.items():
            print ("Total count of CSE starting at {0} reverse strand index {1} is {2}").format(side, k,v)
    else:
        pass        
                 

def write_out_forward(file_name, some_record, some_list):
    with open(file_name, 'wa+') as handle:
        for record in some_record:
            if record.id in some_list:
                SeqIO.write(record, handle, "fastq")
            else:
                pass


def compute_edit_distance_forward(CSE, sequence):
    z = distance(CSE, sequence)
    return z


def compute_edit_distance_reverse(REV_CSE, sequence):
    q = distance(REV_CSE, sequence)
    return q


def find_sequence_forward(sequence):
    if re.findall("ATTGGAGTCCT{e<=1}",sequence):
        return True


def find_sequence_reverse(sequence):
    if re.findall("AGGACTCCAAT{e<=1}",sequence):
        return True



if __name__ == "__main__":
    main()
