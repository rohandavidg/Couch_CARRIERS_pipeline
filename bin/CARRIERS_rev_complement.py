#!/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bcbio/anaconda/bin/python


"""
script to do a rev complement of a seq
"""

from Bio import SeqIO
import argparse
from Bio.SeqRecord import SeqRecord
import regex as re
from collections import defaultdict
from Bio.Seq import Seq


CSE = 'ATTGGAGTCCT'
REV_CSE = 'AGGACTCCAAT'

def main():
    args = parse_args()
    run(args.mayo_primer_panel)


def run(mayo_primer_panel):
    rev_complement = do_complement(mayo_primer_panel)
    print ">CSE"
    print CSE


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-l', dest='mayo_primer_panel',
                        help="Mayo primer file in tsv format", required=True)
    args = parser.parse_args()
    return args


def do_complement(primer_file):
    chrom_seq_dict = {}
    chrom_rev_seq_dict = {}
    with open(primer_file) as fin:
        for raw_line in fin:
            line = raw_line.strip()
            value = line.split('\t')
            key = ">" + value[0] + ":" + value[1] + "-" + value[2]
            print key
            sequence = Seq(value[4])
            reverse_complement = sequence.reverse_complement()
            print sequence
#            print reverse_complement


if __name__ == "__main__":
    main()
