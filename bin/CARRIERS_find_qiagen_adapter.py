#!/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bcbio/anaconda/bin/python

"""
script to find the adapters from the softclips of
the primary bam files
"""

import pysam
import os
import argparse
from collections import defaultdict

def main():
    args = parse_args()
    run(args.bam)


def run(bam):
    r1_rev_st_fd_clip, r1_rev_st_rv_clip, r1_for_st_fd_clip, r1_for_st_rev_clip,\
        r2_rev_st_for_clip, r2_rev_st_rev_clip, r2_for_st_for_clip, r2_for_st_rev_clip = parse_bam(bam)
    show_results = write_out_results(r1_rev_st_fd_clip, "read1_reverse_strand_softclips_at_the_begining_of_cigar") 
    show_results = write_out_results(r1_rev_st_rv_clip, "read1_reverse_strand_softclips_at_the_end_of_cigar")
    show_results = write_out_results(r1_for_st_fd_clip, "read1_forward_strand_softclips_at_the_begining_of_cigar") 
    show_results = write_out_results(r1_for_st_rev_clip, "read1_forward_strand_softclips_at_the_end_of_cigar")
    show_results = write_out_results(r2_rev_st_for_clip, "read2_reverse_strand_softclips_at_the_begining_of_cigar")
    show_results = write_out_results(r2_rev_st_rev_clip, "read2_reverse_strand_softclips_at_the_end_of_cigar")
    show_results = write_out_results(r2_for_st_for_clip, "read2_forward_strand_softclips_at_the_begining_of_cigar")
    show_results = write_out_results(r2_for_st_rev_clip, "read2_forward_strand_softclips_at_the_end_of_cigar")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-b', dest='bam',
                        help="path to bam", required=True)
    args = parser.parse_args()
    return args


def parse_bam(bam):
    samfile = pysam.AlignmentFile(bam, "rb")
    read1_reverse_strand_forward_clip = set()
    read1_reverse_strand_reverse_clip = set()
    read1_forward_strand_forward_clip = set()
    read1_forward_strand_reverse_clip = set()
    read2_reverse_strand_forward_clip = set()
    read2_reverse_strand_reverse_clip = set()
    read2_forward_strand_forward_clip = set()
    read2_forward_strand_reverse_clip = set()
    for read in samfile.fetch():
        if read.is_paired:
            if read.is_read1:
                if read.is_reverse:
                    read1_rev_sd_for_clip, read1_rev_sd_rev_clip = find_common_softclips(read.cigar, read.seq)
                    try:
                        read1_reverse_strand_forward_clip.add(read1_rev_sd_for_clip[0])
                        read1_reverse_strand_reverse_clip.add(read1_rev_sd_rev_clip[0])
                    except IndexError:
                        pass
                else:
                    read1_for_sd_for_clip, read1_for_sd_rev_clip = find_common_softclips(read.cigar, read.seq) 
                    try:
                        read1_forward_strand_forward_clip.add(read1_for_sd_for_clip[0])
                        read1_forward_strand_reverse_clip.add(read1_for_sd_rev_clip[0])
                    except IndexError:
                        pass
            else:
                if read.is_reverse:
                    read2_rev_sd_for_clip, read2_rev_sd_rev_clip = find_common_softclips(read.cigar, read.seq)
                    try:
                        read2_reverse_strand_forward_clip.add(read2_rev_sd_for_clip[0])
                        read2_reverse_strand_reverse_clip.add(read2_rev_sd_rev_clip[0])
                    except IndexError:
                        pass 
                else:
                    read2_for_sd_for_clip, read2_for_sd_rev_clip = find_common_softclips(read.cigar, read.seq)           
                    try:
                        read2_forward_strand_forward_clip.add(read2_for_sd_for_clip[0])
                        read2_forward_strand_reverse_clip.add(read2_for_sd_rev_clip[0])
                    except IndexError:
                        pass
    samfile.close() 
    return read1_reverse_strand_forward_clip, read1_reverse_strand_reverse_clip, read1_forward_strand_forward_clip, read1_forward_strand_reverse_clip, read2_reverse_strand_forward_clip, read2_reverse_strand_reverse_clip, read2_forward_strand_forward_clip, read2_forward_strand_reverse_clip


def  find_common_softclips(cigar, seq):
    forward_clip_read = []
    reverse_clip_read = []
    check_softclip = find_softclips(cigar)
    if check_softclip:
        for k, v in check_softclip.items():
            if k == 'front':
                front_clip = check_softclip['front'][0]
                forward_clip_read.append(seq[:front_clip])
            else:
                back_clip = check_softclip['back'][0]
                reverse_clip_read.append(seq[-back_clip:])                    
    return forward_clip_read, reverse_clip_read

                                
def find_softclips(some_tuple):
    tuple_size = len(some_tuple)
    clip_index_dict = defaultdict(list)
    for i, x in enumerate(some_tuple):
        if i == 0:
            if x[0] == 4:
                index = x[1]
                clip_index_dict['front'].append(index)
            else:
                pass
        if i == 1:
            if x[0] == 4:
                index = x[1]
                clip_index_dict['back'].append(index)
            else:
                pass
        if i == 2:
            if x[0] == 4:
                index = x[1]
                clip_index_dict['back'].append(index)
            else:
                pass
        if i == 3:
            if x[0] == 4:
                index = x[1]
                clip_index_dict['back'].append(index)                
            else:
                pass
        if i == 4:
            if x[0] == 4:
                index = x[1]
                clip_index_dict['back'].append(index)
            else:
                pass
    return clip_index_dict

def write_out_results(some_set, set_name):
    with open(set_name, 'wa') as fout:
        for i in some_set:
            fout.write(i + "\n")


if __name__ == "__main__":
    main()
