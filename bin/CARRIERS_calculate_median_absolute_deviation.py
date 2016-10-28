#!/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bcbio/anaconda/bin/python


import sys
import numpy as np
from numpy import mean, absolute
import pandas as pd
import os


def main():
    args = parse.args()
    run(args.input_file)


def run(input_file):
    reads_mapped_dict, sample_reads_mapped_dict, reads_mapped_list = parse_file(input_file)
#    sample_info_dict = parse_info(info_file)
    compute_mad = get_score(reads_mapped_list)
#    out_file = join_files(sample_reads_mapped_dict, sample_info_dict)
    compute_mad = create_dataframe()


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', dest='test_name',
                        help="ordered service name",
                        required=True)
    args = parser.parse_args()
    return args


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
            ReadsInCaptureRegion = line[5].split(" ")[0]
            ReadsInCaptureRegion_per = line[5].split(" ")[1]
            sample_reads_mapped_dict[sample] = [total_reads, TotalMappedReads, total_reads_per,RealignedReads,
                                                RealignedReads_per, ReadsInCaptureRegion, ReadsInCaptureRegion_per]
            reads_mapped_dict[sample] = reads_maped
    return reads_mapped_dict, sample_reads_mapped_dict, reads_mapped_list


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

def mad(data):
    reads_mapped_int = map(int, data)
    reads_mapped_array = np.asarray(reads_mapped_int)
    arr = np.ma.array(reads_mapped_array).compressed()
    med = np.median(reads_mapped_array)
    per_med = np.abs((reads_mapped_array - med))
    return np.median(np.abs(reads_mapped_array - med))


def get_score(value):
    show_score = mad(value)
    print show_score


def join_files(sample_reads_mapped_dict, sample_info_dict):
    with open( 'outfile.tsv', 'wb') as fout:
        fout.write("sample\tTotal_Reads\tTotal_Reads_Mapped\tPercent_Total_ReadsMapped\tRealigned_reads\tpercent_realigned_reads\treads_mapped_to_target\tpercent_reads_mapped_to_target\tvolume\tindex\tI5_Index_ID\tsource\tSample_ID\tDNA_type\tmethod\n")
        for key, value in sample_reads_mapped_dict.items():
            if sample_info_dict[key]:
                out = value + sample_info_dict[key]
                fout.write(key + "\t" + "\t".join(out) + "\n")
                

def create_dataframe():
    outfile =  os.path.abspath('outfile.tsv')
    df = pd.DataFrame.from_csv(outfile, sep='\t',index_col=None)
    df['reads_mapped_to_target'] =  df['reads_mapped_to_target'].str.replace(",","")
    df['reads_mapped_to_target'] = df['reads_mapped_to_target'].astype(int)
    df['Percent_Total_ReadsMapped'] = df['Percent_Total_ReadsMapped'].str.strip("()")
    df['percent_realigned_reads'] = df['percent_realigned_reads'].str.strip("()")
    df['percent_reads_mapped_to_target'] = df['percent_reads_mapped_to_target'].str.strip("()")
    df['absoute_median'] = abs( df['reads_mapped_to_target'] -  df['reads_mapped_to_target'].median())
    df.to_csv('carriers_coverage_absolute_median_results.tsv', index=False, sep='\t', encoding='utf-8')
    print df['reads_mapped_to_target'].median()


if __name__ == "__main__":
    main()
