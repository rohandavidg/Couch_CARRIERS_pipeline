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

FTP_path = "ftp://research-archive.mayo.edu/delivery/Couch_Fergus_coucf/"

def main():
    args = parse_args()
    run(args.run_info_list, args.database_csv_template, 
        args.mad_table, args.primer_performance)


def parse_args():
    """
    parsing arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-r', dest='run_info_list', nargs='+',
                        help="GGPS run info file",
                        required=True)
    parser.add_argument('-d',dest='database_csv_template',
                        help='database csv template',
                        required=True)
    parser.add_argument('-m', dest='mad_table', 
                        help='MAD table from align out files in GGPS', 
                        required=True)
    parser.add_argument('-p', dest='primer_performance',
                        help='primer performance output file')
    parser.add_argument('-c', dest='target_coverage',
                        help='target coverage output file')
    args = parser.parse_args()
    return args


def run(run_info_list, database_csv_template,
        mad_table, primer_performance):
    overal_run_info_dict = create_run_info_dict(run_info_list)
#    print len(overal_run_info_dict)
    re = get_run_info_values(overal_run_info_dict, database_csv_template, mad_table, primer_performance)
    


def create_run_info_dict(run_info_list):
    proj_run_info_dict = {}
    for run_info in run_info_list:
#        print run_info
        run_info_dict = parse_run_info(run_info)
        proj_run_info_dict.update(run_info_dict)
    return proj_run_info_dict


def get_run_info_values(proj_run_info_dict, database_csv_template, mad_table, primer_performance):
    with open (database_csv_template) as csvfile, open ('out.tsv', 'wb+') as fout:
        reader = csv.DictReader(csvfile)
        rfd_header = reader.fieldnames
        headers_to_keep = rfd_header[:29]
        out_headers = rfd_header + ['Indexes', 'TotalReads', 'OnTargetReads', 'Median absolute deviation', 'Primer range']
        writer = csv.writer(fout, delimiter='\t')
        writer.writerow(out_headers)
        mad_dict = parse_mad_table(mad_table)
        CARRIERS_primer_range_dict = eval_primer_performance(primer_performance)
        for row in reader:
            try:
                mad_dict[row['CARRIERS ID (carriers_id)']] and proj_run_info_dict[row['CARRIERS ID (carriers_id)']]\
                        and CARRIERS_primer_range_dict[row['CARRIERS ID (carriers_id)']]
                out = [row[i] for i in headers_to_keep]  + [proj_run_info_dict[row['CARRIERS ID (carriers_id)']][-2], 
                                                            proj_run_info_dict[row['CARRIERS ID (carriers_id)']][0], 
                                                            mad_dict[row['CARRIERS ID (carriers_id)']][0],  
                                                            proj_run_info_dict[row['CARRIERS ID (carriers_id)']][-1],
                                                            'V4.0.1', 
                                                            proj_run_info_dict[row['CARRIERS ID (carriers_id)']][1], 
                                                            proj_run_info_dict[row['CARRIERS ID (carriers_id)']][2], 
                                                            '','',
                                                            proj_run_info_dict[row['CARRIERS ID (carriers_id)']][3], 
                                                            proj_run_info_dict[row['CARRIERS ID (carriers_id)']][4], 
                                                            mad_dict[row['CARRIERS ID (carriers_id)']][5],
                                                            mad_dict[row['CARRIERS ID (carriers_id)']][6],
                                                            '', 'YES', mad_dict[row['CARRIERS ID (carriers_id)']][1],
                                                            mad_dict[row['CARRIERS ID (carriers_id)']][2], 
                                                            mad_dict[row['CARRIERS ID (carriers_id)']][3],
                                                            mad_dict[row['CARRIERS ID (carriers_id)']][4]]
#                                                                CARRIERS_primer_range_dict[row['CARRIERS ID (carriers_id)']]]
                fout.write('\t'.join(str(i) for i in out) + '\n')
            except KeyError:
                print row['CARRIERS ID (carriers_id)']


def parse_run_info(run_info):        
    flowcell_id = ''
    delivery_folder = ''
    NGS_portal = ''
    sample_list = []
    run_info_dict = {}
    with open(run_info) as fin:
        for raw_line in fin:
            try:
                (key, val) = raw_line.split("=")
                if key == 'DELIVERY_FOLDER':
                    value = shlex.split(val)
                    flowcell_id += value[0].split('/')[4]
                    delivery_folder += shlex.split(val)[0]
                if key == 'PROJECTNAME':
                    NGS_portal += shlex.split(val)[0]
                if key == 'SAMPLENAMES':
                    sample_names = [i.split(':')for i in shlex.split(val)][0]
                    CARRIERS_sample_names = [i.split('_')[1] for i in sample_names]
                    sample_list.append(CARRIERS_sample_names)
            except ValueError:
                pass

    for sample in sample_list[0]:
        run_info_dict[sample] = ['s_' + sample +'.bam',
                                 delivery_folder + '/' + 'bam/' + 's_' +sample +".bam",
                                 FTP_path + "/" + flowcell_id + "/" + 'secondary/bam/' + 's_' + sample +".bam",
                                 delivery_folder + "/" + 'variants/gVCF/' + 's_' + sample + '.g.vcf.gz',
                                 FTP_path + "/" + flowcell_id + "/" + 'variants/gVCF/' + 's_' + sample + '.g.vcf.gz',
                                 NGS_portal, flowcell_id]
    return run_info_dict
                


def parse_mad_table(mad_table):
    mad_dict = {}
    with open (mad_table) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            if float(row['ReadsInCaptureRegion']) < 500000:
                mad_dict[row['sample']] = [row['lanes'], row['Indexes'], row['TotalReads'],
                                           row['ReadsInCaptureRegion'], row['Absolute_median_deviation'], 
                                           'POOR', 'Total reads mapped below 500000']
            else:
                mad_dict[row['sample']] = [row['lanes'], row['Indexes'], row['TotalReads'], row['ReadsInCaptureRegion'], 
                                           row['Absolute_median_deviation'], 'PASS', '']
    return mad_dict


def eval_primer_performance(primer_file):
    CARRIERS_primer_range_dict = {}
    df = pd.read_csv(primer_file, sep='\t')
    df['target'] = df['chrom'] +':' + df['loc5'].astype(str) + '-' + df['loc3'].astype(str)
    headers = df.columns.tolist()
    new_df = df.drop(df.columns[[0,1,2,3,4,5,6]], axis=1)
    new_headers = new_df.columns.tolist()
    rearrange_headers = new_headers[-1:] + new_headers[:-1]
    filtered_new_df = new_df[rearrange_headers]
    transporse_df = filtered_new_df.transpose()
    transporse_df.columns = transporse_df.iloc[0]
    transporse_df = transporse_df[1:]
    transporse_df = transporse_df.reset_index()
    transporse_df = transporse_df.rename(columns ={'index' : 'sample'})
    transporse_df['sample'] = transporse_df['sample'].str.split('_').str[-1]
    targets = transporse_df.columns.tolist()[1:]
    transporse_df['high'] = transporse_df[targets].max(axis=1)
    transporse_df['low'] = transporse_df[targets].min(axis=1)
    transporse_df['range'] = transporse_df[['high']].div(transporse_df.low, axis=0)
    CARRIERS_primer_range_dict = transporse_df.set_index('sample').to_dict()['range']
    return CARRIERS_primer_range_dict
#    print transporse_df[rearrange_headers[1:]]#.idxmax(axis=0, skipna=True)


#['Sample BAM name (bio_sample_bam_name)', 'NGS Lane Number (bio_ngs_lane)', 'NGS Flowcell (bio_ngs_flowcell)', 'Workflow Version (bio_workflow_vs)', 'BAM Local (bio_bam_local)', 'BAM Remote (FTP) (bio_bam_ftp)', 'VCF Local (bio_vcf_local)', 'VCF Remote (bio_vcf_remote)', 'gVCF Local (bio_gvcf_local)', 'gVCF Remote (bio_gvcf_remote)', 'Quality control pass/fail (bio_qc)', 'Quality control failed reason (bio_qc_reason)', 'Bioinformatics comments (bio_comments)', 'Bioinformatics complete (bio_complete)']
#Index, TotalReads, OnTargetReads, Median absolute deviation,Primer range

if __name__== "__main__":
    main()

