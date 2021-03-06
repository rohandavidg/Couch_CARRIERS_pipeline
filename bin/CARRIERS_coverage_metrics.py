#/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bcbio/anaconda/bin/python


"""
creates a dataframe for coverage metrics using bams
"""

import pysam
import os
import collections
import argparse
import sys
import time
import datetime
import subprocess
from crimson import flagstat
import json
from random import randint
#import pyximport
#pyximport.install()
#import _pysam_flagstat


GATK="/data5/bsi/bictools/alignment/gatk/3.5/GenomeAnalysisTK.jar"
JAVA="/usr/local/biotools/java/jdk1.7.0_67/bin/java"
REFERENCE="/data2/bsi/reference/sequence/human/ncbi/37.1/indexed/allchr.fa"
#GENES="/data5/bsi/epibreast/m087494.couch/Couch/Huge_Breast_VCF/CARRIERS/SubProjects/CARRIERS_genes_txt.refgene.sort"
INTERVALS="/data5/bsi/epibreast/m087494.couch/Couch/Huge_Breast_VCF/CARRIERS/SubProjects/CARRIERS_targets.sort.exons.bed"
CAPTURE_BED="/data5/bsi/epibreast/m087494.couch/Couch/Huge_Breast_VCF/CARRIERS_PANC.targets.bed"
QSUB="/home/oge/ge2011.11/bin/linux-x64/qsub"
SAMTOOLS="/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bin/samtools"
GENES="/data5/bsi/epibreast/m087494.couch/Couch/Huge_Breast_VCF/CARRIERS/SubProjects/CARRIERS_genes_txt.refgene.sort.withoutheader.chr.ordered"

def main():
    args = parse_args()
    run(args.bam_path, args.bam_path_list, args.lane, args.project)


def run(bam_path, bam_path_list, lane, project):
#    flagstat_out = get_flagstats(bam_path)
#    flagstat_to_json = convert_flagstat(flagstat_out)
    if bam_path:
        all_bam_path, coverage_path = bam_path_ggps(bam_path)
        if project == 'all':
            create_depth_of_coverage = run_dept_coverage_all(all_bam_path, lane)
        else:
            create_depth_of_coverage = run_dept_coverage(all_bam_path, coverage_path, lane)        
    else:
        pass
    if bam_path_list:
        lane_bam_path, coverage_path = get_bam_all(bam_path_list)
        if project == 'all' or project == 'ALL':
            create_depth_of_coverage = run_dept_coverage_all(lane_bam_path, lane)
        else:
            create_depth_of_coverage = run_dept_coverage(lane_bam_path, coverage_path, lane)
    else:
        pass



def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-b', dest='bam_path',
                        help="path to the ggps bam directory")
    parser.add_argument('-B', dest='bam_path_list', help="list of bams")#, type=argparse.FileType('r'))
    parser.add_argument('-l', dest='lane', help='flowcell lane', required=True)
    parser.add_argument('-a', dest='project', help='when running for all bams in project')
    args = parser.parse_args()
    return args


def get_bam_files(bam_path):
    assert os.path.isdir(bam_path)
    for root, dirs, filename in os.walk(bam_path):
        for name in filename:
            if name.endswith(".bam") and not name.endswith("extra.bam"):
                bam = root + "/" + name
                yield bam

 
def get_flagstats(bam_path):
    flagstat_list = []
    bam_gen = get_bam_files(bam_path)
    dir_name = os.path.dirname(bam_path)
    flagstat_dir = dir_name + "/" + "flagstat"
    if os.path.isdir(flagstat_dir):
        pass
    else:
        os.makedirs(flagstat_dir, 0755 )
        for bam in bam_gen:
            sample_bam_name = os.path.basename(bam)
            sample_name = sample_bam_name.split(".")[0]
            flagstat_file = flagstat_dir + "/" + sample_name + ".flagstat"
            flagstat_list.append(flagstat_file)
            cmd = SAMTOOLS + " flagstat " + bam + " > " + flagstat_file
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            sam_cmd, err = p.communicate()
    return flagstat_list


def get_bam_all(bam_path_list):
    bam_list = []
    bam_dir = os.path.dirname(bam_path_list)
    run_dir = "/".join(bam_dir.split("/")[:-1])
    coverage_dir = run_dir + "/" + "coverage"
    with open(bam_path_list, 'r') as bpin:
        for bam in bpin:
            new_bam = bam.strip()
            bam_list.append(new_bam)

    if os.path.isdir(coverage_dir):
        pass
    else:
        os.makedirs(coverage_dir)
    return bam_list, coverage_dir


def bam_path_ggps(bam_path):
    bam_path_lane_list = []
    bam_file = get_bam_files(bam_path)
    for bam in bam_file:
        bam_dir = os.path.dirname(bam)
        #        print bam_dir
        run_dir = "/".join(bam_dir.split("/")[:-1])
        coverage_dir = run_dir + "/" + "coverage"        
#        print coverage_dir
        bam_path_lane_list.append(" -I " + bam.strip())
        if os.path.isdir(coverage_dir):
            pass
        else:
            os.makedirs(coverage_dir)
    return bam_path_lane_list, coverage_dir


def run_dept_coverage(some_bam_list, coverage_dir, lane):
    if lane:
#        print lane
        lane_dir = coverage_dir + "/" + "lane" + lane
#        print lane_dir
        if os.path.isdir(lane_dir):
            pass
        else:
            os.makedirs(lane_dir)
        script_name = lane_dir + "/" + "run_depth_of_coverage" + lane +".sh"
        with open(script_name, 'wa+') as fout:
            fout.write("#/bin/bash" + "\n")
            fout.write('\n')
            fout.write('#$ -q 7-days' + '\n')
            fout.write('#$ -l h_vmem=10G' + '\n')
            fout.write('#$ -pe threaded 16' + '\n')
            fout.write('#$ -M gnanaolivu.rohandavid@mayo.edu' + '\n')
            fout.write('#$ -m ae' + '\n')
            fout.write('#$ -V' + '\n')
            fout.write('#$ -cwd' + '\n')
            fout.write('#$ -N gatk_depth_of_coveage'+lane + '\n')
            job_string = JAVA + " -Xmx44g -jar " + GATK + " -T DepthOfCoverage -L " + INTERVALS + " -R " + REFERENCE + " "
            fout.write('\n')
            output_file = "CARRIERS_depth_of_coverage.txt"
            output_file_path = lane_dir + "/" + output_file
        #            chrom_variant_list.append(output_file_path)
            if some_bam_list:
                fout.write(job_string + " -I " + " -I ".join(str(i) for i in some_bam_list) + " -o " + output_file + " --summaryCoverageThreshold 1 --summaryCoverageThreshold 10 --summaryCoverageThreshold 20 --summaryCoverageThreshold 50 --summaryCoverageThreshold 100 --summaryCoverageThreshold 150 -dt NONE --calculateCoverageOverGenes:REFSEQ " + GENES + " -omitBaseOutput --omitDepthOutputAtEachBase --omitLocusTable")
            else:
                pass
    else:
        print "lane information not given, specify ALL for all lanes"

def run_dept_coverage_all(some_bam_list, lane):
    coverage_dir = os.getcwd()
    if lane:
#        print lane
        lane_dir = coverage_dir + "/" + "coverage/lane" + lane
#        print lane_dir
        if os.path.isdir(lane_dir):
            pass
        else:
            os.makedirs(lane_dir)
        script_name = lane_dir + "/" + "run_depth_of_coverage" + lane +".sh"
        with open(script_name, 'wa+') as fout:
            fout.write("#/bin/bash" + "\n")
            fout.write('\n')
            fout.write('#$ -q 7-days' + '\n')
            fout.write('#$ -l h_vmem=5G' + '\n')
            fout.write('#$ -pe threaded 16' + '\n')
            fout.write('#$ -M gnanaolivu.rohandavid@mayo.edu' + '\n')
            fout.write('#$ -m ae' + '\n')
            fout.write('#$ -V' + '\n')
            fout.write('#$ -cwd' + '\n')
            fout.write('#$ -N gatk_depth_of_coveage'+lane + '\n')
            job_string = JAVA + " -Xmx44g -jar " + GATK + " -T DepthOfCoverage -L " + INTERVALS + " -R " + REFERENCE + " " 
            fout.write('\n')
            output_file = "CARRIERS_depth_of_coverage.txt"
            output_file_path = lane_dir + "/" + output_file
        #            chrom_variant_list.append(output_file_path)
            if some_bam_list:
                fout.write(job_string + " -I " + " -I ".join(str(i) for i in some_bam_list) + " -o " + output_file + " --summaryCoverageThreshold 1 --summaryCoverageThreshold 10 --summaryCoverageThreshold 20 --summaryCoverageThreshold 50 --summaryCoverageThreshold 100 --summaryCoverageThreshold 150 -dt NONE --calculateCoverageOverGenes:REFSEQ " + GENES + " -omitBaseOutput --omitDepthOutputAtEachBase --omitLocusTable")
            else:
                pass
    else:
        print "lane information not given, specify ALL for all lanes"

if __name__ == '__main__':
    main()
