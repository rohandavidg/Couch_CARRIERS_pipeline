#!/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bcbio/anaconda/bin/python

"""
script to subset the gVCF
"""

import os
import argparse
import csv
import logging
import datetime
import collections
import time
from collections import defaultdict
import subprocess


JAVA="/usr/local/biotools/java/jdk1.7.0_67/bin/java"
GATK="/projects/bsi/bictools/apps/alignment/GenomeAnalysisTK/3.2-2/GenomeAnalysisTK.jar"
REFERENCE="/data2/bsi/reference/sequence/human/ncbi/hg19/indexed/allchr.fa"
QSUB="/home/oge/ge2011.11/bin/linux-x64/qsub"



def main():
    args = parse_args()
    run(args.gcvf_path, args.bed_file)


def run(gvcf_path, bed_file):
    bed_file_dict, bin_count = break_bed_file(bed_file)
    sub_bed_file_path = create_bin_dir(gvcf_path, bed_file_dict, bin_count)
    write_qsub_command, chrom_vcf_list = qsub_command(gvcf_path, sub_bed_file_path)
    run_qsub_command = submit_qsub(write_qsub_command, bin_count)
    write_catvariant_command = create_bash_catvariant(chrom_vcf_list, gvcf_path)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-g', dest='gcvf_path',
                        help="path to the gvcf directory", required=True)
    parser.add_argument('-b', dest='bed_file',
                        help='CARRIERS slopped 5 target bed')
    args = parser.parse_args()
    return args


def validate_gvcf(gvcf_path):
    assert os.path.isdir(gvcf_path)
    for root, dirs, filename in os.walk(gvcf_path):
        for name in filename:
            if name.endswith(".g.vcf.gz"):
                gvcf = root + "/" + name
                yield gvcf


def break_bed_file(bed_file):
    bed_dict = defaultdict(list)
    with open(bed_file) as fin:
        header = ['chrom', 'start','stop']
        reader = csv.DictReader(fin, delimiter='\t', fieldnames=header)
        for row in reader:
            bed_dict[row['chrom']].append([row['start'], row['stop']])
    return bed_dict, len(bed_dict)


def create_bin_dir(gvcf_path, bed_dict, dir_count):
    bed_path = []
    bin_dir = "/".join(gvcf_path.split("/")[:-1]) + "/" + "chrom_vcf"
    for key, value in bed_dict.items():
        chrom_bin_dir = bin_dir + "/" + key
        if os.path.exists(chrom_bin_dir):
            filename = chrom_bin_dir + '/' + key + ".bed"
            bed_path.append(filename)
        else:
            os.makedirs(chrom_bin_dir, 0755)
            filename = chrom_bin_dir + '/' + key + ".bed"
            bed_path.append(filename)
            with open(filename, 'wa+') as fout: 
                for pos in value:
                    fout.write(key + "\t" + pos[0] + "\t" + pos[1] + "\n")
    return bed_path


def qsub_command(gvcf_path, bed_path):
    gvcf_file = validate_gvcf(gvcf_path)
    gvcf_path = []
    bash_genotype_path = []
    chrom_variant_list = []
    for gfile in gvcf_file:
        gvcf_path.append(" --variant " + gfile)

    for bed in bed_path:
        dir_name = "/".join(bed.split("/")[:-1])
        job_script_name = dir_name + "/" + "run_genotype_gVCF.sh"
        bash_genotype_path.append(job_script_name)
        chrom = os.path.basename(bed)
        with open(job_script_name, 'wa+') as fout:
            fout.write("#/bin/bash" + "\n")
            fout.write('\n')
            fout.write('#$ -q 1-day' + '\n')
            fout.write('#$ -l h_vmem=5G' + '\n')
            fout.write('#$ -pe threaded 16' + '\n')
            fout.write('#$ -M gnanaolivu.rohandavid@mayo.edu' + '\n')
            fout.write('#$ -m ae' + '\n')
            fout.write('#$ -V' + '\n')
            fout.write('#$ -cwd' + '\n')
            fout.write('#$ -N ' + chrom + "_GenotypeGvcf" + '\n')
            job_string = JAVA + " -Xmx24g -jar " + GATK + " -T GenotypeGVCFs -L " + bed + " -R " + REFERENCE 
            fout.write('\n')
            output_file = chrom + "_variants.vcf.gz"
            output_file_path = dir_name + "/" + output_file
            chrom_variant_list.append(output_file_path)
            fout.write(job_string + " ".join(str(i) for i in gvcf_path) + " -o " + chrom + "_variants.vcf.gz")
    return bash_genotype_path, chrom_variant_list


def submit_qsub(bash_genotype_path, chrom_count):
    assert len(bash_genotype_path) == chrom_count
    for path in bash_genotype_path:
        CWD = os.path.dirname(path)
        cmd = QSUB + " " + path
        print cmd
#        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=CWD)
#        run_job, error = p.communicate()
#        print run_job


def create_bash_catvariant(chrom_variant_list, gvcf_path):
    vcf_list = sorted([])
    variant_dir = "/".join(gvcf_path.split('/')[:-1])
    for vcf in chrom_variant_list:
        vcf_list.append(" -V " + vcf)
    
    catvariant_job_name = variant_dir + "/" + "run_catvariants.sh"
    with open(catvariant_job_name, 'wa+') as fout:
        fout.write("#/bin/bash" + "\n")
        fout.write('\n')
        fout.write('#$ -q 1-day' + '\n')
        fout.write('#$ -l h_vmem=6G' + '\n')
        fout.write('#$ -pe threaded 16' + '\n')
        fout.write('#$ -M gnanaolivu.rohandavid@mayo.edu' + '\n')
        fout.write('#$ -m ae' + '\n')
        fout.write('#$ -V' + '\n')
        fout.write('#$ -cwd' + '\n')
        fout.write('#$ -N CARRIERS_catvariants' + '\n')
        job_string = JAVA + " -cp " + GATK + " org.broadinstitute.gatk.tools.CatVariants -R " + REFERENCE
        fout.write('\n')
        output_file = "variants.vcf.gz"
        output_file_path = variant_dir + "/" + output_file
        chrom_variant_list.append(output_file_path)
        fout.write(job_string + " ".join(str(i) for i in vcf_list) + " --assumeSorted -out variants.vcf.gz")

    

if __name__ == "__main__":
    main()
