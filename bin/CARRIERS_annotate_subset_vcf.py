#!/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bcbio/anaconda/bin/python

"""
script to annotate a vcf file
"""

import argparse
import os
import tabix
import logging
import csv
import collections
import datetime
import time
import vcf
import json

catalog_drill_map = {'HGMD_2015Q2_GRCh37_nodups':'HGMD', 'dbSNP_142_GRCh37p13':'dbSNP', 'Clinvar_20160515_GRCh37':'ClinVar', 'BlatED':'BLAT_ED', 'pfam_domains':'pfam_functional_domain', 'noTCGA_ExAc':'ExAc', 'dbNSFP_v3a_GRCh37':'dbNSFP', 'OMIM_phenotypes':'OMIM'}

def main():
    args = parse_args()
    run(args.subset_vcf_path, args.catalog_path, args.drill_path)


def run(subset_vcf_path, catalog_path, drill_path):
    logger = configure_logger()
    chrom_vcf = get_subset_vcf(subset_vcf_path, logger, catalog_path, drill_path)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-s', dest='subset_vcf_path',
                        help="path to the subset to chrom vcf directory", required=True)
    parser.add_argument('-c', dest='catalog_path', default="/data5/bsi/epibreast/m087494.couch/Couch/Huge_Breast_VCF/catalog_2016_Q1.file",
                        help='bior catalog path')
    parser.add_argument('-d', dest='drill_path', default="/data5/bsi/epibreast/m087494.couch/Couch/Huge_Breast_VCF/drill_2016_Q2.file",
                        help='bior drill path')
    args = parser.parse_args()
    return args


def configure_logger():
    """
    setting up logging
    """
    logger = logging.getLogger('Annotate_vcf')
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(time.strftime("Annotate_vcf-%Y%m%d.log"))
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s'\t'%(name)s'\t'%(levelname)s'\t'%(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger
    

def validate_gvcf(subset_path, logger):
    assert os.path.isdir(subset_path)
    for root, dirs, filename in os.walk(subset_path):
        for name in filename:
            if name.endswith("bedvariants.vcf.gz"):
                gvcf = root + "/" + name
                yield gvcf
#            else:
#                logger.debug('{0} file not found'.format(name))
            

def get_subset_vcf(subset_path, logger, catalog_path, drill_path):
    subset_vcf_gen = validate_gvcf(subset_path, logger)
    for sub_vcf in subset_vcf_gen:
        show_vcf = parse_vcf(sub_vcf, catalog_path, drill_path)
#        list_catalog = parse_catalog_file(catalog_path)


def parse_vcf(some_vcf, catalog_path, drill_path):
    vcf_reader = vcf.Reader(open(some_vcf, 'r'))
    for record in vcf_reader:
        info = parse_catalog_file(catalog_path, record.CHROM, record.POS, drill_path)
#        print info


def parse_catalog_file(catalog_path, chrom, pos, drill_path):
    with open(catalog_path) as fin:
        raw_line = fin.readlines()
        for line in raw_line:
            catalog = line.strip().split("\t")[2]
            db =catalog.split("/")[6]
            drill_dict = parse_drill(drill_path)
            drill_value = [v for k, v in drill_dict.items() if k == db]
            tb = tabix.open(catalog)
            chrom_number = chrom.replace('chr','')
            records = tb.query(chrom_number,int(pos)-1, int(pos))
            for record in records:
                data =  json.loads(record[3])
                for elem in drill_value:
                    elem_list = elem.split(",")
                    for i in elem_list:
                        print data[i]



def parse_drill(drill_path):
    drill_info_dict = {}
    with open(drill_path) as din:
        raw_line = din.readlines()
        for line in raw_line:
            value = line.strip().split("\t")
            catalog = value[0]
            fields = value[1]
            try:
                drill_info_dict[catalog_drill_map[catalog]] = fields
            except KeyError:
                pass
    return drill_info_dict


if __name__ == "__main__":
    main()
