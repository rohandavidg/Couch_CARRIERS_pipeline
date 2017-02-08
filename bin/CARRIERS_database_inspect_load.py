#!/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bcbio/anaconda/bin/python


"""
this script will inspect the tables in the db
and load the data in the database
"""

import os
import pandas as pd
import numpy as np
import os
import sys
import argparse
from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from pandas.io import sql 
import pypyodbc
import odo
from odo import discover
from odo import odo
from odo import resource

database_path = "sqlite:////data5/bsi/epibreast/m087494.couch/Couch/Huge_Breast_VCF/CARRIERS/SubProjects/CARRIERS_PRIMER_COVERAGE_DATABASE/working/"
Table1 = 'primer_performance'
Table2 = 'target_coverage'

def main():
    args = parse_args()
    run(args.database_name, args.primer_performance, args.coverage_performance)


def run(database_name, primer_performance, coverage_performance):
    if database_name:
        database = database_path + database_name +'.db'
        establish_connection = check_database(database)
        primer_performance_db = primer_data_db(primer_performance, database)
        coverage_performance_db = coverage_data_db(coverage_performance, database)
    else:
        print 'argument not given'


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-d', dest='database_name',
                        help="name of the database to upload data",
                        required=True)
    parser.add_argument('-p', dest='primer_performance',
                        help="primer performance file",
                        required=True)
    parser.add_argument('-c', dest='coverage_performance',
                        help="coverage performance file",
                        required=True)    
    args = parser.parse_args()
    return args


def check_database(database):
    connection = create_engine(database).connect()
#    try:
    result = connection.execute('SELECT * from primer_table')
    if result:
        return True
    else:
        return False
#    except OperationalError:
#        print "ERROR IN CONNECTION: check if name of database is correct"


def primer_data_db(primer_performance, database):
    engine = create_engine(database).connect()
    df = pd.read_csv(primer_performance, sep='\t')
    headers = df.columns.tolist()
    df['primer_target'] = df['chrom'] + '_' + df['loc5'].astype(str) + '_' + df['loc3'].astype(str)
    new_df = df.drop(df.columns[[0,1,2,3,4,6]], axis=1)
    new_headers = new_df.columns.tolist()[:-1]
    new_df = pd.melt(new_df, id_vars=['primer_target', 'target_gene'], var_name='CARRIERS_sample', value_name="primer_value")
    new_headers = new_df.columns.tolist()
    format_header = ['CARRIERS_sample', 'primer_target', 'target_gene', 'primer_value']
    format_new_df = new_df[format_header]
    sorted_new_df = format_new_df.sort(['CARRIERS_sample', 'target_gene'])
    sorted_new_df.to_sql(name='primer_table', con=engine, if_exists = 'replace', index=False, chunksize=200000)

def coverage_data_db(coverage_performance, database):
    engine = create_engine(database).connect()
    df = pd.read_csv(coverage_performance, sep='\t')
    headers = df.columns.tolist()
    new_header = ['Target', 'primers'] + [i for i in headers if '%_above_20' in i] +  [i for i in headers if '%_above_50' in i] + [i for i in headers if '%_above_100' in i]
    new_df =  df[new_header]
    melt_df = pd.melt(new_df, id_vars=['Target', 'primers'], var_name='CARRIERS_sample_with_coverage', value_name='coverage_value')
    melt_df['CARRIERS_sample'] = melt_df['CARRIERS_sample_with_coverage'].str.split('_').str[1]
    melt_df['coverage_percentage_above_X'] = melt_df['CARRIERS_sample_with_coverage'].str.split('_').str[-1]
    melt_df = melt_df.rename(columns={'primers': 'target_gene', 'Target':'CARRIERS_target'})
    coverage_headers = ['CARRIERS_sample', 'CARRIERS_target', 'target_gene', 'coverage_percentage_above_X', 'coverage_value']
    coverage_df =  melt_df[coverage_headers]
    sorted_coverage_df = coverage_df.sort(['CARRIERS_sample', 'target_gene'])
    sorted_coverage_df.to_sql(name='coverage_table',  con=engine, if_exists = 'replace', index=False, chunksize=200000)

if __name__ == "__main__":
    main()
