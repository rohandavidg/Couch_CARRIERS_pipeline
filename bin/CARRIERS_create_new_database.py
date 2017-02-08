#!/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bcbio/anaconda/bin/python

"""
this script will create a CARRIERS sqlite database
for primer performance and target coverage
"""

import os
import pandas as pd
import numpy as np
import os
import sys
from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
#from quickwiki.model.meta import Base


Base = declarative_base()
def main():
    database_name = sys.argv[1]
    if database_name:
        print 'Creating database'
        path = 'sqlite:////data5/bsi/epibreast/m087494.couch/Couch/Huge_Breast_VCF/CARRIERS/SubProjects/CARRIERS_PRIMER_COVERAGE_DATABASE/working/'+ database_name + '.db'
        engine = create_engine(path)
        Base.metadata.create_all(engine)
        print 'created empty datbase ' + path
    else:
        print "name of database not given"


class primer_performance(Base):
    __tablename__ = 'primer_table'
    id = Column('id', Integer, primary_key=True, nullable=False)
    CARRIERS_sample = Column(String(20))# ForeignKey('primer_performance_id'))
    primer_target = Column(String(50))
    target_gene = Column(String(100))
    primer_value = Column(Integer)
#['CARRIERS_sample', 'primer_target', 'target_gene', 'primer_value']

class coverage_performace(Base):
     __tablename__ = 'coverage_table'
     id = Column('id', Integer, primary_key=True, nullable=False)
     CARRIERS_sample = Column(String(20))
     CARRIERS_target = Column(String(50))
     target_gene = Column(String(20))
     coverage_percentage_above_X = Column(Integer)
     coverage_value = Column(Integer)

if __name__ == "__main__":
    main()
