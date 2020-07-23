#!/usr/bin/env python

from collections import OrderedDict
import pandas as pd
import argparse
import pprint



__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2020"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"




class MergeGenePop(object):

   def __init__(self, genepop_files, poolfile, reference):

       self.genepop_files = genepop_files
       self.genepop_DF = OrderedDict()
       self.index = 0
       self.poolfile = poolfile
       self.reference = reference
       self.pos_file = open("pos.txt", "w")
       self.pool_file = open("pool_GenePop.txt", "w")

       
   def process_files(self):
       
       for genepop in self.genepop_files:
           fname_ref = genepop.rsplit(".",1)[0]
           with open(genepop) as genepop_fp:
               index_tmp = self.index 
               genotype_data = OrderedDict()
               i = 1
               for line in genepop_fp:
                   if line.startswith("gt ,"):
                       genotype_data[i] = line.replace("gt ,","").strip().split()
                       i += 1
                       
               self.index = self.index + len(genotype_data[i-1])
               self.genepop_DF[fname_ref] =  pd.DataFrame(genotype_data)
               print("{}\t{}\t{}".format(fname_ref, index_tmp, self.index-1), file = self.reference)
               
       pool_DF = pd.concat(self.genepop_DF).T
       pool_DF.to_csv(self.poolfile, sep = " ", header = False, index = False)
       print(pool_DF)
       self.pos_file.close()
       self.pool_file.close()
       

           
        

if  __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Merge genepop files")
    parser.add_argument('genepop_files', nargs='+', help = "list of genepop files, shell expansion wild cards may be used")
    parser.add_argument('-p', '--poolfile', default = "pool.Genepop",  type=argparse.FileType('w'))
    parser.add_argument('-r', '--reference', default = "Genepop.ref",  type=argparse.FileType('w'))
    args = parser.parse_args()
    genepop = MergeGenePop(args.genepop_files, args.poolfile, args.reference)
    genepop.process_files()
