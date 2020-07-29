#!/usr/bin/env python

from itertools import tee
from get_region import  GetRegion
from collections import defaultdict
from collections import OrderedDict
from collections import Counter
import argparse
import pprint
import csv
import sys



__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2020"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"




class SplitSample(GetRegion):

    """ Count SNPs in a popoolation2 *_rc file   """
    
    def __init__(self, rc_file, target_coverage, sample_list, biallelic):
        
        """
            initialise the GetRegion class and set up file parsing

        """
        super().__init__(rc_file, target_coverage)
        self.sample_list = ( fname.split('.')[0]+".snps" for fname in sample_list.read().splitlines())
        self.biallelic = biallelic
        
        
    def get_data(self):
        """
            extract data from the rc file
        
        """
        self.alleles = []
        self.rc_data = []
        self.maa = []
        self.mia = []
        for row in self.csv_reader:
            if row[1].startswith("#"):
                continue
            self.rc_data.append(row[:7])
            self.snps = list(row[8])
            n = len(self.snps)
            self.alleles.append(self.snps)
            self.maa.append([ row[9+i]  for i in  range(0, n)])
            self.mia.append([ row[9 + n +i]  for i in  range(0, n)])


    def write_samples(self):
         """
            write out sample files according to order in sample file list
         
         """

         biallelic = lambda row: row[3] == '2' if self.biallelic else lambda row: True
         public_snp  = lambda maa, mai:  all([int(maa.split("/")[0]), int(mai.split("/")[0])])
         
         for j,fname in enumerate(self.sample_list):
            sample_file = csv.writer(open(fname, "w"), delimiter="\t")

            for i, row  in  enumerate(self.rc_data):
                snp_check = "-"
                maa  = self.maa[i][j]
                mia  = self.mia[i][j]
                if biallelic(row):
                    if public_snp(maa, mia):
                        snp_check = "snp"
                row = row + [self.alleles[i][j]] + [maa] + [mia] + [snp_check]
                sample_file.writerow(row)

                
if  __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A  parser for popoolation *_rc file")
    parser.add_argument('rc_file', type=argparse.FileType('r'))
    parser.add_argument('-t','--target_coverage', type=int, default = 30)
    parser.add_argument('-s','--sample_list', type = argparse.FileType('r'), default = "bam_files")
    parser.add_argument('-b','--biallelic', action = "store_false", default = True, help = "turn off biallelic filtering")
    args = parser.parse_args()
    split_sample = SplitSample(args.rc_file,
                               args.target_coverage,
                               args.sample_list,
                               args.biallelic)
    split_sample.get_data()
    split_sample.write_samples() 
