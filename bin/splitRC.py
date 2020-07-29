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
    
    def __init__(self, rc_file, target_coverage, sample_list):
        
        """
            initialise the GetRegion class and set up file parsing

        """
        super().__init__(rc_file, target_coverage)
        #self.outfile = csv.writer(outfile, delimiter="\t")
        self.sample_list = ( fname.split('.')[0]+".snps" for fname in sample_list.read().splitlines())
        

        
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
            rc_data.append(row[:7])
            snps = list(row[8])
            alleles.append(snps)
            maa.append([ row[9+i]  for i in  range(0, len(snps))])
            mia.append([ row[9+len(snps)+i]  for i in  range(0, len(snps))])
            
    
            # snp_types, snp_char, snp_type = self.SNP_type(zip(row[7],row[9:]))
    #         islice = len(snp_types) * 2
    #         if snp_type == "private":
    #             #called_row  = row[1:7 + islice] + snp_types + [snp_type]
    #             called_row  = row + snp_types + [snp_type]
    #             print("\t".join(called_row))
    #             #self.outfile.writerow(called_row)

            
    # def SNP_type(self, major_freqs):
        
    #     """
           
    #        calls whether the snp is private and public

    #     """

    #     get_coverage = lambda frac : map(int, frac.split('/'))        
    #     mf, major_freqs  = tee(major_freqs)
    #     snp_types = []
    #     snp_char  = []
    #     pvt_count = []
    #     pop = True
        
    #     for allele_rec in major_freqs:
    #          allele_type = None
    #          snp = allele_rec[0]
    #          count, cov = get_coverage(allele_rec[1])
    #          allele_type = str(count)
    #          if (count == 0):
    #              pvt_count.append(2)
                 
    #          elif (count == cov):
    #              pvt_count.append(1)
                 
    #          else:
    #              if pop:
    #                  pvt_count.append(1)
    #                  pop = False

    #          snp_char.append(snp)
    #          snp_types.append(allele_type)
        
    #     snp_type = "private" if sum(pvt_count) >= 2 else "public"

    #     return snp_types, snp_char, snp_type
            

                
if  __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A  parser for popoolation *_rc file")
    parser.add_argument('rc_file', type=argparse.FileType('r'))
    parser.add_argument('-t','--target_coverage', type=int, default = 30)
    parser.add_argument('-s','--sample_list', type = argparse.FileType('r'), default = "bam_files")
    args = parser.parse_args()
    split_sample = SplitSample(args.rc_file, args.target_coverage, args.sample_list)
    split_sample.get_data()
    #count.getSNPs()    
