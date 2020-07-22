#!/usr/bin/env python
from collections import defaultdict
from collections import OrderedDict

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



class GetRegion(object):

    """ A  parser for popoolation2 *_rc file   """

    def __init__(self, rc_fp, target_coverage):
        
        """
            Gets the file object for reading by the argparse
            
        """


        self.target_coverage = target_coverage
        comment  =  next(rc_fp)
        assert comment.startswith("##"), "First line of the *_rc files must be a comment"
        
        self.csv_reader = csv.reader(rc_fp, delimiter='\t')
        
    def get(self):
        
        """
           get region from the file

        """
        
        region_dict = OrderedDict() 
        for row in self.csv_reader:
            region = row[0]
            #get major allels and corresponding freq
            if not self.check_target(zip(row[7],row[9:])):
                sys.stdout.write("!fail "+"\t".join(row) + "\n")
                continue
            
            region_dict.setdefault(region,[]).append(int(row[1]))

        #pprint.pprint(region_dict)    
        chr_start_end = {}
        min_max = lambda val: (min(val), max(val)) 
        for key,values in region_dict.items():
            min_x, max_x = min_max(values)
            #To stdout
            print("{}:{}-{}".format(key, min_x, max_x))

            
    def check_target(self, major_freqs):
        get_coverage = lambda frac : map(int, frac.split('/'))
        for allele_rec in major_freqs:
            count, cov = get_coverage(allele_rec[1])
            if cov <= self.target_coverage:
                return False
            
        return True
            
            
if  __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A  parser for popoolation *_rc file")
    parser.add_argument('rc_file', type=argparse.FileType('r'))
    parser.add_argument('-t','--target_coverage', type=int, default = 30)
    args = parser.parse_args()
    region = GetRegion(args.rc_file, args.target_coverage)
    region.get()












    
