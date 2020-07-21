#!/usr/bin/env python
from collections import defaultdict
from collections import OrderedDict

import argparse
import pprint
import csv

__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2020"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"



class GetRegion(object):

    """ A  parser for popoolation2 *_rc file   """

    def __init__(self, rc_fp):
        
        """
            Gets the file object for reading by the argparse

        """
    
        csv_reader = csv.reader(rc_fp, delimiter='\t')
        _  =  next(csv_reader)
        
        region_dict = OrderedDict() 
        for row in csv_reader:
            region = row[0]
            region_dict.setdefault(region,[]).append(int(row[1]))

        #pprint.pprint(region_dict)    
        chr_start_end = {}
        min_max = lambda val: (min(val), max(val)) 
        for key,values in region_dict.items():
            min_x, max_x = min_max(values)
            #To stdout
            print("{}:{}-{}".format(key, min_x, max_x))
            
            
if  __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A  parser for popoolation *_rc file")
    parser.add_argument('rc_file', type=argparse.FileType('r'))
    args = parser.parse_args()
    GetRegion(args.rc_file)












    
