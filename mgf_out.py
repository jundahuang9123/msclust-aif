#!/usr/bin/env python
"""
Author: Junda Huang
"""

from sys import argv
from matchms.importing import load_from_mgf
from pyteomics import mgf

def mgf_test(mgf):
    file = load_from_mgf("msclust-aif/out.mgf")
    print(file)

def aifcluster_read(file):
    """
    Read aifcluster information from csv file write into dictionary
    Param:
        file(string): file name of the csv file containing aifcluster info
    output:
        aif_dict(dict): dictionary containing aifcluster information
                        {'Precursor ion': {'info' : [infos]}; {'fragment':}}
    """
    aif_dict = {}
    headers_list = []
    frag_list = []
    with open (file, 'r') as aif_file:
        for i, line in enumerate(aif_file):
            items = line.strip('\n').split(',')
            if i == 0:
                headers_list = items[0:]
                
    return aif_dict

def mgf_write(dict):
    mgf.write()

def main(mgf):
    mgf_test(mgf)
    return
    
# main
if __name__ == '__main__':
    main(argv[1])