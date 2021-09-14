#!/usr/bin/env python
"""
Author: Junda Huang

Script to rewrite AIF clustering output into MGF spectra file type.
Then use matchms functions to get matching scores for the aif clustered 
spectra

Command line input: python3 aifcluster_to_mgf.py <aifcluster>.csv <output>.mgf 
"""

# import statements
from sys import argv
from pyteomics import mgf

# functions
def aifcluster_read(file):
    """
    Read aifcluster information from csv file write into dictionary
    Params:
        file(string): file name of the csv file containing aifcluster info
    Returns:
        aif_dict(dict): dictionary containing aifcluster information
                        {'Precursor ion': {'info' : [infos]}; {'fragment':}}
        headers_list(list): list of information of the aif clusters
                            [retentiontime, m/z, ...]
    """
    aif_dict = {}
    mass_dict = {}
    headers_list = []
    with open (file, 'r') as aif_file:
        for i, line in enumerate(aif_file):
            items = line.strip('\n').split(',')
            if i == 0:
                headers_list = items[1:4] + [items[31]]
            elif 'Pre' in items[0]:
                aif_dict[items[0]] = {}
                aif_dict[items[0]]['params'] = {headers_list[0] : items[1],\
                headers_list[2] : \
                    "{:.4f}".format(float(items[3])/1000000), \
                        headers_list[1] : items[2], \
                headers_list[3] : items[31]}
            else:
                mass_dict[items[0]] = items[1:]
    for mass, peaks in mass_dict.items():
        for prec in aif_dict:
            if peaks[30] == aif_dict[prec]['params']['clusterId1']:
                aif_dict[prec][mass] = [i for i in peaks[:29]]
    aif_dict_list = []
    for prec, peaks in aif_dict.items():
        mgf_dict = {}
        mgf_dict['params'] = peaks['params']
        mgf_dict['m/z array'] = []
        mgf_dict['intensity array'] = []
        for mass in peaks:
            if 'params' not in mass:
                mgf_dict['m/z array'].\
                    append("{:.4f}".format(float(peaks[mass][2])/1000000))
                masslist = peaks[mass][4:]
                mgf_dict['intensity array'].append(max(masslist)) 
        aif_dict_list.append(mgf_dict)
    return aif_dict_list, headers_list

def mgf_write(dict, out, k_order):
    """
    write AIF clusters spectra into mgf files
    Params:
        aif_dict(dict): dictionary containing aifcluster information
                        {'Precursor ion': {'info' : [infos]}; {'fragment':}}
        k_order(list): list of information of the aif clusters
                            [retentiontime, m/z, ...]
    Returns:
        mgf_file(string): filename of the output mgf file
    """
    mgf_file = mgf.write(dict, output = out, key_order = k_order)
    return mgf_file

def main(file, outmgf):
    """
    Main function
    Params:
        file(string): file name of the csv file containing aifcluster info
        outmgf(string): filename of the output mgf file
    Returns:
        outmgf(string): filename of the output mgf file
    """
    mgf_dict, h_list = aifcluster_read(file)
    mgf_file = mgf_write(mgf_dict, outmgf, h_list)
    return
    
# main
if __name__ == '__main__':
    # parsing files from command line
    aifcluster_csv = argv[1]
    mgfoutfile = argv[2]

    # main
    main(aifcluster_csv, mgfoutfile)