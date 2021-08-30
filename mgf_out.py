#!/usr/bin/env python
"""
Author: Junda Huang


Command line input: python3 <aifcluster>.csv 
"""

from sys import argv
from matchms.importing import load_from_mgf
from matchms.importing import load_from_msp
from matchms.filtering import default_filters
from matchms.filtering import normalize_intensities
from matchms import calculate_scores
from matchms.similarity import CosineGreedy
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
    mass_dict = {}
    headers_list = []
    with open (file, 'r') as aif_file:
        for i, line in enumerate(aif_file):
            items = line.strip('\n').split(',')
            if i == 0:
                headers_list = items[1:4] + [items[31]]
                #print(headers_list)
            elif 'Pre' in items[0]:
                aif_dict[items[0]] = {}
                aif_dict[items[0]]['params'] = {headers_list[0] : items[1],\
                headers_list[2] : \
                    "{:.2f}".format(float(items[3])/1000000), \
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
                    append("{:.2f}".format(float(peaks[mass][2])/1000000))
                masslist = peaks[mass][4:]
                mgf_dict['intensity array'].append(max(masslist)) 
        aif_dict_list.append(mgf_dict)
    return aif_dict_list, headers_list

def mgf_write(dict, out, k_order):
    mgf_file = mgf.write(dict, output = out, key_order = k_order)
    return mgf_file

def match_ms_score(mgf_file, references_file_mgf, output):
    queries = load_from_mgf(mgf_file)
    references = load_from_msp(references_file_mgf)
    q_spectra = []
    ref_spectra = []
    for q in queries:
        for r in references:
            q = normalize_intensities(default_filters(q))
            r = normalize_intensities(default_filters(r))
            q_spectra.append(q)
            ref_spectra.append(r)
    matches = calculate_scores(references = ref_spectra, queries = q_spectra,\
     similarity_function = CosineGreedy())
    with open (output, 'w') as out:
        for match in matches:
            (reference, query, match) = match
            if reference is not query and match["matches"] >= 1:
                #out.write(f"Reference scan id: {reference.metadata['scans']}")
                #out.write(f"Query scan id: {query.metadata['scans']}")
                out.write(f"Score: {match['score']:.4f}")
                out.write(f"Number of matching peaks: {match['matches']}")
                out.write("----------------------------")
    return 

def spectra_loading():
    return

def main(file, ref):
    mgf_dict = aifcluster_read(file)
    mgf_file = mgf_write(mgf_dict)
    matches = match_ms_score(mgf_file, ref)
    return
    
# main
if __name__ == '__main__':
    # matches = main(argv[1], argv[2])
    aifdictlist, h_list = aifcluster_read(argv[1])
    mgf_file = mgf_write(aifdictlist, 'test2.mgf', h_list)
    match_ms_score('test2.mgf', argv[2], argv[3])
