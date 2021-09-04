#!/usr/bin/env python
"""
Author: Junda Huang

Script to rewrite AIF clustering output into MGF spectra file type.
Then use matchms functions to get matching scores for the aif clustered 
spectra

Command line input: python3 <aifcluster>.csv <output>.mgf 
                    <references library>.msp/mgf output.txt
"""

from sys import argv
from matchms.importing import load_from_mgf
from matchms.importing import load_from_msp
from matchms.filtering import default_filters
from matchms.filtering import normalize_intensities
from matchms import calculate_scores
from matchms.similarity import CosineGreedy
from pyteomics import mgf

#def mgf_test(mgf):
#   file = load_from_mgf("msclust-aif/out.mgf")
# # print(file)

def aifcluster_read(file):
    """
    Read aifcluster information from csv file write into dictionary
    Params:
        file(string): file name of the csv file containing aifcluster info
    Returns:
        aif_dict(dict): dictionary containing aifcluster information
                        {'Precursor ion': {'info' : [infos]}; {'fragment':}}
        headers_list(list): 
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
        k_order(list):
    Returns:
        mgf_file()
    """
    mgf_file = mgf.write(dict, output = out, key_order = k_order)
    return mgf_file

def match_ms_score(mgf_file, references_file_mgf, output):
    """
    Find spectra matches from references library, calculate matches scores
    Params:
    Returns:
    """
    queries = load_from_mgf(mgf_file)
    #queries = load_from_msp('test.msp')
    #queries = load_from_mgf('Qe05947-neg.mgf')
    #references = load_from_msp(references_file_mgf)
    with open (output, 'w') as out:
        for q in queries: 
            q_spectra = [normalize_intensities(default_filters(q))]
            ref_spectra = []
            for r in load_from_msp(references_file_mgf):
                tolerence = float(q.metadata['precursormz'])
                #tolerence = float(q.metadata['pepmass'][0])
                if 'precursormz' in r.metadata.keys():
                    target = float(r.metadata['precursormz'].replace(',','.'))
                    if abs(target - tolerence) <= 0.5:
                        r = normalize_intensities(default_filters(r))
                        ref_spectra.append(r)
            if ref_spectra != []:
                for match in calculate_scores(ref_spectra, q_spectra, \
                    CosineGreedy()):
                    (reference, query, match) = match
                    if reference is not query and match["matches"] >= 5:
                        out.write(f"Reference precursormz:\
                            {reference.metadata['precursormz']}\n")
                    #print(f"Reference precursormz:\
                    #    {reference.metadata['precursormz']}\n")
                        out.write(f"Query scan id:\
                            {query.metadata['scan nr']}\n")
                        #print(query.metadata)
                        #out.write(f"Query scan id:\
                        #    {query.metadata['scans']}\n")
                    #print(f"Query scan id:\
                     #   {query.metadata['scans']}\n")
                        out.write(f"Query mass:\
                            {query.metadata['precursormz']}\n")
                    #    out.write(f"Query mass:\
                    #        {query.metadata['pepmass']}\n")    
                        out.write(f"Score: {match['score']:.4f}\n")
                        out.write(f"Number of matching peaks:\
                            {match['matches']}\n")
                    #print(f"Query mass:\
                    #    {query.metadata['pepmass']}\n")    
                    #print(f"Score: {match['score']:.4f}\n")
                    #print(f"Number of matching peaks:\
                    #    {match['matches']}\n")
                        out.write("----------------------------\n")
                    #print("----------------------------\n")
    return 

def spectra_loading():
    return

def main(file, outmgf, ref, txt_out):
    """
    Main function
    Params:
    Returns:
    """
    mgf_dict, h_list = aifcluster_read(file)
    mgf_file = mgf_write(mgf_dict, outmgf, h_list)
    match_ms_score(outmgf, ref, txt_out)
    return
    
# main
if __name__ == '__main__':
    # parsing files from command line
    aifcluster_csv = argv[1]
    mgfoutfile = argv[2]
    ref_lib = argv[3]
    output_txt = argv[4]

    # main
    main(aifcluster_csv, mgfoutfile, ref_lib, output_txt)