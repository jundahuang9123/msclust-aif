#!/usr/bin/env python
"""
Author: Junda Huang

Script using matchms scoring functions to get matching scores 
for imput spectra in mgf/msp again mgf/msp spetra references library
spectra

Command line input: python3 matchms_scoring.py <query>.mgf 
                    <referrence_library>.msp/mgf <output>.txt
"""

# import statements
from sys import argv
from matchms.importing import load_from_mgf
from matchms.importing import load_from_msp
from matchms.filtering import default_filters
from matchms.filtering import normalize_intensities
from matchms import calculate_scores
from matchms.similarity import CosineGreedy

# functions

def match_ms_score(mgf_file, references_file, output):
    """
    Find spectra matches from references library, calculate matches scores 
    and write into output file
    Params:
        mgf_file(string): file name of the query mgf file
        references_file(string): file name of the reference library file (msp)
        output(string): file name of the output text file
    Returns:
        output(string): file name of the output text file
    """
    queries = load_from_mgf(mgf_file)
    #queries = load_from_msp('test.msp')
    #queries = load_from_mgf('Qe05947-neg.mgf')
    #references = load_from_msp(references_file_mgf)
    with open (output, 'w') as out:
        for q in queries: 
            q_spectra = [normalize_intensities(default_filters(q))]
            ref_spectra = []
            for r in load_from_msp(references_file):
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
                    if reference is not query and match["matches"] >= 1:
                        out.write(f"Reference precursormz:\
                            {reference.metadata['precursormz']}\n")
                        out.write(f"Reference Name:\
                            {reference.metadata['name']}\n")
                        out.write(f"Reference Formula:\
                            {reference.metadata['formula']}\n")
                        out.write(f"Query scan id:\
                            {query.metadata['scan nr']}\n")
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
    '''
    Load matching spectra...
    '''
    return

def main(query, ref, output):
    '''
    The main function
    Params:
        mgf_file(string): file name of the query mgf file
        references_file(string): file name of the reference library file (msp)
        output(string): file name of the output text file
    Returns:
        output(string): file name of the output text file
    '''
    match_ms_score(query, ref, output)
    #spectra_loading()
    return

# main
if __name__ == '__main__':
    # parsing files from command line
    query = argv[1]
    references_library = argv[2]
    output_file = argv[3]

    # main
    main(query, references_library, output_file)