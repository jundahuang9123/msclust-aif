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
import numpy as np
from matchms.importing import load_from_mgf
from matchms.importing import load_from_msp
from matchms.filtering import default_filters
from matchms.filtering import normalize_intensities
from matchms import calculate_scores, Spectrum, Scores
from matchms.similarity import CosineGreedy
from matchms.Spikes import Spikes

# functions

def matchms_score(q, q_spectra, ref_spectra, reference):
    """
    """
    if '.msp' in reference:
        for r in load_from_msp(reference):
            #tolerence = float(q.metadata['precursormz'])
            tolerence = float(q.metadata['pepmass'][0])
            if 'precursormz' in r.metadata.keys():
                target = float(r.metadata['precursormz'].replace(',','.'))
                if abs(target - tolerence) <= 0.5:
                    r = normalize_intensities(default_filters(r))
                    ref_spectra.append(r)
    elif '.mgf' in reference:
        for r in load_from_mgf(reference):
            tolerence = float(q.metadata['precursormz'])
            #tolerence = float(q.metadata['pepmass'][0])
            if 'pepmass' in r.metadata.keys():
                target = float(r.metadata['pepmass'][0])
                #target = float(r.metadata['precursormz'])
                if abs(target - tolerence) <= 0.01:
                #if abs(target - tolerence) <= 5 * tolerence * (10**(-6)):
                    r = normalize_intensities(default_filters(r))
                    ref_spectra.append(r)
    if ref_spectra != []:
        matches = calculate_scores(ref_spectra, q_spectra, CosineGreedy())
        return matches
    else:
        return None

def un_normalize(spectrum):
    """
    """
    s = spectrum.clone()
    max_in = np.max(s.peaks.intensities)
    mz, intensity = s.peaks
    un_norm = intensity * max_in
    s.peaks = Spikes(mz = mz, intensities = un_norm)
    return s


def matchms_to_file(mgf_file, references_file, output):
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
            matches = matchms_score(q, q_spectra, ref_spectra, \
                references_file)
            if matches != None:
                for match in matches:
                    (reference, query, match) = match
                    if reference is not query and match["matches"] >= 1 and \
                        abs((int(reference.metadata['rtinseconds'])/60) \
                            -(int(query.metadata['retention'])/1000000))<=0.4:
                        #out.write(f"Reference precursormz:\
                        #    {reference.metadata['precursormz']}\n")
                        out.write(f"Reference precursormz:\
                            {reference.metadata['pepmass']}\n")
                        out.write(f"Reference rettime:\
                            {int(reference.metadata['rtinseconds'])/60}\n")
                        #out.write(f"Reference Name:\
                        #    {reference.metadata['name']}\n")
                        if 'formula' in reference.metadata.keys():
                            out.write(f"Reference Formula:\
                                {reference.metadata['formula']}\n")
                        out.write(f"Query rettime:\
                            {int(query.metadata['retention'])/1000000}\n")
                        #out.write(f"Query rettime:\
                        #    {int(query.metadata['rtinseconds'])/60}\n")
                        out.write(f"Score: {match['score']:.4f}\n")
                        out.write(f"Number of matching peaks:\
                            {match['matches']}\n")
                        out.write("----------------------------\n")
                        if match['score'] >= 0.95 and match["matches"] >= 5:
                            spec_r = un_normalize(reference).plot()
                            spec_q = un_normalize(q_spectra[0]).plot()
                            spec_r.savefig\
                                ('msclust-aif/newspecmatch/05012022/outputlibmsms/lib{}.png'.\
                                format(reference.metadata['pepmass'][0]))
                            spec_q.savefig\
                                ('msclust-aif/newspecmatch/05012022/outputlibmsms/aif{}.png'.\
                                format(query.metadata['precursormz']))
                        '''
                        spec_r = un_normalize(reference).plot()
                        spec_q = un_normalize(q_spectra[0]).plot()
                        spec_r.savefig('output/{}.png'.\
                            format(reference.metadata['pepmass']))
                        spec_q.savefig('output/{}.png'.\
                            format(query.metadata['precursormz']))'''
            """
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
                        out.write(f"Score: {match['score']:.4f}\n")
                        out.write(f"Number of matching peaks:\
                            {match['matches']}\n")
                        out.write("----------------------------\n")"""
    return matches

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
    matches = matchms_to_file(query, ref, output)
    return matches

# main
if __name__ == '__main__':
    # parsing files from command line
    query = argv[1]
    references_library = argv[2]
    output_file = argv[3]

    # main
    scores = main(query, references_library, output_file)