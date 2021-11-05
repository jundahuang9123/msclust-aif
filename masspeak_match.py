#!/usr/bin/env python
"""
Author: Junda Huang

Script using matchms scoring functions to get match peaks for spectra 
comparison, regarding abundances/intensity

Command line input: python3 masspeak_match.py <query>.mgf 
                    <referrence_library>.mgf <output>.txt
"""

# import statement
from sys import argv
from matchms.importing import load_from_mgf
from matchms.filtering import default_filters

def matchms_to_file(q_file, r_file, output):
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
    with open (output, 'w') as file:
        queries = load_from_mgf(q_file)
        for q in queries: 
            q = default_filters(q)
            ref_spectra = []
            for r in load_from_mgf(r_file):
                tolerence = float(q.metadata['precursormz'])
                target = float(r.metadata['pepmass'][0])
                if abs(target - tolerence) <= 0.01 and \
                        abs((int(r.metadata['rtinseconds'])/60) \
                            -(int(q.metadata['retention'])/1000000)) < 0.4:
                    ref_spectra.append(default_filters(r))
            if ref_spectra != []:
                for r in ref_spectra:
                    match_peak = 0
                    r_total = len(r.peaks.mz)
                    q_total = len(q.peaks.mz)
                    for r_peak in r.peaks.mz:
                        for q_peak in q.peaks.mz:
                            if abs(float(r_peak) - float(q_peak)) <= 0.01:
                                match_peak += 1
                    if match_peak != 0:
                        file.write("AIF:{} {:2f}\n".format\
                            (q.metadata['precursormz'], \
                                int(q.metadata['retention'])/1000000))
                        file.write("MSMS:{} {:2f}\n".format\
                            (r.metadata['pepmass'][0], \
                                int(r.metadata['rtinseconds'])/60))
                        file.write(f"Match peaks:{match_peak}\n")
                        file.write(f"{q_total - match_peak} unique AIF \
                            {r_total - match_peak} unique MSMS\n")
                        file.write('---------------\n')
    return

# main
if __name__ == '__main__':
    # parsing files from command line
    query = argv[1]
    references_library = argv[2]
    output_file = argv[3]

    matchms_to_file(query, references_library, output_file)
