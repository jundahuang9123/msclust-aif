#!/usr/bin/env python
"""
Author: Junda Huang


Command line input: python3 matchms_scoring.py <query>.mgf 
                    <referrence_library>.msp/mgf <output>.txt
"""

# import statements
from sys import argv
import numpy as np
from matchms.importing import load_from_msp
from matchms import calculate_scores, Spectrum, Scores

# functions
def spectra_loading(msp, mz, name):
    '''
    Load matching spectra...

    '''
    for ref in load_from_msp(msp):
        if 'precursormz' in ref.metadata.keys() \
            and 'name' in ref.metadata.keys():
            if ref.metadata['precursormz'] == mz and \
                ref.metadata['name'] == name:
                spec_r = ref.plot()
                spec_r.savefig('Library{},{}.png'.\
                    format(ref.metadata['precursormz'], \
                        ref.metadata['name']))
    return

# main
if __name__ == '__main__':
    # parsing files from command line
    msp = argv[1]
    mz = argv[2]
    name =argv[3]

    spectra_loading(msp, mz, name)