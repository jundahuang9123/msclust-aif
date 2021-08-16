#!/usr/bin/env python
"""
Author: Junda Huang
"""

from matchms.importing import load_from_mgf
from pyteomics import mgf

def mgf_test(mgf):
    file = load_from_mgf("msclust-aif/out.mgf")
    print(file)

def mgf_write(dict):
    mgf.write()
