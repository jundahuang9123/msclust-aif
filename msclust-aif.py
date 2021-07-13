#!/usr/bin/env python
"""
Author: Junda Huang

Command line input: python3 msclust-aif.py fullscanmsclust.csv 
                    aifmsclustready.csv fullscansim.csv fullscanmic.csv
                    [aifmsclust.csv] retention_time_tolerance 
                    correlation_threshold

"""

# import statements
from sys import argv
from statistics import stdev, mean
from scipy.stats import pearsonr

# functions
def peaks_parse(csv):
    """
    fullscan data parse, peak, ret.time, mass, cluster etc. into a dictionary

    input:
        file(string): input file name of the fullscan file
    output:
        #fullscan_dict(dictionary): dictionary of expression of 
                a list of tuples file
                {gene1: [(treatment1, expression), (treatment2, expression)]}
        #list_treat(list): list of treatments of desire    
    """
    peaks_dict = {}
    file_header = []
    with open (csv, 'r') as full_file:
        for line, i in enumerate(full_file):
            items = line.split(',').rstrip()
            if i == 0 and items[0] == 'peak#':
                file_header = items[0:]
            elif i == 0:
                file_header = ['peak-aif'] + items[0:]
            elif file_header[0] == 'peak-aif':
                peaks_dict[items[0]] = [(file_header[j], items[j]) \
                    for j in range(1, len(items))]
            else:
                peaks_dict[items[0]] = [(file_header[j], items[j]) \
                    for j in range(len(items))]
    return peaks_dict

def sim_mic_parse(csv):
    """
    """
    sim_mic_dict = {}
    return sim_mic_dict

def masspattern_match(list1, list2, threshold):
    """
    Match expression pattern with metabolites pattern

    input:
        list1(list): list of precursor abundances
        list2(list): list of aif fragments abundances
        threshold(float): threshold_correlation (positive number<=1)
    output:
        match(boolean): if the two lists' pattern match
    """
    if len(list1) != len(list2):
        raise ValueError\
            ("precursor/fragment sample sizes do not match")
    elif threshold > 1 or threshold < 0:
        raise ValueError ("threshold has to be postive number < 1")
    else: 
        match = True
        correlation, p = pearsonr(list1, list2)
        if correlation <= threshold and correlation >= -threshold:
            match = False
    return match, correlation

# main
if __name__ == '__main__':

# 1. parse files
    fullscan = peaks_parse(argv[1])
    aif = peaks_parse(argv[2])
    mic = sim_mic_parse(argv[4])
    sim = sim_mic_parse(argv[3])