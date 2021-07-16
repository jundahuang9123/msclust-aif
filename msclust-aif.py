#!/usr/bin/env python
"""
Author: Junda Huang

Command line input: python3 msclust-aif.py fullscanmsclust.csv 
                    aifmsclustready.csv fullscansim.csv fullscanmic.csv
                    [aifmsclust.csv] retention_time_tolerance 
                    correlation_threshold

"""

# import statements
import csv
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
        
    """
    peaks_dict = {}
    file_header = []
    sample_list = []
    with open (csv, 'r') as full_file:
        for i, line in enumerate(full_file):
            items = line.strip('\n').split(',')
            if i == 0 and items[0] == 'peak#':
                file_header = items[0:]
            elif i == 0:
                file_header = ['peak-aif'] + items[0:]
            elif file_header[0] == 'peak-aif':
                items.insert(0, str(i))
                peaks_dict[items[0]] = {}
                for j in range(len(items)):
                    peaks_dict[items[0]][file_header[j]] = items[j]
            else:
                peaks_dict[items[0]] = {}
                for j in range(len(items)):
                    peaks_dict[items[0]][file_header[j]] = items[j]
    for sample in file_header:
        if 'Qe' in sample:
            sample_list.append(sample)
    return peaks_dict, sample_list, file_header

def sim_mic_parse(csv):
    """
    """
    prec_dict = {}
    with open (csv, 'r') as sim_mic_file:
        for i, line in enumerate(sim_mic_file):
            items = line.rstrip('\n').split(',')
            if i == 0:
                header = items[0:]
            else:
                prec_dict[items[0]] = {}
                for j in range(len(items)):
                    prec_dict[items[0]][header[j]] = items[j]
    return prec_dict

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
        if correlation <= threshold and correlation != 'nan' : #and correlation >= -threshold:
            match = False
    return match, correlation

def retention_control(ret1, ret2, retention_time_tolerance):
    """
    """
    if abs(ret1 - ret2) <= retention_time_tolerance:
        control = True
    else:
        control = False
    return control

def aif_cluster(peaks_dict, sample_list, prec_dict, \
    retention_time_tolerance, correlation_threshold):
    """
    """
    for frag, frag_info in peaks_dict.items():
        frag_info['clusterId1'] = ''
        frag_info['membership1'] = ''
    for prec, prec_info in prec_dict.items():
        for frag, frag_info in peaks_dict.items():
            list1 = [] 
            list2 = []
            if retention_control(int(frag_info['retention']), \
                int(prec_info['retention']), float(retention_time_tolerance)):
                for sample in sample_list:
                    list1.append(float(prec_info[sample]))
                    list2.append(float(frag_info[sample]))
                    match, correlation = masspattern_match\
                        (list1, list2, correlation_threshold)
                    if match:
                        frag_info['clusterId1'] = prec_info['clusterId']
                        frag_info['membership1'] = abs(correlation)
                        #print(frag_info['clusterId1'], \
                        #    frag_info['membership1'])
                    else: 
                    #elif match and frag_info['cluster1'] != '':
                    #    frag_info['clusterId2'] = prec_info['clusterId']
                    #    frag_info['membership2'] = correlation
                    #elif match and frag_info['cluster2'] != '':
                    #    frag_info['clusterId3'] = prec_info['clusterId']
                    #    frag_info['membership3'] = correlation
                        continue
            else:
                continue
    return peaks_dict

def output_write(dict, header, out):
    """
    """
    with open(out, 'w', newline='') as file:
        fieldnames = [i for i in dict['8089'].keys()]
        writer = csv.DictWriter(file, fieldnames = fieldnames)
        writer.writeheader()
        for k, v in dict.items():
            #rowdict = dict(zip(*({i : j} for i, j in zip(fieldnames, v))))
            rowdict = {}
            for i in range(len(fieldnames)):
                #print(len(fieldnames), len(v))
                if fieldnames[i] == 'peak-aif':
                    rowdict[fieldnames[i]] = k
                else:
                    rowdict[fieldnames[i]] = v[fieldnames[i]]
            writer.writerow(rowdict)
    return file

# main
if __name__ == '__main__':

# 1. parse files
    #fullscan, sample_list, header = peaks_parse(argv[1])
    aif, sample_list, header = peaks_parse(argv[2])
    mic = sim_mic_parse(argv[4])
    #sim = sim_mic_parse(argv[3])
    retention_time_tolerance = int(argv[-3])
    correlation_threshold = float(argv[-2])
    out = argv[-1]

# 2. cluster aif peaks to fullcan clusters
    out_dict = aif_cluster(aif, sample_list, mic, \
        retention_time_tolerance, correlation_threshold)

# 3: print output in csv
    aif_cluster_csv = output_write(out_dict, header, out)