#!/usr/bin/env python
"""
Author: Junda Huang

Script to cluster AIF fragments to precursor ions calculated by MSClust

Command line input: python3 msclust-aif.py fullscanmsclust.csv 
                    aifmsclustready.csv fullscansim.csv fullscanmic.csv
                    [aifmsclust.csv] retention_time_tolerance 
                    correlation_threshold Correlation_threshold_confidence
                    output

"""

# import statements
import csv
from sys import argv
from statistics import stdev, mean
from numpy.lib.type_check import nan_to_num
from scipy.stats import pearsonr

# functions
def peaks_parse(csv):
    """
    fullscan data parse, peak, ret.time, mass, cluster etc. into a dictionary
    params:
        csv(string): input file name of the fullscan file
    return:
        peaks_dict(dictionary):
        sample_list(list):
        file_header(list):
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
    sim or mic file data parse
    peak, ret.time, mass, cluster etc. into a dictionary
    params:
        csv(string): input file name of the sim or mic file
    return:
        prec_dict(dictionary): sim or mic value as the precursor ion and parse 
                                their info into dictionary
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

def masspattern_match(list1, list2, threshold, threshold_confidence):
    """
    Match expression pattern with metabolites pattern

    input:
        list1(list): list of precursor abundances
        list2(list): list of aif fragments abundances
        threshold(float): threshold_correlation (positive number<=1)
    output:
        match(boolean): if the two lists' pattern match
    """
    match = True
    if len(list1) != len(list2):
        raise ValueError\
            ("precursor/fragment sample sizes do not match")
    elif threshold > 1 or threshold < 0:
        raise ValueError ("threshold has to be postive number < 1")
    elif len(list1) < 2:
        match = False 
        correlation = 0
        return match, correlation
    else: 
        correlation, p = pearsonr(list1, list2)
        if str(correlation) == 'nan':
            correlation = nan_to_num(correlation)
        if correlation <= threshold or p > threshold_confidence:
            match = False
    return match, correlation

def retention_control(ret1, ret2, retention_time_tolerance):
    """
    check if two retention time are within the tolerance level
    params:
        ret1(int): retention time to be compared in mili seconds
        ret2(int): retention time to be compared in mili seconds
        retention_time_tolerance(float): limit set for retention time 
                                            differences window
    return:
        control(boolean): within tolerance(True) or not(False)
    """
    if abs(ret1 - ret2) <= retention_time_tolerance:
        control = True
    else:
        control = False
    return control

def aif_cluster(peaks_dict, sample_list, prec_dict, retention_time_tolerance,\
    correlation_threshold, correlation_threshold_confidence):
    """
    """
    for frag, frag_info in peaks_dict.items():
        frag_info['clusterId1'] = ''
        frag_info['membership1'] = ''
    precfinal = []
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
                    (list1, list2, correlation_threshold, \
                        correlation_threshold_confidence)
                if match:
                    precfinal.append(prec)
                    frag_info['clusterId1'] = prec_info['clusterId1']
                    frag_info['membership1'] = abs(correlation)
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
            precfinal = list(set(precfinal))
    for prec in precfinal:
        prec_new = 'Pre' + str(prec_dict[prec]['clusterId1'])
        peaks_dict[prec_new] = prec_dict[prec]
    return peaks_dict

def clusters_compare(peaks_dict, cluster_dict):
    """
    """
    
    return

def output_write(dict, header, out):
    """
    """
    with open(out, 'w', newline='') as file:
        fieldnames = [i for i in dict['8089'].keys()]
        writer = csv.DictWriter(file, fieldnames = fieldnames)
        writer.writeheader()
        for k, v in dict.items():
            rowdict = {}
            for i in range(len(fieldnames)):
                if fieldnames[i] == 'peak-aif':
                    rowdict[fieldnames[i]] = k
                else:
                    rowdict[fieldnames[i]] = v[fieldnames[i]]
            writer.writerow(rowdict)
    return file

def main():
    """
    The main function
    """
    peaks_parse(csv)
    sim_mic_parse(csv)
    aif_cluster(peaks_dict, sample_list, prec_dict, retention_time_tolerance,\
    correlation_threshold, correlation_threshold_confidence)
    output_write(dict, header, out)
    return

# main
if __name__ == '__main__':

# 1. parse files
    #fullscan, sample_list, header = peaks_parse(argv[1])
    aif, sample_list, header = peaks_parse(argv[2])
    #mic = sim_mic_parse(argv[4])
    sim = sim_mic_parse(argv[3])
    retention_time_tolerance = int(argv[-4])
    correlation_threshold = float(argv[-3])
    correlation_threshold_confidence = float(argv[-2])
    out = argv[-1]

# 2. cluster aif peaks to fullcan clusters
    out_dict = aif_cluster(aif, sample_list, sim, retention_time_tolerance, \
        correlation_threshold, correlation_threshold_confidence)

# 3: print output in csv
    aif_cluster_csv = output_write(out_dict, header, out)