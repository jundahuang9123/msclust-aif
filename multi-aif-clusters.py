#!/usr/bin/env python
"""
Author: Junda Huang
"""
# functions
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
    dif = abs(ret1 - ret2)
    if abs(ret1 - ret2) <= retention_time_tolerance:
        control = True
    else:
        control = False
    return control, dif

def aif_cluster(peaks_dict, sample_list, prec_dict, retention_time_tolerance,\
    correlation_threshold, correlation_threshold_confidence):
    """
    """
    for frag, frag_info in peaks_dict.items():
        frag_info['clusterId1'] = []
        frag_info['membership1'] = []
        frag_info['ratio']
    precfinal = []
    for prec, prec_info in prec_dict.items():
        for frag, frag_info in peaks_dict.items():
            list1 = [] 
            list2 = []
            control, dif = retention_control(int(frag_info['retention']), \
                int(prec_info['retention']), float(retention_time_tolerance))
            if control:
                for sample in sample_list:
                    list1.append(float(prec_info[sample]))
                    list2.append(float(frag_info[sample]))
                match, correlation = masspattern_match\
                    (list1, list2, correlation_threshold, \
                        correlation_threshold_confidence)
                if match:
                    precfinal.append(prec)
                    frag_info['clusterId1'].append(prec_info['clusterId1'])
                    frag_info['membership1'].append(abs(correlation))
                    frag_info['ratio'].append(dif)
            else:
                continue
            precfinal = list(set(precfinal))
    for prec in precfinal:
        prec_new = 'Pre' + str(prec_dict[prec]['clusterId1'])
        peaks_dict[prec_new] = prec_dict[prec]
    return peaks_dict