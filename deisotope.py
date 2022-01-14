#!/usr/bin/env python
"""
Author: Junda Huang

Script to reduce isotope from mgf spectra

Command line input: python3 spec.mgf deiso.mgf
"""

# import statement
from sys import argv

# functions
def parser(mgf):
    specs = {}
    with open(mgf, 'r') as f:
        for i, line in enumerate(f):
            items = line.strip('\n').split(' ')
            if items[0] == 'BEGIN':
                specs[i] = []
                specs[i].append(line)
                key = i
            else:
                specs[key].append(line)

    return specs

def deisotop(dic):
    for key, spec in dic.items():
        peaks = []
        pop = []
        for line in spec:
            if ((line.split(' '))[0][0]).isdigit():
                peaks.append(line)
        for peak1 in peaks:
            mz1 = float(peak1.split(' ')[0])
            intens1 = float(peak1.split(' ')[1])
            for peak2 in peaks:
                if abs(mz1 - float(peak2.split(' ')[0])) <= 1.007 \
                    and intens1 > float(peak2.split(' ')[1]):
                    pop.append(peak2)
        pop = set(pop)
        for i in pop:
            spec.remove(i)
    return dic

def write_mgf(dic, out):
    with open(out, 'w') as f:
        for k, v in dic.items():
            for j in v:
                f.write(j)
    return



# main
if __name__ == '__main__':
    # parsing files from command line
    originmgf = argv[1]
    newmgf = argv[2]
    dic = deisotop(parser(originmgf))
    write_mgf(dic, newmgf)
