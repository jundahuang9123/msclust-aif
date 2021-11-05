#!/usr/bin/env python
"""
Author: Junda Huang

Script reducing small peaks m/z under 90 in msms spectra

Command line input: python3 msms_processing.py msms.mgf less_msms.mgf
"""

# import statements
import os
from sys import argv

# functions

def parse_file(file):
    """
    """
    new_file = []
    with open(file, 'r') as f:
        for line in f:
            #print(line)
            if len(line) > 1:
                if line[2] == ".":
                    if int(line[:2]) <= 90:
                        continue
                    else:
                        new_file.append(line)
                else:
                    new_file.append(line)
            else:
                new_file.append(line)
    return new_file

def write_newmgf(file, new_file):
    """
    """
    new = 'newMSMSprosessed/{}'.format(file)
    with open(new,'w') as f:
        for line in new_file:
            f.write(line)
    return

def main(dir):
    """
    """
    for filename in os.listdir(dir):
        file = os.path.join(dir, filename)
        new = parse_file(file)
        write_newmgf(filename, new)
    return
# main
if __name__ == '__main__':
    # parsing files directory from command line
    dir = argv[1]
    main(dir)
   