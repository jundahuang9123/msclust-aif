#!/usr/bin/env python
"""
Author: Junda Huang

Script to convert result from test file to xlsx for further analysis

Command line input: python3 txt_to_csv.py <output>.txt <output>.csv
"""

# import statements
from sys import argv
import pandas

# functions
def text_reader(file):
    """
    """
    read_file = pandas.read_csv(r'{}'.format(file))
    for i in read_file:
        print(i, read_file[i])
    return read_file

def csv_out(file, csv):
    """
    """
    file.to_csv(r'{}'.format(csv), index = None)
    return

def main(file, csv):
    """
    """
    read_file = text_reader(file)
    csv_out(read_file, csv)
    return

# main
if __name__ == '__main__':
    file = argv[1]
    csv = argv[2]

    main(file, csv)
