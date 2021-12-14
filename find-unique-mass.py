#!/usr/bin/env python

from sys import argv

def parser(file):
    linelist = []
    with open (file, 'r') as f:
        for line in f:
            line = line.strip('\n')
            if line[0:5] == 'Query':
                linelist.append(line)
            else:
                continue
    return linelist

def main(file1, file2):
    l1 = parser(file1)
    l2 = parser(file2)
    a = set(l1)
    b = set(l2)
    len1 = len(a)
    len2 = len(b)
    return len1, len2

if __name__ == '__main__':
    # parsing files from command line
    f1 = argv[1]
    f2 = argv[2]
    l1, l2 = main(f1, f2)
    print(l1, l2)