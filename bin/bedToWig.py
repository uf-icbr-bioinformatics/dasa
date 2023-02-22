#!/usr/bin/env python

import sys
import csv

def main(filename):
    sys.stdout.write("track type=wiggle_0\n")
    with open(filename, "r") as f:
        c = csv.reader(f, delimiter='\t')
        for row in c:
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            sys.stdout.write("variableStep chrom={} span={}\n".format(chrom, end-start))
            sys.stdout.write("{}\t1\n".format(start))

main(sys.argv[1])

