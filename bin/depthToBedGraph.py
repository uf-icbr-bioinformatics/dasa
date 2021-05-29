#!/usr/bin/env python

import sys, csv
import subprocess as sp

PYVER = sys.version_info.major

def Next(c):
    if PYVER == 2:
        return c.next()
    else:
        return c.__next__()

class Converter(object):
    blockChrom = ""
    blockStart = 0
    blockEnd = 0
    blockValue = 0
    scale = 1.0
    chrom = None
    window = None
    total = 0.0                 # For window mode
    bamfile = None
    outfile = "/dev/stdout"
    
    def __init__(self, args):
        prev = ""
        for a in args:
            if prev == "-s":
                self.scale = float(a)
                prev = ""
            elif prev == "-c":
                self.chrom = a
                prev = ""
            elif prev == "-o":
                self.outfile = a
                prev = ""
            elif prev == "-w":
                self.window = int(a)
                prev = ""
            elif a in ["-s", "-c", "-o", "-w"]:
                prev = a
            else:
                self.bamfile = a

    def initialize(self, line):
        self.blockChrom = line[0]
        self.blockStart = int(line[1])
        self.blockEnd = self.blockStart
        self.blockValue = int(line[2])
        self.total = self.blockValue

    def writeBlock(self, out):
        out.write("{}\t{}\t{}\t{}\n".format(self.blockChrom, self.blockStart-1, self.blockEnd-1, self.blockValue * self.scale))        

    def run(self):
        if self.chrom:
            cmdline = "samtools depth -r {} {}".format(self.chrom, self.bamfile)
        else:
            cmdline = "samtools depth {}".format(self.bamfile)
        proc = sp.Popen(cmdline, shell=True, stdout=sp.PIPE)
        with open(self.outfile, "w") as out:
            c = csv.reader(proc.stdout, delimiter='\t')
            line =  Next(c)
            self.initialize(line)
            if self.window:
                self.run_window(c, out)
            else:
                for line in c:
                    if line[0] != self.blockChrom:
                        self.writeBlock(out)
                        self.initialize(line)
                    elif int(line[2]) == self.blockValue:
                        self.blockEnd = int(line[1])
                    else:
                        self.writeBlock(out)
                        self.initialize(line)
                self.writeBlock(out)

    def writeWindow(self, out):
        sz = self.blockEnd - self.blockStart + 1
        out.write("{}\t{}\t{}\t{}\n".format(self.blockChrom, self.blockStart-1, self.blockEnd-1, self.total * self.scale / sz))

    def run_window(self, c, out):
        for line in c:
            if line[0] != self.blockChrom:
                self.writeWindow(out)
                self.initialize(line)
            elif int(line[1]) - self.blockStart < self.window:
                self.blockEnd = int(line[1])
                self.total += int(line[2])
            else:
                self.writeWindow(out)
                self.initialize(line)
        self.writeWindow(out)

# Usage: depthToBedGraph in.bam -o out.bedGraph -c chrom -s scale

if __name__ == "__main__":
    C = Converter(sys.argv[1:])
    C.run()
