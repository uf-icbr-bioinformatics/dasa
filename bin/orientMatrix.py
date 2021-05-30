#!/usr/bin/env python

import sys, csv, gzip
import json
import numpy as np

class MatFlipper(object):
    infile = ""
    genesfile = ""
    outfile = ""

    botstrand = []
    _coords = []
    _matrix = []

    def parseArgs(self, args):
        self.genesfile = args[0]
        self.infile    = args[1]
        self.outfile   = args[2]

    def readGenes(self):
        with open(self.genesfile, "r") as f:
            c = csv.reader(f, delimiter='\t')
            for line in c:
                if len(line) > 3 and line[3] == '-':
                    self.botstrand.append("{}:{}-{}".format(line[0], line[1], line[2]))
        sys.stderr.write("{} regions on bottom strand.\n".format(len(self.botstrand)))

    def readMatrix(self):
        nr = 0
        data = []
        with gzip.open(self.infile, "r") as f:
            l = f.readline().strip()
            if l[0] != "@":
                sys.stderr.write("Malformed input matrix!\n")
                return False
            self._header = json.loads(l[1:])
            c = csv.reader(f, delimiter='\t')
            for line in c:
                self._coords.append(line[:6])
                row = [float(x) for x in line[6:]]
                if line[3] in self.botstrand:
                    row = row[::-1]
                    nr += 1
                data.append(row)
        self._matrix = np.vstack(data)
        sys.stderr.write("Matrix read, {} rows, {} reversed.\n".format(len(self._matrix), nr))
        return True

    def writeMatrix(self):
        with gzip.open(self.outfile, "wt") as out:
            out.write("@")
            json.dump(self._header, out)
            out.write("\n")
            for i in range(len(self._coords)):
                out.write("\t".join(self._coords[i]) + "\t")
                out.write("\t".join([str(x) for x in self._matrix[i]]))
                out.write("\n")
        sys.stderr.write("New matrix written to {}\n".format(self.outfile))

if __name__ == "__main__":
    args = sys.argv[1:]
    MF = MatFlipper()
    MF.parseArgs(args)
    MF.readGenes()
    MF.readMatrix()
    MF.writeMatrix()
