#!/usr/bin/env python

import os
import sys
import csv
import glob
import time
import math
import os.path

from BIutils import BIfile, BIexperiment, BImisc, BIcsv, BIcolors, BItext

"""
Strategy:

1. Run MACS on each sample - generate list of peaks for each sample (done in advance).
2. Merge peaks for samples in each condition. Find common / unique peaks.
   Compute number of reads and open region for each sample.
3. Determine normalization factor for each sample.
4. Quantify common peaks in each sample applying normalization factors.
5. Write matrix for DESeq2.
"""

## Utils

def parseCoords(s):
    """Parse a string of the form chr:start-end and return the three components as a tuple."""
    p1 = s.find(":")
    p2 = s.find("-")
    return (s[:p1], s[p1+1:p2], s[p2+1:])

def parseClassification(filename):
    """Parse the output of the classify script."""
    d = {}
    with open(filename, "r") as f:
        for line in f:
            count = line[:7].strip()
            record = line[8:].split("_")[0]
            d[int(record)] = int(count)
    return d

### Writer (helper class that writes scripts)

class Writer(object):

    def writePeakMerger(self, scriptname, outfile, infiles):
        with BImisc.ShellScript(scriptname) as out:
            out.write("""

sort -m -k1,1 -k2,2n -k3,3n {} | mergeBed -i - > {}.pre
            """.format(" ".join([ s.pathname("bedfile") for s in infiles ]), outfile))
        return scriptname

    def writeFRIPQsub(self):
        scriptname = "_scripts/frip.sh"
        with BImisc.ShellScript(scriptname) as out:
            out.write("""
module load bedtools
coverageBed -sorted -a $1 -b $2 | awk '{sum+=$7} END {print sum}' - > $3
""")
        return scriptname

    def writeQsubAnnotate(self, mgr):
        name = "_scripts/annot.qsub"
        with BImisc.ShellScript(name) as out:
            out.write("""
module load dibig_tools
IN=$1                           # file of peak regions
OUT=$2                          # genes associated with peak regions

genes.py classify -t -ca -X -db {} -o $OUT -d {} -c 4:log2\(FC\) @$IN
""".format(mgr.genesdb, mgr.distance))
        return name

    def writeQsubClassify(self, mgr):
        name = "_scripts/classify.qsub"
        with BImisc.ShellScript(name) as out:
            out.write("""
IN=$1
OUT=$2
DB="{}"

module load bedtools
bedtools intersect -a $IN -b $DB -wa -wb | cut -f 9 | sort | uniq -c > $OUT
""".format(mgr.hmmdb))
        return name

### Class that generates a WashU hub configuration file

class WashUEntry():
    name = ""
    filename = ""
    color = ""
    group = ""
    ttype = "bigwig"

    def __init__(self, name, filename, color="#000000", group=1, ttype="bigwig"):
        self.name = name
        self.filename = filename
        self.color = color
        self.colornegative = '"#000000"'
        self.group = group
        self.ttype = ttype

    def write(self, hub, s):
        """Write this entry to stream s."""
        mode = "show"
        horiz = ""
        #fixedscale = ""

        if self.ttype == "bigbed":
            self.colorpositive = '"#0000FF"'
            mode = "thin" 
            horiz = '{"value": 0, "color": "#000000"}'
        elif self.ttype == "updown":
            self.ttype = "bigwig"
            horiz = '{"value": 0, "color": "#000000"}'
            self.color = '"#FF0000"'
            self.colornegative = '"#0000FF"'
            #fixedscale = """\n  "fixedscale": {"min": -1, "max": 1},"""
        s.write("""{{
  "type": "{}",
  "name": "{}",
  "height": 50,
  "group": {},
  "colorpositive": {},
  "colornegative": {},
  "horizontalLines": [{}],
  "mode": "{}",
  "url": "{}{}"
}}""".format(self.ttype, self.name, self.group, self.color, self.colornegative, horiz, mode, hub.baseurl, self.filename))

class WashUHub():
    entries = []
    colgroups = []
    palette = None
    baseurl = ""

    def __init__(self, filename=None):
        self.entries = []
        if filename:
            self.parseFromConf(filename)
            
    def parseFromConf(self, filename):
        group = 1
        colgroup = 1
        samecol = 1

        with open(filename, "r") as f:
            for line in f:
                line = line.strip()
                words = line.split("\t")
                if len(words) < 2:
                    continue
                w0 = words[0]
                w1 = words[1]
                if w0 == "group":
                    group = int(w1)
                elif w0 == "url":
                    self.baseurl = w1
                elif w0 == "samecol":
                    samecol = int(w1)
                else:
                    if len(words) > 2:
                        ttype = words[2]
                    else:
                        ttype  = "bigwig"
                    e = WashUEntry(w0, w1, group=group, ttype=ttype)
                    self.entries.append(e)
                    self.colgroups.append(colgroup)
                    samecol += -1
                    if samecol == 0:
                        colgroup += 1
                        samecol = 1

        if len(self.entries) > 27:
            self.palette = BIcolors.Palette(size=64)
        else:
            self.palette = BIcolors.Palette()
        c = 0
        color = ""
        for (e, cg) in zip(self.entries, self.colgroups):
            if cg != c:
                color = self.palette.nextColor()
                c = cg
            e.color = color

    def add(self, e):
        """Add entry E to this WashU hub."""
        self.entries.append(e)

    def write(self, s):
        """Write this hub to stream s."""
        comma = False
        s.write("[\n")
        for e in self.entries:
            if comma:
                s.write(",\n")
            else:
                comma = True
            e.write(self, s)
        s.write("]\n")

class Alignment(object):
    nreads = 0
    npeaks = 0
    totalopen = 0
    openfact = 1.0
    nrip = 0

    def getFactor(self, mgr):
        # fact = 1.0
        # if mgr.readsnorm:
        #     fact = fact * self.readsfact
        # if mgr.opennorm:
        #     fact = fact / self.openfact
        # return fact 

        #return 1.0 * self.totalopen / self.nreads
        return self.openfact

class Sample(BIfile.Filer, Alignment):
    name = ""
    parent = None
    readsfact = 1.0
    openfact = 1.0
    factor = 1.0

    def __init__(self, name):
        super(Sample, self).__init__("")
        self.name = name
        self.addFile("peaksfile", name + "_peaks.xls")
        self.addFile("bamfile", name + ".bam")
        self.addFile("bwfile", name + ".bw")
        self.addFile("bbfile", name + ".bb")
        # The following two files go into the condition directory
        self.directory = self.parent
        self.addFile("bedfile", name + ".macspeaks.bed")
        self.addFile("frip", name + ".frip")

    def FRIPcmdline(self, name):
        return " ".join([name, self.pathname("bedfile"), self.pathname("bamfile"), self.pathname("frip")])

    def readFRIP(self):
        with open(self.pathname("frip"), "r") as f:
            s = f.read()
            self.nrip = int(s)

class Condition(BIfile.Filer, Alignment):
    name = ""
    samples = []

    def __init__(self, name):
        super(Condition, self).__init__("")
        self.name = name
        self.directory = name
        self.samples = []
        self.addFile("bedfile", name + ".mergedpeaks.bed")
        self.addFile("bamfile", name + ".bam")
        self.addFile("bwfile", name + ".bw", dir=".") # This goes in the contrast
        self.addFile("bbfile", name + ".bb")

    def setAlignment(self):
        self.nreads = sum([s.nreads for s in self.samples])

    # def getFactor(self, mgr):
    #     fact = 0.0
    #     for s in self.samples:
    #         fact = fact + s.getFactor(mgr)
    #     return fact / len(self.samples)

    def mergePeaks(self, W, n):
        return W.writePeakMerger("_scripts/merge-{}.sh".format(n), self.pathname("bedfile"), self.samples)

    def convertBED(self):
        npeaks = 0
        totopen = 0
        regnum = 1
        with open(self.pathname("bedfile"), "w") as out:
            with open(self.pathname("bedfile") + ".pre", "r") as f:
                c = csv.reader(f, delimiter='\t')
                for line in c:
                    npeaks += 1
                    totopen += (int(line[2]) - int(line[1]))
                    out.write("{}\t{}\t{}\treg{}\t{}\t+\n".format(line[0], line[1], line[2], regnum, regnum))
                    regnum += 1
        return (npeaks, totopen)

class Contrast(BIfile.Filer):
    test =  None
    ctrl = None
    label = ""
    countsfiles = {}
    nup = 0
    ndown = 0

    def __init__(self, test, ctrl):
        super(Contrast, self).__init__("")
        self.countsfiles = {}
        self.test = test
        self.ctrl = ctrl
        self.label = "{}.vs.{}".format(test.name, ctrl.name)
        self.directory = self.label
        self.addFile("bedfile", self.label + ".commonpeaks.bed")
        self.addFile("counts", self.label + ".count.csv")
        self.addFile("testsizes", self.label + ".test-sizes.csv")
        self.addFile("ctrlsizes", self.label + ".ctrl-sizes.csv")
        self.addFile("sizes", self.label + ".sizes.csv")
        self.addFile("avgsizes", self.label + ".avgsizes.csv")
        self.addFile("matrix", self.label + ".matrix.csv")
        self.addFile("diff", self.label + ".diff.csv")
        self.addFile("diffx", self.label + ".diffpeaks.xlsx")
        self.addFile("diffbdg", self.label + ".diff.bedGraph")
        self.addFile("diffbw", self.label + ".diff.bw")
        self.addFile("sig", self.label + ".sig.csv")
        self.addFile("uptest", self.label + ".test-up.csv")
        self.addFile("upctrl", self.label + ".ctrl-up.csv")
        self.addFile("testgenes", self.label + ".test-genes.csv")
        self.addFile("ctrlgenes", self.label + ".ctrl-genes.csv")
        self.addFile("testregs", self.label + ".test-class.csv")
        self.addFile("ctrlregs", self.label + ".ctrl-class.csv")
        self.addFile("regcounts", self.label + ".regcounts.csv")
        self.addFile("annx", self.label + ".annot.xlsx")
        self.addFile("sizesplot", self.label + ".sizes")
        self.addFile("tssplot", self.label + ".tss")
        self.addFile("testplot", self.label + ".testpeaks")
        self.addFile("ctrlplot", self.label + ".ctrlpeaks")
        for smp in self.test.samples:
            self.countsfiles[smp.name] = self.addFile(smp.name + "_counts", "{}-counts.csv".format(smp.name))
        for smp in self.ctrl.samples:
            self.countsfiles[smp.name] = self.addFile(smp.name + "_counts", "{}-counts.csv".format(smp.name))

    def hubfiles(self, mgr):
        keys = ["testgenes", "ctrlgenes", "uptest", "upctrl"]
        if mgr.hmmdb:
            keys.append("regcounts")
        return [ self.pathname(k) for k in keys ]

    def plotfiles(self):
        keys = ["tssplot", "testplot", "ctrlplot", "sizesplot"]
        return [ self.pathname(k) + ".png" for k in keys ]

    def commonPeaks(self, n):
        """Write a script that computes the intersection of the peaks for the two conditions in this contrast."""
        scriptfile = "_scripts/intersect-{}.sh".format(n)
        with BImisc.ShellScript(scriptfile) as out:
            out.write("""

intersectBed -a {} -b {} -u | cut -f 1-3 | sort -k1,1 -k2,2n -k3,3n | uniq > {}
""".format(self.test.pathname("bedfile"), self.ctrl.pathname("bedfile"), self.pathname("bedfile")))
        return scriptfile

    def mergedPeaks(self, n):
        """Writes a script that computes the union (merge) of the peaks for the two conditions in this contrast."""
        scriptfile = "_scripts/cmerge-{}.sh".format(n)
        with BImisc.ShellScript(scriptfile) as out:
            out.write("""

sort -m -k1,1 -k2,2n -k3,3n {} {} | mergeBed -i - > {}
""".format(self.test.pathname("bedfile"), self.ctrl.pathname("bedfile"), self.pathname("bedfile")))
        return scriptfile

    def quantifyPeaks(self, scriptname):
        cmdlines = []
        for smp in self.test.samples:
            cmdlines.append(" ".join([scriptname, self.pathname("bedfile"), smp.pathname("bamfile"), self.countsfiles[smp.name].pathname()]))
        for smp in self.ctrl.samples:
            cmdlines.append(" ".join([scriptname, self.pathname("bedfile"), smp.pathname("bamfile"), self.countsfiles[smp.name].pathname()]))
        return cmdlines

    def getFactors(self, mgr):
        """Returns a list of two values, representing the normalization factors for the test and control conditions respectively (based on totopen)."""
        ft = self.test.totopen
        fc = self.ctrl.totopen
        mr = max(ft, fc)
        return [1.0 * ft / mr, 1.0 * fc / mr]

    def makeMatrix(self, mgr):
        samples = self.test.samples + self.ctrl.samples
        nsamples = len(samples)
        if mgr.opennorm:
            factors = [ smp.openfact for smp in samples ]
        else:
            factors = [ 1.0 for smp in samples ]
        streams = [ open(self.countsfiles[smp.name].pathname(), "r") for smp in samples ]
        with open(self.pathname("matrix"), "w") as out:
            out.write("\t".join([smp.name for smp in samples]) + "\n")
            readers = [ csv.reader(s, delimiter='\t') for s in streams ]
            try:
                for l1 in readers[0]:
                    # d = int(l1[2]) - int(l1[1]) # Do we need to normalize by peak size?
                    v = factors[0] * float(l1[3])
                    out.write("{}:{}-{}\t{}".format(l1[0], l1[1], l1[2], int(10.0 * v)))
                    for i in range(1, nsamples):
                        li = readers[i].next()
                        v = factors[i] * float(li[3])
                        out.write("\t{}".format(int(10.0 * v)))
                    out.write("\n")
            finally:
                for r in streams:
                    r.close()

    def writeDEGscript(self, scriptname):
        labels = ["1" for smp in self.test.samples] + ["2" for smp in self.ctrl.samples]
        with open(scriptname, "w") as out:
            out.write("""## Script to run DESeq2 on table of counts data.

library("DESeq2")

#args <- commandArgs(TRUE)
#datafile = args[1]
datafile = "{}"
counts = as.matrix(read.csv(datafile, sep='\t', row.names=1))
levels = c("{}", "{}")
labels = c({})
sampleTable = data.frame(condition=factor(levels[labels]))
rownames(sampleTable) = colnames(counts)
dds = DESeqDataSetFromMatrix(countData = counts, colData = sampleTable, design = ~ condition)
dds = estimateSizeFactors(dds)
keep = rowSums(counts(dds)) >= 5
dds = dds[keep,]
dds = DESeq(dds)
res = results(dds, contrast=c("condition", "{}", "{}"))
write.table(res, file="{}", sep='\t')
""".format(self.pathname("matrix"), self.test.name, self.ctrl.name,
           ",".join(labels), self.test.name, self.ctrl.name,
           self.pathname("diff")))

    def extractSignificant(self, mgr):
        up = []
        dn = []
        with open(self.pathname("sig"), "w") as out:
            with open(self.pathname("diffbdg"), "w") as bed:
                for line in BIcsv.CSVreader(self.pathname("diff"), skip=1):
                    try:
                        fc = float(line[2])
                        p  = float(line[6])
                    except ValueError: # Some P-values are NA
                        continue
                    (chrom, start, end) = parseCoords(line[0].strip('"'))
                    bed.write("{}\t{}\t{}\t{}\n".format(chrom, start, end, fc))
                    if abs(fc) >= mgr.log2fc and p <= mgr.pval:
                        if fc > 0:
                            self.nup += 1
                            up.append((chrom, start, end, fc, p))
                        else:
                            self.ndown += 1
                            dn.append((chrom, start, end, -fc, p))
                        out.write("{}\t{}\t{}\t{}\t{}\n".format(chrom, start, end, fc, p))

        up.sort(key=lambda s: s[3], reverse=True)
        dn.sort(key=lambda s: s[3], reverse=True)
        with open(self.pathname("uptest"), "w") as out:
            for rec in up:
                out.write("\t".join([str(s) for s in rec]) + "\n")

        with open(self.pathname("upctrl"), "w") as out:
            for rec in dn:
                out.write("\t".join([str(s) for s in rec]) + "\n")

    def annotateDifferential(self, script):
        cmdlines = []
        if BImisc.missingOrStale(self.pathname("testgenes"), self.pathname("uptest")):
            cmdlines.append(" ".join([script, self.pathname("uptest"), self.pathname("testgenes")]))
        if BImisc.missingOrStale(self.pathname("ctrlgenes"), self.pathname("upctrl")):
            cmdlines.append(" ".join([script, self.pathname("upctrl"), self.pathname("ctrlgenes")]))
        return cmdlines

    def classifyChromRegions(self, script):
        cmdlines = []
        if BImisc.missingOrStale(self.pathname("testregs"), self.pathname("uptest")):
            cmdlines.append(" ".join([script, self.pathname("uptest"), self.pathname("testregs")]))
        if BImisc.missingOrStale(self.pathname("ctrlregs"), self.pathname("upctrl")):
            cmdlines.append(" ".join([script, self.pathname("upctrl"), self.pathname("ctrlregs")]))
        return cmdlines

    def addCounts(self, data):
        pointers = [["Active Promoter", [1]],
                    ["Weak Promoter", [2]],
                    ["Poised Promoter", [3]],
                    ["Strong Enhancer", [4, 5]],
                    ["Weak Enhancer", [6, 7]],
                    ["Insulator", [8]],
                    ["Transcribed", [9, 10, 11]],
                    ["Repressed", [12]],
                    ["Heterochromatin/lo", [13]],
                    ["Repetitive/CNV", [14, 15]]]
        result = []
        for rec in pointers:
            sum = 0
            for idx in rec[1]:
                sum += data[idx]
            result.append([rec[0], sum])
        return result

    def parseClassifications(self):
        testdata = parseClassification(self.pathname("testregs"))
        testsum = float(sum(testdata.values()))
        testcounts = self.addCounts(testdata)
        ctrldata = parseClassification(self.pathname("ctrlregs"))
        ctrlsum = float(sum(ctrldata.values()))
        ctrlcounts = self.addCounts(ctrldata)

        with open(self.pathname("regcounts"), "w") as out:
            out.write("#Region\tPeaks in {}\t% in {}\tPeaks in {}\t% in {}\n".format(self.test.name, self.test.name, self.ctrl.name, self.ctrl.name))
            for (rec1, rec2) in zip(testcounts, ctrlcounts):
                out.write("{}\t{}\t{:.2f}\t{}\t{:.2f}\n".format(
                    rec1[0], rec1[1], 100.0 * rec1[1] / testsum,
                    rec2[1], 100.0 * rec2[1] / ctrlsum))

    def countClassified(self):
        a = 0
        b = 0
        with open(self.pathname("regcounts"), "r") as f:
            f.readline()
            c = csv.reader(f, delimiter='\t')
            for line in c:
                a += int(line[1])
                b += int(line[3])
        return (a, b)

### Manager
    
class Manager():
    samplesfile = None
    contrastfile = None
    steps = "12345678"
    experiment = None
    conditions = {}
    condnames = []
    samples = {}                # Sample names to samples
    samplenames = []
    contrasts = []
    dryrun = False

    # Params
    log2fc = 1                  # -fc
    pval = 0.01                 # -pval
    readsnorm = True            # -norm [R][O]
    opennorm = True
    mode = "I"                  # How to compare peaks in contrasts (M = merge, I = intersect)
    logfile = "dasa.log"        # -log
    genomever = "hg38"
    genesdb = "/ufrc/data/reference/icbr/GRCh38/Homo_sapiens.GRCh38.95.db"
    generegions = "/ufrc/data/reference/icbr/GRCh38/all-genes-top.csv"
    distance = 2000                         # -d
    hmmdb = None # "/ufrc/data/reference/icbr/GRCh38/wgEncodeBroadHmmHmecHMM.hg38.txt"

    # Hub
    hubname = "hub"
    huburl = "http://licht:licht@lichtlab.cancer.ufl.edu/reports/"
    chromsizes = "/apps/dibig_tools/manorm/hg38.chrom.sizes"
#"/ufrc/data/reference/icbr/mm10/ucsc/mm10.chrom.sizes"

    # Internal
    _meths = {}
    _opts = {}

    def __init__(self):
        self.conditions = {}
        self.condnames = []
        self.samples = {}
        self.samplenames = []
        self.contrasts = []
        self._meths = { "1": self.Step1,
                        "2": self.Step2,
                        "3": self.Step3,
                        "4": self.Step4,
                        "5": self.Step5,
                        "6": self.Step6,
                        "7": self.Step7,
                        "8": self.Step8 }
        self.steps = "".join(sorted(self._meths.keys()))
        self.initOptions()

    def error(self, msg, *args):
        sys.stdout.write("Error: " + msg.format(*args) + "\n")
        sys.exit(1)

    def validate(self):
        if not (0.0 <= self.pval <= 1.0):
            self.error("P-value (-p) should be between 0 and 1.")
        self.mode = self.mode.upper()
        if not self.mode in "MI":
            self.error("Mode (-m) should be one of M and I.")
        return True

    def initOptions(self):
        options = [(["-h", "--help"], self.usage),
                   (["-l", "--log"], lambda a: setattr(self, "logfile", a)),
                   (["-s", "--steps"], lambda a: setattr(self, "steps", a)),
                   (["-n", "--norm"], self.setNorm),
                   (["-m", "--mode"], lambda a: setattr(self, "mode", a)),
                   (["-f", "--log2fc"], lambda a: setattr(self, "log2fc", float(a))),
                   (["-p", "--pval"], lambda a: setattr(self, "pval", float(a))),
                   (["-gv", "--genomever"], lambda a: setattr(self, "genomever", a)),
                   (["-gdb", "--genesdb"], lambda a: setattr(self, "genesdb", a)),
                   (["-gr", "--generegions"], lambda a: setattr(self, "generegions", a)),
                   (["-gd", "--genedist"], lambda a: setattr(self, "distance", int(a))),
                   (["-mdb", "--hmmdb"], lambda a: setattr(self, "hmmdb", a)),
                   (["-hc", "--chromsizes"], lambda a: setattr(self, "chromsizes", a)),
                   (["-hn", "--hubname"], lambda a: setattr(self, "hubname", a)),
                   (["-hu", "--huburl"], lambda a: setattr(self, "huburl", a)),

               ]
        for opt in options:
            for k in opt[0]:
                self._opts[k] = opt[1]

    def setNorm(self, a):
        a = a.upper()
        self.readsnorm = ("R" in a)
        self.opennorm = ("O" in a)

    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev:
                f = self._opts[prev]
                f(a)
                prev = ""
            elif a in self._opts.keys():
                prev = a
            elif a == "-x":
                self.dryrun = True
            elif self.samplesfile is None:
                self.samplesfile = a
            else:
                self.contrastsfile = a
        if prev in ["-h", "--help"]:
            return self.usage()
        return True

    def usage(self, what=None):
        sys.stdout.write("""dasa.py - Differential ATAC-Seq Analysis

Usage: dasa.py [options] conditions contrasts

`Conditions' is a tab-delimited file with two columns containing a condition name and
a comma-delimited list of its samples, respectively. `Contrasts' is a tab-delimited
file with two columns, each line represent a contrast between the sample in the first column
(test) and the one in the second column (control). The program assumes that for each sample
S there will be a MACS output file called S.xls and a BAM file called S.bam with the 
associated .bai index in the current directory.

Options:
  -s S    | Execute the listed steps (default: {})
  -x      | Dry run - generate report, don't submit any jobs.
  -norm N | Perform read-number normalization if N contains 'R', open-regions
            normalization if N contains "O" (default: RO)
  -fc F   | Fold change threshold for significant peaks (default: {})
  -pval P | P-value threshold for significant peaks (default: {})
  -log L  | Write execution log to file L.

Steps:

  1. Compute number of peaks, scaling factors
  2. Compute unique/common peaks for each contrast
  3. Compute peak coverage
  4. Determine differential peaks
  5. Add gene annotations
  6. Generate track files
  7. Generate plots
  8. Generate final reports

""".format(self.steps, self.log2fc, self.pval))
        sys.exit(0)

    def openLog(self):
        self.logstream = open(self.logfile, "w")
        self.logstream.write(" ".join(sys.argv) + "\n")

    def closeLog(self):
        self.logstream.close()

    def log(self, msg, *args):
        s = msg.format(*args) + "\n"
        sys.stderr.write(s)
        self.logstream.write(s)

    def run(self):
        self.validate()
        self.openLog()
        try:
            self.initialize()
            W = Writer()
            BImisc.shell("rm -f *.done; mkdir -p _scripts {}".format(" ".join([ c.name for c in self.conditions.values() ] +
                                                                              [ c.label for c in self.contrasts ])))
            for j in self.steps:
                if j in self._meths:
                    m = self._meths[j]
                    m(W)
        finally:
            self.closeLog()
        sys.stderr.write(BItext.RED("=== Done. ===\n"))

    # def addSample(self, name):
    #     self.samplenames.append(name)
    #     s = Sample(name)
    #     self.samples[name] = s
    #     return s

    # def addCondition(self, name):
    #     c = Condition(name)
    #     self.conditions[name] = c
    #     return c
        
    def initialize(self):
        self.experiment = BIexperiment.Experiment()
        self.experiment.initConditionsFromFile(self.samplesfile)
        self.experiment.initContrastsFromFile(self.contrastsfile)
        self.condnames = self.experiment.conditions
        self.samplenames = self.experiment.samples

        for sname in self.samplenames:
            self.samples[sname] = Sample(sname)
        for cname in self.condnames:
            cond = Condition(cname)
            for sn in self.experiment.condsamples[cname]:
                self.samples[sn].parent = cname
                cond.samples.append(self.samples[sn])
            self.conditions[cname] = cond
        for [testname, ctrlname] in self.experiment.contrasts:
            self.contrasts.append(Contrast(self.conditions[testname], self.conditions[ctrlname]))

    def submit(self, cmdline, prefix="cov", mem="20G", ntasks=1, cpus=1, generic=True):
        w = "generic.qsub" if generic else ""
        subline = "submit -done {}.@.done -o --mem={},--ntasks={},--cpus-per-task={} {} {}".format(prefix, mem, ntasks, cpus, w, cmdline)
        sys.stderr.write("[Submitting: {}]\n".format(subline))
        return BImisc.shell(subline)

    def waitFor(self, nwanted, prefix="cov"):
        patt = prefix + ".*.done"
        sys.stderr.write("[Waiting for {} jobs to complete.]\n".format(nwanted))
        while True:
            found = glob.glob(patt)
            n = len(found)
            if n >= nwanted:
                for f in found:
                    os.remove(f)
                return
            else:
                time.sleep(30)

    def convertPeaks(self, peaksfile, bedfile):
        """Convert a MACS output file `peaksfile' to a BED file."""
        regnum = 1
        with open(bedfile, "w") as out:
            with open(peaksfile, "r") as f:

                # skip header
                for line in f:
                    if line != '\n' and line[0] != '#':
                        break

                tot = 0
                chrom = ""
                start = 0
                end = 0
                c = csv.reader(f, delimiter='\t')
                for line in c:
                    bchrom = line[0]
                    if "_" in bchrom: # get rid of weird chromosomes
                        continue

                    # New chromosome?
                    if bchrom != chrom:
                        if end > 0:
                            out.write("{}\t{}\t{}\treg{}\t{}\t+\n".format(chrom, start, end, regnum, regnum))
                            regnum += 1
                        chrom = bchrom
                        start = 0
                        end = 0

                    # Unwanted chromosome?
                    if bchrom == 'chrM' or "random" in bchrom:
                        start = 0
                        end = 0
                        continue

                    # Good line
                    bstart = int(line[1])
                    bend   = int(line[2])
                    if start <= bstart <= end:
                        # Extend current region
                        end = bend
                    else:
                        # Start new region
                        tot += (end - start)
                        if end > 0:
                            out.write("{}\t{}\t{}\treg{}\t{}\t+\n".format(chrom, start, end, regnum, regnum))
                            regnum += 1
                        start = bstart
                        end = bend

                out.write("{}\t{}\t{}\treg{}\t{}\t+\n".format(chrom, start, end, regnum, regnum))
                tot += (end - start)
        return (tot, regnum)

    def countReads(self, bamfile):
        """Counts the number of aligned reads in a BAM file."""
        nreads = 0
        data = BImisc.shell("""samtools idxstats {} | grep -v ^\*""".format(bamfile))
        for line in data.split("\n"):
            fields = line.split("\t")
            if len(fields) >= 3:
                nreads += int(fields[2])
        return nreads

    def WashUlink(self, json, text=None):
        if not text:
            text = os.path.split(json)[1]
        return '<A target="_blank" href="http://epigenomegateway.wustl.edu/browser/?genome={}&datahub={}/{}/{}">{}</A>'.format(
            self.genomever, self.huburl, self.hubname, json, text)

    def Step1(self, W):
        sys.stderr.write(BItext.RED("=== Step 1: Initial peak analysis ===\n"))
        samples = sorted(self.samples.values())

        self.log("Parsing MACS peaks...")
        for smp in samples:
            self.log("  {} => {}", smp.file("peaksfile"), smp.file("bedfile"))
            (smp.totalopen, smp.npeaks) = self.convertPeaks(smp.pathname("peaksfile"), smp.pathname("bedfile"))

        self.log("\nComputing number of reads...")
        for smp in samples:
            smp.nreads = self.countReads(smp.file("bamfile"))

        maxopen = max([smp.totalopen for smp in samples])
        maxreads = max([smp.nreads for smp in samples])
        for smp in samples:
            smp.openfact = 1.0 * smp.totalopen / maxopen
            smp.readsfact = 1.0 * maxreads / smp.nreads
            smp.factor = smp.getFactor(self)

        self.log("\nTotal open bases:")
        for smp in samples:
            self.log("  {}\t{}\t{}", smp.name, smp.totalopen, smp.openfact)

        self.log("\nNumber of reads:")
        for smp in samples:
            self.log("  {}\t{}\t{}", smp.name, smp.nreads, smp.readsfact)

        self.log("\nScaling factors:")
        for smp in samples:
            self.log("  {}\t{}", smp.name, smp.factor)
        self.log("\n")

    def Step2(self, W):
        sys.stderr.write(BItext.RED("=== Step 2: Peak merging ===\n"))
        n = 0
        njobs = 0
        self.log("\nCondition\tNumber of reads\tNumber of peaks\tOpen bases")
        for cond in self.conditions.values():
            cond.setAlignment()
            script = cond.mergePeaks(W, n)
            BImisc.shell(script)
            (cond.npeaks, cond.totopen) = cond.convertBED()
            cond.npeaks = cond.file("bedfile").nlines()
            cond.totopen = int(BImisc.shell("awk '{sum += ($3-$2)} END {print sum}' " + cond.pathname("bedfile")))
            self.log("{}\t{}\t{}\t{}".format(cond.name, cond.nreads, cond.npeaks, cond.totopen))
            for smp in cond.samples:
                self.log("  {}\t{}\t{:.1f}%".format(smp.name, smp.npeaks, 100.0 * smp.npeaks / cond.npeaks))
            self.log("")
            n += 1
            if not self.dryrun and BImisc.missingOrStale(cond.pathname("bamfile")):
                self.submit("samtools merge {} {}".format(cond.pathname("bamfile"), " ".join([smp.pathname("bamfile") for smp in cond.samples])), mem="2G")
                njobs += 1
        self.waitFor(njobs)
        njobs = 0
        for cond in self.conditions.values():
            if not self.dryrun and BImisc.missingOrStale(cond.pathname("bamfile") + ".bai", cond.pathname("bamfile")):
                self.submit("samtools index {}".format(cond.pathname("bamfile")), mem="1G")
                njobs += 1
        self.waitFor(njobs)

        n = 0
        for contr in self.contrasts:
            if self.mode == "I":
                script = contr.commonPeaks(n)
            else:
                script = contr.mergePeaks(n)
            BImisc.shell(script)
            n += 1
        self.log("\nCommon peaks in each contrast:")
        for contr in self.contrasts:
            self.log("  {}\t{}", contr.label, contr.file("bedfile").nlines())
        self.log("\nAverage size of common peaks in each contrast:")
        for contr in self.contrasts:
            # Maybe these shells should become a single script?
            BImisc.shell("intersectBed -a {} -b {} -wa |awk -v OFS='\t' '{{print $1,$2,$3,$3-$2}}' - > {}".format(
                contr.test.pathname("bedfile"), contr.ctrl.pathname("bedfile"), contr.pathname("testsizes")))
            BImisc.shell("intersectBed -a {} -b {} -wa |awk -v OFS='\t' '{{print $1,$2,$3,$3-$2}}' - > {}".format(
                contr.ctrl.pathname("bedfile"), contr.test.pathname("bedfile"), contr.pathname("ctrlsizes")))
            BImisc.shell("paste {} {} > {}".format(contr.pathname("testsizes"), contr.pathname("ctrlsizes"), contr.pathname("sizes")))
            BImisc.shell("awk -v OFS='\t' '{{sum1+=$4; sum2+=$8}} END {{print sum1/NR,sum2/NR}}' {} > {}".format(contr.pathname("sizes"), contr.pathname("avgsizes")))
            with open(contr.pathname("avgsizes"), "r") as f:
                line = f.readline().strip().split("\t")
                self.log("  {}\t{}\t{}", contr.label, line[0], line[1])

    def Step3(self, W):
        sys.stderr.write(BItext.RED("=== Step 3: Peak quantification ===\n"))
        samples = self.samples.values()
        njobs = 0

        # Submit jobs for FRIP computation
        scriptname = W.writeFRIPQsub()
        for smp in samples:
            cmdline = smp.FRIPcmdline(scriptname)
            if not self.dryrun:
                self.submit(cmdline, mem="10G")
                njobs += 1

        # Submit jobs for quantification
        scriptname = "_scripts/quantify.sh"
        with BImisc.ShellScript(scriptname) as out:
            out.write("""
samtools bedcov $1 $2 > $3
""")                                      # .format(P.overlap))
        for contr in self.contrasts:
            cmdlines = contr.quantifyPeaks(scriptname)
            for cmd in cmdlines:
                if not self.dryrun:
                    self.submit(cmd)
                    njobs += 1

        # Wait for everything to be finished
        self.waitFor(njobs)

        self.log("\nFraction of reads in peaks:")
        for smp in samples:
            smp.readFRIP()
            self.log("  {}\t{}\t{}\t{:.2f}".format(smp.name, smp.nrip, smp.nreads, 100.0 * smp.nrip / smp.nreads))

        for contr in self.contrasts:
            contr.makeMatrix(self)
        self.log("\n")

    def Step4(self, W):
        sys.stderr.write(BItext.RED("=== Step 4: Differential analysis ===\n"))
        njobs = 0
        for contr in self.contrasts:
            scriptname = "_scripts/rundeseq-{}.R".format(njobs)
            contr.writeDEGscript(scriptname)
            if not self.dryrun:
                self.submit("Rscript " + scriptname + " module:R/3.5.1")
                njobs += 1
        self.waitFor(njobs)
        self.log("\nSignificantly different peaks in each contrast:")
        for contr in self.contrasts:
            contr.extractSignificant(self)
            self.log("  {}\t{} up \t{} down", contr.label, contr.nup, contr.ndown)
        self.log("\n")

    def Step5(self, W):
        sys.stderr.write(BItext.RED("=== Step 5: Annotation of differential peaks ===\n"))
        """Annotate common and unique peaks."""
        self.log("*** Annotating differential peaks ***")
        name = W.writeQsubAnnotate(self)
        if self.hmmdb:
            name2 = W.writeQsubClassify(self)
        njobs = 0
        if not self.dryrun:
            for contr in self.contrasts:
                cmdlines = contr.annotateDifferential(name)
                for cmd in cmdlines:
                    self.submit(cmd, mem="5G")
                    njobs += 1
                if self.hmmdb:
                    cmdlines = contr.classifyChromRegions(name2)
                    for cmd in cmdlines:
                        self.submit(cmd, mem="5G")
                        njobs += 1
        self.waitFor(njobs)
        if self.hmmdb:
            for contr in self.contrasts:
                contr.parseClassifications()
        self.log("\n")

    def Step6(self, W):
        sys.stderr.write(BItext.RED("=== Step 6: Generation of genome browser track files ===\n"))
        """Generate bigwig and bigbed files."""
        BImisc.shell("mkdir -p {} {}/peaks {} {}".format(self.hubname, self.hubname, " ".join([ self.hubname + "/" + contr.label for contr in self.contrasts]), " ".join([ self.hubname + "/" + cond.name for cond in self.conditions.values()])))

        n = 0
        for smp in self.samples.values():
            if not self.dryrun:
                dest = self.hubname + "/peaks/" + smp.pathname("bwfile")
                if BImisc.missingOrStale(dest, smp.pathname("bamfile")):
                    self.submit("bamtowig.qsub {} {} deep=Y scale={}".format(smp.file("bamfile"), dest, 100.0 * smp.openfact), generic=False) # we only normalize by openfactor, because deeptools already does CPM normalization
                    n += 1
        for cond in self.conditions.values():
            if not self.dryrun:
                # dest = self.hubname + "/" + cond.pathname("bwfile")
                # if BImisc.missingOrStale(dest, cond.pathname("bamfile")):
                #     self.submit("bamtowig.qsub {} {} deep=Y scale={}".format(cond.pathname("bamfile"), dest, 100.0 * cond.getFactor(self)), generic=False)
                #     n += 1
                dest = self.hubname + "/" + cond.pathname("bbfile")
                if BImisc.missingOrStale(dest, cond.pathname("bedfile")):
                    BImisc.shell("bedToBigBed -type=bed3+3 -tab {} {} {}".format(cond.file("bedfile"), self.chromsizes, dest))

        for contr in self.contrasts:
            f = contr.getFactors(self)  # factors based on totopen
            dest = "{}/{}/{}".format(self.hubname, contr.label, contr.test.pathname("bwfile"))
            if BImisc.missingOrStale(dest, contr.test.pathname("bamfile")):
                self.submit("bamtowig.qsub {} {} deep=Y scale={}".format(contr.test.pathname("bamfile"), dest, 100.0 * f[0]), generic=False)
                n += 1
            dest = "{}/{}/{}".format(self.hubname, contr.label, contr.ctrl.pathname("bwfile"))
            if BImisc.missingOrStale(dest, contr.ctrl.pathname("bamfile")):
                self.submit("bamtowig.qsub {} {} deep=Y scale={}".format(contr.ctrl.pathname("bamfile"), dest, 100.0 * f[1]), generic=False)
                n += 1
            dest =  "{}/{}".format(self.hubname, contr.pathname("diffbw"))
            if BImisc.missingOrStale(dest, contr.pathname("diffbdg")):
                BImisc.shell("bedGraphToBigWig {} {} {}".format(contr.pathname("diffbdg"), self.chromsizes, dest))

        self.waitFor(n)
        for smp in self.samples.values():
            dest = self.hubname + "/peaks/" + smp.pathname("bbfile")
            if BImisc.missingOrStale(dest, smp.pathname("bedfile")):
                BImisc.shell("bedToBigBed -type=bed3+3 -tab {} {} {}".format(smp.file("bedfile"), self.chromsizes, dest))

        # Generate JSON configuration file for hub
        for contr in self.contrasts:
            BImisc.shell("cp {} {}/{}/".format(" ".join(contr.hubfiles(self)), self.hubname, contr.label))

        # Write track configuration file and index.html
        with open(self.hubname + "/peaks.conf", "w") as out:
            out.write("url\t{}/{}/\n".format(self.huburl, self.hubname))
            for (smpname, smp) in self.samples.iteritems():
                out.write("\nsamecol\t2\n")
                out.write("group\t1\n")
                out.write("{}\t{}\n".format(smp.name, "peaks/" + smp.pathname('bwfile')))
                out.write("group\t2\n")
                out.write("{}_peaks\t{}\tbigbed\n".format(smp.name, "peaks/" + smp.pathname('bbfile')))
        H = WashUHub(self.hubname + "/peaks.conf")
        with open(self.hubname + "/peaks.json", "w") as out:
            H.write(out)

        # Now write tracks for contrasts
        for contr in self.contrasts:
            cond1 = contr.test
            cond2 = contr.ctrl
            contr.hubconf = "{}/{}.conf".format(self.hubname, contr.label)
            contr.hubjson = "{}/{}.json".format(self.hubname, contr.label)

            with open(contr.hubconf, "w") as out:
                out.write("url\t{}/{}/\n".format(self.huburl, self.hubname))
                out.write("group\t1\n")
                out.write("{}\t{}/{}\n".format(cond1.name, contr.label, cond1.pathname('bwfile')))
                out.write("{}_peaks\t{}\tbigbed\n".format(cond1.name, cond1.pathname('bbfile')))
                out.write("{}\t{}/{}\n".format(cond2.name, contr.label, cond2.pathname('bwfile')))
                out.write("{}_peaks\t{}\tbigbed\n".format(cond2.name, cond2.pathname('bbfile')))
                out.write("group\t2\n")
                out.write("{}\t{}\tupdown\n".format(contr.label, contr.pathname("diffbw")))
            H2 = WashUHub(contr.hubconf)
            sys.stderr.write("{} => {}\n".format(contr.hubconf, contr.hubjson))
            with open(contr.hubjson, "w") as out:
                H2.write(out)

        self.log("\n")

    def Step7(self, W):
        sys.stderr.write(BItext.RED("=== Step 7: Plotting ===\n"))
        """Generate plots."""
        nplot = 0
        for contr in self.contrasts:
            if not self.dryrun:
                dest = contr.pathname("sizesplot") + ".png"
                if BImisc.missingOrStale(dest, contr.pathname("sizes")):
                    self.submit("density_scatterplot.py -cx 8 -cy 4 -l10 -yl {} -xl {} {} {}".format
                                (contr.test.name, contr.ctrl.name, contr.pathname("sizes"), dest),
                                mem="2G")
                    nplot += 1
                testbw = "{}/{}/{}".format(self.hubname, contr.label, contr.test.pathname('bwfile'))
                ctrlbw = "{}/{}/{}".format(self.hubname, contr.label, contr.ctrl.pathname('bwfile'))
                self.submit("volcano.sh tss {} {} {} {} ".format(contr.pathname("tssplot"), self.generegions, testbw, ctrlbw), cpus=8)
                self.submit("volcano.sh regions {} {} {} {} ".format(contr.pathname("testplot"), contr.pathname("uptest"), testbw, ctrlbw), cpus=8)
                self.submit("volcano.sh -s 2 regions {} {} {} {} ".format(contr.pathname("ctrlplot"), contr.pathname("upctrl"), testbw, ctrlbw), cpus=8)
                nplot += 3
        self.waitFor(nplot)

        for contr in self.contrasts:
            BImisc.shell("cp {} {}/{}/".format(" ".join(contr.plotfiles()), self.hubname, contr.label))

        self.log("\n")

    def Step8(self, W, dry=False):
        """Generate track hub."""
        
        # for contr in self.contrasts:
        #     shell("cp {} {}/{}/".format(" ".join(contr.hubfiles()), P.hubname, contr.dirname))

        # # Write track configuration file and index.html
        # with open(P.hubname + "/peaks.conf", "w") as out:
        #     out.write("url\t{}/{}/\n".format(P.huburl, P.hubname))
        #     for smpname in self.samplenames:
        #         smp = self.samples[smpname]
        #         out.write("\nsamecol\t2\n")
        #         out.write("group\t1\n")
        #         out.write("{}\t{}\n".format(smp.name, "peaks/" + smp.pathname('bwfile')))
        #         out.write("group\t2\n")
        #         out.write("{}_peaks\t{}\tbigbed\n".format(smp.name, "peaks/" + smp.pathname('bbfile')))
        # shell("{} {} > {}".format(P.mkhub, P.hubname + "/peaks.conf", P.hubname + "/peaks.json"))

        # # Now write tracks for contrasts
        # for contr in self.contrasts:
        #     smp1 = self.samples[contr.sample1]
        #     smp2 = self.samples[contr.sample2]
        #     contr.hubconf = "{}/{}.vs.{}.conf".format(P.hubname, contr.sample1, contr.sample2)
        #     contr.hubjson = "{}/{}.vs.{}.json".format(P.hubname, contr.sample1, contr.sample2)

        #     with open(contr.hubconf, "w") as out:
        #         out.write("url\t{}/{}/\n".format(P.huburl, P.hubname))
        #         out.write("group\t1\n")
        #         out.write("{}\t{}\n".format(smp1.name, "peaks/" + smp1.pathname('bwfile')))
        #         out.write("{}_peaks\t{}\tbigbed\n".format(smp1.name, "peaks/" + smp1.pathname('bbfile')))
        #         out.write("{}\t{}\n".format(smp2.name, "peaks/" + smp2.pathname('bwfile')))
        #         out.write("{}_peaks\t{}\tbigbed\n".format(smp2.name, "peaks/" + smp2.pathname('bbfile')))
        #         out.write("group\t2\n")
        #         src = contr.unique1diff
        #         dst = P.hubname + "/" + contr.dirname + contr.sample1 + ".diffPeaks.bw"
        #         nl = nlines(src)
        #         sys.stderr.write("{}: {} lines\n".format(src, nl))
        #         if nl > 0:
        #             shell("{}/bedGraphToBigWig {} {} {}".format(P.ucscdir, src, P.genome, dst))
        #             out.write(contr.sample1 + "_u\t" + contr.dirname + contr.sample1 + ".diffPeaks.bw" + "\n")
        #         src = contr.unique2diff
        #         dst = P.hubname + "/" + contr.dirname + contr.sample2 + ".diffPeaks.bw"
        #         nl = nlines(src)
        #         sys.stderr.write("{}: {} lines\n".format(src, nl))
        #         if nl > 0:
        #             shell("{}/bedGraphToBigWig {} {} {}".format(P.ucscdir, src, P.genome, dst))
        #             out.write(contr.sample2 + "_u\t" + contr.dirname + contr.sample2 + ".diffPeaks.bw" + "\n")
        #         src = contr.commondiff
        #         dst = "{}/{}/{}.vs.{}.diffPeaks.bw".format(P.hubname, contr.dirname, contr.sample1, contr.sample2)
        #         nl = nlines(src)
        #         sys.stderr.write("{}: {} lines\n".format(src, nl))
        #         if nl > 0:
        #             shell("{}/bedGraphToBigWig {} {} {}".format(P.ucscdir, src, P.genome, dst))
        #             out.write("Common\t" + contr.dirname + "{}.vs.{}.diffPeaks.bw".format(contr.sample1, contr.sample2))
        #     shell("{} {} > {}".format(P.mkhub, contr.hubconf, contr.hubjson))

        sys.stderr.write(BItext.RED("=== Step 8: Generating final report ===\n"))
        self.log("Writing index.html")
        with open(self.hubname + "/index.html", "w") as out:
            out.write("""<!DOCTYPE html>
<html>
<head>
<style>
TABLE {{width: 80%; border-collapse: collapse; border: 2px solid black;}}
TD, TH {{border: 1px solid grey}}
</style>
<body>
<center><h1>{} - Differential ATAC-Seq Analysis</h1>
<table>
<caption>Table 1. Number of aligned reads and open bases in each sample.</caption>
<tr><th>Sample</th><th>Aligned reads</th><th>Open bases</th><th>Factor</th><th>Peaks</th><th>FRIP</th></tr>
""".format(self.hubname))

            for smpname in self.samplenames:
                smp = self.samples[smpname]
                smp.readFRIP()
                out.write("""<tr><th>{}</th><td align='right'>{:,}</td><td align='right'>{:,}</td><td align='right'>{:.3f}</td><td align='right'>{:,}</td><td align='right'>{:.1f}%</td></tr>
                """.format(smp.name, smp.nreads, smp.totalopen, smp.getFactor(self), smp.npeaks, 100.0 * smp.nrip / smp.nreads))
            out.write("""</table>
<br>
Genome browser track: {}
<br><br>
""".format(self.WashUlink("peaks.json", "Peak calling ({} samples)".format(len(self.samples)))))

            out.write("""<table>
<caption>Table 2. Number of peaks and average peak sizes in each contrast.</caption>
<tr><th>Test</th><th>Ctrl</th><th>Number of peaks</th><th>Avg size in Test</th><th>Avg size in Ctrl</th></tr>
""")
            for contr in self.contrasts:
                npeaks = contr.file("bedfile").nlines()
                with open(contr.pathname("avgsizes"), "r") as f:
                    (tsize, csize) = f.readline().strip().split("\t")
                out.write("<TR><TH>{}</TH><TH>{}</TH><TD align='right'>{:,}</TD><TD align='right'>{:.1f}bp</TD><TD align='right'>{:.1f}bp</TD></TR>\n".format(contr.test.name, contr.ctrl.name, npeaks, float(tsize), float(csize)))
            out.write("""</table>\n<br>\n""")

            out.write("""<table>
<caption>Table 3. Differential peaks in each contrast.</caption>
<tr><th>Test</th><th>Ctrl</th><th>Significant</th><th>Up in Test</th><th>Up in Ctrl</th><th>Full results</th><th>WashU hub</th></tr>
""")
            for contr in self.contrasts:
                dest = self.hubname + "/" + contr.pathname("diffx")
                if BImisc.missingOrStale(dest, contr.pathname("sig")):
                    BImisc.shell("csvtoxls.py -q {} {} -name Significant {} -name Up-{} {} -name Up-{}".format(
                        dest, contr.pathname("sig"), contr.pathname("uptest"), contr.test.name, contr.pathname("upctrl"), contr.ctrl.name))
                nsig = contr.file("sig").nlines()
                nup = contr.file("uptest").nlines()
                ndown = contr.file("upctrl").nlines()
                json = contr.label + ".json"
                out.write("<TR><TH>{}</TH><TH>{}</TH><TD align='right'>{:,}</TD><TD align='right'>{:,}</TD><TD align='right'>{:,}</TD><TD align='center'>{}</TD><TD align='center'>{}</TD></TR>\n".format(
                    contr.test.name, contr.ctrl.name, nsig, nup, ndown, BImisc.linkify(contr.pathname("diffx"), None), 
                    self.WashUlink(json)))
            out.write("""</table>\n<br>\n""")

            out.write("""<table>
<caption>Table 4. Annotation of differential peaks in each contrast.</caption>
<tr><th>Test</th><th>Ctrl</th><th>Genes - Test</th><th>Genes - Ctrl</th><th>Classified - Test</TH><TH>Classified - Ctrl</TH><th>Full results</th></tr>
""")
            for contr in self.contrasts:
                dest = self.hubname + "/" + contr.pathname("annx")
                if BImisc.missingOrStale(dest, contr.pathname("testgenes")):
                    cmdline = "csvtoxls.py -q {} {} -name Genes-{} {} -name Genes-{}".format(
                        dest, contr.pathname("testgenes"), contr.test.name, contr.pathname("ctrlgenes"), contr.ctrl.name)
                    if self.hmmdb:
                        cmdline += "{} -name Classification".format(contr.addFile("regcounts"))
                    BImisc.shell(cmdline)
                if self.hmmdb:
                    (tregs, cregs) = contr.countClassified()
                else:
                    tregs = 0
                    cregs = 0

                tgenes = contr.file("testgenes").nlines()
                cgenes = contr.file("ctrlgenes").nlines()
                out.write("<TR><TH>{}</TH><TH>{}</TH><TD align='right'>{:,}</TD><TD align='right'>{:,}</TD><TD align='right'>{:,}</TD><TD align='right'>{:,}</TD><TD align='center'>{}</TD></TR>\n".format(
                    contr.test.name, contr.ctrl.name, tgenes, cgenes, tregs, cregs, BImisc.linkify(contr.pathname("annx"), None)))
            out.write("""</table>\n<br>\n""")

            out.write("""<table>
<caption>Table 5. Plots based on contrast data.</caption>
<tr><th>Test</th><th>Ctrl</th><th>Peak sizes</th><th>TSS</th><th>Test Peaks</TH><TH>Ctrl Peaks</TH></tr>
""")
            for contr in self.contrasts:
                out.write("<TR><TH>{}</TH><TH>{}</TH><TD align='center'>{}</TD><TD align='center'>{}</TD><TD align='center'>{}</TD><TD align='center'>{}</TD></TR>".format(
                    contr.test.name, contr.ctrl.name, 
                    BImisc.linkify(contr.pathname("sizesplot") + ".png", None), 
                    BImisc.linkify(contr.pathname("tssplot") + ".png", None),
                    BImisc.linkify(contr.pathname("testplot") + ".png", None), 
                    BImisc.linkify(contr.pathname("ctrlplot") + ".png", None)))
            out.write("""</table>\n<br>\n""")

            out.write("""
</center>
</body>
</html>
""")
        # Finally, zip everything
        # print "zipping: zip -r {}.zip {}".format(P.hubname, P.hubname)
        self.log("Creating zip file")
        BImisc.shell("zip -v -r {}.zip {}/".format(self.hubname, self.hubname))
        
### DUMMY
#     def dummy(self):
#         with open("dummy", "w") as out:
#             out.write("""<table>
# <caption>Table 2. Number of unique and common peaks in each contrast.</caption>
# <tr><th rowspan='2'>Test</th><th rowspan='2'>Ctrl</th><th colspan='2'>Unique in Test</th><th>Common</th><th colspan='2'>Unique in Ctrl</th></tr>
# <tr><th>Number of peaks</th><th>Avg peak size</th><th>Number of peaks</th><th>Number of peaks</th><th>Avg peak size</th></tr>
# """)
#             for contr in self.contrasts:
#                 out.write("""<tr><th>{}</th><th>{}</th><td align='right'>{:,}</td><td align='right'>{:.1f} bp</td><td align='right'>{:,}</td><td align='right'>{:,}</td><td align='right'>{:.1f} bp</td></tr>
#                 """.format(contr.sample1, contr.sample2, contr.nunique1, contr.avgsize1, contr.ncommon, contr.nunique2, contr.avgsize2))
#             out.write("""</table>
# <br>
# <table>
# <caption>Table3. Number of significantly increased/decreased peaks in each contrast.</caption>
# <tr><th>Test</th><th>Ctrl</th><th>Unique in Test</th><th>Common</th><th>Unique in Ctrl</th></tr>
# """)
#             for contr in self.contrasts:
#                 cnt = contr.counts
#                 out.write("""<tr><th>{}</th><th>{}</th><td align='right'>{:,} / {:,}</td><td align='right'>{:,} / {:,}</td><td align='right'>{:,} / {:,}</td></tr>
# """.format(contr.sample1, contr.sample2, cnt['u1up'], cnt['u1dn'], cnt['coup'], cnt['codn'], cnt['u2up'], cnt['u2dn']))
#             out.write("""</table>
# <br>
# <table>
# <caption>Table 4. Genes associated with significantly different peaks.</caption>
# <tr><th rowspan='2'>Test</th><th rowspan='2'>Ctrl</th><th colspan='2'>Unique in Test</th><th colspan='2'>Common</th><th colspan='2'>Unique in Ctrl</th></tr>
# <tr><th>Peaks</th><th>Genes</th><th>Peaks</th><th>Genes</th><th>Peaks</th><th>Genes</th></tr>
# """)
#             for contr in self.contrasts:
#                 out.write("""<tr><th>{}</th><th>{}</th>
#   <td><A href="{}">Peaks</a></td><td><A href="{}">up</A> ({}), <A href="{}">down</A> ({})</td>
#   <td><A href="{}">Peaks</a></td><td><A href="{}">up</A> ({}), <A href="{}">down</A> ({})</td>
#   <td><A href="{}">Peaks</a></td><td><A href="{}">up</A> ({}), <A href="{}">down</A> ({})</td>
# </tr>
# """.format(contr.sample1, contr.sample2, 
#            contr.unique1genesx, contr.unique1up, nlines(contr.unique1up), contr.unique1dn, nlines(contr.unique1dn),
#            contr.commongenesx,  contr.commonup, nlines(contr.commonup), contr.commondn, nlines(contr.commondn),
#            contr.unique2genesx, contr.unique2up, nlines(contr.unique2up), contr.unique2dn, nlines(contr.unique2dn)))
#             out.write("""</table>
# <br>
# <table><caption>Table 5. Genome browser tracks for differential peak analysis.</caption>
# <tr><th>Test</th><th>Ctrl</th><th>Tracks</th></tr>
# """)
#             for contr in self.contrasts:
#                 out.write("""<tr><th>{}</th><th>{}</th>
#   <td><A target="_blank" href="http://epigenomegateway.wustl.edu/browser/?genome={}&datahub={}/{}">{} vs {} differential peaks</A></td></tr>
# """.format(contr.sample1, contr.sample2, P.genomever, P.huburl, contr.hubjson, contr.sample1, contr.sample2))

#             out.write("""</table>
# <br>
# <p>Table 6. Heatmaps and volcano plots.</p>
# """)
#             for contr in self.contrasts:
#                 out.write("""<table>
# <tr><th colspan='2'>{} vs {}</th></tr>
# <tr><th>Common peak sizes</th><td>{}</td></tr>
# <tr><th>TSS volcano plot</th><td>{}</td></tr>
# <tr><th>{} unique peaks volcano plot</th><td>{}</td></tr>
# <tr><th>{} unique peaks volcano plot</th><td>{}</td></tr>
# <tr><th>Peak classification</th><td>{}</td></tr>
# </table><br><br>
# """.format(contr.sample1, contr.sample2, 
#            linkify(contr.commonhmap, "Common peak sizes"), 
#            linkify(contr.tssvolcano + ".png", "TSS volcano plot"),
#            contr.sample1,
#            linkify(contr.peaks1volcano + ".png", "Test peaks volcano plot"),
#            contr.sample2,
#            linkify(contr.peaks2volcano + ".png", "Ctrl peaks volcano plot"),
#            linkify(contr.allClassifyx, "Peak classification")))
#             out.write("""
# </center>
# </body>
# </html>
# """)

#         # Finally, zip everything
#         # print "zipping: zip -r {}.zip {}".format(P.hubname, P.hubname)
#         shell("zip -r {}.zip {}/".format(P.hubname, P.hubname))


### Top-level

def verifyEnv():
    try:
        what = "samtools"
        BImisc.shell("samtools --version")
        what = "bedtools"
        BImisc.shell("bedtools --version")
        return True
    except:
        sys.stderr.write("This program requires {}, but it was not found in $PATH.\n".format(what))
        return False

if __name__ == "__main__":
    if verifyEnv():
        M = Manager()
        if M.parseArgs(sys.argv[1:]):
            M.run()

