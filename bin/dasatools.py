#!/usr/bin/env python

import sys
import csv
import shlex
import os.path
import xlsxwriter
from collections import defaultdict

PYVER = sys.version_info.major

# Utils

def CSVreader(stream):
    return csv.reader(stream, delimiter='\t')

def Next(reader):
    if PYVER == 2:
        return reader.next()
    else:
        return reader.__next__()

def convertPeaks(peaksfile, bedfile):
    """Convert a MACS output file `peaksfile' to a BED file. Also works if the input is already in BED format."""
    regnum = 1
    with open(bedfile, "w") as out:
        with open(peaksfile, "r") as f:
            tot = 0
            chrom = ""
            start = 0
            end = 0
            c = CSVreader(f)
            for line in c:
                
                if len(line) == 0 or line[0][0] == '#' or line[0] == 'chr':
                    continue

                bchrom = line[0]
                if "_" in bchrom: # get rid of weird chromosomes
                    continue

                # New chromosome?
                if bchrom != chrom:
                    if end > 0:
                        out.write("{}\t{}\t{}\treg{}\t{}\t+\n".format(chrom, start, end, regnum, regnum))
                        #out.write("{}\t{}\t{}\tid:{}\n".format(chrom, start, end, regnum))
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
                        #out.write("{}\t{}\t{}\tid:{}\n".format(chrom, start, end, regnum))
                        out.write("{}\t{}\t{}\treg{}\t{}\t+\n".format(chrom, start, end, regnum, regnum))
                        regnum += 1
                    start = bstart
                    end = bend

            out.write("{}\t{}\t{}\treg{}\t{}\t+\n".format(chrom, start, end, regnum, regnum))
            #out.write("{}\t{}\t{}\tid:{}\n".format(chrom, start, end, regnum))
            tot += (end - start)
    return (tot, regnum)

def computeStats(args):
    maxreads = 0
    maxopen = 0

    data = []
    condnreads = defaultdict(int)
    outfile = args[0]
    cndfile = args[1]

    for filename in args[2:]:
        with open(filename, "r") as f:
            row = f.readline().strip().split("\t")
        smpdata = { 'sample': row[1], 'cond': row[0], 'totalopen': int(row[5]), 'npeaks': int(row[6]), 'nreads': int(row[7]) }
        condnreads[smpdata['cond']] += smpdata['nreads']
        if smpdata['nreads'] > maxreads:
            maxreads = smpdata['nreads']
            maxopen = smpdata['totalopen']
        data.append(smpdata)
    data.sort(key=lambda r: r['cond'] + ":" + r['sample'])
    fields = ['sample', 'cond', 'totalopen', 'npeaks', 'nreads', 'readsfact', 'openfact', 'factor']

    # Write stats for each sample
    with open(outfile, "w") as out:
        out.write("\t".join(fields) + "\n")
        for smpdata in data:
            smpdata['readsfact'] = 1.0 * maxreads / smpdata['nreads']
            smpdata['openfact']  = 1.0 * maxopen / smpdata['totalopen']
            smpdata['factor'] = smpdata['readsfact'] / smpdata['openfact']
            out.write("\t".join([str(smpdata[fld]) for fld in fields]) + "\n")

    # Write stats for each condition (only nreads for now)
    with open(cndfile, "w") as out:
        out.write("cond\tnreads\n")
        for cond in condnreads.keys():
            out.write("{}\t{}\n".format(cond, condnreads[cond]))

def writeMatrixOld(names, countsarg, factorsarg):
    countfiles = countsarg[1:].split(",")
    factors = [float (f) for f in factorsarg.split(",")]
    nsamples = len(countfiles)

    sys.stdout.write(names.replace(",", "\t") + "\n")
    streams = [ open(cf, "r") for cf in countfiles]
    readers = [ CSVreader(s) for s in streams ]
    try:
        for l1 in readers[0]:
            v = float(l1[3]) / factors[0]
            sys.stdout.write("{}:{}-{}\t{}".format(l1[0], l1[1], l1[2], int(10.0 * v)))
            for i in range(1, nsamples):
                li = Next(readers[i])
                v = float(li[3]) / factors[0]
                sys.stdout.write("\t{}".format(int(10.0 * v)))
            sys.stdout.write("\n")
    finally:
        for r in streams:
            r.close()

def readFactors(factorsfile):
    factors = {}
    with open(factorsfile, "r") as f:
        c = CSVreader(f)
        for line in c:
            factors[line[0]] = float(line[3])
    return factors

def writeMatrix(testCond, conditions, samples, factorsfile, countfiles):
    factors = readFactors(factorsfile)
    conditions = conditions.split(",")
    samples = samples.split(",")
    testnames = []
    testfiles = []
    ctrlnames = []
    ctrlfiles = []
    nsamples = len(conditions)
    for i in range(nsamples):
        cond = conditions[i]
        smp = samples[i]
        count = countfiles[i]
        if cond == testCond:
            testnames.append(smp)
            testfiles.append(count)
        else:
            ctrlnames.append(smp)
            ctrlfiles.append(count)

    ls = ["1"]*len(testnames) + ["2"]*len(ctrlnames)
    with open("labels.txt", "w") as out:
        out.write(",".join(ls) + "\n")

    names = testnames + ctrlnames
    files = testfiles + ctrlfiles

    sys.stdout.write("\t".join(names) + "\n")
    streams = [ open(cf, "r") for cf in files]
    readers = [ CSVreader(s) for s in streams ]
    facts   = [ factors[s] for s in names ]
    try:
        for l1 in readers[0]:
            v = float(l1[3]) / facts[0]
            sys.stdout.write("{}:{}-{}\t{}".format(l1[0], l1[1], l1[2], int(10.0 * v)))
            for i in range(1, nsamples):
                li = Next(readers[i])
                v = float(li[3]) / facts[0]
                sys.stdout.write("\t{}".format(int(10.0 * v)))
            sys.stdout.write("\n")
    finally:
        for r in streams:
            r.close()

def parseCoords(s):
    """Parse a string of the form chr:start-end and return the three components as a tuple."""
    p1 = s.find(":")
    p2 = s.find("-")
    return (s[:p1], s[p1+1:p2], s[p2+1:])

def extractSignificant(diffpeaks, log2fc, pval):
    nin = 0
    up = []
    dn = []
    header = "#Chrom\tStart\tEnd\tlog2(FC)\tP-value\tMean1\tMean2\n"
    with open(diffpeaks, "r") as f, open("alldiffpeaks.csv", "w") as dout, open("sigpeaks.csv", "w") as out, open("diffpeaks.bedGraph", "w") as diffbed, open("sigpeaks.bedGraph", "w") as sigbed:
        out.write(header)
        dout.write(header)
        f.readline()
        c = CSVreader(f)
        for line in c:
            try:
                fc = float(line[2])
                p  = float(line[6])
                nin += 1
            except ValueError: # Some P-values are NA
                continue
            mean1 = float(line[1])
            mean2 = mean1 * fc**2
            (chrom, start, end) = parseCoords(line[0].strip('"'))
            diffbed.write("{}\t{}\t{}\t{:.2f}\n".format(chrom, start, end, fc))
            dout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, start, end, fc, p, mean1, mean2))
            if abs(fc) >= log2fc and p <= pval:
                if fc > 0:
                    up.append((chrom, start, end, fc, p, mean1, mean2))
                else:
                    dn.append((chrom, start, end, -fc, p, mean1, mean2))
                out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, start, end, fc, p, mean1, mean2))
                sigbed.write("{}\t{}\t{}\t{:.2f}\n".format(chrom, start, end, fc))

        sys.stderr.write("{} in, {} up, {} down\n".format(nin, len(up), len(dn)))
        up.sort(key=lambda s: s[3], reverse=True)
        dn.sort(key=lambda s: s[3], reverse=True)
        with open("test-up.csv", "w") as out:
            out.write(header)
            for rec in up:
                out.write("\t".join([str(s) for s in rec]) + "\n")

        with open("ctrl-up.csv", "w") as out:
            out.write(header)
            for rec in dn:
                out.write("\t".join([str(s) for s in rec]) + "\n")

def extractSignificantGenes(diffpeaks, log2fc, pval):
    nin = 0
    up = []
    dn = []
    with open(diffpeaks, "r") as f, open("sigpeaks.csv", "w") as out:
        out.write("#Gene\tlog2(FC)\tP-value\n")
        f.readline()
        c = CSVreader(f)
        for line in c:
            gname = line[0]
            try:
                fc = float(line[2])
                p  = float(line[6])
                nin += 1
            except ValueError: # Some P-values are NA
                continue
            if abs(fc) >= log2fc and p <= pval:
                if fc > 0:
                    up.append((gname, fc, p))
                else:
                    dn.append((gname, -fc, p))
                out.write("{}\t{}\t{}\n".format(gname, fc, p))

        sys.stderr.write("{} in, {} up, {} down\n".format(nin, len(up), len(dn)))
        up.sort(key=lambda s: s[1], reverse=True)
        dn.sort(key=lambda s: s[1], reverse=True)
        with open("test-up.csv", "w") as out:
            out.write("#Gene\tlog2(FC)\tP-value\n")
            for rec in up:
                out.write("\t".join([str(s) for s in rec]) + "\n")

        with open("ctrl-up.csv", "w") as out:
            out.write("#Gene\tlog2(FC)\tP-value\n")
            for rec in dn:
                out.write("\t".join([str(s) for s in rec]) + "\n")

def combineChromSizes(filenames):
    sizes = {}
    for filename in filenames:
        with open(filename, "r") as f:
            c = CSVreader(f)
            for line in c:
                chrom = line[0]
                size = int(line[1])
                if chrom == "*":
                    continue
                if chrom in sizes:
                    sizes[chrom] = max(sizes[chrom], size)
                else:
                    sizes[chrom] = size
    for chrom in sorted(sizes.keys()):
        sys.stdout.write("{}\t{}\n".format(chrom, sizes[chrom]))

def computeFactors(statsfile, nreadsfile):
    data = {}
    with open(statsfile, "r") as f:
        c = CSVreader(f)
        for line in c:
            smp = line[0]
            totopen = int(line[1])
            data[smp] = {"totopen": totopen}
    with open(nreadsfile, "r") as f:
        c = CSVreader(f)
        for line in c:
            smp = line[0]
            nreads = int(line[1])
            data[smp]["nreads"] = nreads

    maxnreads = 0
    maxopen = 0
    for smp in data.keys():
        smpdata = data[smp]
        if smpdata["nreads"] > maxnreads:
            maxnreads = smpdata["nreads"]
            maxopen = smpdata["totopen"]
    for smp in data.keys():
        smpdata = data[smp]
        smpdata["readsfact"] = 1.0 * maxnreads / smpdata["nreads"]
        smpdata["openfact"] = 1.0 * maxopen / smpdata["totopen"]
        smpdata["factor"] = smpdata["readsfact"] / smpdata["openfact"]
        sys.stdout.write("{}\t{}\t{}\t{}\n".format(smp, smpdata["readsfact"], smpdata["openfact"], smpdata["factor"]))

def computeFrip(fripsfile, nreadsfile):
    frips = {}
    with open(fripsfile, "r") as f:
        c = CSVreader(f)
        for line in c:
            frips[line[0]] = int(line[1])
    with open(nreadsfile, "r") as f:
        c = CSVreader(f)
        for line in c:
            smp = line[0]
            nreads = int(line[1])
            sys.stdout.write("{}\t{}\t{}\t{}\n".format(smp, 1.0 * frips[smp] / nreads, frips[smp], nreads))

## Code to read a nextflow config file and return it as a dict of dicts.

def readNextflowConfig(filename):
    dic = {}
    with open(filename, "r") as f:
        text = f.read().replace("\n", " ")
    words = shlex.split(text)

    state = 0
    key = ""
    key2 = ""

    for w in words:
        if w == '=':
            continue

        if state == 0:          # we're at top level, current work is top-level key
            key = w
            dic[key] = {}
            state = 1
        elif state == 1:
            if w == '{':
                state = 2
            else:
                sys.stderr.write("Error parsing `{}': '{{' expected, found '{}'.\n".format(filename, w))
                return False
        elif state == 2:
            if w == '}':
                state = 0
            else:
                key2 = w
                state = 3
        elif state == 3:
            if w == '}':
                sys.stderr.write("Error parsing `{}': word expected, found '}'.\n".format(filename))
                return False
            else:
                dic[key][key2] = w
                state = 2
    return dic

## Code to write WashU hubs

def readIDs(filename, col=0, unique=False):
    with open(filename, "r") as f:
        data = f.read().rstrip("\n").split("\n")
    ids = [ row.split("\t")[col] for row in data ]
    if unique:
        uids = []
        for id in ids:
            if id not in uids:
                uids.append(id)
        return uids
    return ids

def readTable(filename):
    """Read a tab-delimited file, returning a list of lists."""
    with open(filename, "r") as f:
        rows = f.read().rstrip("\n").split("\n")
    return [r.split("\t") for r in rows]

def toDict(table):
    d = {}
    for row in table:
        k = row[0]
        if k:
            d[k] = row[1:]
    return d

def hubentry(params):
    return """{{
  "type": "bigwig",
  "name": "{smp}",
  "height": 50,
  "group": 1,
  "colorpositive": "#000000",
  "colornegative": "#000000",
  "horizontalLines": [],
  "mode": "show",
  "url": "{url}/{smp}/{smp}.bw"
}},
{{
  "type": "bigbed",
  "name": "{smp}_peaks",
  "height": 50,
  "group": 2,
  "colorpositive": "#000000",
  "colornegative": "#000000",
  "horizontalLines": [{{"value": 0, "color": "#000000"}}],
  "mode": "thin",
  "url": "{url}/{smp}/{smp}.bb"
}},
""".format(**params)

def writeHubs(samplesfile, contrastsfile, outdir, url):
    samples = readIDs(samplesfile, col=1)
    conditions = readIDs(samplesfile, unique=True)
    trackurl = url

    # Hub for samples
    sys.stderr.write("Writing peaks.json\n")
    with open(outdir + "/peaks.json", "w") as out:
        out.write("[\n")
        for smp in samples:
            out.write(hubentry({"smp": smp, "url": trackurl}))
        out.write("]\n")
        
    # Hub for conditions
    sys.stderr.write("        condpeaks.json\n")
    with open(outdir + "/condpeaks.json", "w") as out:
        out.write("[\n")
        for cond in conditions:
            out.write(hubentry({"smp": cond, "url": trackurl}))
        out.write("]\n")

    # Hubs for contrasts
    contr = readTable(contrastsfile)
    for pair in contr:
        label = "{}.vs.{}".format(pair[0], pair[1])
        sys.stderr.write("        {}.json\n".format(label))
        with open("{}/{}.json".format(outdir, label), "w") as out:
            out.write("[\n")
            out.write(hubentry({"smp": pair[0], "url": trackurl}))
            out.write(hubentry({"smp": pair[1], "url": trackurl}))
            #out.write(hubentry({"smp": label, "url": trackurl}))
            out.write("""{{
  "type": "bigwig",
  "name": "{lab}",
  "height": 50,
  "group": 3,
  "colorpositive": "#FF0000",
  "colornegative": "#0000FF",
  "horizontalLines": [{{"value": 0, "color": "#000000"}}],
  "mode": "show",
  "url": "{url}/{lab}/{lab}.sig.bw"
}},
""".format(**{"lab": label, "url": trackurl}))
            out.write("""{{
  "type": "bigwig",
  "name": "{lab}",
  "height": 50,
  "group": 3,
  "colorpositive": "#FF0000",
  "colornegative": "#0000FF",
  "horizontalLines": [{{"value": 0, "color": "#000000"}}],
  "mode": "show",
  "url": "{url}/{lab}/{lab}.diff.bw"
}}
""".format(**{"lab": label, "url": trackurl}))
            out.write("]\n")

## Write final html report and README

def getTag(s):
    p = s.find(">")
    if p > 0:
        return s[4:p-2]
    else:
        return ""

def linkify(path, name, target=""):
    tg = " target='{}'".format(target) if target else ""
    return "<A href='{}/{}'{}>{}</A>".format(path, name, tg, name)

def writeReport(samplesfile, contrastsfile, outdir, template, configfile): # template, title, baseurl, flags="T"):
    params = readNextflowConfig(configfile)["params"]
    template = params["reportTemplate"] if "reportTemplate" in params else template
    title = params["reportTitle"]
    params["washurl"] = "http://epigenomegateway.wustl.edu/browser/?genome={}&datahub={}/{}".format(params["hubOrganism"], params["hubURL"], params["hubName"])
    outdir = params["reportDir"]
    samples = readTable(samplesfile)
    conditions = readIDs(samplesfile, unique=True)
    contrasts = readTable(contrastsfile)
    indexfile = "index.html"

    with open(indexfile, "w") as out:
        with open(template, "r") as f:
            for line in f:
                if line.startswith("<!--"):
                    code = getTag(line)
                    if code == "Title":
                        out.write(title)
                    elif code == "Table1":
                        writeTable1(out, samples)
                    elif code == "Table2":
                        writeTable2(out, params)
                    elif code == "Table3":
                        writeTable3(out, contrasts)
                    elif code == "Table4":
                        writeTable4(out, contrasts, params)
                    elif code == "Table5":
                        writeTable5(out, contrasts, params)
                    elif code == "Table6":
                        writeTable6(out, contrasts, params)
                    elif code == "Table7":
                        writeTable7(out, contrasts, params)
                    elif code == "Table8":
                        writeTable8(out, contrasts, params)
                    else:
                        out.write(line)
                else:
                    out.write(line)

    with open("README", "w") as out:
        hubname = params["hubName"]
        out.write("""DASA Results Directory

The file "index.html" is the main output file. It should be opened with a web browser.
Plots contained in the report are atored in the plots/ subdirectory.

The data/ directory contains the results of differential analysis. There is one
subdirectory for each contrast (e.g. KO.vs.WT) containing the following files:

  KO.vs.WT-common-peaks.bed       Coordinates of common peaks that were compared
  KO.vs.WT-diffpeaks.csv          Results of DESeq2 on all common peaks
  KO.vs.WT-sigpeaks.csv           Peaks showing significant difference
  KO.vs-WT-test-up.csv            Peaks that are significantly higher in test condition (e.g. KO1)
  KO.vs-WT-ctrl-up.csv            Peaks that are significantly higher in control condition (e.g. WT)

The {} directory contains the bigWig and bigBed files to display these results in the
WashU epigenetic browser. You should move the contents of this directory to a web server
so that they can be accessed through the following URL:

  {}

If configuration is correct, the links in Sections 2 and 4 of the HTML report will work properly.
""".format(params["hubName"], params["hubURL"]))

def writeTable1(out, samples):
    factors = toDict(readTable("data/sample-factors.txt"))
    nreads = toDict(readTable("data/all-sample-reads.txt"))
    stats = toDict(readTable("data/all-sample-stats.txt"))
    frips = toDict(readTable("data/all-frips.txt"))
    out.write("<tr><th>Sample</th><th>Aligned reads</th><th>Open bases</th><th>Factor</th><th>Peaks</th><th>FRIP</th></tr>")
    for row in samples:
        smp = row[1]
        out.write("<tr><th>{}</th><td align='right'>{:,}</td><td align='right'>{:,} bp</td><td align='right'>{:.3f}</td><td align='right'>{:,}</td><td align='right'>{:.1f}%</td></tr>\n".format(smp, int(nreads[smp][0]), int(stats[smp][0]), float(factors[smp][2]), int(stats[smp][1]), 100.0*float(frips[smp][0])))

def writeTable2(out, params):
    out.write("""<TR>
  <TD align='center' width='50%'><A href='{}/peaks.json' target='washu'>SAMPLES</A></TD>
  <TD align='center'><A href='{}/condpeaks.json' target='washu'>CONDITIONS</A></TD>
</TR>""".format(params["washurl"], params["washurl"]))

def writeTable3(out, contrasts):
    contrdata = toDict(readTable("data/all-contr-stats.txt"))
    out.write("<tr><th>Test</th><th>Ctrl</th><th>Number of peaks</th><th>Avg size in Test</th><th>Avg size in Ctrl</th></tr>")
    for contr in contrasts:
        label = contr[0] + ".vs." + contr[1]
        ctrdata = contrdata[label]
        out.write("<tr><th>{}</th><th>{}</th><td align='right'>{:,}</td><td align='right'>{:.1f} bp</td><td align='right'>{:.1f} bp</td></tr>".format(
            contr[0], contr[1], int(ctrdata[0]), float(ctrdata[1]), float(ctrdata[2])))

def writeTable4(out, contrasts, params):
    contrdata = toDict(readTable("data/all-contr-counts.txt"))
    out.write("<tr><th>Test</th><th>Ctrl</th><th>Significant</th><th>Up in Test</th><th>Up in Ctrl</th><th>Peaks</th><th>WashU hub</th></tr>")
    for contr in contrasts:
        label = contr[0] + ".vs." + contr[1]
        ctrdata = contrdata[label]
        peaksetlink = """<BR><A target="_geneset" href="/tools/rs/index.cgi?cmd=add&name={}&source={}&org={}&gf={}">Create RegionSet</A>""".format(params["hubName"] + "-" + label, params["hubName"], params["hubOrganism"], 
                                                                                                                                               params["hubURL"] + "/data/" + label + "/" + label + "-sigpeaks.csv")

        out.write("<tr><th>{}</th><th>{}</th><td align='right'>{:,}</td><td align='right'>{:,}</td><td align='right'>{:,}</td><td align='center'>Sig: {}<BR>All: {}{}</td><td align='center'>{}</td></tr>".format(
            contr[0], contr[1], int(ctrdata[0]), int(ctrdata[1]), int(ctrdata[2]),
            linkify("data/" + label + "/", label + "-sigpeaks.xlsx"),
            linkify("data/" + label + "/", label + "-diffpeaks.xlsx"),
            peaksetlink if "sets" in params else "",
            linkify(params["washurl"], label + ".json", target="washu")))
    
def writeTable5(out, contrasts, params):
    contrdata = toDict(readTable("data/all-genediff-counts.txt"))
    out.write("<tr><th rowspan='2'>Test</th><th rowspan='2'>Ctrl</th><th rowspan='2'>Significant</th><th colspan='2'>Accessibility</th><th rowspan='2'>Genes</th></tr>\n")
    out.write("<tr><th>Increased</th><th>Decreased</th></tr>\n")
    for contr in contrasts:
        label = contr[0] + ".vs." + contr[1]
        ctrdata = contrdata[label]
        out.write("<tr><th>{}</th><th>{}</th><td align='right'>{:,}</td><td align='right'>{:,}</td><td align='right'>{:,}</td><td align='center'>{}</td></tr>\n".format(
            contr[0], contr[1], int(ctrdata[0]), int(ctrdata[1]), int(ctrdata[2]), linkify("data/" + label + "/", label + ".TSS.xlsx")))
    
def writeTable6(out, contrasts, params):
    contrdata = toDict(readTable("data/all-genebodydiff-counts.txt"))
    out.write("<tr><th rowspan='2'>Test</th><th rowspan='2'>Ctrl</th><th rowspan='2'>Significant</th><th colspan='2'>Accessibility</th><th rowspan='2'>Genes</th></tr>\n")
    out.write("<tr><th>Increased</th><th>Decreased</th></tr>\n")
    for contr in contrasts:
        label = contr[0] + ".vs." + contr[1]
        ctrdata = contrdata[label]
        out.write("<tr><th>{}</th><th>{}</th><td align='right'>{:,}</td><td align='right'>{:,}</td><td align='right'>{:,}</td><td align='center'>{}</td></tr>\n".format(
            contr[0], contr[1], int(ctrdata[0]), int(ctrdata[1]), int(ctrdata[2]), linkify("data/" + label + "/", label + ".genes.xlsx")))

def writeTable7(out, contrasts, params):
    datafile = "data/all-enhancersdiff-counts.txt"
    if os.path.isfile(datafile):
        contrdata = toDict(readTable("data/all-enhancersdiff-counts.txt"))
        if contrdata:
            out.write("<tr><th rowspan='2'>Test</th><th rowspan='2'>Ctrl</th><th rowspan='2'>Significant</th><th colspan='2'>Accessibility</th><th rowspan='2'>Genes</th></tr>\n")
            out.write("<tr><th>Increased</th><th>Decreased</th></tr>\n")
            for contr in contrasts:
                label = contr[0] + ".vs." + contr[1]
                ctrdata = contrdata[label]
                out.write("<tr><th>{}</th><th>{}</th><td align='right'>{:,}</td><td align='right'>{:,}</td><td align='right'>{:,}</td><td align='center'>{}</td></tr>\n".format(
                    contr[0], contr[1], int(ctrdata[0]), int(ctrdata[1]), int(ctrdata[2]), linkify("data/" + label + "/", label + ".enhancers.xlsx")))

def writeTable8(out, contrasts, flags):
    out.write("<tr><th>Test</th><th>Ctrl</th><th>Peak sizes</th><th>TSS</th><th>Test Peaks</TH><TH>Ctrl Peaks</TH></tr>")
    for contr in contrasts:
        label = contr[0] + ".vs." + contr[1]
        out.write("<tr><th>{}</th><th>{}</th><td align='center'>{}</td><td align='center'>{}</td><td align='center'>{}</td><td align='center'>{}</td></tr>".format(
            contr[0], contr[1], 
            linkify("plots/" + label, label + ".scatterplot.png"),
            linkify("plots/" + label, label + ".TSS.png"),
            linkify("plots/" + label, label + ".testpeaks.png"),
            linkify("plots/" + label, label + ".ctrlpeaks.png")))

# Gene matrix

def geneMatrix(samplesfile, factorsfile):
    samples = readTable(samplesfile)
    samplenames = [smp[1] for smp in samples]
    factors = readFactors(factorsfile)

    conds = []
    labels = []
    idx = 0
    cond = ""
    for smp in samples:
        if smp[0] != cond:
            conds.append('"' + smp[0] + '"')
            cond = smp[0]
            idx += 1
        labels.append(idx)
    with open("levels.txt", "w") as out:
        out.write(",".join(conds) + "\n")
    with open("labels.txt", "w") as out:
        out.write(",".join([str(x) for x in labels]) + "\n")
        
    sys.stdout.write("#Gene\t" + "\t".join(samplenames) + "\n")
    c = CSVreader(sys.stdin)
    prevkey = ""
    for row in c:
        key = row[0] + ":" + row[1] + "-" + row[2]
        if key == prevkey:
            continue
        sys.stdout.write("\t".join(row[4:]) + "\n")
        prevkey = key

# Convert tab-delimited files to excel

def toExcel(filenames):
    dest = filenames[0]
    srcs = filenames[1:]
    sys.stderr.write("Writing {}...\n".format(dest))
    workbook = xlsxwriter.Workbook(dest, {'strings_to_numbers': True})
    bold = workbook.add_format({'bold': 1})
    for src in srcs:
        if ":" in src:
            parts = src.split(":")
            src = parts[0]
            sname = parts[1].replace("_", " ")
        else:
            sname = os.path.splitext(os.path.split(src)[1])[0]
        sys.stderr.write("+ {}\n".format(src))
        ws = workbook.add_worksheet(sname)
        r = c = 0
        with open(src, "r") as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                c = 0
                for col in row:
                    if r == 0:
                        ws.write(r, c, col, bold)
                    else:
                        ws.write(r, c, col)
                    c += 1
                r += 1
    workbook.close()

def makeRegions(infile, outfile, code, upstream, downstream):
    """Adjust regions contained in `infile' writing them to `outfile'. If code is 't', 
output region is from 'upstream' bases before TSS to 'downstream' bases after it. If
code is 'b', region is from 'upstream' before TSS to 'downstream' bases after TES.
Strand (in column 4) is taken into account."""
    with open(outfile, "w") as out:
        with open(infile, "r") as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if row[0][0] == '#':
                    continue
                start = int(row[1])
                end   = int(row[2])
                if code == 'b':
                    if row[3] == "+":
                        nstart = start - upstream
                        nend   = end + downstream
                    else:
                        nstart = start - downstream
                        nend   = end + upstream
                elif code == 't':
                    if row[3] == "+":
                        nstart = start - upstream
                        nend   = start + downstream
                    else:
                        nstart = end - downstream
                        nend   = end + upstream
                out.write("{}\t{}\t{}\t{}\t{}\n".format(row[0], max(nstart, 0), nend, row[3], row[4]))
                

def main(cmd, args):
    if cmd == "convert":
        smp = args[0]
        (totopen, npeaks) = convertPeaks(args[1], args[2])
        sys.stdout.write("{}\t{}\t{}\n".format(smp, totopen, npeaks))
    elif cmd == "stats":
        computeStats(args)
    elif cmd == "matrix":
        writeMatrix(args[0], args[1], args[2], args[3], args[4:])
    elif cmd == "sig":
        extractSignificant(args[0], float(args[1]), float(args[2]))
    elif cmd == "gsig":
        extractSignificantGenes(args[0], float(args[1]), float(args[2]))
    elif cmd == "sizes":
        combineChromSizes(args)
    elif cmd == "factors":
        computeFactors(args[0], args[1])
    elif cmd == "gmatrix":
        geneMatrix(args[0], args[1])
    elif cmd == "frip":
        computeFrip(args[0], args[1])
    elif cmd == "hubs":
        writeHubs(args[0], args[1], args[2], args[3])
    elif cmd == "report":
        writeReport(args[0], args[1], args[2], args[3], args[4])
    elif cmd == "xlsx":
        toExcel(args)
    elif cmd == "regions":
        makeRegions(args[0], args[1], args[2], int(args[3]), int(args[4]))
    else:
        sys.stderr.write("Missing dasatools.py command!\n")
        sys.exit(1)

if __name__ == "__main__":
    args = sys.argv[1:]
    main(args[0], args[1:])
