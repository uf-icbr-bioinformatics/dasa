#!/usr/bin/env python

import sys
import csv
import os.path
import xlsxwriter
from collections import defaultdict

PYVER = sys.version_info.major

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
    with open(diffpeaks, "r") as f, open("sigpeaks.csv", "w") as out, open("sigpeaks.bedGraph", "w") as bed:
        out.write("#Chrom\tStart\tEnd\tlog2(FC)\tP-value\n")
        f.readline()
        c = CSVreader(f)
        for line in c:
            try:
                fc = float(line[2])
                p  = float(line[6])
                nin += 1
            except ValueError: # Some P-values are NA
                continue
            (chrom, start, end) = parseCoords(line[0].strip('"'))
            bed.write("{}\t{}\t{}\t{}\n".format(chrom, start, end, fc))
            if abs(fc) >= log2fc and p <= pval:
                if fc > 0:
                    up.append((chrom, start, end, fc, p))
                else:
                    dn.append((chrom, start, end, -fc, p))
                out.write("{}\t{}\t{}\t{}\t{}\n".format(chrom, start, end, fc, p))

        sys.stderr.write("{} in, {} up, {} down\n".format(nin, len(up), len(dn)))
        up.sort(key=lambda s: s[3], reverse=True)
        dn.sort(key=lambda s: s[3], reverse=True)
        with open("test-up.csv", "w") as out:
            out.write("#Chrom\tStart\tEnd\tlog2(FC)\tP-value\n")
            for rec in up:
                out.write("\t".join([str(s) for s in rec]) + "\n")

        with open("ctrl-up.csv", "w") as out:
            out.write("#Chrom\tStart\tEnd\tlog2(FC)\tP-value\n")
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
        d[row[0]] = row[1:]
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
    with open(outdir + "/peaks.json", "w") as out:
        out.write("[\n")
        for smp in samples:
            out.write(hubentry({"smp": smp, "url": trackurl}))
        out.write("]\n")
        
    # Hub for conditions
    with open(outdir + "/condpeaks.json", "w") as out:
        out.write("[\n")
        for cond in conditions:
            out.write(hubentry({"smp": cond, "url": trackurl}))
        out.write("]\n")

    # Hubs for contrasts
    contr = readTable(contrastsfile)
    for pair in contr:
        label = "{}.vs.{}".format(pair[0], pair[1])
        with open("{}/{}.json".format(outdir, label), "w") as out:
            out.write("[\n")
            out.write(hubentry({"smp": pair[0], "url": trackurl}))
            out.write(hubentry({"smp": pair[1], "url": trackurl}))
            out.write(hubentry({"smp": label, "url": trackurl}))
            out.write("""{{
  "type": "bigwig",
  "name": "{smp}",
  "height": 50,
  "group": 3,
  "colorpositive": "#FF0000",
  "colornegative": "#0000FF",
  "horizontalLines": [{{"value": 0, "color": "#000000"}}],
  "mode": "show",
  "url": "{url}/{smp}/{smp}.diff.bw"
}}
""".format(**{"smp": label, "url": trackurl}))
            out.write("]\n")

## Write final html report and README

def getTag(s):
    p = s.find(">")
    if p > 0:
        return s[4:p-2]
    else:
        return ""

def linkify(path, name):
    return "<A href='{}/{}'>{}</A>".format(path, name, name)

def writeReport(samplesfile, contrastsfile, outdir, template, title, baseurl, flags="T"):
    samples = readTable(samplesfile)
    conditions = readIDs(samplesfile, unique=True)
    contrasts = readTable(contrastsfile)
    indexfile = outdir + "/index.html"

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
                        writeTable2(out, baseurl)
                    elif code == "Table3":
                        writeTable3(out, contrasts)
                    elif code == "Table4":
                        writeTable4(out, contrasts, baseurl)
                    elif code == "Table5":
                        writeTable5(out, contrasts, flags)
                    elif code == "Table6":
                        writeTable6(out, contrasts, flags)
                    else:
                        out.write(line)
                else:
                    out.write(line)

    readmefile = outdir + "/README"

    with open(readmefile, "w") as out:
        hubname = baseurl.rstrip("/").split("/")[-1]
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
""".format(hubname, baseurl))

def writeTable1(out, samples):
    factors = toDict(readTable("sample-factors.txt"))
    nreads = toDict(readTable("all-sample-reads.txt"))
    stats = toDict(readTable("all-sample-stats.txt"))
    frips = toDict(readTable("all-frips.txt"))
    out.write("<tr><th>Sample</th><th>Aligned reads</th><th>Open bases</th><th>Factor</th><th>Peaks</th><th>FRIP</th></tr>")
    for row in samples:
        smp = row[1]
        out.write("<tr><th>{}</th><td align='right'>{:,}</td><td align='right'>{:,} bp</td><td align='right'>{:.3f}</td><td align='right'>{:,}</td><td align='right'>{:.1f}%</td></tr>\n".format(smp, int(nreads[smp][0]), int(stats[smp][0]), float(factors[smp][2]), int(stats[smp][1]), 100.0*float(frips[smp][0])))

def writeTable2(out, baseurl):
    out.write("""<TR>
  <TD align='center' width='50%'><A href='{}/peaks.json' target='washu'>SAMPLES</A></TD>
  <TD align='center'><A href='{}/condpeaks.json' target='washu'>CONDITIONS</A></TD>
</TR>""".format(baseurl, baseurl))

def writeTable3(out, contrasts):
    contrdata = toDict(readTable("all-contr-stats.txt"))
    out.write("<tr><th>Test</th><th>Ctrl</th><th>Number of peaks</th><th>Avg size in Test</th><th>Avg size in Ctrl</th></tr>")
    for contr in contrasts:
        label = contr[0] + ".vs." + contr[1]
        ctrdata = contrdata[label]
        out.write("<tr><th>{}</th><th>{}</th><td align='right'>{:,}</td><td align='right'>{:.1f} bp</td><td align='right'>{:.1f} bp</td></tr>".format(
            contr[0], contr[1], int(ctrdata[0]), float(ctrdata[1]), float(ctrdata[2])))

def writeTable4(out, contrasts, baseurl):
    contrdata = toDict(readTable("all-contr-counts.txt"))
    out.write("<tr><th>Test</th><th>Ctrl</th><th>Significant</th><th>Up in Test</th><th>Up in Ctrl</th><th>Full results</th><th>WashU hub</th></tr>")
    for contr in contrasts:
        label = contr[0] + ".vs." + contr[1]
        ctrdata = contrdata[label]
        out.write("<tr><th>{}</th><th>{}</th><td align='right'>{:,}</td><td align='right'>{:,}</td><td align='right'>{:,}</td><td align='center'>{}</td><td align='center'>{}</td></tr>".format(
            contr[0], contr[1], int(ctrdata[0]), int(ctrdata[1]), int(ctrdata[2]),
            linkify("data/" + label + "/", label + "-diffpeaks.csv"),
            linkify(baseurl, label + ".json")))
    
def writeTable5(out, contrasts, flags):
    contrdata = toDict(readTable("all-genediff-counts.txt"))
    out.write("<tr><th rowspan='2'>Test</th><th rowspan='2'>Ctrl</th><th rowspan='2'>Significant</th><th colspan='2'>Accessibility</th><th rowspan='2'>Genes</th></tr>\n")
    out.write("<tr><th>Increased</th><th>Decreased</th></tr>\n")
    for contr in contrasts:
        label = contr[0] + ".vs." + contr[1]
        ctrdata = contrdata[label]
        out.write("<tr><th>{}</th><th>{}</th><td align='right'>{:,}</td><td align='right'>{:,}</td><td align='right'>{:,}</td><td align='center'>{} / {}</td></tr>\n".format(
            contr[0], contr[1], int(ctrdata[0]), int(ctrdata[1]), int(ctrdata[2]), linkify("data/" + label + "/", label + ".genes.test.csv"), linkify("data/" + label + "/", label + ".genes.ctrl.csv")))
    
def writeTable6(out, contrasts, flags):
    out.write("<tr><th>Test</th><th>Ctrl</th><th>Peak sizes</th><th>TSS</th><th>Test Peaks</TH><TH>Ctrl Peaks</TH></tr>")
    for contr in contrasts:
        label = contr[0] + ".vs." + contr[1]
        out.write("<tr><th>{}</th><th>{}</th><td align='center'>{}</td><td align='center'>{}</td><td align='center'>{}</td><td align='center'>{}</td></tr>".format(
            contr[0], contr[1], 
            linkify("plots/" + label, label + ".scatterplot.png"),
            linkify("plots/" + label, label + ".TSS.png") if "T" in flags else "-",
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
        writeReport(args[0], args[1], args[2], args[3], args[4], args[5])
    elif cmd == "xlsx":
        toExcel(args)
    else:
        sys.stderr.write("Missing dasatools.py command!\n")
        sys.exit(1)

if __name__ == "__main__":
    args = sys.argv[1:]
    main(args[0], args[1:])
