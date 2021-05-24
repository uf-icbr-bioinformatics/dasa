#!/usr/bin/env python

# coding: utf-8
import sys
import numpy as np
import pandas as pd
import matplotlib
from matplotlib import gridspec
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde, linregress
from mpl_toolkits.axes_grid1 import make_axes_locatable

class DensityPlot():
    infile = None
    outfile = None
    log = False
    hasHeader = None
    skipRows = None
    xlabel = None
    ylabel = None
    xsize = 10
    ysize = 8
    pointSize = 50
    imgformat = "png"
    cx = 4
    cy = 3
    regprob = 0.0
    boxplots = False
    outliers = False
    fixed_scale = False
    max_x = None
    max_y = None

    def parseArgs(self, args):
        if "-h" in args or "--help" in args:
            return usage()
        prev = ""
        for a in args:
            if prev == "-cx":
                self.cx = int(a) - 1
                prev = ""
            elif prev == "-cy":
                self.cy = int(a) - 1
                prev = ""
            elif prev == "-c":
                self.outliers = float(a)
                prev = ""
            elif prev == "-s":
                self.skipRows = int(a)
                prev = ""
            elif prev == "-p":
                self.pointSize = int(a)
                prev = ""
            elif prev == "-xl":
                self.xlabel = a
                prev = ""
            elif prev == "-yl":
                self.ylabel = a
                prev = ""
            elif prev == "-xs":
                self.xsize = float(a)
                prev = ""
            elif prev == "-ys":
                self.ysize = float(a)
                prev = ""
            elif prev == "-f":
                self.imgformat = a
                prev = ""
            elif prev == "-R":
                self.regprob = float(a)
                prev = ""
            elif prev == "-m":
                parts = a.split(",")
                self.max_x = float(parts[0])
                self.max_y = float(parts[1])
                self.fixed_scale = True
                prev = ""
            elif a in ["-c", "-cx", "-cy", "-s", "-p", "-xl", "-yl", "-f", "-xs", "-ys", "-R", "-m"]:
                prev = a
            elif a == "-l":
                self.log = np.log(2)
            elif a == "-l10":
                self.log = np.log(10)
            elif a == "-t":
                self.hasHeader = 0
            elif a == "-r":
                self.regprob = 1.0
            elif a == "-b":
                self.boxplots = True
            elif self.infile is None:
                self.infile = a
            elif self.outfile is None:
                self.outfile = a
        if self.infile and self.outfile:
            return True
        else:
            return usage()

    def run(self):
        regx = []
        regy = []

        df = pd.read_table(self.infile, header=self.hasHeader, skiprows=self.skipRows)
        sys.stderr.write("Data file read.\n")

        dx = df[self.cx]
        dy = df[self.cy]

        if self.outliers:
            sys.stderr.write("Removing top/bottom {}% outliers:\n".format(self.outliers))
            xlimit1 = np.percentile(dx, self.outliers)
            ylimit1 = np.percentile(dy, self.outliers)
            xlimit2 = np.percentile(dx, 100-self.outliers)
            ylimit2 = np.percentile(dy, 100-self.outliers)

            newx = []
            newy = []
            for i in range(len(dx)):
                if (xlimit1 <= dx[i] <= xlimit2) and (ylimit1 <= dy[i] <= ylimit2):
                    newx.append(dx[i])
                    newy.append(dy[i])
            dx = newx
            dy = newy

        if self.log:
            sys.stderr.write("Applying log transformation.\n")
            dy = [ np.log(y+0.5) / self.log for y in dy ]
            dx = [ np.log(x+0.5) / self.log for x in dx ]
            min_x = min(dx)-0.1
            min_y = min(dy)-0.1

        sys.stderr.write("Computing densities.\n")
        xy = np.vstack([dx, dy])
        z = gaussian_kde(xy)(xy)

        sys.stderr.write("Plotting.\n")
        idx = z.argsort()
        x, y, z = np.array(dx)[idx], np.array(dy)[idx], np.array(z)[idx]

        if self.boxplots:
            # create figure (note: will get automatically used as current "main" figure)
            plt.figure(figsize=(self.xsize, self.ysize))
            # create scatterplot subplot with an aspect ratio of 1.0 since we want it to be a square
            ax = plt.subplot(aspect=1.0)
            # create a divider based on the scatterplot subplot
            div1 = make_axes_locatable(ax)
            # append a subplot to the right of the divider to be the first boxplot
            bp1 = div1.append_axes("right", size="35%", pad=0.4, sharey=ax)
            # append another subplot to be the second boxplot
            bp2 = div1.append_axes("right", size="35%", pad=0.4, sharey=ax)
            # hide the y tick labels for the boxplots
            for bp in (bp1, bp2):
                plt.setp(bp.get_yticklabels(), visible=False)
        else:
            fig, ax = plt.subplots(figsize=(self.xsize, self.ysize))
        ax.scatter(x, y, c=z, s=self.pointSize, edgecolor='')

        if self.log:
            plt.xlim(xmin=min_x)
            plt.ylim(ymin=min_y)
        if self.fixed_scale:
            plt.xlim(xmax=self.max_x)
            plt.ylim(ymax=self.max_y)
        diag_line, = ax.plot(ax.get_xlim(), ax.get_ylim(), ls="-", c=".3")

        # Do we want boxplots?
        if self.boxplots:
            bp1.boxplot(dx)
            bp2.boxplot(dy)

        # Do we want a regression line?
        if self.regprob:
            sys.stderr.write("Computing regression.\n")
            for i in range(len(dx)):
                if np.random.random() < self.regprob:
                    regx.append(dx[i])
                    regy.append(dy[i] - dx[i])
            regx = np.array(regx)
            regy = np.array(regy)
            slope, intercept, r_value, p_value, std_err = linregress(regx, regy)
            slope += 1
            sys.stderr.write("Slope: {}\nP-value: {}\n".format(slope, p_value))
            #print (slope, intercept, r_value, p_value, std_err)
            yline = regx * slope + intercept
            plt.plot(regx, yline, color="red")

        if self.xlabel:
            ax.set_xlabel(self.xlabel)
            if self.boxplots:
                bp1.set_title(self.xlabel)
        if self.ylabel:
            ax.set_ylabel(self.ylabel)
            if self.boxplots:
                bp2.set_title(self.ylabel)
        plt.savefig(self.outfile, format=self.imgformat)

def usage():
    sys.stdout.write("""density_scatterplot.py - Draw density scatterplots of paired data.

Usage: density_scatterplot.py [options] datafile imgfile

Read data from two columns of file `datafile' and draw a density heatmap 
of their scatterplot to `imgfile'. 

Options related to input data:

  -cx C | Use column C for X axis coordinates (default: {}).
  -cy C | Use column C for Y axis coordinates (default: {}).
  -s  S | Skip S rows from top of input file (default: {}).
  -c  C | Outlier removal. Clip the top and bottom C% of the values.
  -l    | If supplied, log-transform data (base 2).
  -l10  | If supplied, log-transform data (base 10).
  -r    | If supplied, draw regression line.
  -R R  | Same as -r, but subsample points with probability R.
  -b    | If supplied, draw box plots.

Graphical options:

  -xs S | Set X dimension of image to S inches (default: {}).
  -ys S | Set Y dimension of image to S inches (default: {}).
  -xl L | Set X axis label to L.
  -yl L | Set Y axis label to L.
  -p P  | Set dot size to P (default: {}).
  -f F  | Set output image format to F (default: {}).

""".format(DensityPlot.cx + 1, DensityPlot.cy + 1, DensityPlot.skipRows, DensityPlot.xsize, DensityPlot.ysize, DensityPlot.pointSize, DensityPlot.imgformat))
    return False

if __name__ == "__main__":
    DP = DensityPlot()
    if DP.parseArgs(sys.argv[1:]):
        DP.run()
