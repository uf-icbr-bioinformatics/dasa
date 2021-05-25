# DASA
Differential ATAC-Seq Analysis

DASA is a complete tool for differential ATAC-Seq analysis. Starting from the output
of a peak caller like MACS2, it combines peaks from replicates of the same condition,
and it finds a set of comparable peaks between the two conditions. These peaks are then
quantified, using appropriate normalizations, and compared using DESeq2. The result
is a set of differential peaks, with an associated fold change and P-value.

DASA generates a full report of its analysis, including links to Excel files, plots, 
and hubs for the WashU epigenome browser. 

## Requirements

DASA currently needs the following programs to be installed and in PATH:

- Python (with numpy, scipy, pandas)
- samtools
- bedtools
- R (with the DESeq2 package)
- bedGraphToBigWig, bedToBigBed

A containerized version of DASA, that eliminates the need for dependencies, will be released in the near future.

## Usage

DASA is implemented with NextFlow. To run it, execute:

```
nextflow run /path/to/dasa.nf
```

When the pipeline is complete, it will create a directory called `Report` (this can be changed with the
--reportDir option) containing a full HTML report, output data files, and tracks for the WashU Epigenome Browswer.
See the README file in that directory for more details.

## Configuration

Parameters can be specified on the command line or in a `nextflow.config` file. The following parameters
are currently recognize:

Name | Description
-----|------------
samples    | File describing samples and experimental conditions. See description in the Input section below. Required.
contrasts  | File listing desired comparisons. See description in the Input section below. Required.
log2fc     | Threshold on the log2(fold change) for differential analysis. Default: 1.0.
pvalue     | Threshold on FDR-corrected P-value for differential analysis. Default: 0.01.
merge_mode | Determines which area of partially-overlapping peaks to use in differential analysis. Can be either "I" (intersecting part is used) or "U" (the union of the two peaks is used). Default: I.
reportDir  | Name of final report directory. Default: "Report".
reportName | Title of final report. Default: "ATAC-Seq Differential Analysis".
reportTemplate | Filename of HTML template for report. Default: builtin template.
hubName     | Directory name of WashU hub.
hubOrganism | Organism code for genome browser.
hubURL      | URL at which genome browser tracks will be published.

## Input

The Samples file is a tab-delimited file describing experimental design. It should have four columns:

```
  condition   sample   peaks_file   bam_file
```

Each experimental condition should be represented by at least two samples (biological replicates).
The peaks file is assumed to be in MACS format, and bam_file is the BAM file that MACS called
peaks on. Peak files and BAM files should be located in the launch directory.

The Contrasts file is a tab-delimited file with two columns, indicating the pairs of conditions to be
compared with each other. Use the --example option to see an example of these two files.

For example, let's assume we performed and ATAC-Seq experiment on three different experimental conditions
(WT, KO1, KO2) each one having three biological replicates (WT_1, WT_2, etc). We wish
to compare peaks in KO1 and KO2 versus WT, and also KO1 vs KO2. We assume that peak
files are called S_peaks.xls and BAM files are called S.bam, where S is the sample name.

The SAMPLES file should loook like this:

```
  WT    WT_1    WT_1_peaks.xls    WT_1.bam
  WT    WT_2    WT_2_peaks.xls    WT_2.bam
  WT    WT_3    WT_3_peaks.xls    WT_3.bam
  KO1   KO1_1   KO1_1_peaks.xls   KO1_1.bam
  KO1   KO1_2   KO1_2_peaks.xls   KO1_2.bam
  KO1   KO1_3   KO1_3_peaks.xls   KO1_3.bam
  KO2   KO2_1   KO2_1_peaks.xls   KO2_1.bam
  KO2   KO2_2   KO2_2_peaks.xls   KO2_2.bam
  KO2   KO2_3   KO2_3_peaks.xls   KO2_3.bam
```

The CONTRASTS file should be:

```
  KO1    WT
  KO2    WT
  KO1    KO2
```

