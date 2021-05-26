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
The peaks file can be in BED format (with at least three columns: chromosome, start, end) or in MACS format, 
and bam_file is the BAM file that peaks were called on. Peak files and BAM files should be located in the launch directory.

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

## Operation

DASA will perform the following steps:

1. For each sample, compute total number of reads (from the BAM file), total number of peaks, and total size of all peaks ("open" region).
1. Compute normalization factors for all samples. The normalization factor for each sample is the combination of two values, one based on the number of reads, the second on the total open region. This accounts for the fact that peak height is directly proportional to the total number of reads, and inversely proportional to the total open region.
2. For each contrast, determine the set of peaks that can be compared. For every pair of partially-overlapping peaks in the two conditions being compared, we take either the common part (if merge_mode is "I") or the union of the two peaks ("U").
3. The "size" of each common peaks is computed, and scaled using the appropriate factor from Step 2. This produces a matrix with one row for each common peak, and one column for each sample in the conditions being compared, containing peak sizes.
4. Differential analysis is performed on the matrix using DESeq2. The results are filtered using the supplied log2fc and pvalue parameters, producing a list of differential peaks.
5. Finally, generate plots, genome browser tracks, and the final report.

## Graphical output

The pipeline will generate four plots for each contrast:

1. A scatterplot of peak sizes in the two conditions being compared. This is useful to determine if the average size of peaks changes in response to an experimental variable.
2. A tornado plot showing the average profile of all peaks that are higher in the test condition, and of the same peaks in the control condition.
3. A tornado plot showing the average profile of all peaks that are higher in the control condition, and of the same peaks in the test condition.
4. If the `tssfile` parameter is supplied: a tornado plot showing the average profile of the test and control condition around Transcription Start Sites.

## To Do

1. Containerization.
2. Conversion of output data files to Excel format.
3. More fine-grained control over normalization factors.
4. Better handling of `-resume`.
5. Speed up scatterplot generation.
