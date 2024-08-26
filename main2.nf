#!/usr/bin/env nextflow

/*
========================================================================================
                         dasa.nf
========================================================================================

Differential ATAC-Seq Analysis Pipeline.
#### Homepage / Documentation
https://github.com/uf-icbr-bioinformatics/dasa
#### Authors
Alberto Riva (ariva@ufl.edu), ICBR Bioinformatics Core, University of Florida.
----------------------------------------------------------------------------------------
*/

include { InitializeSamples; TSSRegions; GeneRegions; CountReads; MergeBEDs; ComputeFactors; PreCondFactors; MergeBAMs } from './processes.nf'

/* Parameters */

params.h = false
params.help = false
params.example = false
params.samples = "SAMPLES"
params.contrasts = "CONTRASTS"
params.merge_mode = "U"
params.log2fc = 1.0
params.pvalue = 0.01
params.log2fc_genes = 1.0
params.pvalue_genes = 0.01
params.reportDir = "Report"
params.deeptools = true
params.washu = true
params.genesfile = false
params.enhancersfile = false

/* For WashU hub creation */
params.hub = true // Set to false to disable hub generation
params.hubURL = "https://lichtlab.cancer.ufl.edu/reports/"
params.hubName = "hub"
params.hubOrganism = "hg38"

/* HTML report */
params.reportTitle = "ATAC-Seq Differential Analyisis"
params.reportTemplate = "${workflow.projectDir}/templates/report-template.html"
log.info """Report template: ${params.reportTemplate}"""

if (workflow.containerEngine) {
   log.info """Running in ${workflow.containerEngine} container ${workflow.container}."""
}

/* Internal variables */

samplesfile = file(params.samples, checkIfExists: true)
contrastsfile = file(params.contrasts, checkIfExists: true)
outdir = "${workflow.launchDir}/${params.reportDir}/"
log.info """Genes file: ${params.genesfile}"""

/* Functions */

def helpMessage() {
    log.info"""
Usage: nextflow run dasa.nf --samples SAMPLES --contrasts CONTRASTS [other options]

SAMPLES is a tab-delimited file describing experimental design. It should have four columns:

  condition   sample   peaks_file   bam_file

Each experimental condition should be represented by at least two samples (biological replicates).
The peaks file is assumed to be in MACS format or in BED format, and bam_file is the BAM file that 
the peaks were called on. Peak files and BAM files should be located in the launch directory. 

CONTRASTS is a tab-delimited file with two columns, indicating the pairs of conditions to be 
compared with each other. Use the --example option to see an example of these two files.

Analysis options:

  --log2fc          Threshold on log2(fold change) for differential analysis. Default: 1.0.

  --pvalue          Threshold on FDR-corrected P-value for differential analysis. Default: 0.01.

  --log2fc_genes    Threshold on log2(fold change) for differential TSS accessibility analysis. Default: 1.0.

  --pvalue_genes    Threshold on FDR-corrected P-value for differential TSS accessibility analysis. Default: 0.01.

  --merge_mode      Determines which area of partially-overlapping peaks to use in differential 
                    analysis. Can be either "I" (intersecting part is used) or "U" (the union of
                    the two peaks is used). Default: U.

  --genesfile       A file containing locations of genes. Should be a tab-delimited file with five columns:
                    chromosome, start, end, strand, name.

  --enhancersfile   A file containing location of enhancers. Should be a tab-delimited file with at least four columns:
  		    chromosome, start, end, linked gene.		    

Report options:

  --reportDir       Name of final report directory. Default: "Report".

  --reportTitle     Title of final report. Default: "ATAC-Seq Differential Analysis"

  --reportTemplate  Filename of HTML template for report.

WashU genome browser options (see README file in report directory for details):

  --hubName         Directory name of this hub.

  --hubOrganism     Organism code for genome browser.

  --hubURL          URL at which genome browser tracks will be published.
  
(c) 2021, A.Riva, University of Florida
"""
}

def example() {
    log.info"""
In this example, ATAC-Seq was performed on three different experimental conditions
(WT, KO1, KO2) each one having three biological replicates (WT_1, WT_2, etc). We wish
to compare peaks in KO1 and KO2 versus WT, and also KO1 vs KO2. We assume that peak
files are called S_peaks.xls and BAM files are called S.bam, where S is the sample name.

The SAMPLES file should loook like this:

  WT    WT_1    WT_1_peaks.xls    WT_1.bam
  WT    WT_2    WT_2_peaks.xls    WT_2.bam
  WT    WT_3    WT_3_peaks.xls    WT_3.bam
  KO1   KO1_1   KO1_1_peaks.xls   KO1_1.bam
  KO1   KO1_2   KO1_2_peaks.xls   KO1_2.bam
  KO1   KO1_3   KO1_3_peaks.xls   KO1_3.bam
  KO2   KO2_1   KO2_1_peaks.xls   KO2_1.bam
  KO2   KO2_2   KO2_2_peaks.xls   KO2_2.bam
  KO2   KO2_3   KO2_3_peaks.xls   KO2_3.bam

The CONTRASTS file should be:

  KO1    WT
  KO2    WT
  KO1    KO2

"""
}

if (params.help || params.h) {
   helpMessage()
   exit 0
}

if (params.example) {
   example()
   exit 0
}

workflow {
  samples = channel.fromPath(samplesfile)
	           .splitCsv(header:false, sep:'\t')
	           .map( row -> tuple(row[1], row[0], "${workflow.launchDir}/" + row[2], "${workflow.launchDir}/" + row[3]) )

  TSSRegions()
  GeneRegions()
  smpdata = InitializeSamples(samples)
  counts  = CountReads(smpdata.sample_files)

  totreads_per_condition = counts.cond_totreads_ch.groupTuple().map( { [ it[0], it[1].sum { it as Integer } ] } )

  MergeBEDs(smpdata.cond_beds)

  combined_stats = smpdata.combine_convert_stats.toSortedList( { a,b -> a[0] <=> b[0] } ).map( { it -> it.collect( {x -> x[1]} ).join(" ") } )
  combined_nreads = counts.combined_nreads.toSortedList( { a,b -> a[0] <=> b[0] } ).map( { it -> it.collect( {x -> x[1]} ).join(" ") } )

  ComputeFactors(combined_stats, combined_nreads)
  MergeBAMs(smpdata.sample_files.groupTuple(by: 1).map( row -> tuple(row[1], row[2].join(" ") ) ))
}
