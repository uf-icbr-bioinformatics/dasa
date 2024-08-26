#!/bin/bash

SAMPLES=$1
CONTRASTS=$2
shift 2

if [[ -z "$SAMPLES" ]];
then
cat <<EOF

Usage: dasa.sh SAMPLES CONTRASTS [nextflow-args]

Any nextflow-args, if supplied, are passed to the nextflow run command.

Pipeline parameters that can be specified in nextflow.config:

SAMPLES is a tab-delimited file describing experimental design. It should have four columns:

  condition   sample   peaks_file   bam_file

Each experimental condition should be represented by at least two samples (biological replicates).
The peaks file is assumed to be in MACS format or in BED format, and bam_file is the BAM file that 
the peaks were called on. Peak files and BAM files should be located in the launch directory. 

CONTRASTS is a tab-delimited file with two columns, indicating the pairs of conditions to be 
compared with each other. Use the --example option to see an example of these two files.

Analysis options (can be specified in nextflow.config or passed on the command line):

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
  
EOF

  exit 1
fi

set -e

if [[ ! -f nextflow.config ]];
then 
  echo "No nextflow.config file!"
  exit 1
fi

nextflow run /apps/dibig_tools/dasa/main.nf \
  -resume -N ${USER}@ufl.edu -with-report dasa-report.html \
  -with-singularity /apps/dibig_tools/dasa/container/dasa.img \
  --samples $SAMPLES \
  --contrasts $CONTRASTS \
  $*

if [[ -f "mkreport.sh" ]];
then
  echo Generating report...
  bash mkreport.sh
  echo done.
fi
