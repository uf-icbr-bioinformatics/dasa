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
params.tss_range = 2000
params.gene_body_upstream = 1000
params.gene_body_downstream = 1000

/* For WashU hub creation */
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

/* Input channels */

/* For every sample, output all info in the SAMPLES file */
Channel
	.fromPath(samplesfile)
	.splitCsv(header:false, sep:'\t')
	.map( row -> tuple(row[1], row[0], "${workflow.launchDir}/" + row[2], "${workflow.launchDir}/" + row[3]) )
	.into { samples_ch; cond_samples_ch; contr_samples_ch; contr_samples_ch_2 } /* smp, cond, peaks, bam */

/* For every condition, output the info in the SAMPLES file for all its samples. */
Channel
	.fromPath(samplesfile)
	.splitCsv(header:false, sep:'\t')
	.map( row -> tuple(row[0], row[1]) )
	.groupTuple()
	.set { conds_ch } /* cond, [smp1, smp2, ...] */

/* Output label, test, and ctrl condition for each contrast. */
Channel
	.fromPath(contrastsfile)
	.splitCsv(header:false, sep:'\t')
	.map( row -> tuple(row[0] + ".vs." + row[1], row[0], row[1]) )
	.into { contrasts_ch; contrasts_ch_2; contrasts_ch_3; contrasts_ch4; contr_for_tornado_ch; contr_for_scatterplot_ch; contr_for_genediff; contr_for_genebodydiff; contr_for_enhancersdiff } /* TEST.vs.CTRL, TEST, CTRL */

Channel
	.from contrasts_ch4
	.combine(contr_samples_ch_2) /* label, test, ctrl, smp, cond, peaks, bam */
	.filter( { it[4] == it[1] || it[4] == it[2] } )
	.map(row -> tuple(row[0], row[4], row[3]))

/* Processes */

/* Check that all input files exist, convert peaks to BED */
process InitializeSamples {
	executor "local"

	input: 
	tuple smp, condName, xlsfile, bamfile from samples_ch

	output:
	tuple smp, condName, bamfile into smp_count_reads
	tuple condName, file("${smp}.bed") into merge_bed_ch
	tuple smp, file("${smp}-convert-stats.txt") into combine_convert_stats /* convert-stats contains: totopen, npeaks */
	tuple smp, file("${smp}.bed"), bamfile into frip_ch, sample_bed_to_bb
	tuple smp, bamfile into sample_bam_to_bw_ch

	script:
	"""
	# Check that BAM files exist
	if [[ ! -e "$bamfile" ]];
	then
	  echo "Error: file $bamfile for condition $condName does not exist!"
	  exit 1
	fi
	if [[ ! -e "$xlsfile" ]];
	then
	  echo "Error: file $xlsfile for condition $condName does not exist!"
	  exit 1
	fi

	bedfile=${smp}.bed
	dasatools.py convert ${smp} $xlsfile \$bedfile >> ${smp}-convert-stats.txt
	"""
}

/* ** Generate file of TSS and gene body regions ** */

process TSSRegions {
	executor "local"

	output:
	file("tss-regions.txt") into tss_plots, tss_matrix

	script:
	"""
	#!/bin/bash

	dasatools.py regions ${params.genesfile} tss-regions.txt t ${params.tss_range} ${params.tss_range}
	"""
}

process GeneRegions {
	executor "local"

	output:
	path("gene-regions.txt") into gene_matrix

	script:
	"""
	#!/bin/bash

	dasatools.py regions ${params.genesfile} gene-regions.txt b ${params.gene_body_upstream} ${params.gene_body_downstream}
	"""
}


/* ** Determine number of reads for each BAM file ** */
process CountReads {
	executor "local"

	input:
	tuple smp, condName, bamfile from smp_count_reads

	output:
	stdout into max_nreads_ch
	tuple condName, stdout into cond_totreads_ch
	tuple smp, file("${smp}-nreads.txt") into combine_nreads

	script:
	"""
	NREADS=\$(samtools idxstats $bamfile | awk '{sum += \$3} END {print sum}')
	echo -e "${smp}\t\$NREADS" > ${smp}-nreads.txt
	echo \$NREADS
	"""
}

/* ** Determine total number of reads for each condition ** */
Channel
	.from cond_totreads_ch
	.groupTuple()
	.map( { [ it[0], it[1].sum { it as Integer } ] } )
	.set { compute_cond_factors_ch } /* cond, totalreads */

/* ** Determine the maximum number of reads ** */
Channel
	.from max_nreads_ch
	.map( { it as Integer } )
	.max()

/* ** Merge BED files for each condition ** */
Channel
	.from merge_bed_ch
	.groupTuple()
	.map( row -> tuple(row[0], row[1].join(" ")) )
	.set { merge_bed }

process MergeBEDs {
	executor "local"

	input:
	tuple cond, bedfiles from merge_bed

	output:
	tuple cond, file("${cond}.bed") into common_peaks_ch, post_merge_beds
	tuple cond, file("${cond}-bed-stats.txt") into cond_bed_stats_ch

	script:
	"""
	out=${cond}.bed
	sort -m -k1,1 -k2,2n -k3,3n $bedfiles | mergeBed -i - > \$out
	NPEAKS=\$(grep -c ^ \$out)
	TOTOPEN=\$(awk '{sum += (\$3-\$2)} END {print sum}' \$out)
	echo -e "$cond\t\$NPEAKS\t\$TOTOPEN" > ${cond}-bed-stats.txt
	"""
}

/* ** Compute normalization factors for samples ** */
Channel
	.from combine_convert_stats
	.toSortedList( { a,b -> a[0] <=> b[0] } )
	.map( { it -> it.collect( {x -> x[1]} ).join(" ") } )
	.set { combined_stats }

Channel
	.from combine_nreads
	.toSortedList( { a,b -> a[0] <=> b[0] } )
	.map( { it -> it.collect( {x -> x[1]} ).join(" ") } )
	.set { combined_nreads }

process ComputeFactors {
	executor "local"

	input:
	val allstats from combined_stats
	val allreads from combined_nreads

	output:
	file("sample-factors.txt") into factors_ch, factors_ch2, factors_ch3, factors_ch4
	file("all-sample-reads.txt") into nreads_for_frip_ch
	tuple file("sample-factors.txt"), file("all-sample-stats.txt"), file("all-sample-reads.txt") into report_stats

	publishDir "$outdir/data", mode: "copy"

	script:
	"""
	cat $allstats > all-sample-stats.txt
	cat $allreads > all-sample-reads.txt
	dasatools.py factors all-sample-stats.txt all-sample-reads.txt > sample-factors.txt
	"""
}

/* ** Compute normalization factors for each condition ** */
Channel
	.from compute_cond_factors_ch /* cond, totalreads */
	.join(cond_bed_stats_ch) /* cond, totalreads, statsfile */
	.set { cond_factors }

process PreCondFactors {
	executor "local"

	input:
	tuple cond, totalreads, statsfile from cond_factors

	output:
	tuple file("cond-nreads.txt"), statsfile into cond_factors_2

	script:
	"""
	echo -e "$cond\t$totalreads" > cond-nreads.txt
	"""
}

Channel
	.from cond_factors_2
	.reduce(tuple("", "")) { a, b -> tuple(a[0] + " " + b[0], a[1] + " " + b[1]) } 
	.set { cond_factors_3 }

process CondFactors {
	executor "local"

	input:
	tuple readsfiles, statsfiles from cond_factors_3

	output:
	file("cond-factors.txt") into condfactors_ch

	script:
	"""
	cat $readsfiles > all-cond-reads.txt
	cat $statsfiles > all-cond-stats.txt
	dasatools.py factors all-cond-stats.txt all-cond-reads.txt > cond-factors.txt
	"""
}

/* ** Merge BAM files for each condition ** */
Channel
	.from cond_samples_ch
	.groupTuple(by: 1)
	.map( row -> tuple(row[1], row[3].join(" ") ) )
	.set { merge_bams_ch }

process MergeBAMs {
	time "12h"
	memory "5G"

	input:
	tuple cond, bamfiles from merge_bams_ch

	output:
	tuple cond, file("${cond}.bam") into cond_bam_to_bw_ch, post_merge_bams

	script:
	"""
	samtools merge ${cond}.bam $bamfiles
	samtools index ${cond}.bam
	"""
}

/* ** Compute common peaks for each contrast ** */
Channel
	.from contrasts_ch
	.combine(common_peaks_ch)
	.filter( { it[3] == it[1] || it[3] == it[2] } )
	.groupTuple()
	.map( row -> tuple(row[0], row[4][0], row[4][1]) )
	.view()
	.set { common_peaks }

process CommonPeaks {
	executor "local"

	input:
	tuple contr, testbed, ctrlbed from common_peaks

	output:
	tuple contr, file("${contr}-common-peaks.bed") into quantify_peaks_ch
	tuple contr, file("sizes.txt") into size_plot_ch
	file("${contr}-peaks-stats.txt") into combine_cond_stats_ch

	publishDir "$outdir/data/$contr/", mode: "copy"

	script:
	"""
	if [[ "${params.merge_mode}" == "I" ]];
	then
	  intersectBed -a $testbed -b $ctrlbed -u | cut -f 1-3 | sort -k1,1 -k2,2n -k3,3n | uniq > ${contr}-common-peaks.bed
	else
	  sort -m -k1,1 -k2,2n -k3,3n $testbed $ctrlbed | mergeBed -i - > ${contr}-common-peaks.bed
	fi
	intersectBed -a $testbed -b $ctrlbed -wa |awk -v OFS='\t' '{print \$1,\$2,\$3,\$3-\$2}' > test-sizes.txt
	intersectBed -a $ctrlbed -b $testbed -wa |awk -v OFS='\t' '{print \$1,\$2,\$3,\$3-\$2}' > ctrl-sizes.txt
	paste test-sizes.txt ctrl-sizes.txt > sizes.txt
	TESTAVG=\$(awk -v OFS='\t' '{sum1+=\$4} END {print sum1/NR}' sizes.txt)
	CTRLAVG=\$(awk -v OFS='\t' '{sum2+=\$8} END {print sum2/NR}' sizes.txt)
	NPEAKS=\$(grep -c ^ ${contr}-common-peaks.bed)
	echo -e "${contr}\t\${NPEAKS}\t\${TESTAVG}\t\${CTRLAVG}" > ${contr}-peaks-stats.txt
	"""
}

/* ** Quantify common peaks in each contrast ** */

/* Things get hairy here. For each contrast:
   - take test and control condition
   - take the bam files for all samples in test condition, and all samples in control condition
   - call bedcov on commonpeaksbed against each bam file
   - collect resulting output files, pass them to dasatools matrix together with factors
*/

Channel
	.from contrasts_ch_2
	.combine(contr_samples_ch) /* label, testcond, ctrlcond, smp, cond, peaks, bam */
	.filter( { it[4] == it[1] || it[4] == it[2] } )
	.map( row -> tuple(row[0], row[3], row[4], row[6]) ) /* label, smp, cond, bam */
	.combine(quantify_peaks_ch) /* label, smp, cond, bam, label, bed */
	.filter( { it[0] == it[4] } )
	.set { quantify_peaks }

process QuantifyPeaks {
	time "5h"

	input:
	tuple contr, smp, cond, bamfile, dummy, bedfile from quantify_peaks

	output:
	tuple contr, cond, smp, file("counts.txt") into post_quantify

	script:
	"""
	samtools bedcov $bedfile $bamfile > counts.txt
	"""
}

/* ** Build the counts matrix for each contrast ** */
Channel
	.from contrasts_ch_3
	.combine(post_quantify) /* label, testcond, ctrlcond, label, cond, smp, counts */
	.filter( { it[0] == it[3] && (it[1] == it[4] || it[2] == it[4]) } )
	.map( row -> tuple(row[0], row[1], row[2], row[4], row[5], row[6]) )
	.groupTuple()
	.map( row -> tuple(row[0], row[1][0], row[2][0], row[3].join(","), row[4].join(","), row[5].join(" ")) )
	.set { build_matrix_ch }

process BuildMatrix {
	executor "local"

	input:
	tuple contr, testcond, ctrlcond, conditions, samples, countfiles from build_matrix_ch
	val factors from factors_ch

	output:
	tuple contr, testcond, ctrlcond, file("${contr}.matrix.csv"), file("labels.txt") into diffpeak_analysis_ch

	script:
	"""
	dasatools.py matrix $testcond $conditions $samples $factors $countfiles > ${contr}.matrix.csv
	"""
}

/* ** Perform differential analysis with DESeq2 ** */
process DiffPeakAnalysis {
	memory "10G"
	time "2h"

	input:
	tuple contr, testcond, ctrlcond, matrixfile, labelsfile from diffpeak_analysis_ch

	output:
	tuple contr, file("diffpeaks.csv") into extract_sig_ch

	publishDir "$outdir/data/$contr/", mode: "copy", pattern: "diffpeaks.csv", saveAs: { filename -> "${contr}-${filename}" }

	script:
	labels = file(labelsfile).text
	"""
#!/usr/bin/env Rscript

library("DESeq2")

datafile = "$matrixfile"
counts = as.matrix(read.csv(datafile, sep='\t', row.names=1))
counts[counts>2147483647]=2147483647
levels = c("$testcond", "$ctrlcond")
labels = c($labels)
sampleTable = data.frame(condition=factor(levels[labels]))
rownames(sampleTable) = colnames(counts)
dds = DESeqDataSetFromMatrix(countData = counts, colData = sampleTable, design = ~ condition)
dds = estimateSizeFactors(dds)
keep = rowMeans(counts(dds)) >= 5
dds = dds[keep,]
dds = DESeq(dds)
res = results(dds, contrast=c("condition", "$testcond", "$ctrlcond"))
write.table(res, file="diffpeaks.csv", sep='\t')
        """
}

process ExtractSignificant {
	executor "local"

	input:
	tuple label, diffpeaks from extract_sig_ch

	output:
	tuple label, file("diffpeaks.bedGraph"), file("sigpeaks.bedGraph") into bdg_to_bw_ch
	tuple label, file("test-up.csv"), file("ctrl-up.csv") into regions_to_tornado_ch
	file("sigpeaks.csv")
	file("sigpeaks.xlsx")
	file("diffpeaks.xlsx")
	file("contr-counts.txt") into combine_contr_counts_ch

	publishDir "$outdir/data/$label/", mode: "copy", pattern: "*.csv",  saveAs: { filename -> "${label}-${filename}" }
	publishDir "$outdir/data/$label/", mode: "copy", pattern: "*.xlsx", saveAs: { filename -> "${label}-${filename}" }

	"""
	#!/bin/bash
	# generate alldiffpeaks.csv sigpeaks.csv, sigpeaks.bedGraph, test-up.csv, ctrl-up.csv
	dasatools.py sig $diffpeaks ${params.log2fc} ${params.pvalue}
	N1=\$(grep -c ^ sigpeaks.csv)
	N2=\$(grep -c ^ test-up.csv)
	N3=\$(grep -c ^ ctrl-up.csv)
	N1=\$((N1-1))
	N2=\$((N2-1))
	N3=\$((N3-1))
	echo -e "$label\t\$N1\t\$N2\t\$N3" > contr-counts.txt
	dasatools.py xlsx sigpeaks.xlsx test-up.csv ctrl-up.csv
	dasatools.py xlsx diffpeaks.xlsx alldiffpeaks.csv
	"""
}

/* ** Compute FRIP score for each sample ** */

process FRIP {
	memory "16G"
	time "2h"

	input:
	tuple smp, bedfile, bamfile from frip_ch

	output:
	file("frip.txt") into combine_frip_ch

	script:
	"""
	FRIP=\$(coverageBed -sorted -a $bedfile -b $bamfile | awk '{sum += \$7} END {print sum}')
	echo -e "${smp}\t\$FRIP" > frip.txt
	"""
}

Channel
	.from combine_frip_ch
	.reduce("") { a, b -> a + " " + b }
	.set { combined_frip }

process CollectFRIP {
	executor "local"

	input:
	val fripsfiles from combined_frip
	val nreadsfile from nreads_for_frip_ch

	output:
	file("all-frips.txt") into report_frips

	publishDir "$outdir/data", mode: "copy"

	script:
	"""
	cat $fripsfiles > combined-frips.txt
	dasatools.py frip combined-frips.txt  $nreadsfile > all-frips.txt
	"""
}

/* ** Convert BAMs for each sample to BW ** */

process SampleBAMtoBW {
	memory "10G"
	time "5h"
	cpus 8

	input:
	tuple smp, bamfile from sample_bam_to_bw_ch
	val factors from factors_ch2

	output:
	file("${smp}.bw") into smp_bw_for_hub_ch

	publishDir "$outdir/${params.hubName}/$smp/", mode: "copy", pattern: "*.bw"

	script:
	"""
	bw=${smp}.bw
	fact=\$(grep ${smp} $factors | cut -f 4)
	bamBigWig.sh $bamfile \$bw \$fact
	"""
}

/* ** Convert BAM for each condition to BW ** */
process CondBAMtoBW {
	memory "10G"
	time "5h"
	cpus 8

	input:
	tuple cond, bamfile from cond_bam_to_bw_ch
	val condfactors from condfactors_ch

	output:
	tuple cond, file("${cond}.bw") into combine_for_tornado_ch
	file("${cond}.bw") into cond_bw_for_hub_ch

	publishDir "$outdir/${params.hubName}/$cond/", mode: "copy", pattern: "*.bw"

	script:
	"""
	bw=${cond}.bw
	fact=\$(grep ${cond} $condfactors | cut -f 4)
	bamBigWig.sh $bamfile \$bw \$fact
	"""
}

/* ** Convert BED files for samples and conditions into bigBed ** */

process SampleBEDtoBB {
	memory "2G"
	time "1h"

	input:
	tuple smp, bedfile, bamfile from sample_bed_to_bb

	output:
	file("${smp}.pk.bw") into smp_bb_for_hub_ch

	publishDir "$outdir/${params.hubName}/$smp/", mode: "copy", pattern: "*.pk.bw"

	script:
	"""
	samtools idxstats $bamfile | cut -f 1,2 | grep -v _ | grep -v ERCC > chrom.sizes
	#bedToBigBed -type=bed3+3 -tab $bedfile chrom.sizes ${smp}.bb
	bedToWig.py $bedfile chrom.sizes > ${smp}.wig
	wigToBigWig ${smp}.wig chrom.sizes ${smp}.pk.bw
	"""
}

Channel
	.from post_merge_beds
	.join(post_merge_bams)
	.set { cond_bed_to_bb }

process CondBEDtoBB {
	memory "2G"
	time "1h"

	input:
	tuple cond, bedfile, bamfile from cond_bed_to_bb

	output:
	file("${cond}.bb") into cond_bb_for_hub_ch
	file("chrom.sizes") into merge_chrom_sizes_ch

	publishDir "$outdir/${params.hubName}/$cond/", mode: "copy", pattern: "*.bb"

	script:
	"""
	samtools idxstats $bamfile | cut -f 1,2 | grep -v _ | grep -v ERCC > chrom.sizes
	dasatools.py convert $cond $bedfile temp.bed
	bedToBigBed -type=bed3+3 -tab temp.bed chrom.sizes ${cond}.bb
	"""
}

/* ** Convert differential bedgraph to BW ** */

Channel
	.from merge_chrom_sizes_ch
	.reduce("") { a, b -> a + " " + b }
	.set { all_chromsizes }

process CombineChromSizes {
	executor "local"
	
	input:
	val chromsizes from all_chromsizes

	output:
	file("chrom.sizes") into combined_chromsizes

	script:
	"""
	dasatools.py sizes $chromsizes > chrom.sizes
	"""
}

process BDGtoBW {
	memory "2G"
	time "1h"

	input:
	tuple label, diffbdg, sigbdg from bdg_to_bw_ch
	file chromsizes from combined_chromsizes

	output:
	file("${label}.diff.bw") into contr_bw_for_hub_ch
	file("${label}.sig.bw") into sig_bw_for_hub_ch

	publishDir "$outdir/${params.hubName}/$label/", mode: "copy", pattern: "*.bw"

	script:
	"""
	bedGraphToBigWig $diffbdg $chromsizes ${label}.diff.bw
	bedGraphToBigWig $sigbdg $chromsizes ${label}.sig.bw
	"""
}

/* ** Create WashuU hubs ** */

/* We need: 
  - per-sample bw and bb
  - per-condition bw and bb
  - per-contrast bw
*/
/* In the following channels we only take the last item to ensure work is done -
   we deduce the filenames from the inputs. */

Channel
	.from smp_bw_for_hub_ch
	.last()
	.set { smp_bw_for_hub }

Channel
	.from smp_bb_for_hub_ch
	.last()
	.set { smp_bb_for_hub }

Channel
	.from cond_bw_for_hub_ch
	.last()
	.set { cond_bw_for_hub }

Channel
	.from cond_bb_for_hub_ch
	.last()
	.set { cond_bb_for_hub }

Channel
	.from contr_bw_for_hub_ch
	.last()
	.set { contr_bw_for_hub }

Channel
	.from sig_bw_for_hub_ch
	.last()
	.set { sig_bw_for_hub }

process MakeHubs {
	executor "local"

	input:
	val d1 from smp_bw_for_hub
	val d2 from smp_bb_for_hub
	val d3 from cond_bw_for_hub
	val d4 from cond_bb_for_hub
	val d5 from contr_bw_for_hub
	val d6 from sig_bw_for_hub

	script:
	"""
	dasatools.py hubs ${workflow.launchDir}/${params.samples} ${workflow.launchDir}/${params.contrasts} \
          ${workflow.launchDir}/${params.reportDir}/${params.hubName} ${params.hubURL}/${params.hubName}
	"""
}

/* ** Draw tornado plots using deeptools. ** */
Channel
	.from contr_for_tornado_ch                            /* label, test, ctrl */
	.combine(combine_for_tornado_ch)                      /* label, test, ctrl, cond, bw */
	.filter( { it[3] == it[1] || it[3] == it[2] } )
	.map( row -> tuple(row[0], row[3], row[4]) )          /* label, cond, bw */
	.groupTuple()                                         /* label, [cond], [bw] */
	.join(regions_to_tornado_ch)                          /* label, [cond], [bw], upregions, dnregions */
	.into { regions_tornado_ch; tss_tornado_ch }

process RegionsTornado {
	time "2h"
	memory "2G"

	input:
	tuple label, samples, bigwigs, upregions, dnregions from regions_tornado_ch

	output:
	file("*.png")

	publishDir "$outdir/plots/$label/", mode: "copy", pattern: "*.png"

	script:
	"""
#module purge
#module load deeptools

function draw_tornado() {
  NAME=\$1
  REGIONS=\$2
  BW1=\$3
  BW2=\$4
  D=1000

  computeMatrix reference-point -p max -R \$REGIONS -S \$BW1 \$BW2 -o \${NAME}.mat.gz \
    --referencePoint center -b \$D -a \$D --skipZeros --missingDataAsZero
  plotHeatmap -m \${NAME}.mat.gz -o \${NAME}.png --heatmapWidth 6 \
    --samplesLabel ${samples[0]} ${samples[1]} \
    --xAxisLabel "distance (bp)" --refPointLabel "Peak" --yAxisLabel "Regions" \
    --regionsLabel "ATAC" \
    --colorList white,red --sortUsingSamples 1 --sortUsing mean
}
nup=\$(grep -c ^ $upregions)
ndn=\$(grep -c ^ $dnregions)
if [[ \$nup -gt 1 ]];
then
  draw_tornado ${label}.testpeaks $upregions ${bigwigs[0]} ${bigwigs[1]}
fi
if [[ \$ndn -gt 1 ]];
then
  draw_tornado ${label}.ctrlpeaks $dnregions ${bigwigs[0]} ${bigwigs[1]}
fi
	"""
}

process TSStornado {
	time "6h"
	memory "4G"

	input:
	tuple label, samples, bigwigs, upregions, dnregions from tss_tornado_ch
	file(tssfile) from tss_plots

	output:
	file("*.png")

	publishDir "$outdir/plots/$label/", mode: "copy", pattern: "*.png"

	script:
	"""
#module purge
#module load deeptools

function draw_tornado() {
  NAME=\$1
  REGIONS=\$2
  BW1=\$3
  BW2=\$4
  D=${params.tss_range}

  computeMatrix reference-point -p max -R \$REGIONS -S \$BW1 \$BW2 -o \${NAME}.mat0.gz \
    --referencePoint center \
    -b \$D -a \$D --skipZeros --missingDataAsZero
  orientMatrix.py \$REGIONS \${NAME}.mat0.gz \${NAME}.mat.gz
  plotHeatmap -m \${NAME}.mat.gz -o \${NAME}.png --heatmapWidth 6 \
    --samplesLabel ${samples[0]} ${samples[1]} \
    --xAxisLabel "distance (bp)" --refPointLabel "TSS" --yAxisLabel "Genes" \
    --regionsLabel "ATAC" \
    --colorList white,red --sortUsingSamples 1 --sortUsing mean
}

draw_tornado ${label}.TSS ${tssfile} ${bigwigs[0]} ${bigwigs[1]}
	"""
}

/* ** Size scatterplots ** */

Channel
	.from size_plot_ch /* label, sizes */
	.join(contr_for_scatterplot_ch) /* label, sizes, test, ctrl */
	.set { size_plot }

process SizeScatterplot {
	time "12h"

	input:
	tuple label, sizes, testCond, ctrlCond from size_plot

	output:
	file("${label}.scatterplot.png") into dummy_scatterplot

	publishDir "$outdir/plots/$label/", mode: "copy", pattern: "*.png"

	script:
	"""
	density_scatterplot.py -b -cx 8 -cy 4 -l10 -yl "$testCond" -xl "$ctrlCond" $sizes ${label}.scatterplot.png
	"""
}

/* ** Differential gene accessibility analysis - at TSS. ** */

process GeneMatrix {
	time "4h"
	cpus 20

	input:
	val samplefactors from factors_ch3
	file(tssfile) from tss_matrix

	output:
	tuple file("gene-matrix.txt"), file("levels.txt"), file("labels.txt") into genediff_ch

	"""
#!/bin/bash

BAMS=""
for bam in \$(cut -f 4 ${samplesfile}); do
  BAMS="\$BAMS ${workflow.launchDir}/\$bam"
done

split -l 2000 ${tssfile} genes-

for glist in genes-*;
do
  bedtools multicov -bams \$BAMS -bed \${glist} > \${glist}.matrix.txt &
done
wait
cat *.matrix.txt | dasatools.py gmatrix ${samplesfile} ${samplefactors} > gene-matrix.txt
"""
}

Channel
	.from genediff_ch
	.combine(contr_for_genediff)
	.set { genediff }

process GeneDiff {
	time "2h"
	executor "local"

	input:
	tuple genematrix, levelsfile, labelsfile, contr, testcond, ctrlcond from genediff

	output:
	tuple contr, file("genediff.csv") into extract_siggene_ch

	publishDir "$outdir/data/$contr/", mode: "copy", pattern: "genediff.csv", saveAs: { filename -> "${contr}-${filename}" }

	script:
	levels = file(levelsfile).text
	labels = file(labelsfile).text

	"""
#!/usr/bin/env Rscript

library("DESeq2")

datafile = "$genematrix"
message("Starting read csv")
counts = as.matrix(read.csv(datafile, sep='\t', row.names=1))
counts[counts>2147483647]=2147483647
message("Read csv done")
levels = c($levels)
labels = c($labels)
#message(1)
sampleTable = data.frame(condition=factor(levels[labels]))
rownames(sampleTable) = colnames(counts)
#message(2)
dds = DESeqDataSetFromMatrix(countData = counts, colData = sampleTable, design = ~ condition)
#message(2.5)
dds = estimateSizeFactors(dds)
#message(3)
keep = rowMeans(counts(dds)) >= 5
dds = dds[keep,]
dds = DESeq(dds)
#message(4)
res = results(dds, contrast=c("condition", "$testcond", "$ctrlcond"))
#message(5)
write.table(res, file="genediff.csv", sep='\t')
        """
}

process ExtractSignificantGenes {
	executor "local"

	input:
	tuple contr, genediff from extract_siggene_ch

	output:
	file("genediff-counts.txt") into combine_genediff_counts_ch
	file("${contr}.TSS.xlsx")

	publishDir "$outdir/data/$contr/", mode: "copy", pattern: "*.xlsx"

	"""
	#!/bin/bash
	dasatools.py gsig $genediff ${params.log2fc_genes} ${params.pvalue_genes}
	N1=\$(grep -c ^ sigpeaks.csv)
	N2=\$(grep -c ^ test-up.csv)
	N3=\$(grep -c ^ ctrl-up.csv)
	N1=\$((N1-1))
	N2=\$((N2-1))
	N3=\$((N3-1))
	echo -e "$contr\t\$N1\t\$N2\t\$N3" > genediff-counts.txt
	dasatools.py xlsx ${contr}.TSS.xlsx test-up.csv:Increased_in_test ctrl-up.csv:Increased_in_ctrl ${genediff}:All_genes
	"""
}

Channel
	.from combine_genediff_counts_ch
	.reduce("") { a, b -> a + " " + b }
	.set { combined_genediff_counts }

/* ** Differential gene accessibility analysis - at gene bodies. ** */

Channel
	.from gene_matrix
	.splitText(by: 1000, file: true)
	.set { genebody_chunks }

process GeneBodyMatrixChunk {
	time "6h"
	memory "2G"

	input:
	val gbchunk from genebody_chunks

	output:
	file("gene-counts.mat") into make_gene_matrix

	"""
#!/bin/bash

BAMS=""
for bam in \$(cut -f 4 ${samplesfile}); do
  BAMS="\$BAMS ${workflow.launchDir}/\$bam"
done

bedtools multicov -bams \$BAMS -bed ${gbchunk} > gene-counts.mat
"""
}

process GeneBodyMatrix {
	time "1h"

	input:
	val samplefactors from factors_ch4
	file '*.mat' from make_gene_matrix.collect()

	output:
	tuple file("genebody-counts.fullmat"), file("levels.txt"), file("labels.txt") into genebodydiff_ch

	"""
#!/bin/bash
cat *.mat | dasatools.py gmatrix ${samplesfile} ${samplefactors} > genebody-counts.fullmat
"""
}	

/*
process GeneBodyMatrix {
	time "12h"
	cpus 30

	input:
	val samplefactors from factors_ch4
	val genebodyfile from gene_matrix

	output:
	tuple file("genebody-matrix.txt"), file("levels.txt"), file("labels.txt") into genebodydiff_ch

	"""
#!/bin/bash

BAMS=""
for bam in \$(cut -f 4 ${samplesfile}); do
  BAMS="\$BAMS ${workflow.launchDir}/\$bam"
done

split -l 2000 ${genebodyfile} genes-

for glist in genes-*;
do
  bedtools multicov -bams \$BAMS -bed \${glist} > \${glist}.matrix.txt &
done
wait
cat *.matrix.txt | dasatools.py gmatrix ${samplesfile} ${samplefactors} > genebody-matrix.txt
"""
}
*/

Channel
	.from genebodydiff_ch
	.combine(contr_for_genebodydiff)
	.set { genebodydiff }

process GeneBodyDiff {
	time "2h"
	executor "local"

	input:
	tuple genematrix, levelsfile, labelsfile, contr, testcond, ctrlcond from genebodydiff

	output:
	tuple contr, file("genebodydiff.csv") into extract_siggenebody_ch

	publishDir "$outdir/data/$contr/", mode: "copy", pattern: "genebodydiff.csv", saveAs: { filename -> "${contr}-${filename}" }

	script:
	levels = file(levelsfile).text
	labels = file(labelsfile).text

	"""
#!/usr/bin/env Rscript

library("DESeq2")

datafile = "$genematrix"
message("Starting read csv")
counts = as.matrix(read.csv(datafile, sep='\t', row.names=1))
counts[counts>2147483647]=2147483647
message("Read csv done")
levels = c($levels)
labels = c($labels)
#message(1)
sampleTable = data.frame(condition=factor(levels[labels]))
rownames(sampleTable) = colnames(counts)
#message(2)
dds = DESeqDataSetFromMatrix(countData = counts, colData = sampleTable, design = ~ condition)
#message(2.5)
dds = estimateSizeFactors(dds)
#message(3)
keep = rowMeans(counts(dds)) >= 5
dds = dds[keep,]
dds = DESeq(dds)
#message(4)
res = results(dds, contrast=c("condition", "$testcond", "$ctrlcond"))
#message(5)
write.table(res, file="genebodydiff.csv", sep='\t')
        """
}

process ExtractSignificantGeneBodies {
	executor "local"

	input:
	tuple contr, genediff from extract_siggenebody_ch

	output:
	file("genebodydiff-counts.txt") into combine_genebodydiff_counts_ch
	file("${contr}.genes.xlsx")

	publishDir "$outdir/data/$contr/", mode: "copy", pattern: "*.xlsx"

	"""
	#!/bin/bash
	dasatools.py gsig $genediff ${params.log2fc_genes} ${params.pvalue_genes}
	N1=\$(grep -c ^ sigpeaks.csv)
	N2=\$(grep -c ^ test-up.csv)
	N3=\$(grep -c ^ ctrl-up.csv)
	N1=\$((N1-1))
	N2=\$((N2-1))
	N3=\$((N3-1))
	echo -e "$contr\t\$N1\t\$N2\t\$N3" > genebodydiff-counts.txt
	dasatools.py xlsx ${contr}.genes.xlsx test-up.csv:Increased_in_test ctrl-up.csv:Increased_in_ctrl ${genediff}:All_genes
	"""
}

Channel
	.from combine_genebodydiff_counts_ch
	.reduce("") { a, b -> a + " " + b }
	.set { combined_genebodydiff_counts }

/* ** Differential accessibility - at enhancers ** */

process EnhancerMatrix {
	time "6h"
	memory "2G"

	when:
	params.enhancersfile

	output:
	tuple file("enhancers-matrix.txt"), file("levels.txt"), file("labels.txt") into enhancersdiff_ch

	script:
	"""
#!/bin/bash

BAMS=""
for bam in \$(cut -f 4 ${samplesfile}); do
  BAMS="\$BAMS ${workflow.launchDir}/\$bam"
done

# Put enhancers file in same format as genes file
awk -F"\t" 'BEGIN {OFS=FS} { \$5=\$1"_"\$2"_"\$3"_"\$4; \$4="+";  print }' ${params.enhancersfile} > enhancers.txt
bedtools multicov -bams \$BAMS -bed enhancers.txt > enhancer-counts.mat
cat enhancer-counts.mat | dasatools.py gmatrix ${samplesfile} /dev/null > enhancers-matrix.txt
"""
}

Channel
	.from enhancersdiff_ch
	.combine(contr_for_enhancersdiff)
	.set { enhancersdiff }

process EnhancersDiff {
	time "2h"
	executor "local"

	input:
	tuple  enhmatrix, levelsfile, labelsfile, contr, testcond, ctrlcond from enhancersdiff

	output:
	tuple contr, file("enhancersdiff.csv") into extract_enhancers_ch

	publishDir "$outdir/data/$contr/", mode: "copy", pattern: "enhancersdiff.csv", saveAs: { filename -> "${contr}-${filename}" }

	script:
	levels = file(levelsfile).text
	labels = file(labelsfile).text
	"""
#!/usr/bin/env Rscript

library("DESeq2")

datafile = "$enhmatrix"
message("Starting read csv")
counts = as.matrix(read.csv(datafile, sep='\t', row.names=1))
counts[counts>2147483647]=2147483647
message("Read csv done")
levels = c($levels)
labels = c($labels)
#message(1)
sampleTable = data.frame(condition=factor(levels[labels]))
rownames(sampleTable) = colnames(counts)
#message(2)
dds = DESeqDataSetFromMatrix(countData = counts, colData = sampleTable, design = ~ condition)
#message(2.5)
dds = estimateSizeFactors(dds)
#message(3)
keep = rowMeans(counts(dds)) >= 5
dds = dds[keep,]
dds = DESeq(dds)
#message(4)
res = results(dds, contrast=c("condition", "$testcond", "$ctrlcond"))
#message(5)
write.table(res, file="enhancersdiff.csv", sep='\t')
        """
}

process ExtractSignificantEnhancers {
	executor "local"

	input:
	tuple contr, enhdiff from extract_enhancers_ch

	output:
	file("enhancersdiff-counts.txt") into combine_enhancersdiff_counts_ch
	file("${contr}.enhancers.xlsx")

	publishDir "$outdir/data/$contr/", mode: "copy", pattern: "*.xlsx"

	"""
	#!/bin/bash

	function convert() {
	  echo -e "Chrom\tStart\tEnd\tGene\tLog2(FC)\tP-value" > \$2
	  tail -n +2 \$1 | tr '_' '\t' >> \$2
	}

	function convert2() {
	  echo -e "Chrom\tStart\tEnd\tGene\tLog2(FC)\tP-value" > \$2
	  cut -f 1,3,7 \$1 | tail -n +2 | tr -d '"' | tr '_' '\t' >> \$2
	}

	dasatools.py gsig $enhdiff ${params.log2fc_genes} ${params.pvalue_genes}
	N1=\$(grep -c ^ sigpeaks.csv)
	N2=\$(grep -c ^ test-up.csv)
	N3=\$(grep -c ^ ctrl-up.csv)
	N1=\$((N1-1))
	N2=\$((N2-1))
	N3=\$((N3-1))
	echo -e "$contr\t\$N1\t\$N2\t\$N3" > enhancersdiff-counts.txt

	convert test-up.csv test-up.c.csv
	convert ctrl-up.csv ctrl-up.c.csv
	convert2 ${enhdiff} all.c.csv

	dasatools.py xlsx ${contr}.enhancers.xlsx test-up.c.csv:Increased_in_test ctrl-up.c.csv:Increased_in_ctrl all.c.csv:All_genes
	"""
}

Channel
	.from combine_enhancersdiff_counts_ch
	.reduce("") { a, b -> a + " " + b }
	.set { combined_enhancersdiff_counts }


/* ** Write final HTML report. ** */
Channel
	.from combine_cond_stats_ch
	.reduce("") { a, b -> a + " " + b }
	.set { combined_cond_stats }

Channel
	.from combine_contr_counts_ch
	.reduce("") { a, b -> a + " " + b }
	.set { combined_contr_counts }

process Report {
	executor "local"

	input:
	val frips from report_frips
	tuple factors, allstats, allreads from report_stats
	val condstats from combined_cond_stats
	val contrcounts from combined_contr_counts
	val diffgenecounts from combined_genediff_counts
	val diffgenebodycounts from combined_genebodydiff_counts
	val diffenhancerscounts from combined_enhancersdiff_counts

	output:
	file("all-contr-stats.txt")
	file("all-contr-counts.txt")
	file("all-genediff-counts.txt")
	file("all-genebodydiff-counts.txt")
	file("all-enhancersdiff-counts.txt")

	publishDir "$outdir/data/", mode: "copy"


	script:
	"""
	cat $condstats > all-contr-stats.txt
	cat $contrcounts > all-contr-counts.txt
	cat $diffgenecounts > all-genediff-counts.txt
	cat $diffgenebodycounts > all-genebodydiff-counts.txt
	cat $diffenhancerscounts > all-enhancersdiff-counts.txt

	#cp $frips $factors $allstats $allreads .
	URL="http://epigenomegateway.wustl.edu/browser/?genome=${params.hubOrganism}&datahub=${params.hubURL}/${params.hubName}/"
	cat > ${workflow.launchDir}/mkreport.sh <<EOF
	#!/bin/bash
	mkdir -p ${params.reportDir}
	pushd ${params.reportDir}
	${workflow.projectDir}/bin/dasatools.py report ${workflow.launchDir}/${params.samples} ${workflow.launchDir}/${params.contrasts} $outdir ${params.reportTemplate} ../nextflow.config 
	rm -f METADATA
	echo -e "DATE\t\$(date +%Y-%m-%dT%H:%M)" > METADATA
	echo -e "NAME\t${params.reportDir}" >> METADATA
	echo -e "PIPELINE\tDA" >> METADATA
	popd
	EOF
	"""
}

