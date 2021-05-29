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
params.merge_mode = "I"
params.log2fc = 1.0
params.pvalue = 0.01
params.reportDir = "Report"
params.deeptools = true
params.washu = true
params.tssfile = false

/* For WashU hub creation */
params.hubURL = "http://lichtlab.cancer.ufl.edu/reports/"
params.hubName = "hub"
params.hubOrganism = "hg38"

/* HTML report */
params.reportName = "ATAC-Seq Differential Analyisis"
if (workflow.containerEngine) {
   log.info """Running in ${workflow.containerEngine} container ${workflow.container}.\n"""
   params.reportTemplate = "/usr/local/share/dasa/report-template.html"
} else {
   params.reportTemplate = "${workflow.projectDir}/bin/report-template.html"
}
log.info """Report template: ${params.reportTemplate}"""

/* Internal variables */

samplesfile = file(params.samples, checkIfExists: true)
contrastsfile = file(params.contrasts, checkIfExists: true)
outdir = "${workflow.launchDir}/${params.reportDir}/"

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

  --merge_mode      Determines which area of partially-overlapping peaks to use in differential 
                    analysis. Can be either "I" (intersecting part is used) or "U" (the union of
                    the two peaks is used). Default: I.

  --tssfile         A file containing locations of Transcription Start Sites (optional).

Report options:

  --reportDir       Name of final report directory. Default: "Report".

  --reportName      Title of final report. Default: "ATAC-Seq Differential Analysis"

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
	.into { contrasts_ch; contrasts_ch_2; contrasts_ch_3; contrasts_ch4; contr_for_tornado_ch; contr_for_scatterplot_ch } /* TEST.vs.CTRL, TEST, CTRL */

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
	file("${smp}-convert-stats.txt") into combine_convert_stats /* convert-stats contains: totopen, npeaks */
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

/* ** Determine number of reads for each BAM file ** */
process CountReads {
	executor "local"

	input:
	tuple smp, condName, bamfile from smp_count_reads

	output:
	stdout into max_nreads_ch
	tuple condName, stdout into cond_totreads_ch
	file("${smp}-nreads.txt") into combine_nreads

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
	.reduce("") { a, b -> a + " " + b }
	.set { combined_stats }

Channel
	.from combine_nreads
	.reduce("") { a, b -> a + " " + b }
	.set { combined_nreads }

process ComputeFactors {
	executor "local"

	input:
	val allstats from combined_stats
	val allreads from combined_nreads

	output:
	file("sample-factors.txt") into factors_ch, factors_ch2
	file("all-sample-reads.txt") into nreads_for_frip_ch
	tuple file("sample-factors.txt"), file("all-sample-stats.txt"), file("all-sample-reads.txt") into report_stats

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
levels = c("$testcond", "$ctrlcond")
labels = c($labels)
sampleTable = data.frame(condition=factor(levels[labels]))
rownames(sampleTable) = colnames(counts)
dds = DESeqDataSetFromMatrix(countData = counts, colData = sampleTable, design = ~ condition)
dds = estimateSizeFactors(dds)
keep = rowSums(counts(dds)) >= 5
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
	tuple label, file("sigpeaks.bedGraph") into bdg_to_bw_ch
	tuple label, file("test-up.csv"), file("ctrl-up.csv") into regions_to_tornado_ch
	file("sigpeaks.csv") into dummy_sig
	file("contr-counts.txt") into combine_contr_counts_ch

	publishDir "$outdir/data/$label/", mode: "copy", pattern: "*.csv", saveAs: { filename -> "${label}-${filename}" }

	"""
	# generate sigpeaks.csv, sigpeaks.bedGraph, test-up.csv, ctrl-up.csv
	dasatools.py sig $diffpeaks ${params.log2fc} ${params.pvalue}
	N1=\$(grep -c ^ sigpeaks.csv)
	N2=\$(grep -c ^ test-up.csv)
	N3=\$(grep -c ^ ctrl-up.csv)
	echo -e "$label\t\$N1\t\$N2\t\$N3" > contr-counts.txt
	"""
}

/* ** Compute FRIP score for each sample ** */

process FRIP {
	memory "2G"

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
	file("${smp}.bb") into smp_bb_for_hub_ch

	publishDir "$outdir/${params.hubName}/$smp/", mode: "copy", pattern: "*.bb"

	script:
	"""
	samtools idxstats $bamfile | cut -f 1,2 | grep -v _ | grep -v ERCC > chrom.sizes
	bedToBigBed -type=bed3+3 -tab $bedfile chrom.sizes ${smp}.bb
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
	bedToBigBed -type=bed3 -tab $bedfile chrom.sizes ${cond}.bb
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
	tuple label, diffbdg from bdg_to_bw_ch
	file chromsizes from combined_chromsizes

	output:
	file("${label}.bw") into contr_bw_for_hub_ch

	publishDir "$outdir/${params.hubName}/$label/", mode: "copy", pattern: "*.bw"

	script:
	"""
	bedGraphToBigWig $diffbdg $chromsizes ${label}.bw
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

process MakeHubs {
	executor "local"

	input:
	val d1 from smp_bw_for_hub
	val d2 from smp_bb_for_hub
	val d3 from cond_bw_for_hub
	val d4 from cond_bb_for_hub
	val d5 from contr_bw_for_hub

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
	file("*.png") into dummy_png

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
    -b \$D -a \$D --skipZeros --missingDataAsZero
  plotHeatmap -m \${NAME}.mat.gz -o \${NAME}.png --heatmapWidth 6 \
    --samplesLabel ${samples[0]} ${samples[1]} \
    --xAxisLabel "distance (bp)" --refPointLabel "Peak" --yAxisLabel "Regions" \
    --regionsLabel "ATAC" \
    --colorList white,red --sortUsingSamples 1 --sortUsing mean
}
draw_tornado ${label}.testpeaks $upregions ${bigwigs[0]} ${bigwigs[1]}
draw_tornado ${label}.ctrlpeaks $dnregions ${bigwigs[0]} ${bigwigs[1]}
	"""
}

process TSStornado {
	time "2h"
	memory "2G"

	input:
	tuple label, samples, bigwigs, upregions, dnregions from tss_tornado_ch

	output:
	file("*.png") into dummy_png_2

	publishDir "$outdir/plots/$label/", mode: "copy", pattern: "*.png"

	when:
	params.tssfile

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

draw_tornado ${label}.TSS ${params.tssfile} ${bigwigs[0]} ${bigwigs[1]}
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

	script:
	"""
	cat $condstats > all-contr-stats.txt
	cat $contrcounts > all-contr-counts.txt
	cp $frips $factors $allstats $allreads .
	URL="http://epigenomegateway.wustl.edu/browser/?genome=${params.hubOrganism}&datahub=${params.hubURL}/${params.hubName}/"
	dasatools.py report ${workflow.launchDir}/${params.samples} ${workflow.launchDir}/${params.contrasts} \
          $outdir ${params.reportTemplate} "${params.reportName}" \$URL
	"""
}

/*
workflow.onComplete {
  cmdline = ["${workflow.projectDir}/bin/dasatools.py", "report", "$outdir"]
  cmdline.execute()
}
*/