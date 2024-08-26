params.tss_range = 2000
params.gene_body_upstream = 1000
params.gene_body_downstream = 1000

/* Check that all input files exist, convert peaks to BED */
process InitializeSamples {
	executor "local"

	input: 
	tuple val(smp), val(condName), val(xlsfile), val(bamfile)

	output:
	tuple val(smp), val(condName), val(bamfile), emit: sample_files      /* into smp_count_reads */
	tuple val(condName), path("${smp}.bed"), emit: cond_beds          /* into merge_bed_ch */
	tuple val(smp), path("${smp}-convert-stats.txt"), emit: combine_convert_stats /* convert-stats contains: totopen, npeaks */
	tuple val(smp), path("${smp}.bed"), val(bamfile) /* into frip_ch, sample_bed_to_bb */
	tuple val(smp), val(bamfile)                     /* into sample_bam_to_bw_ch */

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
	path("tss-regions.txt") /* into tss_plots, tss_matrix */

	script:
	"""
	#!/bin/bash

	dasatools.py regions ${params.genesfile} tss-regions.txt t ${params.tss_range} ${params.tss_range}
	"""
}

process GeneRegions {
	executor "local"

	output:
	path("gene-regions.txt") /* into gene_matrix */

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
	tuple val(smp), val(condName), val(bamfile)

	output:
	stdout                                    /* into max_nreads_ch */
	tuple val(condName), stdout, emit: cond_totreads_ch               /* into cond_totreads_ch */
	tuple val(smp), file("${smp}-nreads.txt"), emit: combined_nreads

	script:
	"""
	NREADS=\$(samtools idxstats $bamfile | awk '{sum += \$3} END {print sum}')
	echo -e "${smp}\t\$NREADS" > ${smp}-nreads.txt
	echo \$NREADS
	"""
}

process MergeBEDs {
	executor "local"

	input:
	tuple val(cond), val(bedfiles) /* from merge_bed */

	output:
	tuple val(cond), path("${cond}.bed")           /* into common_peaks_ch, post_merge_beds */
	tuple val(cond), path("${cond}-bed-stats.txt") /* into cond_bed_stats_ch */

	script:
	"""
	out=${cond}.bed
	sort -m -k1,1 -k2,2n -k3,3n $bedfiles | mergeBed -i - > \$out
	NPEAKS=\$(grep -c ^ \$out)
	TOTOPEN=\$(awk '{sum += (\$3-\$2)} END {print sum}' \$out)
	echo -e "$cond\t\$NPEAKS\t\$TOTOPEN" > ${cond}-bed-stats.txt
	"""
}

process ComputeFactors {
	executor "local"

	input:
	val allstats /* from combined_stats */
	val allreads /* from combined_nreads */

	output:
	path("sample-factors.txt")      /* into factors_ch, factors_ch2, factors_ch3, factors_ch4 */
	path("all-sample-reads.txt")    /* into nreads_for_frip_ch */
	tuple file("sample-factors.txt"), file("all-sample-stats.txt"), file("all-sample-reads.txt") /* into report_stats */

	publishDir "$outdir/data", mode: "copy"

	script:
	"""
	cat $allstats > all-sample-stats.txt
	cat $allreads > all-sample-reads.txt
	dasatools.py factors all-sample-stats.txt all-sample-reads.txt > sample-factors.txt
	"""
}

process PreCondFactors {
	executor "local"

	input:
	tuple cond, totalreads, statsfile /* from cond_factors */

	output:
	tuple file("cond-nreads.txt"), statsfile /* into cond_factors_2 */

	script:
	"""
	echo -e "$cond\t$totalreads" > cond-nreads.txt
	"""
}

process MergeBAMs {
	time "12h"
	memory "5G"

	input:
	tuple val(cond), val(bamfiles) /* from merge_bams_ch */

	output:
	tuple val(cond), path("${cond}.bam") /* into cond_bam_to_bw_ch, post_merge_bams */

	script:
	"""
	samtools merge ${cond}.bam $bamfiles
	samtools index ${cond}.bam
	"""
}

