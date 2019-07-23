# dasa
Differential ATAC-Seq analysis

DASA is a complete tool for differential ATAC-Seq analysis. Starting from the output
of a peak caller like MACS2, it combines peaks from replicates of the same condition,
and it finds a set of comparable peaks between the two conditions. These peaks are then
quantified, using appropriate normalizations, and compared using DESeq2. The result
is a set of differential peaks, with an associated fold change and P-value.

DASA generates a full report of its analysis, including links to Excel files, plots, 
and hubs for the WashU epigenome browser. 

NOTE: this is a pre-release of DASA. More documentation and information on its dependencies
is coming soon. We will also provide DASA as a Singularity container for ease of
installation.

