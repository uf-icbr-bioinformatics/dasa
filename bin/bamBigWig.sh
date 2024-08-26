#!/bin/bash

BAM=$1
BW=$2
factor=$3

TMPD=$(mktemp -d -p .)
mkdir -p $TMPD

samtools idxstats $BAM | cut -f 1,2 | grep -v \* | grep -v ERCC > $TMPD/chrom.sizes

BDGS=""

echo $(date --rfc-3339=seconds) Converting chroms to bedGraph...

for chrom in $(cut -f 1 $TMPD/chrom.sizes);
do 
  BDG=$TMPD/$chrom.bedGraph
  BDGS="$BDGS $BDG"
  depthToBedGraph.py $BAM -c $chrom -o $BDG -s $scale &
done

wait

echo $(date --rfc-3339=seconds) Concatenating bedGraphs...

cat $BDGS > $TMPD/all.bedGraph

echo $(date --rfc-3339=seconds) Converting to bigWig...

bedGraphToBigWig $TMPD/all.bedGraph $TMPD/chrom.sizes $BW

if [ "$?" == "0" ];
then
  rm -fr $TMPD
fi

echo $(date --rfc-3339=seconds) Done.
