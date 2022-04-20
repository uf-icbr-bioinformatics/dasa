#!/bin/bash

SAMPLES=$1
CONTRASTS=$2
CONFIG=$3

if [[ "$SAMPLES" == "" ]];
then
  echo "Usage: dasa.sh SAMPLES CONTRASTS [CONFIG]"
  echo 
  echo "If CONFIG file is specified, it should be in Nexflow config"
  echo "format (param = value) and can contain the following parameters:"
  echo 
  echo "tssfile = path to file containing TSS coordinates for all genes"
  echo "hubURL = URL where this genome browser hub will be hosted"
  echo "hubName = name of this genome browser hub"
  echo "hubOrganism = organism identified for this hub (e.g. hg38, mm10)"
  exit 1
fi

if [[ "$CONFIG" != "" ]];
then
  conf="-c $CONFIG"
else
  conf=""
fi

nextflow run $conf -resume \
  -with-singularity /apps/dibig_tools/dasa/container/dasa.img \
  /apps/dibig_tools/dasa/dasa.nf \
  --samples $SAMPLES \
  --contrasts $CONTRASTS

#  --hubName $NAME \
#  --hubOrganism hg38 \
#  --hubURL http://lichtlab.cancer.ufl.edu/reports/${HUBDIR}/${NAME}/ \
#  --tssfile /apps/dibig_tools/dasa/data/GRCh38-all-tss.csv
