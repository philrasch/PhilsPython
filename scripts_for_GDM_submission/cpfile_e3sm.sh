#!/bin/bash
#set -x
infile='/tmp/flname_e3sm'
outdir='/e3sm_prod/phil/yyy/e3sm_prod/haruki/data_for_protocol_paper/UKESM1/'
outdir='/e3sm_prod/haruki/data_for_protocol_paper/UKESM1/'
outdir='/e3sm_prod/haruki/data_for_protocol_paper/CESM2/'
outdir='/e3sm_prod/haruki/data_for_protocol_paper/E3SMv2/'
flist=`cat $infile| sort| uniq | grep -v xxx `
for file in $flist
do
    command="cp -i  $file $outdir"
    echo $command
    eval $command
done
