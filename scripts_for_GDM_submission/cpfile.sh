#!/bin/bash
#set -x
infile='/tmp/flname'
outdir='/e3sm_prod/phil/yyy/e3sm_prod/haruki/data_for_protocol_paper/UKESM1/'
outdir='/e3sm_prod/haruki/data_for_protocol_paper/UKESM1/'
flist=`cat $infile| sort| uniq | grep -v xxx `
for file in $flist
do
    command="cp -i  $file $outdir"
    echo $command
    eval $command
done
