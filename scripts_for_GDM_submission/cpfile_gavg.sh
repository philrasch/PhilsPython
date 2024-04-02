#!/bin/bash
#set -x
infile='/tmp/flname_gavg'
outdirc='/e3sm_prod/haruki/data_for_protocol_paper/CESM2/coupled/'
outdire='/e3sm_prod/haruki/data_for_protocol_paper/E3SMv2/coupled/'
flist=`cat $infile| sort| uniq | grep -v xxx | grep e3sm-reshaped`
for file in $flist
do
    command="cp -i  $file $outdire"
    echo $command
    eval $command
done
flist=`cat $infile| sort| uniq | grep -v xxx | grep cesm2`
for file in $flist
do
    command="cp -i  $file $outdirc"
    echo $command
    eval $command
done
