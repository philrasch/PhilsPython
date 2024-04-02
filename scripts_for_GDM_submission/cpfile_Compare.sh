#!/bin/bash
#set -x
infile='/tmp/flname'
outdirc='/e3sm_prod/haruki/data_for_protocol_paper/CESM2/coupled/'
outdire='/e3sm_prod/haruki/data_for_protocol_paper/E3SMv2/coupled/'
outdiru='/e3sm_prod/haruki/data_for_protocol_paper/UKESM1/'
outdiro='/e3sm_prod/haruki/data_for_protocol_paper/'

flist=`cat $infile| sort| uniq | grep -v xxx | grep PJR`
echo other
for file in $flist
do
    command="cp -i  $file $outdiro"
    echo $command
    eval $command
done

echo e3sm
flist=`cat $infile| sort| uniq | grep -v xxx | grep e3sm-reshaped`
for file in $flist
do
    command="cp -i  $file $outdire"
    echo $command
    eval $command
done

echo cesm
flist=`cat $infile| sort| uniq | grep -v xxx | grep cesm2`
for file in $flist
do
    command="cp -i  $file $outdirc"
    echo $command
    eval $command
done

echo ukesm
flist=`cat $infile| sort| uniq | grep -v xxx | grep UKESM1`
for file in $flist
do
    command="cp -i  $file $outdiru"
    echo $command
    eval $command
done
