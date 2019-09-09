#!/bin/bash
input_name=$1
ws_name=$2
mu_name=$3
data_name=$4
fix_vars=$5
cal_opt=$6 ## limit,pvalue
data_opt=$7 ## obs,exp
out_dir=$8
strategy=$9
fixOther=${10}
#limit inputs/2015_Graviton_histfactory_EKEI_v6.root combWS xs combDatabinned mG:500,GkM:0.01
source /afs/cern.ch/user/x/xju/work/h4l/h4lcode/workspaces/setup.sh
which gcc
which root
which get_stats
echo $cal_opt 
echo $data_opt
echo $fixOther

get_stats $input_name $ws_name $mu_name $data_name ModelConfig $strategy "$fix_vars" "$cal_opt" "$data_opt" $fixOther >& fit.log

if [ ! -d ${out_dir} ];then
    mkdir -vp ${out_dir}
fi
echo "save outputs to ${out_dir}"
out_name=`echo $fix_vars | sed 's/:/_/g' | sed 's/,/_/g'`
echo "out name: $out_name"
cp fit.log ${out_dir}/${out_name}.log
cp stats_results.txt ${out_dir}/${out_name}.txt
