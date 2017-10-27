#!/bin/bash
out_dir=$1
input_name=$2
poi_name=$3
poi_value="$4"
mass=$5
ntoys=$6
seed=$7
ws_name=$8
data_name=$9

source /afs/cern.ch/user/x/xju/work/h4l/h4lcode/workspaces/setup.sh
which gcc
which root
which Expected_qmu

out_name="expected_q_seed${seed}.root"
Expected_qmu $input_name $out_name $poi_name "$poi_value" $mass $ws_name $data_name $ntoys $seed

if [ ! -d ${out_dir} ];then
    mkdir -vp ${out_dir}
fi
echo "saved to ${out_dir}"
cp $out_name ${out_dir}/
