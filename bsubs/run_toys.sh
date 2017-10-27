#!/bin/bash
out_dir=$1
action=$2
input_name=$3
mass=$4
poi_value=$5
options=$6

source /afs/cern.ch/user/x/xju/work/h4l/h4lcode/workspaces/setup.sh
which gcc
which root
which mainToys

mainToys $action $input_name $mass "$poi_value" $options

if [ ! -d ${out_dir} ];then
    mkdir -vp ${out_dir}
fi
echo "saved to ${out_dir}"
cp *root ${out_dir}/
