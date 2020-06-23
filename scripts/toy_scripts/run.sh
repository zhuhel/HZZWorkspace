

ws=/usatlas/groups/bnl_local/hzhu/HZZAna/statanalysis/HZZtoys/workspace_combination_cut/combined_afterPara_mH1500.root
#ws=/usatlas/groups/bnl_local/hzhu/HZZAna/statanalysis/HZZtoys/llvv_ws/ws_NWA_fullrunII_allSyst_ZZEWCorr_mH1500.root
setup=/usatlas/groups/bnl_local/hzhu/HZZAna/CMakeWS/setup.sh

#wsname=combined
wsname=combWS
dataname=combData
poi=mu_ggF

#asylimits=/usatlas/groups/bnl_local/hzhu/HZZAna/statanalysis/HZZtoys/llvv_ws/llvv_ggF_limits_test.txt
asylimits=/usatlas/groups/bnl_local/hzhu/HZZAna/statanalysis/HZZtoys/workspace_combination_cut/comb_ggF_limits_test.txt
job_type=expected

python submit_toys_limit.py ${ws} ${setup} -w ${wsname} -d ${dataname} -p ${poi} -n 10 -j 500 -a ${asylimits} -t ${job_type}
