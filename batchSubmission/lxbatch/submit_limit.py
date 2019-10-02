#!/usr/bin/env python

import os
import sys
import commands

workdir = os.getcwd()

#input_ws = "/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2016/Workspaces/HighMass/Prod_v04/2016_07_11/combined04_LWA_b.root"
#input_ws = "/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2016/Workspaces/HighMass/Prod_v04/2016_07_11/combined04_NWA_b.root"
input_ws = "/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2016/Workspaces/HighMass/Prod_v04/2016_07_15/HM_NWA_v05.root"
#input_ws = "/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2016/Workspaces/HighMass/Prod_v04/2016_07_15/HM_LWA_v05.root"
ana_name = "VBF" #VBF, ggF
frac_width = 0.0001 # need to change for LWA
submit = True
fix_other_poi = "0:0"

mG_low = 200
mG_hi = 1000
mG_step = 5

exe = "/afs/cern.ch/user/x/xju/work/h4l/h4lcode/workspaces/bsubs/run_limit.sh"
data_opt = "obs,exp" #obs, exp
cal_opt = "pvalue,limit" # limit,pvalue
ws_name = "combined"
mu_name = "XS_"+ana_name
out_name = workdir
do_LWA = "LWA" in input_ws
width_name = "gamma" ## natural width

data_name = "obsData"
if not do_LWA:
    out_name += "/NWA_"+ana_name+"/"
else:
    out_name += "/LWA_"+ana_name+"_"+width_name+"_"+str(frac_width)+"/"

goodjobs = []
badjobs = []
print out_name
for mG in range(mG_low, mG_hi+mG_step, mG_step):
    if not do_LWA:
        fix_vars = "mH:"+str(mG)
    else:
        fix_vars = "mH:"+str(mG)+","+width_name+":"+str(mG*frac_width)

    run_cmd = exe+" "+input_ws+" "+ws_name+" "+mu_name+" "+\
            data_name+" "+fix_vars+" "+cal_opt+" "+data_opt+" "+\
            out_name+" 1 "+fix_other_poi
    if not submit: print run_cmd
    #-G u_zp -q 8nh for atlas sources
    #-G ATLASWISC_GEN -q wisc for wisconsin sources
    bsubs_cmd = "bsub -q wisc -R 'pool>4000' -C 0 -o" + \
            workdir+ "/output "+ run_cmd
    if submit:
        status,output=commands.getstatusoutput(bsubs_cmd)
    else:
        continue
    if status != 0:
        print output
        badjobs.append(0)
    else:
        goodjobs.append(1)

print "Good jobs: "+ str(len(goodjobs))+", "+str(len(badjobs))+" failed!"
