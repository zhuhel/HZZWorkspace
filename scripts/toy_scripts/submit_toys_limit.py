#!/usr/bin/env python

import os
import sys
import commands
import datetime
import shutil

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '.')))

from optparse import OptionParser
workdir = os.getcwd()

exe = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'run_toys.sh')

def get_mass(input_str):
    if not ',' in input_str:
      return input_str

    items = input_str.split(',')
    mass_str = 0
    for item in items:
        if "mH" in item:
            mass_str = item
            break

    if ":" in mass_str:
        return mass_str.split(':')[1]
    elif "=" in mass_str:
        return mass_str.split('=')[1]
    else:
        return mass_str

def get_limit_dic(file_name):
    limit_dic = {}
    with open(file_name, 'r') as f:
        for line in f:
            if line[0] == "#":
                continue
            items = line.split()
            mass_ = int(get_mass(items[0]))
            start_i = 3
            w_xs = 1.
            one_entry = [
                float(items[start_i])*w_xs, # -2
                float(items[start_i+1])*w_xs, # -1
                float(items[start_i+2])*w_xs, # median
                float(items[start_i+3])*w_xs, # +1
                float(items[start_i+4])*w_xs, # +2
                float(items[start_i+5])*w_xs, # obs
            ]
            limit_dic[mass_] = one_entry
    return limit_dic

class submit_toys:
    def __init__(self, poi_name, ws_name, data_name, mass_dic, input_pattern, job_type):
        self.poi_name = poi_name
        self.ws_name = ws_name
        self.data_name = data_name

        self.mass_dic = mass_dic
        self.input_p = input_pattern

        self.job_type = job_type #"toys", "expected", "observed" 
        self.njobs = 300
        self.ntoys = 10
   
        self.subfile="condorSubmit.tmp.sub"
        self.exefile="condorRunFile.tmp.sh"
        self.subinputfile="condorInputFile.tmp.txt"
       
        self.islxplus=('lxplus' in os.environ["HOSTNAME"])
        self.isbnl   =('bnl' in os.environ["HOSTNAME"])

    def backup_file(self, afile):
        """
        backup file based on its time stamp
        """
 
        # get the timestamp
        st = os.stat(afile)
        ts = datetime.datetime.fromtimestamp(st.st_mtime).strftime('%Y-%m-%d-%H-%M')
        nfile = afile+"_"+ts
        print("Info=> backup file {} -> {}".format(afile, nfile))
        shutil.move( afile, nfile )

    def write_subfile(self, queue="longlunch"):

        if os.path.isfile(self.subfile): self.backup_file(self.subfile)

        jfb = ""
        #+JobFlavour   = \"microcentury\" \n\
        #+JobFlavour   = \"longlunch\" \n\
        #+JobFlavour   = \"workday\" \n\
        #+JobFlavour   = \"tomorrow\" \n\
        #+JobFlavour   = \"testmatch\" \n\
        if self.islxplus:
            jfb="+JobFlavour   = \"{}\"".format(queue)
	elif self.isbnl:
	    jfb="+Experiment     = \"atlas\" \naccounting_group = group_atlas.bnl"

        submitfile=open(self.subfile,"w")
        submitfile.write("""Universe      = vanilla \n\
GetEnv          = True \n\
executable    = {0} \n\
output        = log/log.$(ClusterId).$(ProcId).out \n\
error         = log/log.$(ClusterId).$(ProcId).err \n\
log           = log/log.$(ClusterId).$(ProcId).log \n\
{1}
\n\
queue arguments from {2}""".format(self.exefile, jfb, self.subinputfile))
        submitfile.close()

    def write_exefile(self, setup_script):

        if os.path.isfile(self.exefile): self.backup_file(self.exefile)

        runfile_name=self.exefile

        runfile=open(runfile_name,"w")
        runfile.write("""#!/bin/bash \n\
ls *\n\
echo \"Running in: $PWD  @ $HOSTNAME\"
echo \"Time : \" $(date -u)
  
uname -a
ulimit -S -s 20000
ulimit -Sc 0
ulimit -Hc 0
ulimit -c 0
ulimit -d unlimited
ulimit -f unlimited
ulimit -l unlimited
ulimit -n unlimited
ulimit -s unlimited
ulimit -t unlimited

echo \"Current folder is\"
pwd
ls -l

out_dir=$1
action=$2
input_name=$3
mass=$4
poi_value=$5
ws_name=$6
data_name=$7
poi_name=$8
ntoys=$9
seed=${{10}}

echo \"Setting up environment...\"
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${{ATLAS_LOCAL_ROOT_BASE}}/user/atlasLocalSetup.sh

source {0}

which gcc
which root
which mainToys

mainToys $action $input_name $mass "$poi_value" -w=${{ws_name}} -d=${{data_name}} -p=${{poi_name}} -n=${{ntoys}} -s=${{seed}}

if [ ! -d $out_dir ];then
    mkdir -vp $out_dir
fi
echo "saved to $out_dir"
cp *root $out_dir/
  
ls *\n\
""".format(setup_script))
        runfile.close()

        os.system("chmod +x "+runfile_name)

    def write_subinputfile(self):

        if os.path.isfile(self.subinputfile): self.backup_file(self.subinputfile)

        inputfile=open(self.subinputfile,"w")

        out_name =  os.getcwd()
        for mass, poi_values in self.mass_dic.iteritems():
            print "mass: ", mass
            input_ws = self.input_p.format(mass)

            print "Input:", input_ws
            if self.job_type == "toys":
                # throw toys for the NLL_SB and NLL_bonly
                for poi in poi_values:
                    for job in range(self.njobs):
                        seed = job*self.ntoys + 1
                        run_cmd = " ".join([out_name, self.job_type,
                                            input_ws, str(mass), str(poi), self.ws_name, self.data_name, self.poi_name, str(self.ntoys), str(seed)])
                        inputfile.write("{}\n".format( run_cmd ))
            elif self.job_type == "expected":
                poi_input = ",".join(map(lambda x: str(x), poi_values))
                for job in range(self.njobs):
                    seed = job*self.ntoys + 1
                    run_cmd = " ".join([out_name, self.job_type,
                                        input_ws, str(mass), str(poi_input), self.ws_name, self.data_name, self.poi_name, str(self.ntoys), str(seed)])
                    inputfile.write("{}\n".format( run_cmd ))
            elif self.job_type == "observed":
                poi_input = ",".join(map(lambda x: str(x), poi_values))
                run_cmd = " ".join([out_name, self.job_type,
                                    input_ws, str(mass), str(poi_input), self.ws_name, self.data_name, self.poi_name, str(1), str(1)])
                inputfile.write("{}\n".format( run_cmd ))

            else:
                pass

        inputfile.close()

    def submit_condor(self):

        os.system("mkdir -p log")
        os.system("condor_submit "+self.subfile)


if __name__ == "__main__":
    usage = "%prog action ws_pattern "
    version="%prog 1.0"
    parser = OptionParser(usage=usage, description="submit jobs for toys", version=version)
    parser.add_option('-w', dest="ws_name", help="workspace name", default="combined")
    parser.add_option('-d', dest="data_name", help="dataset name", default="obsData")
    parser.add_option('-p', dest="poi_name", help="POI name", default="mu_ggF")
    parser.add_option('-n', dest="n_toys", help="number of toys", type='int', default=1)
    parser.add_option('-j', dest="n_jobs", help="number of jobs", type='int', default=1)
    parser.add_option('-a', dest="asymp", help="inputs from asymptotic", default="limit_asym.txt")
    parser.add_option('-t', dest="job_type", help="job type: toys, expected, observed", default="toys")

    (options,args) = parser.parse_args()
    if len(args) < 2:
        parser.print_help()
        exit(1)

    job_type = options.job_type
    poi_name = options.poi_name
    ws_name = options.ws_name
    data_name = options.data_name
    if os.path.isfile(options.asymp):
        mass_dic  = get_limit_dic(options.asymp)
        formated_mass_dic = {}
        for x,y in mass_dic.iteritems():
            new_y = map(lambda x: round(x, 4), y)
            #new_y = new_y[:-1]
            formated_mass_dic[x] = new_y
    else:
        formated_mass_dic = {
            2000: [0.0369, 0.022]
        }

    ## workspace file
    input_p = args[0]
    ## HZZws setup script
    setup_script = args[1] 

    job = submit_toys(poi_name, ws_name, data_name, formated_mass_dic, input_p, job_type)
    job.njobs = options.n_jobs
    job.ntoys = options.n_toys
    job.write_subfile()
    job.write_exefile(setup_script)
    job.write_subinputfile()
    job.submit_condor()
