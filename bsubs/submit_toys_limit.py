#!/usr/bin/env python

import os
import sys
import commands


sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '.')))
from submit import BsubHandle

sys.path.insert(0, '/afs/cern.ch/work/x/xju/code/check_workspace')
from check_workspace.helper import get_limit_dic

from optparse import OptionParser
workdir = os.getcwd()

exe = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'run_toys.sh')

class submit_toys:
    def __init__(self, poi_name, ws_name, data_name, mass_dic, input_pattern):
        self.poi_name = poi_name
        self.ws_name = ws_name
        self.data_name = data_name

        self.mass_dic = mass_dic
        self.input_p = input_pattern

        self.job_type = "toys" ## "expected", "observed" 
        self.njobs = 300
        self.ntoys = 10

        self.bsub_handle = BsubHandle()
        self.bsub_handle.no_submit = True

    def submit(self):
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
                        opt = '\"-w={} -d={} -p={} -n={} -s={}\"'.format(
                            self.ws_name, self.data_name,
                            self.poi_name, self.ntoys, seed
                        )
                        run_cmd = " ".join([exe, out_name, self.job_type,
                                            input_ws, str(mass), str(poi), opt])
                        self.bsub_handle.submit(run_cmd)
            elif self.job_type == "expected":
                poi_input = '"'+",".join(map(lambda x: str(x), poi_values))+'"'
                for job in range(self.njobs):
                    seed = job*self.ntoys + 1
                    opt = '\"-w={} -d={} -p={} -n={} -s={}\"'.format(
                        self.ws_name, self.data_name,
                        self.poi_name, self.ntoys, seed
                    )
                    run_cmd = " ".join([exe, out_name, self.job_type,
                                        input_ws, str(mass), str(poi_input), opt])
                    self.bsub_handle.submit(run_cmd)
            elif self.job_type == "observed":
                poi_input = '"'+",".join(map(lambda x: str(x), poi_values))+'"'
                opt = '\"-w={} -d={} -p={}\"'.format(
                    self.ws_name, self.data_name,
                    self.poi_name
                )
                run_cmd = " ".join([exe, out_name, self.job_type,
                                    input_ws, str(mass), str(poi_input), opt])
                self.bsub_handle.submit(run_cmd)
            else:
                pass

        self.bsub_handle.print_summary()

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

    (options,args) = parser.parse_args()
    if len(args) < 2:
        parser.print_help()
        exit(1)

    poi_name = options.poi_name
    ws_name = options.ws_name
    data_name = options.data_name
    if False:
        mass_dic  = get_limit_dic(options.asymp)
        formated_mass_dic = {}
        for x,y in mass_dic.iteritems():
            new_y = map(lambda x: round(x, 4), y)
            new_y = new_y[:-1]
            formated_mass_dic[x] = new_y
    else:
        formated_mass_dic = {
            2000: [0.0369, 0.022]
        }

    #print formated_mass_dic
    #input_p = "/afs/cern.ch/user/x/xju/work/h4l/highmass/workspaces/Graviton/HZZllvv_Graviton{0}_combined_HZZllvv_Graviton{0}_model_forLimits.root"
    input_p = args[1]
    job = submit_toys(poi_name, ws_name, data_name, formated_mass_dic, input_p)
    job.njobs = options.n_jobs
    job.ntoys = options.n_toys
    job.job_type = args[0]
    job.bsub_handle.no_submit = False
    job.submit()
