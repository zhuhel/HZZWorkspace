#!/usr/bin/env python
import subprocess
import glob
import os
from  sets import Set
import shutil
import re
import sys
from optparse import OptionParser

import ROOT

import high_mass
import low_mass
import monoH
import common
import dump_yields


class PrepareWS:
    def __init__(self, info):
        print "preparing ",info[0],"workspace"
        self.wsdir="workspaces/"
        if not os.path.isdir(self.wsdir):
            os.mkdir(self.wsdir)
        self.ws_type, self.config = info 
        self.sample_dict = common.samples
        if (self.config == monoH):
            print "You are doing monoH!"
            self.sample_dict = common.samples_monoH
        
        #########################
        # smooth histograms
        #########################
        #"/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2015/MiniTrees/Test_v01/mc/Nominal/"
        self.rho = 1.0

        self.obs_name = "m4l"
        self.no_tau = True
        self.out_prefix = "test"


        # key: mc_channel_number [str]
        # value: mass value [str]
        self.channel_dic = {}
        self.channel_dic_name = "channel_dic.txt"
        
        #########################
        # systematic uncertainties
        #########################
        #Test_v01, 
        self.nusiance_dir =\
            "/afs/cern.ch/user/x/xizhao/work/public/newOutput_v02/"+ self.ws_type +"/"
        self.np_list_name = "nuisance.txt"
        self.np_norm = Set([])
        self.np_shape = Set([])
        
        #########################
        # workspace input
        #########################
        self.ws_input_dir = "./"
        self.statistic_thr = -0.02

        #########################
        # ZZ shape
        #########################
        self.zz_config_name = "zz_theory.ini"
        self.qqZZ_pdfname = "pdf_VV_highMass"
        self.qqZZ_qcdname = "QCDscale_VV"
        self.ggZZ_pdfname = self.qqZZ_pdfname
        self.ggZZ_qcdname = "QCDscale_ggVV"

        #########################
        # yields related
        #########################
        self.yield_output = ""

    def get_list(self):
        """ select input minitrees for each sample using pre-define patterns
        and saved to a text file, which will be used in smoothing
        """
        for key in self.config.samples:
            value = self.sample_dict[key]
            #print key
            out_name = key + ".list"
            file_info = ""
            if value[1]:
                for im in range(self.config.mass_points):
                    mass = str(self.config.get_mass(im))
                    try:
                        file_path = self.config.samples_sig[mass][key]
                        common.check_file(file_path)
                        if not os.path.isfile(file_path):
                            print file_path," does not exist!"
                        file_info += file_path + ' '+mass+'\n'
                    except KeyError:
                        pass
            else:
                file_path = self.config.samples_bkg[key]
                common.check_file(file_path)
                file_info = file_path + '\n'
            with open(out_name, 'w') as f:
                f.write(file_info)
            print out_name," written"

    def get_smooth_config(self, smooth_config_name):
        """ generate a config file for smoothing """
        out = "[main]\n"
        out += 'categories = '+','.join(iter(self.config.categories)) + '\n'
        out += 'filedir = ./\n'
        out += 'outdir = '+self.ws_input_dir+'\n'
        out += 'outputname = '+self.out_prefix+'\n'
        for key in iter(self.config.samples):
            out += key+' = '+key+".list, "+str(self.rho) + '\n'
        
        signals = []
        bkgs = []
        for key in self.config.samples:
            value = self.sample_dict[key]
            if value[1]:
                signals.append(key)
            else:
                bkgs.append(key)

        out += "signals = " + ','.join(signals) + '\n'
        out += "backgrounds = " + ','.join(bkgs) + '\n'
        out += "observables = " + self.obs_name + '\n'
        out += "treename = tree_incl_all\n"
        out += "branch = " + self.config.branch + '\n'
        for category, cuts in self.config.categories.iteritems():
            out += "["+category+"]\n"
            out += "cut = " + cuts + '\n'

        with open(smooth_config_name, 'w') as f:
            f.write(out)

        print smooth_config_name," has been created"


    def get_ws_config(self, higgs_mass):
        out = "[main]\n"
        out += "data = "+self.config.data + '\n'
        out += "fileDir = "+self.ws_input_dir + '\n'
        out += 'NPlist = '+self.np_list_name+'\n'
        out += "observable = "+self.obs_name+", "+self.config.obs_binning+"\n"
        out += 'categories = '+','.join(sorted(iter(self.config.categories))) + '\n'
        sample_list = self.config.samples
        if higgs_mass > 400 and self.config == high_mass:
            sample_list = self.config.samples_gt_400
        out += 'mcsets = '+','.join(sample_list)+'\n'
        
        ### normalization info ###
        yield_table = "yields_13TeV_"+str(higgs_mass)+".txt"
        if not os.path.isfile(yield_table):
            dump = dump_yields.DumpYield(self.config)
            with open(yield_table, 'w') as f:
                f.write(dump.get_table(str(higgs_mass)))
        out += 'normalization = '+yield_table+'\n'

        ### sample information ####
        out += '[samples]\n'
        #for key, value in sorted(self.samples.iteritems()):
        for key in sample_list:
            value = self.sample_dict[key]
            out += key+' = '
            hist_input_name = self.out_prefix+'_'+key
            shape_input_name = key
            norm_input_name = "norm_"+key
            if value[1]:
                hist_input_name += '_'+str(higgs_mass)
                shape_input_name = key+str(higgs_mass)
                norm_input_name += str(higgs_mass)
            hist_input_name += '.root'
            if self.ws_type != "":
                shape_input_name += "_"+self.ws_type+'_Shape.root'
                norm_input_name += "_"+self.ws_type+'.txt'
            else:
                shape_input_name += '_Shape.root'
                norm_input_name += '.txt'
            self.check_file(self.ws_input_dir+hist_input_name)
            self.check_file(self.ws_input_dir+shape_input_name)
            self.check_file(self.ws_input_dir+norm_input_name)

            out += hist_input_name+', '+shape_input_name+', '+\
                norm_input_name+', '+value[2]
            if self.statistic_thr > 0:
                out += ', '+str(self.statistic_thr) 
            out += '\n'
        out_name = self.ws_input_dir+"config/wsconfig_"+str(higgs_mass)+".ini"
        if not os.path.isdir(os.path.dirname(out_name)):
            os.mkdir(os.path.dirname(out_name))
        with open(out_name, 'w') as f :
            f.write(out)
        print out_name," written"
        return out_name

    def rename(self, name):
        """ rename shape sys, changing mc_channel_number to a string"""
        if "qqZZ_Powheg" in name:
            return name.replace("qqZZ_Powheg", "qqZZ")
        suffix = name.split('.')[1]
        real_name = name.split('.')[0]
        newname = ""
        if len(self.channel_dic) < 1:
            if os.path.isfile(self.ws_input_dir+self.channel_dic_name):
                self.read_chan()
            else:
                self.get_list()
        for index, token in enumerate(real_name.split('_')):
            
            if index > 0:
                newname += '_'
            try:
                app = self.channel_dic[token]
            except KeyError:
                app = token
            newname += app
        newname += '.'+suffix
        return newname
        #print name, newname

    def read_chan(self):
        with open(self.ws_input_dir+self.channel_dic_name, 'r') as f:
            for line in f:
                key, value = line[:-1].split(' ')
                self.channel_dic[key] = value

    def get_zz_shape_config(self):
        out = "[main]\n"
        out += 'categories = '+','.join(iter(self.config.categories)) + '\n'
        out += self.obs_name+' = '+self.config.obs_binning+'\n'
        out += '[ggZZ]\n'
        out += 'pdfName = '+self.ggZZ_pdfname+'\n'
        out += 'qcdName = '+self.ggZZ_qcdname+'\n'
        out += '[qqZZ]\n'
        out += 'pdfName = '+self.qqZZ_pdfname+'\n'
        out += 'qcdName = '+self.qqZZ_qcdname+'\n'
        out_name = self.ws_input_dir+self.zz_config_name
        with open(out_name, 'w') as f:
            f.write(out)
        print out_name," written"

    def check_file(self, name):
        if not os.path.isfile(name):
            print "[ERROR]: ", name, " does not exist"

    def get_Analytic(self):
        out = "[main]\n"
        out += "fileDir = "+self.ws_input_dir + '\n'
        out += 'NPlist = '+self.np_list_name+'\n'
        out += "observable = "+self.obs_name+", "+self.config.obs_binning+"\n"
        out += 'categories = '+','.join(sorted(iter(self.config.categories))) + '\n'
        sample_list = self.config.samples_para
        out += 'mcsets = '+','.join(sample_list)+'\n'
        
        ### normalization info ###
        yield_table = "yields_13TeV.txt"
        dump = dump_yields.DumpYield(self.config)
        with open(yield_table, 'w') as f:
            f.write(dump.get_table(str(125)))
        out += 'normalization = '+yield_table+'\n'

        ### sample information ####
        out += '[samples]\n'
        #for key, value in sorted(self.samples.iteritems()):
        for key in sample_list:
            value = self.sample_dict[key]
            out += key+' = '
            if value[1]: # keys based
                config_file = key+"_config.ini"
                out +=  config_file+', '+value[2] + '\n'
                self.write_para_config(key, config_file)
            else:  #hist based
                hist_input_name = self.out_prefix+'_'+key
                hist_input_name += '.root'
                shape_input_name = key
                norm_input_name = "norm_"+key
                if self.ws_type != "":
                    shape_input_name += "_"+self.ws_type+'_Shape.root'
                    norm_input_name += "_"+self.ws_type+'.txt'
                else:
                    shape_input_name += '_Shape.root'
                    norm_input_name += '.txt'

                self.check_file(self.ws_input_dir+hist_input_name)
                self.check_file(self.ws_input_dir+shape_input_name)
                self.check_file(self.ws_input_dir+norm_input_name)
                out += hist_input_name+', '+shape_input_name+', '+\
                    norm_input_name+', '+value[2]
                if self.statistic_thr > 0:
                    out += ', '+str(self.statistic_thr) 
                out += '\n'

        out_name = self.ws_input_dir+"wsconfig.ini"
        with open(out_name, 'w') as f :
            f.write(out)
        print out_name," written"
        return out_name

    def write_para_config(self, key, config_file):
        if not key in self.config.sig_samples:
            return 
        out = "[Init]\n"
        good_masses = []
        mass_info = ""
        dump2 = dump_yields.DumpYield(self.config)
        for im in range(self.config.mass_points):
            mass = str(self.config.get_mass(im))
            try: 
                file_pa = self.config.samples_sig[mass][key]
            except:
                print "cannot find",mass,key
                continue
            good_masses.append(mass)
            mass_info += "["+mass+"]\n"
            mass_info += 'minitree = ' + os.path.basename(file_pa) + '\n'
            mass_info += 'mH = ' + mass+ '\n'
            mass_info += 'shape = mean_'+key+mass+"_"+self.config.name+".txt\n"
            mass_info += 'norm = norm_'+key+mass+"_" + self.config.name + '.txt\n'
            text_info, list_info = dump2.get_yield(file_pa)
            mass_info += 'yield = '+','.join(list_info) + '\n'

        out += 'mcsets = ' + ','.join(good_masses)
        out += '\n'
        out += 'minitree_path = ' +common.minitree_dir+"\n"
        out += mass_info
        with open(config_file, 'w') as f:
            f.write(out)

    def make_ws(self, wsconfig_name, mh):
        cmds = ["mainCombiner", wsconfig_name,
                self.wsdir+"combined_"+mh+".root"]
        p = subprocess.Popen(cmds, stdout=subprocess.PIPE).communicate()[0]
        print p

    def get_np_list(self):
        all_norm_files = glob.glob(self.ws_input_dir+"norm_*.txt")
        for norm_file in all_norm_files:
            if not os.path.isfile(norm_file): continue
            with open(norm_file, 'r') as f:
                for line in f:
                    if line[0] == '[':
                        continue
                    self.np_norm.add(line.split('=')[0].strip())
        print len(self.np_norm)," norm systematics found"

        all_shape_files = glob.glob(self.ws_input_dir+"*.root")
        for shape_file in all_shape_files:
            handle = ROOT.TFile(shape_file, "read")
            for key in handle.GetListOfKeys():
                tokens = key.GetName().split('-')
                if len(tokens) < 2 or "Nominal" in key.GetName():
                    continue
                np_name = key.GetName().split('-')[1]
                self.np_shape.add(np_name)
        print len(self.np_shape), " shape systematics found"

        self.np_shape.add(self.qqZZ_pdfname)
        self.np_shape.add(self.ggZZ_pdfname)
        self.np_shape.add(self.qqZZ_qcdname)
        self.np_shape.add(self.ggZZ_qcdname)
        combined_np = self.np_shape.union(self.np_norm)
        print "total: ", len(combined_np)
        print sorted(combined_np)
        out_list = self.ws_input_dir+self.np_list_name
        with open(out_list, 'w') as f:
            f.write('\n'.join(sorted(combined_np)))
        print "write to: ",out_list
            

def run_pvalue(config):
    out =""
    for im in range(config.mass_points):
        mh = config.get_mass(im)
        out += str(mh) + "\n"
        ws_name = wsdir+"combined_"+str(mh)+".root"
        cmds = [common.pvalue_cmd, ws_name, "mu_BSM"]
        p = subprocess.Popen(cmds, stdout=subprocess.PIPE).communicate()[0]
        out += p
    pattern = re.compile("expected")
    print re.findall("expected p0: [0-9]*.[0-9]*", out)
    with open("pvalue.txt", 'w') as f:
        f.write(out)

if __name__ == '__main__':
    usage = "usage: "+sys.argv[0]+"[options] Low/High/MonoH"
    parser = OptionParser(usage, version="0.1")
    parser.add_option("-l", "--list", action="store_true", dest="get_list",
                      default=False, help="get_list")
    parser.add_option("-s", "--smooth", action="store_true", dest="smooth",
                      default="", help="get_smooth_config")
    parser.add_option("-w","--workspace", action="store_true", dest="ws",
                      default="", help="get_ws_config")
    parser.add_option("-a","--analytic", action="store_true", dest="ana",
                      default="", help="get_Analytic")
    parser.add_option("-m","--mass", dest="mass",
                      default="", help="higgs mass")
    parser.add_option("-r","--run", action="store_true", dest="run",
                      default=False, help="make workspace/smooth")
    parser.add_option("-z","--zzshape", action="store_true", dest="zzshape",
                      default=False, help="get_zz_shape_config")
    parser.add_option("-n","--nplist", action="store_true", dest="nplist",
                      default=False, help="get_zz_shape_config")
    options,args = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        exit(1)

    ws_type = args[0]
    if "lo" in ws_type.lower():
        info = ("Low", low_mass) ## Low mass
    elif "hi" in ws_type.lower():
        info = ("High", high_mass) ## High mass
    elif "mo" in ws_type.lower():
        info = ("MonoH", monoH)
    else:
        print "I don't understand: ", ws_type
        exit(2)

    prepare = PrepareWS(info)
    if options.get_list:
        prepare.get_list()

    if options.zzshape:
        prepare.get_zz_shape_config()

    if options.smooth:
        prepare.get_smooth_config("smooth.ini")

    if options.ana:
        prepare.get_Analytic()

    if options.nplist:
        prepare.get_np_list()

    if options.ws:
        if options.mass != "":
            try:
                wsconfig_name = prepare.get_ws_config(int(options.mass))
            except:
                print options.mass," is not right!"
                exit(3)
            if options.run:
                prepare.make_ws(wsconfig_name, options.mass)
        else:
            for im in range(prepare.config.mass_points):
                mh = prepare.config.get_mass(im)
                if mh > 1000 and prepare.config == high_mass: continue
                wsconfig_name = prepare.get_ws_config(mh)
                if options.run:
                    prepare.make_ws(wsconfig_name, str(mh))

