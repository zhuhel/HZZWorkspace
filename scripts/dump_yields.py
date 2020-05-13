#!/usr/bin/env python
import ROOT
ROOT.gROOT.SetBatch()
if not hasattr(ROOT, "loader"):
    ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/loader.c")
import os

import low_mass
import high_mass
import monoH
import zp2hdm
import common

class DumpYield:
    def __init__(self, config):
        self.config = config
        #print self.config.categories
        self.scale = 1.0
        self.branching_ratio = 1.0
        self.lumi_weight = 1.0

    def get_yield(self, file_name, w_name = "weight"):
        out = ""
        yy = []
        if "342556" in file_name or "qq2ZZ" in file_name:
            if w_name != "": w_name = w_name+"*w_EW"
        if not os.path.isfile(file_name):
            print file_name," cannot be found"
            return (None, None)
        chain = ROOT.loader(file_name, common.tree_name)
        h1 = ROOT.TH1F("h1", "h1", 3, 0, 3000);
        ntotal = chain.GetEntries()
        branch_name = self.config.branch.split(',')[0]
        try:
            min_value = float(self.config.branch.split(',')[2])
            max_value = float(self.config.branch.split(',')[3])
        except ValueError:  
            print "cannot find range of: ", self.config.branch
            min_value = 110
            max_value = 140
            #print file_name, branch_name, max_value

        for category, cuts in sorted(self.config.categories.iteritems()):
            if w_name == "":  newcut = ""
            else: newcut = w_name+"*3.2095/3.316*"
            if self.scale < 0:
                newcut = w_name+"/w_br/w_xs*"
            newcut += cuts[:-1]+" && "+branch_name+">"+str(min_value)+\
                    " && "+branch_name+"<"+str(max_value)+")"
            #print newcut
            chain.Draw(branch_name+">>h1", newcut)
            #print category,h1.GetEntries()
            #result = "{:.7f}".format(h1.Integral()*self.branching_ratio*self.lumi_weight)
            result = "{:.1E}".format(h1.Integral()*self.branching_ratio*self.lumi_weight)
            out += " & "+ result
            yy.append(result)
        out += ' \\\\ \n'
        #print out
        return out,yy
    
    def get_met_yield(self, file_name, met_name):
        yy = []
        if not os.path.isfile(file_name):
            print file_name," cannot be found"
            return (None, None)
        chain = ROOT.TChain(common.tree_name, common.tree_name)
        chain.Add(file_name)
        h1 = ROOT.TH1F("h1", "h1", 3, 0, 3000);
        ntotal = chain.GetEntries()
        branch_name = self.config.branch.split(',')[0]
        for category, cuts in sorted(self.config.categories.iteritems()):
            newcut = "weight*"+cuts.replace("met_et", met_name)
            #print newcut
            chain.Draw(branch_name+">>h1", newcut)
            result = "{:.3E}".format(h1.Integral()*self.branching_ratio*self.lumi_weight)
            yy.append(result)
        return yy

    def print_yield(self,samples):
        out = ""
        for sp_key, sp_value in samples.iteritems():
            if sp_value is None or not os.path.isfile(sp_value): 
                print sp_value,"not found!"
                continue
            out += sp_key
            print sp_value 
            out += self.get_yield(sp_value)[0]
        return out

    def get_table(self, mass):
        out = ""
        ncount = 0
        for key, value in sorted(self.config.categories.iteritems()):
            if ncount == 0: out += key
            else: out += ' & '+ key
            ncount += 1
        out += '\n' 
        try:
            self.scale = self.config.samples_sig_scale
            if self.config == monoH: self.branching_ratio = 1.25e-4
            out += self.print_yield(self.config.samples_sig[mass])
        except:
            pass
        self.scale = self.config.samples_bkg_scale
        self.branching_ratio = 1.0
        out += self.print_yield(self.config.samples_bkg)
        return out

    def get_sys_text(self, in_name, out_name, add_met = False):
        #set the weight name by yourself
        # not use weight, to get rid of statistic uncertainties
        weight_name=""
        base_name = os.path.basename(in_name)
        file_name = common.minitree_dir+base_name
        if not os.path.isfile(file_name):
            print file_name,"not found"
            return
        text_info, nom = self.get_yield(file_name, weight_name)
        out_sys = ""
        sys_list = []

        #### Shape systematics ###
        for shape_sys,shape_name in common.shape_sys_list:
            file_down,file_up = common.get_shape_downup(file_name, shape_sys)
            text_info, down = self.get_yield(file_down, weight_name)
            text_info, up = self.get_yield(file_up, weight_name)
            out_sys += shape_sys
            sys_list_cats = [shape_name]
            #print file_down, file_up, file_name
            for ich in range(len(nom)):
                n = float(nom[ich])
                u = float(up[ich])
                if abs(n) < 1E-9: 
                    out_sys += " & ({:.2f} {:.2f})".format(-999, 999)
                    sys_list_cats.append((-999, 999))
                    continue
                up_var = (u-n)/n
                if not down:
                    out_sys += " & ({:.2f} {:.2f})".format(-1*100*up_var, 100*up_var)
                    sys_list_cats.append((1-up_var, 1+up_var))
                    continue
                d = float(down[ich])
                #print n,u,d
                down_var = (d-n)/n
                if "MET_SoftTrk_Reso" in shape_sys:
                    max_var = max([u, d])
                    down_var = (n - max_var)/n
                    up_var = (max_var - n )/n
                out_sys += " & ({:.2f} {:.2f})".format(100*down_var, 100*up_var)
                sys_list_cats.append((1+down_var, 1+up_var))
            sys_list.append(sys_list_cats)
            out_sys += " \\\\ \n"

        #### Normalization systematics ###
        file_norm = common.minitree_sys+"NormSystematic/"+base_name
        if "342556" in base_name:
            file_norm = "/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2015/MiniTrees/Prod_v03/mc_15b/Systematics/NormSystematic/mc15_13TeV.342556.PowhegPy8EG_CT10nloME_AZNLOCTEQ6L1_ZZllll_mll4_m4l_100_150.root"
            #file_norm = "/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2015/MiniTrees/Prod_v03/mc/Systematics/NormSystematic/mc15_13TeV.342556.PowhegPy8EG_CT10nloME_AZNLOCTEQ6L1_ZZllll_mll4_m4l_100_150.root"
        print file_norm
        for norm_sys, norm_name in common.norm_sys_list:
            weight_nom = "weight"
            weight_down = "weight_"+norm_sys+"__1down"
            weight_up = "weight_"+norm_sys+"__1up"
            text_info, nom = self.get_yield(file_norm, weight_nom)
            text_info, down = self.get_yield(file_norm, weight_down)
            text_info, up = self.get_yield(file_norm, weight_up)
            out_sys += norm_sys 
            sys_list_cats = [norm_name]
            for ich in range(len(nom)):
                n = float(nom[ich])
                u = float(up[ich])
                if abs(n) < 1E-9: 
                    out_sys += " & ({:.2f} {:.2f})".format(-999, 999)
                    sys_list_cats.append((-999, 999))
                    continue
                up_var = (u-n)/n
                if not down:
                    out_sys += " & ({:.2f} {:.2f})".format(100*up_var, 100*up_var)
                    sys_list_cats.append((1-up_var, 1+up_var))
                    continue
                d = float(down[ich])
                #print norm_sys,n,u,d
                down_var = (d-n)/n
                out_sys += " & ({:.2f} {:.2f})".format(100*down_var, 100*up_var)
                sys_list_cats.append((1+down_var, 1+up_var))
            sys_list.append(sys_list_cats)
            out_sys += " \\\\ \n"

        if add_met:
            file_norm = common.minitree_sys+"METSystematic/"+base_name
            print file_norm
            for met_name, met_sys in common.met_sys_list:
                weight_name = "weight"
                nom = self.get_met_yield(file_norm, "met_et")
                out_sys += met_name
                sys_list_cats = [met_sys]
                if "SoftTrk_Reso" in met_name:
                    down = self.get_met_yield(file_norm,
                                              "met_et_MET_SoftTrk_ResoPara")
                    up = self.get_met_yield(file_norm,
                                              "met_et_MET_SoftTrk_ResoPerp")
                else:
                    down = self.get_met_yield(file_norm,met_name+"Down")
                    up = self.get_met_yield(file_norm,met_name+"Up")
                for ich in range(len(nom)):
                    n = float(nom[ich])
                    u = float(up[ich])
                    d = float(down[ich])
                    if abs(n) < 1E-9: 
                        out_sys += " & ({:.2f} {:.2f})".format(-999, 999)
                        sys_list_cats.append((-999, 999))
                        continue
                    up_var = (u-n)/n
                    down_var = (d-n)/n
                    if "SoftTrk_Reso" in met_name:
                        max_var = max([u,d])
                        down_var = (n - max_var)/n
                        up_var = (max_var -n )/n
                    out_sys += " & ({:.2f} {:.2f})".format(100*down_var, 100*up_var)
                    sys_list_cats.append((1+down_var, 1+up_var))
                sys_list.append(sys_list_cats)
                out_sys += " \\\\ \n"

        print out_sys
        out_text = ""
        index  = 1
        for key in sorted(self.config.categories.keys()):
            #print key
            out_text += "["+key+"]\n"
            for sys in sys_list:
                out_text += "{} = {:.5f} {:.5f} \n".format(sys[0], sys[index][0],sys[index][1])
            index += 1
        with open(out_name, 'w') as f:
            f.write(out_text)
        #print out_text

def get_monoH_signal():
    dump = DumpYield(monoH)
    #dump = DumpYield(zp2hdm)
    dump.branching_ratio = 1.25E-4
    dump.lumi_weight = 1
    out = ""
    for idm in range(dump.config.dm_points):
        m_dm = dump.config.get_dm(idm)
        dump.config.m_dm = m_dm
        dump.config.samples_sig = dump.config.get_signal_dict()
        for im in range(dump.config.mass_points):
            mass = str(dump.config.get_mass(im))
            try:
                #file_name = dump.config.samples_sig[mass]["shxx"]
                file_name = dump.config.samples_sig[mass]["zphxx"]
                #file_name = dump.config.samples_sig[mass]["zp2hdm"]
                #print file_name
                if not os.path.isfile(file_name):
                    print "cannot find", file_name
                    continue
            except:
                print m_dm,mass," missing"
                continue

            text, number = dump.get_yield(file_name)
            out += str(m_dm)+ " & " + str(mass)
            out += text
    print out

def call_sys():
    dump = DumpYield(monoH)
    bkgs = ["ttH"]
    #for sample in dump.config.bkg_samples:
    for sample in bkgs:
        print sample
        out_name = "norm_"+sample+"_"+dump.config.name+".txt"
        pass
        dump.get_sys_text(dump.config.samples_bkg[sample], out_name, True)
    
    sample = "zphxx"
    for imass in range(dump.config.mass_points):
        mass = str(dump.config.get_mass(imass))
        out_name = "norm_"+sample+mass+"_"+dump.config.name+".txt"
        pass
        #dump.get_sys_text(dump.config.samples_sig[mass][sample], out_name, True)

    #out_name = "jj_"+sample+mass+"_"+dump.config.name+".txt"
    #dump.get_sys_text(dump.config.samples_sig["200"][sample], out_name)
    #out_name = "norm_ZHlvlv_"+dump.config.name+".txt"
    #dump.get_sys_text(dump.config.samples_bkg["ZHlvlv"], out_name)

def get_ana_yield(config, mass):
    dump = DumpYield(config)
    dump.lumi_weight = 1.0
    return dump.get_table(str(mass))

def get_low_mass():
    dump = DumpYield(low_mass)
    f1="/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2015/MiniTrees/Prod_v03/mc/Nominal/mc15_vvv.root"
    #text, info = dump.get_yield("signal.root")
    text, info = dump.get_yield(f1)
    print text

def get_sys():
    dump = DumpYield(low_mass)
    #f1 = "/afs/cern.ch/user/s/schaffer/h4l_2015/MiniTrees/Test_v04/RD_extra/mc/Systematics/NormSystematic/mc15_13TeV.342556.PowhegPy8EG_CT10nloME_AZNLOCTEQ6L1_ZZllll_mll4_m4l_100_150.root"
    #f1 ="/afs/cern.ch/user/s/schaffer/h4l_2015/MiniTrees/Test_v04/RD_extra/mc/Systematics/NormSystematic/mc15_13TeV.341505.PowhegPythia8EvtGen_CT10_AZNLOCTEQ6L1_ggH125_ZZ4lep_noTau.root"
    f1 = "/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2015/MiniTrees/Prod_v03/mc_15b/Systematics/NormSystematic/mc15_13TeV.342556.PowhegPy8EG_CT10nloME_AZNLOCTEQ6L1_ZZllll_mll4_m4l_100_150.root"
    tt, norm = dump.get_yield(f1, "weight*w_EW")
    print tt
    tt, up = dump.get_yield(f1, "weight_HOEW_syst__1up*w_EW")
    print tt
    tt, down = dump.get_yield(f1, "weight_HOEW_syst__1down*w_EW")
    print tt
    sys_list_cats = []
    for ich in range(len(norm)):
        n = float(norm[ich])
        u = float(up[ich])
        d = float(down[ich])
        down_var = (d-n)/n
        up_var = (u-n)/n
        sys_list_cats.append((1+down_var, 1+up_var))

    out_text = ""
    index = 0
    for key in sorted(low_mass.categories.keys()):
        out_text += "["+key+"]\n"
        out_text += "{} = {:.5f} {:.5f} \n".format("HOEW",
                                                       sys_list_cats[index][0],sys_list_cats[index][1])
        index += 1
    print out_text 

def get_ttbar_zjets_for_monoH():
    tt_llmm = 0.23
    zjets_llmm = 0.65
    
    #ggF_2e2mu_13TeV & ggF_2mu2e_13TeV & ggF_4e_13TeV & ggF_4mu_13TeV
    red_bkg = [0.191 , 0.2388 , 0.2388 , 0.191]
    met_eff_tt = 0.2653
    met_eff_zee = 0.0058
    met_eff_zmm = 0.0063
    ratio_tt = tt_llmm/(tt_llmm+zjets_llmm)
    total_tt = (red_bkg[0]+red_bkg[3])*ratio_tt
    out = ""
    out += 'ttbar: {:.2E} {:.2E}\n'.format(total_tt*met_eff_tt, total_tt*(1-met_eff_tt))
    total_zmm = (red_bkg[0]+red_bkg[3])*(1-ratio_tt)
    total_zee = (red_bkg[1]+red_bkg[2])
    total_z_hi = total_zmm*met_eff_zmm + total_zee*met_eff_zee
    total_z_low = total_zmm*(1-met_eff_zmm) + total_zee*(1-met_eff_zee)
    out += 'zjets: {:.2E} {:.2E}\n'.format(total_z_hi, total_z_low)
    print out

if __name__ == '__main__':
    get_monoH_signal()
    #call_shape_sys()
    #print get_ana_yield(monoH, "200")
    #get_sys()

    #get_low_mass()
    #call_sys()

    #get_ttbar_zjets_for_monoH()
