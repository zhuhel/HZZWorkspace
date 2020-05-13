#!/usr/bin/env python
import glob
import sys
import ROOT
import subprocess
from optparse import OptionParser

el_shape_list = ["ATLAS_EG_RESOLUTION_ALL", "ATLAS_EG_SCALE_ALLCORR"]
mu_shape_list = ["ATLAS_MUONS_ID", "ATLAS_MUONS_MS", "ATLAS_MUONS_SCALE"]
el_norm_list = ["ATLAS_EL_EFF_ID",
                "ATLAS_EL_EFF_RECO",
                "ATLAS_EL_EFF_Iso"]
el_sys_list = el_shape_list + el_norm_list
mu_norm_list = ["ATLAS_MU_EFF_STAT",
                "ATLAS_MU_EFF_SYS",
                "ATLAS_MU_EFF_STAT_LOWPT",
                "ATLAS_MU_EFF_SYS_LOWPT",
                "ATLAS_MU_ISO_STAT",
                "ATLAS_MU_ISO_SYS",
                "ATLAS_MU_TTVA_STAT",
                "ATLAS_MU_TTVA_SYS"
               ]
mu_sys_list = mu_shape_list + mu_norm_list
other_shapes = [
    "ATLAS_Shape_Zjet"
]
theory_list = [
    "QCDscale_VV",
    "QCDscale_ggVV",
    "pdf_gg",
    "pdf_qq",
    #"ATLAS_HOQCD_SCALE",
    #"ATLAS_HOEW",
    #"ATLAS_HOEW_QCD",
    "ATLAS_norm_Zjet_llee",
    "ATLAS_norm_Zjet_llmumu",
    "ATLAS_Shape_Zjet",
    #"ATLAS_PRW_DATASF",
]

def clean_sys(file_name, el_list, mu_list, is_test):
    out = ""
    all_list = el_list + mu_list + theory_list
    with open(file_name, 'r') as f:
        cat_name = ""
        for line in f:
            if line[0] == '[':
                out += line
                cat_name = line[1:-2]
                #print "In ", cat_name
                continue
            else:
                np_name = line[:-1].split('=')[0].strip()
                if np_name in el_list and "4mu" in cat_name:
                    continue
                if np_name in mu_list and "4e" in cat_name:
                    continue
                if np_name in all_list:
                    out += line
    if is_test: out_name = "test_norm.txt"
    else: out_name = file_name
    with open(out_name, 'w') as f:
        f.write(out)

def clean(pattern, func):
    file_lists = glob.glob(pattern)
    for infile in file_lists:
        print infile
        func(infile)

def clean_mean(file_name):
    clean_sys(file_name, el_shape_list, mu_shape_list)

def clean_norm(file_name, is_test):
    clean_sys(file_name, el_sys_list, mu_sys_list, is_test)

def smooth_ratio(h1):
    hold = h1.Clone("hold")
    total_bins = hold.GetNbinsX()
    i = 0
    while i < total_bins:
        ibin = i+1
        bin_val = hold.GetBinContent(ibin)
        if bin_val == 0 or abs(bin_val) < 1e-3: ##too small
            h1.SetBinContent(ibin, 1)
            i += 1
            continue
        elif bin_val > 2: ##too large
            if ibin-2 > 0 and ibin+2 <= total_bins:
                avg_val = bin_val + \
                        hold.GetBinContent(ibin-2) +\
                        hold.GetBinContent(ibin-1) +\
                        hold.GetBinContent(ibin+2) +\
                        hold.GetBinContent(ibin+1)
                for jj in range(5):
                    h1.SetBinContent(ibin-2+jj, avg_val/5.0)
                i += 3
                continue
            elif ibin - 1 > 0 and ibin +1 <= total_bins:
                avg_val = bin_val + \
                        hold.GetBinContent(ibin-1) +\
                        hold.GetBinContent(ibin+1)
                h1.SetBinContent(ibin, avg_val/3.0)
                for jj in range(3):
                    h1.SetBinContent(ibin-1+jj, avg_val/3.0)
                i += 2
                continue
            elif ibin != 1:
                h1.SetBinContent(ibin, h1.SetBinContent(ibin-1))
            else:
                pass
        else: ## in the middle
            pass
        i += 1
    #h1.Rebin(2)
    #h1.Scale(0.5)

def clean_shape_sys(file_name, el_list, mu_list, is_test):
    fin =  ROOT.TFile.Open(file_name, "read")
    hist_lists = []
    all_list = el_list + mu_list + other_shapes 
    for key in fin.GetListOfKeys():
        name_list = key.GetName().split('-')
        if(len(name_list) < 4): continue
        var,np_name,cat_name, vary_name = name_list
        if np_name in el_list and "4mu" in cat_name:
            continue
        if np_name in mu_list and "4e" in cat_name:
            continue
        if np_name in all_list:
            h1 = key.ReadObj()
            #smooth_ratio(h1)
            #if "ATLAS_HOEW" in np_name:
            #    h1.SetName(key.GetName().replace("ATLAS_HOEW","qqZZ_HOEW"))
            hist_lists.append(h1)
    fout = ROOT.TFile.Open("test.root", "recreate") 
    for hist in hist_lists:
        hist.Write()
    fout.Close()
    fin.Close()
    if not is_test: subprocess.call(['mv',"test.root",file_name])

def clean_shape(file_name, is_test):
    clean_shape_sys(file_name, el_shape_list, mu_shape_list, is_test)

def change_EL_EFF(file_name):
    out = ""
    with open(file_name, 'r') as f:
        cat_name = ""
        for line in f:
            if line[0] == '[':
                out += line
                cat_name = line[1:-2]
                continue
            else:
                if "ATLAS_EL_EFF_ID" in line:
                    if "4e" in cat_name:
                        out += "ATLAS_EL_EFF_ID = 0.90 1.1\n"
                    else:
                        out += "ATLAS_EL_EFF_ID = 0.95 1.05\n"
                elif "ATLAS_EL_EFF_RECO" in line:
                    if "4e" in cat_name:
                        out += "ATLAS_EL_EFF_RECO = 0.90 1.1\n"
                    else:
                        out += "ATLAS_EL_EFF_RECO = 0.95 1.05\n"
                else:
                    out += line
    with open(file_name, 'w') as f:
        f.write(out)

def change_EXP(file_name):
    out = ""
    with open(file_name, 'r') as f:
        cat_name = ""
        for line in f:
            if line[0] == '[':
                out += line
                cat_name = line[1:-2]
                continue
            else:
                if "ATLAS_QCDscale_ggVV" in line: 
                    out += line.replace("ATLAS_QCDscale_ggVV",
                                        "QCDscale_ggVV")
                elif "ATLAS_pdf_gg" in line:
                    out += line.replace("ATLAS_pdf_gg",
                                        "pdf_gg")
                else:
                    out += line
    with open(file_name, 'w') as f:
        f.write(out)

def add_sys_per_file(file_name, add_lines):
    out = ""
    with open(file_name, 'r') as f:
        for line in f:
            out += line
            if line[0] == '[':
                out += add_lines

    with open(file_name, 'w') as f:
        f.write(out)

def update_lumi(file_name):
    out = ""
    with open(file_name, 'r') as f:
        for line in f:
            if "ATLAS_lumi" in line:
                out += "ATLAS_lumi = 0.966 1.034\n"
            else:
                out += line
    with open(file_name, 'w') as f:
        f.write(out)

def add_lumi(file_name):
    out = ""
    with open(file_name, 'r') as f:
        for line in f:
            if line[0] == "[":
                out += line
                out += "ATLAS_lumi = 0.966 1.034\n"
            else:
                out += line
    with open(file_name, 'w') as f:
        f.write(out)


def add_sys_monoH():
    lumi_sys = "ATLAS_lumi_2015 = 0.91 1.09\n"
    line_dict = {
        "zphxx": lumi_sys+"QCDscale_qqVV = 0.972 1.021\npdf_qq =  0.833 1.099\n",
        "shxx": lumi_sys+"QCDscale_qqVV = 0.972 1.021\npdf_qq =  0.833 1.099\n",
        "ggH": lumi_sys+"QCDscale_ggVV = 0.921 1.074\npdf_gg = 0.940 1.071\n",
        "VBFH": lumi_sys+"QCDscale_qqVV = 0.993 1.007\npdf_qq = 0.968 1.032\n",
        "WH": lumi_sys+"QCDscale_qqVV = 0.993 1.015\npdf_qq = 0.978 1.022\n",
        "ZH": lumi_sys+"QCDscale_qqVV = 0.924 1.092\npdf_qq = 0.973 1.027\n",
        "qqZZ": lumi_sys+"QCDscale_qqVV = 0.97 1.03\npdf_qq = 0.96 1.04\n",
    }
    for file_in in glob.glob("norm_zphxx*.txt"):
        pass
        #add_sys_per_file(file_in, line_dict["zphxx"])

    for file_in in glob.glob("norm_shxx*.txt"):
        add_sys_per_file(file_in, line_dict["shxx"])

    for key in line_dict.keys():
        if key == "zphxx": continue
        #add_sys_per_file("norm_"+key+"_MonoH.txt", line_dict[key])

    #add_sys_per_file("norm_ZHlvlv_MonoH.txt", line_dict["ZH"])

if __name__ == "__main__":
    usage = "usage: "+sys.argv[0]+"[options] pattern"
    parser = OptionParser(usage, version="0.1")
    parser.add_option("--clean_norm", action="store_true", dest="clean_norm",
                      default=False, help="clean normalization files")
    parser.add_option("--clean_shape", action="store_true", dest="clean_shape",
                      default=False, help="clean shape files")
    parser.add_option("--up_lumi", action="store_true", dest="up_lumi",
                      default=False, help="update luminosity")
    parser.add_option("--add_lumi", action="store_true", dest="add_lumi",
                      default=False, help="add lumi sys")
    parser.add_option("--test", action="store_true", dest="test", default=False,
                      help="test only")
    options,args = parser.parse_args()
    if len(args) < 1:
        parser.print_help()
        exit(1)

    pattern = args[0]
    print "pattern: ", pattern
    for f in glob.glob(pattern):
        print f
        if options.clean_norm:
            clean_norm(f, options.test)
        if options.clean_shape:
            clean_shape(f, options.test)
        if options.up_lumi:
            update_lumi(f)
        if options.add_lumi:
            add_lumi(f)
