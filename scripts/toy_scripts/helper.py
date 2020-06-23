# -*- coding: utf-8 -*-

import ROOT
ROOT.gROOT.SetBatch()

import math
import array


def get_mass(input_str):
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
            if line[0] == '#':
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
            if len(items) > 9:
                one_entry.append(float(items[9]))

            limit_dic[mass_] = one_entry
    return limit_dic

def tabulize_limit(file_name):
    out = ""
    out +="\\begin{table}[!htb]\n \t\\centering\n\t\\caption{}\n"
    out +="\t\\vspace{0.3cm}\n \t\\resizebox*{!}{\\dimexpr\\textheight-2\\baselineskip\\relax}{\n \t\\begin{tabular}{*{7}{c}}\n"
    out += "\tmH & -2~$\sigma$ & -1~$\sigma$ & median & 1~$\sigma$ & 2~$\sigma$ & observed \\\\ \\hline"
    out += '\n'
    indices = 1, 2
    with open(file_name) as f:
        for line in f:
            items = line[:-1].split()
            mass = get_mass(items[0])
            items[0] = mass
            new_list = [i for j,i in enumerate(items) if j not in indices]
            out += '\t'
            out += ' & '.join(new_list)
            out += ' \\\\ \n'
    out += '\\hline\n'
    out += "\t\\end{tabular}\n}\n\\end{table}\n"

    print(out)

def set_rooArgSet(roo_arg, value):
    itr = ROOT.TIter(roo_arg.createIterator())
    var = itr()
    while var:
        var.setVal(value)
        var = itr()

def get_total_yields(file_name):
    fin = ROOT.TFile.Open(file_name)
    ws = fin.Get("combWS")
    mc = ws.obj("ModelConfig")
    pdf = mc.GetPdf()
    ws.var("XS_ggF").setVal(0)
    ws.var("XS_VBF").setVal(0)
    ws.var("mH_llll").setVal(700)
    mc.GetParametersOfInterest().Print("v")
    observables = mc.GetObservables()
    obs = pdf.getObservables(observables)
    obs.Print('v')
    pdf.Print('v')

    nuisances = mc.GetNuisanceParameters()
    set_rooArgSet(nuisances, 0)
    total = pdf.expectedEvents(obs)

    this_pdf = ws.obj('model_ggF_2mu2e_13TeV_llll_deComposed')
    this_obs = this_pdf.getObservables(observables)
    total_4l_2mu2e = this_pdf.expectedEvents(this_obs)

    # get up error
    nuis_itr = ROOT.TIter(nuisances.createIterator())
    var = nuis_itr()
    up_error = 0
    while var:
        var.setVal(1.)
        yields = pdf.expectedEvents(obs)
        up_error += (total - yields)**2
        var.setVal(0.)
        var = nuis_itr()
    total_up = math.sqrt(up_error)

    # get down error
    nuis_itr.Reset()
    var = nuis_itr()
    down_error = 0
    while var:
        var.setVal(-1.)
        yields = pdf.expectedEvents(obs)
        down_error += (total - yields)**2
        var.setVal(0.)
        var = nuis_itr()
    total_down = math.sqrt(down_error)

    fin.Close()
    print("{:.2f} +{:.2f} -{:.2f}".format(total, total_up, total_down))
    return total, total_up, total_down

def compatability(obs, exp, up, down):
    sig_var = ROOT.RooRealVar("sig_var", "signal yields", 226., -500., 500.)
    exp_var = ROOT.RooRealVar("exp_var", "expected yields", exp)

    nuis_exp = ROOT.RooRealVar("nuis_exp", "nuisance parameter", 0., -5., 5.)
    low_values = ROOT.std.vector('double')()
    low_values.push_back(1 - down/exp)
    highValues = ROOT.std.vector('double')()
    highValues.push_back(1 + up/exp)
    exp_sys = ROOT.RooStats.HistFactory.FlexibleInterpVar(
        "exp_sys", "exp_sys",
        ROOT.RooArgList(nuis_exp), 1., low_values, highValues)
    exp_var_sys = ROOT.RooFormulaVar("exp_var_sys", "@0 + @1*@2", ROOT.RooArgList(sig_var, exp_var, exp_sys))
    exp_var_sys.Print()

    obs_var = ROOT.RooRealVar("obs_var", "observed events", obs)
    nll = ROOT.RooFormulaVar("nll", "-1 * @0 * TMath::Log(@1) + @1 + 0.5*@2*@2",
                             ROOT.RooArgList(obs_var, exp_var_sys, nuis_exp))

    # create NLL and minize
    minim = ROOT.RooMinimizer(nll)
    minim.setStrategy(1)
    minim.minimize("Minuit2", ROOT.Math.MinimizerOptions.DefaultMinimizerAlgo())

    minim_nll = nll.getVal()
    sig_list = []
    nll_list = []
    for sig in range(-100, 400, 1):
        sig_list.append(sig)
        sig_var.setVal(sig)
        sig_var.setConstant()
        minim.minimize("Minuit2", ROOT.Math.MinimizerOptions.DefaultMinimizerAlgo())
        nll_list.append(2*(nll.getVal() - minim_nll))

    sig_var.setVal(0.)
    minim.minimize("Minuit2", ROOT.Math.MinimizerOptions.DefaultMinimizerAlgo())
    print("sigma",math.sqrt(2*(nll.getVal() - minim_nll)))
    return math.sqrt(2*(nll.getVal() - minim_nll))

    #gr = ROOT.TGraph(len(sig_list), array.array('f', sig_list), array.array('f', nll_list))
    #c1 = ROOT.TCanvas("c1", "c1", 600, 600)
    #gr.Draw("ALP")
    #gr.GetYaxis().SetRangeUser(0., gr.GetYaxis().GetXmax())
    #c1.SaveAs("test.pdf")

if __name__ == "__main__":
    # compatability(1870.00, 1643.18, 164.05, 163.91)
    # compatability(681.0, 612.6, 36.7, 36.7) # llvv
    # compatability(9.0, 4.6, 0.98, 0.98) # llvv
    compatability(31, 19.5, 8., 8.0)
