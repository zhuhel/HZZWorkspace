#!/usr/bin/env python
"""
get acceptance for high mass analysis from workspace
"""

import ROOT
import sys
def get_acceptance(file_name, mH):
    fin = ROOT.TFile.Open(file_name)
    ws = fin.Get("combined")
    categories = {
        "$2\mu 2e$ ggF-like":"ggF_2mu2e_13TeV",
        "$4\mu $ ggF-like":"ggF_4mu_13TeV",
        "$4e$ ggF-like":"ggF_4e_13TeV",
        "VBF-like":"VBF_incl_13TeV"
    }
    samples = ["ggF", "VBF"]
    masses = [300, 800]
    #masses = [float(mH)]
    out = "category & "+" & ".join(samples)+" \\\\ \n"

    for category in categories.keys():
        out += category
        for mass in masses:
            ws.var("mH").setVal(mass)
            for sample in samples:
                poly_name = "ACpol_"+sample+"_"+categories[category]
                var = ws.obj(poly_name)
                if not var:
                    print poly_name,"does not exist"
                    continue

                val = var.getVal()
                out += " & {:.2f}".format(val*100)

        out += " \\\\ \n"
    print out



if __name__ == "__main__":
    if len(sys.argv) < 3:
        print sys.argv[0]," file_name mH"
        sys.exit(1)

    file_name = sys.argv[1]
    mH = sys.argv[2]
    get_acceptance(file_name, mH)
