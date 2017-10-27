#!/usr/bin/env python
import ROOT
import os
from optparse import OptionParser



def plot(f1_name,
         poi_name="mu", x_title="#mu",
         tree_name="physics", nll_name="NLL",
        ):
    canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
    ROOT.SetAtlasStyleCanvas(canvas, True)
    minMu = ROOT.Double(-999)
    mulow = ROOT.Double(-999)
    muHi =  ROOT.Double(-999)
    file_name, l1_name = f1_name.split(':')

    graph = ROOT.getGraphFromFile(
        file_name,
        tree_name, nll_name, poi_name,
        minMu, mulow, muHi)
    graph.SetLineWidth(2)

    mg = ROOT.TMultiGraph()
    mg.Add(graph, "L")
    mg.Draw("A")
    mg.GetXaxis().SetTitle(x_title)
    mg.GetYaxis().SetTitle("-2ln#Lambda")
    xmin = mg.GetXaxis().GetXmin()
    xmax = mg.GetXaxis().GetXmax()

    onesigma = ROOT.TLine()
    onesigma.SetLineStyle(2)
    onesigma.SetLineWidth(2)
    y_var = 1.0
    x_off_sigma = xmax+0.08
    y_off_sigma = 0.08
    onesigma.DrawLine(xmin, y_var, xmax, y_var)
    ROOT.myText(x_off_sigma, y_var-y_off_sigma, 1, "1#sigma", False)
    y_var = 4.0
    onesigma.DrawLine(xmin, y_var, xmax, y_var)
    ROOT.myText(x_off_sigma, y_var-y_off_sigma, 1, "2#sigma", False)
    y_var = 9.0
    onesigma.DrawLine(xmin, y_var, xmax, y_var)
    ROOT.myText(x_off_sigma, y_var-y_off_sigma, 1, "3#sigma", False)

    x_off_title = 0.185
    lumi_weight = 36.5
    ROOT.myText(x_off_title, 0.85, 1, "#bf{#it{ATLAS}} Internal")
    #ROOT.myText(x_off_title, 0.80, 1, "m_{#chi} = 1 GeV")
    ROOT.myText(x_off_title, 0.75, 1, "13 TeV, {:.1f} fb^{{-1}}".format(lumi_weight))
    legend = ROOT.myLegend(0.6, 0.7, 0.8, 0.9)
    legend.AddEntry(graph, l1_name, "l")
    legend.Draw()

    canvas.SaveAs(file_name.replace("root",'pdf'))

if __name__ == "__main__":
    # plot("test_mass_new.root:observed m_{H}", "mH", "m_{H} [GeV]")
    usage = "%prog filename:title poi_name x_title"
    version="%prog 1.0"
    parser = OptionParser(usage=usage, description="plot 1D NLL scan", version=version)
    (options,args) = parser.parse_args()

    if len(args) < 3:
        print parser.print_help()
        exit(1)
    f1_name = args[0]
    poi_name = args[1]
    x_title = args[2]

    SCRIPT_DIR = os.getenv("HZZWSCODEDIR")+"/scripts/"
    print SCRIPT_DIR
    if not hasattr(ROOT, "getGraphFromFile"):
        ROOT.gROOT.LoadMacro(SCRIPT_DIR+"/draw1DNLL.cxx")

    if not hasattr(ROOT, "loader"):
        ROOT.gROOT.LoadMacro(SCRIPT_DIR+"/loader.c")

    import AtlasStyle
    ROOT.gROOT.SetBatch()

    # call real program
    plot(f1_name, poi_name, x_title)
