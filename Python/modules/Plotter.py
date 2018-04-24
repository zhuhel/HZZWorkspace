import ROOT
import logging
from array import array
import sys, re, Utilities


def plot_NPs(data, str_id, output_dir='.', save_file=None, prune=0.01, publicity="Simulation Internal"):
    """
    Takes a dictionary with [NP][variation]["norm"], makes a pretty plot to print, and saves them to a root file optionally.
    """
    logging.debug("Plotting up/down variation of %s NPs", len(data))
    Utilities.check_and_mkdir(output_dir)
    plot_points = Utilities.data_to_plotpoints(data, prune)

    # Sort by largest overall deviation
    plot_points.sort(key=lambda x: x[2] - x[1], reverse=True)
    n_points = len(plot_points)
    graph = ROOT.TGraphAsymmErrors(n_points, array('d', range(n_points)), array('d', [0.0 for i in range(n_points)]),
                                   array('d', [0.45 for i in range(n_points)]), array('d', [0.45 for i in range(n_points)]),
                                   array('d', [abs(x[1]) for x in plot_points]), array('d', [abs(x[2]) for x in plot_points]))

    graph.SetTitle("")
    graph.SetName("graph_%s" % str_id)
    graph.GetYaxis().SetTitle("Systematic variation [%]")

    canvas = ROOT.TCanvas("canvas_%s" % str_id, str_id, 800, 400)
    canvas.cd()
    
    graph.Draw("A5")
    graph.SetFillColorAlpha(9, 0.8)

    y_min, y_max = graph.GetYaxis().GetXmin(), graph.GetYaxis().GetXmax()
    line = ROOT.TLine(0.0, y_min, 0.0, y_max)    
    line.SetLineColor(16)
    for i in range(n_points):
        line.DrawLine(float(i), y_min*1.025, float(i), plot_points[i][1])

    graph.GetXaxis().SetLabelSize(0)
    text = ROOT.TLatex()
    text.SetTextAlign(33)
    text.SetTextAngle(45)
    text.SetTextSize(0.02)
    if len(plot_points) > 50:
        text.SetTextSize(0.015)
    text.SetTextFont(42)
    for x, n in enumerate(plot_points):
        np = n[0]
        text.DrawLatex(x, y_min*1.05, np)
    

    legend = ROOT.TLatex()
    legend.SetTextSize(0.035)
    legend.SetTextFont(42)
    legend.SetTextAlign(11)
    legend.DrawLatexNDC(0.1, 0.91, "#bf{#it{ATLAS}} " + publicity)
    legend.SetTextAlign(31)
    legend.DrawLatexNDC(0.95, 0.91, str_id)

    canvas.SetGridy()
    canvas.SetBottomMargin(0.35)
    canvas.SetTopMargin(0.1)
    canvas.SetRightMargin(0.05)

    canvas.Print(output_dir + "/" + str_id + ".pdf")
    
    if save_file:
        save_file.cd()
        canvas.Write()


def plot_compare_NPs(filenames, titles, str_id, output_dir='.', save_file=None, prune=0.01, publicity="Simulation Internal"):
    """
    Takes a list of two norm_... filenames, their titles, and an ID to make a pretty plot comparing them, and saves them to a root file optionally.
    """
    logging.debug("Plotting NP comparisons for files %s", filenames)
    Utilities.check_and_mkdir(output_dir)
    if len(filenames) != 2 or len(titles) != 2:
        logging.error("Please ensure you have specified only 2 input files and titles")
        sys.exit(0)

    plot_points = {}
    for filename, title in zip(filenames, titles):
        plot_points[title] = {}
        with open(filename, 'r') as f:
            category = None
            for line in f:
                if line.startswith("["):
                    category = line.strip().replace("[", "").replace("]", "")
                    plot_points[title][category] = []
                else:
                    var_options = [float(line.split(" ")[2]) - 1.0, float(line.split(" ")[3]) - 1.0, 0.0]
                    up_variation = 100*max(var_options)
                    down_variation = 100*min(var_options)
                    name = line.split(" ")[0].replace("ATLAS_", "")
                    plot_points[title][category].append((name, down_variation, up_variation))

    categories = sorted(plot_points[titles[0]].keys())
    for cat in categories:
        plot_points_primary = sorted(plot_points[titles[0]][cat], key=lambda x: x[0])
        plot_points_secondary = sorted(plot_points[titles[1]][cat], key=lambda x: x[0])
        #  The below command sorts both lists the same way, in order of highest total variation for the first filename given
        plot_points_primary, plot_points_secondary = zip(*sorted(zip(plot_points_primary, plot_points_secondary), key=lambda x: x[0][2] - x[0][1], reverse=True))
        
        plot_points_primary = list(plot_points_primary)
        plot_points_secondary = list(plot_points_secondary)
        
        #  Pruning
        prune_index = 0
        for p, s in zip(plot_points_primary, plot_points_secondary):
            if p[2] - p[1] < prune and s[2] - s[1] < prune:
                del plot_points_primary[prune_index]
                del plot_points_secondary[prune_index]
                prune_index -= 1
            prune_index += 1

        n_points = len(plot_points_primary)

        graph_primary = ROOT.TGraphAsymmErrors(n_points, array('d', range(n_points)), array('d', [0.0 for i in range(n_points)]),
                                       array('d', [0.35 for i in range(n_points)]), array('d', [0.0 for i in range(n_points)]),
                                       array('d', [abs(x[1]) for x in plot_points_primary]), array('d', [abs(x[2]) for x in plot_points_primary]))

        graph_primary.SetTitle("")
        graph_primary.SetName("graph_%s_%s_%s" % (str_id, cat, titles[0]))
        graph_primary.GetYaxis().SetTitle("Systematic variation [%]")

        canvas = ROOT.TCanvas("canvas_%s_%s" % (str_id, cat), str_id, 800, 400)
        canvas.cd()
    
        graph_primary.Draw("A5")
        graph_primary.SetFillColorAlpha(9, 0.8)

        graph_primary.GetXaxis().SetLimits(-2.0, float(n_points) + 0.5)

        graph_secondary = ROOT.TGraphAsymmErrors(n_points, array('d', range(n_points)), array('d', [0.0 for i in range(n_points)]),
                                       array('d', [0.0 for i in range(n_points)]), array('d', [0.35 for i in range(n_points)]),
                                       array('d', [abs(x[1]) for x in plot_points_secondary]), array('d', [abs(x[2]) for x in plot_points_secondary]))
        graph_secondary.Draw("5")
        graph_secondary.SetFillColorAlpha(46, 0.8)

        y_min, y_max = graph_primary.GetYaxis().GetXmin(), graph_primary.GetYaxis().GetXmax()
        line = ROOT.TLine(0.0, y_min, 0.0, y_max)    
        line.SetLineColor(16)
        for i in range(n_points):
            line.DrawLine(float(i), y_min*1.025, float(i), plot_points_primary[i][1])

        graph_primary.GetXaxis().SetLabelSize(0)
        text = ROOT.TLatex()
        text.SetTextAlign(33)
        text.SetTextAngle(45)
        text.SetTextSize(0.02)
        if n_points > 50:
            text.SetTextSize(0.015)
        text.SetTextFont(42)
        for x, n in enumerate(plot_points_primary):
            np = n[0]
            text.DrawLatex(x, y_min*1.05, np)
    
        legend = ROOT.TLatex()
        legend.SetTextSize(0.035)
        legend.SetTextFont(42)
        legend.SetTextAlign(11)
        legend.DrawLatexNDC(0.1, 0.91, "#bf{#it{ATLAS}} " + publicity)
        legend.SetTextAlign(31)
        legend.DrawLatexNDC(0.95, 0.91, str_id + "_" + cat)

        canvas.SetGridy()
        canvas.SetBottomMargin(0.35)
        canvas.SetTopMargin(0.1)
        canvas.SetRightMargin(0.05)

        legend = ROOT.TLegend(0.8, 0.75, 0.92, 0.87)
        legend.AddEntry(graph_primary, titles[0], "f")
        legend.AddEntry(graph_secondary, titles[1], "f")
        legend.Draw()

        canvas.Print(output_dir + "/" + str_id + "_" + cat + "_" + "_".join(titles) + "_compare.pdf")
    
        if save_file:
            save_file.cd()
            canvas.Write()
