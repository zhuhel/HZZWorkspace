import ROOT
import logging
from array import array


def plot_NPs(data, str_id, save_file=None):
    """
    Takes a dictionary with [NP][variation]["norm"], makes a pretty plot to print, and saves them to file optionally.
    """
    logging.debug("Plotting up/down variation of %s NPs", len(data))
    plot_points = []
    for np, var in data.iteritems():
        if np == 'Nominal':
            continue
        # Store the percentage deviation
        up_variation = 100*(max([var['up']['norm'] - 1.0, var['down']['norm'] - 1.0, 0.0]))
        down_variation = 100*(min([var['up']['norm'] - 1.0, var['down']['norm'] - 1.0, 0.0]))
        plot_points.append((np.replace("ATLAS_", ""), down_variation, up_variation))
    
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
    graph.SetFillColorAlpha(9, 0.62)

    y_min, y_max = graph.GetYaxis().GetXmin(), graph.GetYaxis().GetXmax()
    line = ROOT.TLine(0.0, y_min, 0.0, y_max)    
    line.SetLineStyle(3)
    line.SetLineColor(14)
    for i in range(n_points):      
        line.DrawLine(float(i), y_min, float(i), 0.0)

    graph.GetXaxis().SetLabelSize(0)
    text = ROOT.TLatex()
    text.SetTextAlign(33)
    text.SetTextAngle(45)
    text.SetTextSize(0.02)
    text.SetTextFont(42)
    for x, n in enumerate(plot_points):
        np = n[0] + " "
        text.DrawLatex(x, y_min, np)
    
    graph.Draw("5")

    legend = ROOT.TLatex()
    legend.SetTextSize(0.035)
    legend.SetTextFont(42)
    legend.SetTextAlign(11)
    legend.DrawLatexNDC(0.1, 0.91, "#bf{#it{ATLAS}} Simulation Internal")
    legend.SetTextAlign(31)
    legend.DrawLatexNDC(0.9, 0.91, str_id)

    canvas.SetGridy()
    canvas.SetBottomMargin(0.35)
    canvas.SetTopMargin(0.1)

    canvas.Print(str_id + ".pdf")
    
    if save_file:
        save_file.cd()
        canvas.Write()
