#!/usr/bin/env python

from argparse import ArgumentParser
import Utilities
import logging.config
import sys
import os
from array import array

parser = ArgumentParser()
parser.add_argument(dest='config_file', help="The input config file", type=str)
parser.add_argument('-o', '--output_dir', help="The output results directory", type=str, default='.')
args = parser.parse_args()

script_loc = os.path.dirname(os.path.realpath(__file__))
logging.config.fileConfig(script_loc + '/logging.ini')

try:
    logging.info("Loading ROOT module...")
    import ROOT
    dummy = ROOT.RooRealVar()
    del dummy  # this is just to print the RooFit text first, to neaten up the output
    ROOT.gROOT.SetBatch(True)  # to prevent TCanvases showing up
except ImportError:
    logging.error("Could not import ROOT module. Make sure your root is configured to work with Python.")
    sys.exit()


def get_hist(tree_name, file_names, observables, hist_name, weight_name, cuts, smooth=False, bins=None):
    logging.info("Making chain with name %s", tree_name)
    chain = ROOT.TChain(tree_name)
    for f in file_names:
        logging.info("Adding file %s", f)
        chain.Add(f)
    if chain.GetEntries() == 0:
        logging.warn("Chain %s has no events", tree_name)
    else:
        logging.info("Chain %s has %i events", tree_name, chain.GetEntries())

    tree_obs = []
    for obs in observables:
        obs_params = obs.split(",")
        if len(obs_params) != 4:
            logging.warn("Wrong number of parameters in observable %s. Skipping...", obs)
            return None
        obs_name = obs_params[0].split(":")[0]
        variable = ROOT.RooRealVar(obs_name, obs_name, float(obs_params[2]), float(obs_params[3]))
        variable.setBins(int(obs_params[1]))
        tree_obs.append(variable)
        
    if len(tree_obs) < 1:
        logging.error("Could not find observables. Exiting...")
        return None
    elif len(tree_obs) > 2:
        raise NotImplementedError("No support for more than 2 observables")

    return_hist = None

    if not smooth:
        x = tree_obs[0]
        if len(tree_obs) == 1:
            return_hist = ROOT.TH1F(hist_name, hist_name, x.getBinning().numBins(), x.getMin(), x.getMax())
            logging.debug("Drawing %s with weight %s*(%s)", x.GetName(), weight_name, cuts)
            chain.Draw("{0}>>{1}".format(x.GetName(), hist_name), "{0}*({1})".format(weight_name, cuts))
            if bins:
                return_hist.SetName('unbinned_hist')
                binning = array('f', [float(b) for b in bins.split('/')])
                return_hist = return_hist.Rebin(len(binning) - 1, hist_name, binning)

        elif len(tree_obs) == 2:
            y = tree_obs[1]
            return_hist = ROOT.TH2F(hist_name, hist_name, x.getBinning().numBins(), x.getMin(), x.getMax(), 
                                    y.getBinning().numBins(), y.getMin(), y.getMax())
            logging.debug("Drawing %s and %s (2D) with weight %s*(%s)", x.GetName(), y.GetName(), weight_name, cuts)
            chain.Draw("{0}:{1}>>{2}".format(x.GetName(), y.GetName(), hist_name), "{0}*({1})".format(weight_name, cuts))
    else:
        raise NotImplementedError("Smoothing not implemented yet...")

    logging.info("Events after cut: %i; Integral: %s", return_hist.GetEntries(), return_hist.Integral())
    return return_hist


logging.info("Calculating systematics from config file %s", args.config_file)

top_config = Utilities.config_dict(args.config_file, tokenize=[('samples', ','), ('categories', ','), ('observables', ';')])
logging.debug(Utilities.config_display(top_config))
np_config = Utilities.config_dict(top_config['main']['NPlist'])
logging.debug(Utilities.config_display(np_config))

if 'outdir' in top_config['main']:
    logging.warn("Output directory %s found in config file, replacing the command line argument of %s",
                 top_config['main']['outdir'], args.output_dir)

if not Utilities.check_config(top_config, ['categories', 'samples', 'path', 'treename', 'NPlist', 'weightName']):
    logging.error("Config file does not contain required information. Exiting...")
    sys.exit()

if not Utilities.check_np(np_config):
    logging.error("NP list file does not contain even 'up' and 'down' variations. Exiting...")
    sys.exit()

data = {}

for sample in top_config['main']['samples']:
    logging.info("Starting calculation for %s sample\n%s", sample, '-'*(39 + len(sample)))    
    data[sample] = {}

    try:
        output_text_file = open("norm_%s.txt" % sample, 'w')
    except:
        logging.error("Failed to open file norm_%s.txt for writing results", sample)

    for category in top_config['main']['categories']:
        data[sample][category] = {}
        output_text_file.write("[%s]\n" % category)

        observable_fullname = '_'.join([o.split(',')[0].split(':')[0] for o in top_config[category]['observables']])
        logging.info("Category: %s", category)
        if 'smooth' not in top_config[category]:
            smooth = False
            logging.debug("Using binned histograms, since smoothing information not found in config file")
        else:
            smooth = top_config[category]['smooth']

        if 'bins' in top_config[category]:
            bins = top_config[category]['bins']
        else:
            bins = None

        # Start with the Nominal configuration
        hist_name = '-'.join([observable_fullname, sample, category])
        data[sample][category]["Nominal"] = {}
        logging.debug("Histogram name: %s", hist_name)

        logging.info("Starting with nominal configuration")
        try:
            file_names = [top_config['main']['path'] + "/Nominal/" + s for s in top_config['samples'][sample].split(",")]
        except KeyError:
            logging.error("Could not find %s sample filename(s) in config file. Exiting...", sample)
            sys.exit()

        data[sample][category]["Nominal"]["hist"] = get_hist(top_config['main']['treename'], file_names, top_config[category]['observables'], 
                                                             hist_name, top_config['main']['weightName'], top_config[category]['cuts'], smooth, bins)

        # Shape-like systematics
        for syst_name in np_config['shapeLike']:
            syst_title = np_config['shapeLike'][syst_name]
            file_names = [top_config['main']['path'] + "/{0}/{1}/".format(top_config['main']['sysDir'], syst_name) + s for s in top_config['samples'][sample].split(",")]
            if syst_title not in data[sample][category]:
                data[sample][category][syst_title] = {"up": {}, "down": {}}

            variation = "up"
            if "__1down" in syst_name:
                variation = "down"

            hist_name = '-'.join([observable_fullname, sample, category, syst_title, variation])
            
            if "JER" in syst_name:
                if variation == "down":
                    continue
                else:
                    data[sample][category][syst_title][variation]['hist'] = get_hist(top_config['main']['treename'], file_names, top_config[category]['observables'],
                                                                                     hist_name, top_config['main']['weightName'], top_config[category]['cuts'], smooth, bins)
                    data[sample][category][syst_title][variation]['norm'] = data[sample][category][syst_title][variation]['hist'].Integral()/data[sample][category]['Nominal']['hist'].Integral()
                    data[sample][category][syst_title][variation]['mean'] = data[sample][category][syst_title][variation]['hist'].GetMean()/data[sample][category]['Nominal']['hist'].GetMean()
                    data[sample][category][syst_title][variation]['sigma'] = data[sample][category][syst_title][variation]['hist'].GetRMS()/data[sample][category]['Nominal']['hist'].GetRMS()
                    # Symmetrise the JER
                    data[sample][category][syst_title]['down']['norm'] = 2 - data[sample][category][syst_title]['up']['norm']
                    data[sample][category][syst_title]['down']['mean'] = 2 - data[sample][category][syst_title]['up']['mean']
                    data[sample][category][syst_title]['down']['sigma'] = 2 - data[sample][category][syst_title]['up']['sigma']

            else:        
                data[sample][category][syst_title][variation]['hist'] = get_hist(top_config['main']['treename'], file_names, top_config[category]['observables'],
                                                                                 hist_name, top_config['main']['weightName'], top_config[category]['cuts'], smooth, bins)
                data[sample][category][syst_title][variation]['norm'] = data[sample][category][syst_title][variation]['hist'].Integral()/data[sample][category]['Nominal']['hist'].Integral()
                data[sample][category][syst_title][variation]['mean'] = data[sample][category][syst_title][variation]['hist'].GetMean()/data[sample][category]['Nominal']['hist'].GetMean()
                data[sample][category][syst_title][variation]['sigma'] = data[sample][category][syst_title][variation]['hist'].GetRMS()/data[sample][category]['Nominal']['hist'].GetRMS()
        
        # Norm-like systematics
        for syst_name in np_config['normLike']:
            syst_title = np_config['normLike'][syst_name]
            file_names = [top_config['main']['path'] + "/{0}/NormSystematic/".format(top_config['main']['sysDir']) + s for s in top_config['samples'][sample].split(",")]
            if syst_title not in data[sample][category]:
                data[sample][category][syst_title] = {"up": {}, "down": {}}

            variation = "up"
            if "__1down" in syst_name:
                variation = "down"

            hist_name = '-'.join([observable_fullname, sample, category, syst_title, variation])
            data[sample][category][syst_title][variation]['hist'] = get_hist(top_config['main']['treename'], file_names, top_config[category]['observables'],
                                                                             hist_name, "*".join([top_config['main']['weightName'], syst_name]), top_config[category]['cuts'], "", bins)
            data[sample][category][syst_title][variation]['norm'] = data[sample][category][syst_title][variation]['hist'].Integral()/data[sample][category]['Nominal']['hist'].Integral()
            data[sample][category][syst_title][variation]['mean'] = data[sample][category][syst_title][variation]['hist'].GetMean()/data[sample][category]['Nominal']['hist'].GetMean()
            data[sample][category][syst_title][variation]['sigma'] = data[sample][category][syst_title][variation]['hist'].GetRMS()/data[sample][category]['Nominal']['hist'].GetRMS()
            

        for syst_title in sorted(data[sample][category]):
            if syst_title == "Nominal":
                continue
            output_text_file.write("%s = %s %s\n" % (syst_title, data[sample][category][syst_title]['down']['norm'], data[sample][category][syst_title]['up']['norm']))

    output_text_file.close()
    
