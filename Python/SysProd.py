#!/usr/bin/env python

from argparse import ArgumentParser
import Utilities
import logging.config
import sys
import os

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
except ImportError:
    logging.error("Could not import ROOT module. Make sure your root is configured to work with Python.")
    sys.exit()


def get_hist(tree_name, file_name, observables):
    logging.info("Making chain with name %s", tree_name)
    chain = ROOT.TChain(tree_name)
    logging.info("Adding file %s", file_name)  # could make this into a loop over files
    chain.Add(file_name)
    if chain.GetEntries() == 0:
        logging.warn("Chain %s has no events", tree_name)
    else:
        logging.info("Chain %s has %i events", tree_name, chain.GetEntries())

    tree_obs = ROOT.RooArgSet()
    for obs in observables:
        obs_params = obs.split(",")
        if len(obs_params) == 4:
            obs_name = obs_params[0].split(":")[0]
            variable = ROOT.RooRealVar(obs_name, obs_name, float(obs_params[2]), float(obs_params[3]))
            variable.setBins(int(obs_params[1]))
            tree_obs.add(variable)
        else:
            logging.warn("Wrong number of parameters in observable %s. Skipping...", obs)
    
    return 5


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

for sample in top_config['main']['samples']:
    logging.info("Starting calculation for %s sample", sample)

    for category in top_config['main']['categories']:
        observable_fullname = '_'.join([o.split(',')[0].split(':')[0] for o in top_config[category]['observables']])
        logging.info("Category: %s", category)
        if "smooth" not in top_config[category]:
            smooth = False
            logging.debug("Using binned histograms, since smoothing information not found in config file")
        else:
            smooth = top_config[category]["smooth"]

        hist_name = '-'.join([observable_fullname, sample, category])
        logging.debug("Histogram name: %s", hist_name)

        try:
            file_name = top_config['main']['path'] + "/Nominal/" + top_config['samples'][sample]
        except KeyError:
            logging.error("Could not find %s sample filename in config file. Exiting...", sample)
            sys.exit()

        hist = get_hist(top_config['main']['treename'], file_name, top_config[category]['observables'])
