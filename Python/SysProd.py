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
    import ROOT
except ImportError:
    logging.error("Could not import ROOT module. Make sure your root is configured to work with Python.")
    sys.exit()


def get_hist(h_name, sample, category, file_name, observables, cuts, weight_label, smooth):
    t_cut = ROOT.TChain()

    return_hist = None
    return return_hist


logging.info("Calculating systematics from config file %s", args.config_file)

top_config = Utilities.config_dict(args.config_file, tokenize=['samples', 'categories', 'observables'])
logging.debug(Utilities.config_display(top_config))
np_config = Utilities.config_dict(top_config['main']['NPlist'])
logging.debug(Utilities.config_display(np_config))

if 'outdir' in top_config['main']:
    logging.warn("Output directory %s found in config file, replacing the command line argument of %s",
                 top_config['main']['outdir'], args.output_dir)

if not Utilities.check_config(top_config, ['categories', 'samples', 'path', 'treename', 'NPlist']):
    logging.error("Config file does not contain required information. Exiting...")
    sys.exit()

if not Utilities.check_np(np_config):
    logging.error("NP list file does not contain even 'up' and 'down' variations. Exiting...")
    sys.exit()

for sample in top_config['main']['samples']:
    logging.info("Starting calculation for %s sample", sample)

    for category in top_config['main']['categories']:
        logging.info("Category: %s", category)
        if "smooth" not in top_config[category]:
            smooth = False
            logging.info("Using binned histograms, since smoothing information not found in config file")
        else:
            smooth = top_config[category]["smooth"]

        hist_name = '-'.join([top_config[category]['observables'][0].split(':')[0], sample, category])
        logging.debug("Histogram name: %s", hist_name)