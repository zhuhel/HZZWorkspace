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
except ImportError:
    logging.error("Could not import ROOT module. Make sure your root is configured to work with Python.")
    sys.exit()


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
        logging.info("Category: %s", category)
        if "smooth" not in top_config[category]:
            smooth = False
            logging.debug("Using binned histograms, since smoothing information not found in config file")
        else:
            smooth = top_config[category]["smooth"]

        hist_name = '-'.join([top_config[category]['observables'][0].split(':')[0], sample, category])
        logging.debug("Histogram name: %s", hist_name)

        try:
            file_name = top_config['main']['path'] + "/Nominal/" + top_config['samples'][sample]
        except KeyError:
            logging.error("Could not find %s sample filename in config file. Exiting...", sample)
            sys.exit()
        logging.info("Making chain with name %s", top_config['main']['treename'])
        chain = ROOT.TChain(top_config['main']['treename'])
        logging.info("Adding file %s", file_name)  # could make this into a loop over files
        chain.Add(file_name)
        if chain.GetEntries() == 0:
            logging.warn("Chain %s has no events", top_config['main']['treename'])
        else:
            logging.info("Chain %s has %i events", top_config['main']['treename'], chain.GetEntries())
