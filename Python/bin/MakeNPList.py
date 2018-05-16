#!/usr/bin/env python

import Utilities, sys, os, logging.config
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument(dest='config_file', help="The input config file", type=str)
args = parser.parse_args()

script_loc = os.path.dirname(os.path.realpath(__file__))
logging.config.fileConfig(script_loc + '/../configuration/logging_screen.ini')

logging.info("Using config file %s to make nuisance parameter list", args.config_file)

top_config = Utilities.config_dict(args.config_file, tokenize=[('samples', ',')])
try:
    np_list_filename = top_config['main']['NPlist']
except KeyError:
    logging.error("Please specify 'NPlist' in the config file to use as output filename and try again")
    sys.exit()

logging.info("Using output filename: %s", np_list_filename)

