#!/usr/bin/env python

from argparse import ArgumentParser
import Utilities
import logging.config

parser = ArgumentParser()
parser.add_argument(dest='config_file', help="The input config file", type=str)
args = parser.parse_args()

logging.config.fileConfig('logging.ini')

logging.info("Calculating systematics from config file %s", args.config_file)

top_config = Utilities.config_dict(args.config_file, tokenize=['categories', 'observables'])
logging.debug(Utilities.config_display(top_config))
np_config = Utilities.config_dict(top_config['main']['NPlist'])
logging.debug(Utilities.config_display(np_config))

