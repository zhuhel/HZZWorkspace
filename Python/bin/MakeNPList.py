#!/usr/bin/env python

import Utilities
from argparse import ArgumentParser
import logging.config

parser = ArgumentParser()
parser.add_argument(dest='config_file', help="The input config file", type=str)
args = parser.parse_args()

script_loc = os.path.dirname(os.path.realpath(__file__))
logging.config.fileConfig(script_loc + '/../configuration/logging_screen.ini')
