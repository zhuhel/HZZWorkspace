#!/usr/bin/env python

from ConfigParser import SafeConfigParser
from argparse import ArgumentParser
import Utilities
from datetime import datetime
import pprint

parser = ArgumentParser()
parser.add_argument(dest='config_file', help="The input config file", type=str)
parser.add_argument('--log', help="The output log file name", type=str,
                    default=datetime.now().strftime("log_%Y-%m-%d_%H:%M:%S.txt"))
args = parser.parse_args()

top_config = Utilities.configDict(args.config_file, tokenize=['categories', 'observables'])
np_config = Utilities.configDict(top_config['main']['NPlist'])


