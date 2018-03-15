#!/usr/bin/env python

from ConfigParser import SafeConfigParser
from argparse import ArgumentParser
import Utilities
from datetime import datetime

parser = ArgumentParser()
parser.add_argument(dest='config_file', help="The input config file", type=str)
parser.add_argument('--log', help="The output log file name", type=str,
                    default=datetime.now().strftime("log_%Y-%m-%d_%H:%M:%S.txt"))
args = parser.parse_args()

log = Utilities.Logger
#log.info('jknasd')

top_config = SafeConfigParser()
top_config.readfp(open(args.config_file))

print args.log

for s in top_config.sections():
    for i in top_config.items(s):
        print i[0], i[1]
