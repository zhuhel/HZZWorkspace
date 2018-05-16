#!/usr/bin/env python

import Utilities, sys, os, logging.config
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument(dest='config_file', help="The input config file", type=str)
args = parser.parse_args()

script_loc = os.path.dirname(os.path.realpath(__file__))
logging.config.fileConfig(script_loc + '/../configuration/logging_screen.ini')

logging.info("Loading ROOT")
try:
    import ROOT
except ImportError:
    logging.error("ROOT is needed to open the sample files. Please set up Python with ROOT")
    sys.exit()

logging.info("Using config file %s to make nuisance parameter list", args.config_file)

top_config = Utilities.config_dict(args.config_file, tokenize=[('samples', ',')])
try:
    np_list_filename = top_config['main']['NPlist']
except KeyError:
    logging.error("Please specify 'NPlist' in the config file to use as output filename and try again")
    sys.exit()

try:
    sample_dir = top_config['main']['path']
    sample_names = top_config['main']['samples']
    sys_dir = top_config['main']['sysDir']
    tree_name = top_config['main']['treename']
except KeyError:
    logging.error("Please specify some samples, their path, and/or the name of their tree in the config file")

chosen_sample = None
for s in sample_names:
    try:
        sample_filename = top_config['samples'][s]
    except KeyError:
        logging.warn("No file name given for sample %s", s)
        continue
    if os.path.isfile("{}/{}/NormSystematic/{}".format(sample_dir, sys_dir, sample_filename)):
        chosen_sample = s
        break
    else:
        logging.warn("{}/{}/NormSystematic/{} does not exist".format(sample_dir, sys_dir, sample_filename))

if not chosen_sample:
    logging.error("Could not find a sample to use for NPs. Check the config file and try again")
    sys.exit()
else:
    logging.info("Using sample %s to find weight names", chosen_sample)

logging.info("Using output filename: %s", np_list_filename)
if os.path.isfile(np_list_filename):
    check = raw_input("File %s already exists... Continue? (y/n) " % np_list_filename).lower()
    if check != 'y':
        logging.info("Exiting")
        sys.exit()
    else:
        logging.info("Overwriting existing file")

output_file = open(np_list_filename, 'w')
output_file.write("[shapeLike]\n")

from glob import glob
shapelike_list = [d.split("/")[-1] for d in glob("{}/{}/*".format(sample_dir, sys_dir))]
try:
    shapelike_list.remove("NormSystematic")
except ValueError:
    pass
for np in sorted(shapelike_list):
    output_file.write("{} = ATLAS_{}\n".format(np, np.split("__")[0]))

output_file.write("[normLike]\n")
sample_root_file = ROOT.TFile("{}/{}/NormSystematic/{}".format(sample_dir, sys_dir, top_config['samples'][chosen_sample]), "READ")
sample_tree = sample_root_file.Get(tree_name)
normlike_list = [b.GetName() for b in sample_tree.GetListOfBranches() if b.GetName().startswith("weight_var")]
for np in sorted(normlike_list):
    output_file.write("{} = {}\n".format(np, np.split("__")[0].replace("weight_var", "ATLAS")))

sample_root_file.Close()
output_file.close()
