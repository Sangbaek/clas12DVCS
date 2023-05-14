from glob import glob
import argparse
import subprocess
import itertools
import numpy as np


parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-s","--step", help="what to do?", default = None, type = int)
args   = parser.parse_args()

configs = {
  "bh/inb"  : [6043],
  "bh/outb" : [6066],
  "vgg/inb" : [6044],
  "vgg/outb": [6067],
  "pi0/inb" : [6047],
  "pi0/outb": [6069]
}

default_dir    = "/volatile/clas12/sangbaek/jan2023/sbatch_files"

#step 1.
#merge simulations.
#filter files into a hipo files that have

if args.step == 1:
  this_step_dir  = "collecting_dsts/sim/"
  for config, entry in (itertools.product(configs.items(), [1,2,3,4])):
    mode, runs = config
    for run in runs:
      subprocess.run(['cp','{}/{}/.run_{}'.format(default_dir, this_step_dir, entry), '{}/{}/{}_{}'.format(default_dir, this_step_dir, run, entry)])
      subprocess.run(['sed', '-i', 's/run/{}/g'.format(run), '{}/{}/{}_{}'.format(default_dir, this_step_dir, run, entry)])

#step 2
#convert that into root files
if args.step == 2:
  this_step_dir  = "conversion/"
  entries        = np.linspace(0, 19, 20, dtype = int)
  for config, entry in (itertools.product(configs.items(), entries)):
    mode, runs = config
    for run in runs:
      subprocess.run(['cp','{}/{}/{}/.run_0'.format(default_dir, this_step_dir, mode, entry), '{}/{}/{}/{}_{}'.format(default_dir, this_step_dir, mode, run, entry)])
      subprocess.run(['sed', '-i', 's/run/{}/g;s/_0/_{}/g'.format(run, entry), '{}/{}/{}/{}_{}'.format(default_dir, this_step_dir, mode, run, entry)])
