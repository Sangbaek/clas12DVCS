from glob import glob
import argparse
import subprocess
import itertools


parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
args   = parser.parse_args()
parser.add_argument("-s","--step", help="what to do?", default = None)

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

this_step_dir  = "collecting_dsts/sim/"
for config, i in (itertools.product(configs.items(), [1,2,3,4])):
  mode, runs = config
  for run in runs:
    subprocess.run(['cp','{}/{}/.run_{}'.format(default_dir, this_step_dir, i), '{}/{}/{}_{}'.format(default_dir, this_step_dir, run, i)])
    subprocess.run(['sed', '-i', 's/run/{}/g'.format(run), '{}/{}/{}_{}'.format(default_dir, this_step_dir, run, i)])
