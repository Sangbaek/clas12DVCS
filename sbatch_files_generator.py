from glob import glob
import argparse
import subprocess


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

default_dir = "/volaite/clas12/sangbaek/jan2023/sbatch_files"

#step 1.
#merge simulations.
#filter files into a hipo files that have

mode        = "collecting_dsts"
for key, val in configs.items():
  for i in [1,2,3,4]:
    subprocess.run(['cp','{}/{}/.run_{}'.format(default_dir, mode, i), '{}/{}/{}_{}'.format(default_dir, mode, val, i)])
    subprocess.run(['sed', '-i', 's/run/{}/g'.format(run), '{}/{}/{}_{}'.format(default_dir, mode, val, i)])
