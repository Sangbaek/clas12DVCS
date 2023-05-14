from glob import glob
import argparse
import subprocess
import itertools
import numpy as np


parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-s","--step", help="what to do?", default = 0, type = int)
parser.add_argument("-d","--data", help="sim or exp", default = 'sim', type = str)
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
step           = args.step
this_step_dirs = {1: "collecting_dsts/sim/", 3: "filtering/", 4: "conversion/", 5: "final/", 6: "Gen/", 7: "Gen2/"}

if args.data == "sim":
  entries        = np.linspace(0, 19, 20, dtype = int)

  #step 1.
  #merge simulations.
  if step == 1:
    this_step_dir = this_step_dirs[step]
    for config, entry in (itertools.product(configs.items(), [1,2,3,4])):
      mode, runs = config
      for run in runs:
        subprocess.run(['cp','{}/{}/.run_{}'.format(default_dir, this_step_dir, entry), '{}/{}/{}_{}'.format(default_dir, this_step_dir, run, entry)])
        subprocess.run(['sed', '-i', 's/run/{}/g'.format(run), '{}/{}/{}_{}'.format(default_dir, this_step_dir, run, entry)])

  #step 2.
  #cp the simulation files to the temp directory
  elif step == 2:
    for config, entry in (itertools.product(configs.items(), entries)):
      mode, runs = config
      for run in runs:
        subprocess.run(['mv','/volatile/clas12/sangbaek/jan2023/to_be_filtered/{}_{}.hipo'.format(run, entry), '/volatile/clas12/sangbaek/regular_backup/temp/{}/{}_{}.hipo'.format(mode, run, entry)])


  #step 3, 4, 5: filtering to hipo files that have epg, 4: convert to root files, 5: convert to pkl files.
  #convert that into root files
  else:
    this_step_dir = this_step_dirs[step]
    for config, entry in (itertools.product(configs.items(), entries)):
      mode, runs = config
      for run in runs:
        subprocess.run(['cp','{}/{}/{}/.run_0'.format(default_dir, this_step_dir, mode, entry), '{}/{}/{}/{}_{}'.format(default_dir, this_step_dir, mode, run, entry)])
        subprocess.run(['sed', '-i', 's/run/{}/g;s/_0/_{}/g'.format(run, entry), '{}/{}/{}/{}_{}'.format(default_dir, this_step_dir, mode, run, entry)])

elif args.data == "exp":
  #step 3
  #filter files into a hipo files that have epg
  if step == 3:
    #DVCS wagon inb
    this_step_dir = this_step_dirs[step]
    file_locations = {
      "epg/inb": "/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v1/dst/train/DVCSWagon/",
      "epg/outb": "/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass1/v1/dst/train/DVCSWagon/",
      "epgg/inb": "/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v1/dst/train/DVPi0P/",
      "epgg/outb": "/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass1/v1/dst/train/DVPi0P/",
    }
    for mode in file_locations:
      file_location = file_locations[mode]
      files = glob("{}/*".format(file_location))
      for file in files:
        run = file.split(file_location)[-1][-9:-5]
        subprocess.run(['cp','{}/{}/exp/{}/.run'.format(default_dir, this_step_dir, mode), '{}/{}/exp/{}/{}'.format(default_dir, this_step_dir, mode, run)])

