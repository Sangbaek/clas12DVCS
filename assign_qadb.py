import pandas as pd
import numpy as np
import argparse
from copy import copy

parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-p","--polarity", help="polarity: inb or outb", default = "inb", type = str)
parser.add_argument("-r","--run", help="run number", default = 5032, type = int)
args  = parser.parse_args()

run      = args.run
polarity = args.polarity


qaTree     = pd.read_json("/volatile/clas12/sangbaek/clasqaDB/qadb/qa.rga_{}ending/qaTree.json".format(polarity)).T
chargeTree = pd.read_json("/volatile/clas12/sangbaek/clasqaDB/qadb/qa.rga_{}ending/chargeTree.json".format(polarity)).T
filename   = "/volatile/clas12/sangbaek/jan2023/convPkl/exp/epg/{}/{}.pkl".format(polarity, run)
file       = pd.read_pickle(filename)

for i in range(339):
  filenum    = 5*i
  qadb       = qaTree.loc[qaTree.index == run, filenum]
  if not isinstance(qadb.values[0], dict):
    break    
  chargedb   = chargeTree.loc[qaTree.index == run, filenum]
  evnumMin   = qadb.values[0]['evnumMin']
  evnumMax   = qadb.values[0]['evnumMax']
  defect     = qadb.values[0]['defect']
  chargemin  = chargedb.values[0]['fcChargeMin']
  chargemax  = chargedb.values[0]['fcChargeMax']
  file.loc[(file.EventNum>=evnumMin) & (file.EventNum<=evnumMax), "filenum"] = filenum
  file.loc[(file.EventNum>=evnumMin) & (file.EventNum<=evnumMax), "defect"]  = defect
  file.loc[(file.EventNum>=evnumMin) & (file.EventNum<=evnumMax), "chargemin"] = chargemin
  file.loc[(file.EventNum>=evnumMin) & (file.EventNum<=evnumMax), "chargemax"] = chargemax

file.to_pickle(filename)