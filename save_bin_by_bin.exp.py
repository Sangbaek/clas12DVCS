import numpy as np
import pandas as pd
from utils.const import *
from utils.physics import *
from utils.fiducial import *
import os

import argparse
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
warnings.simplefilter(action='ignore', category=pd.errors.SettingWithCopyWarning)


def main(mode):

	suffix = schema_suffices[mode]

	if not os.path.exists("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_inb/{}/excl_level_2/pkl_{}".format("dvcs", suffix)):
		return

	print(suffix)

	binnum, bin_volume = np.loadtxt('volume_list.csv', skiprows = 1, delimiter = ',').T
	bin_volume = {int(binnum[i]): bin_volume[i] for i in range(len(binnum))}

	#Exp - fall 2018 inbending ep->epg
	#for directory in ["dvcs", "pi0"]:
	for directory in [ "pi0"]:
		for polarity in ["inb", "outb"]:
			if polarity == "inb":
				qadb_range = 339
				runlist    = runlist_inb
			else:
				qadb_range = 346
				runlist    = runlist_outb
			pol      = "{}ending".format(polarity)
			qaTree     = pd.read_json("/volatile/clas12/sangbaek/clasqaDB/qadb/qa.rga_{}ending/qaTree.json".format(polarity)).T
			chargeTree = pd.read_json("/volatile/clas12/sangbaek/clasqaDB/qadb/qa.rga_{}ending/chargeTree.json".format(polarity)).T
			df_merged = []
			for runnum in runlist:
				print(runnum)
				try:
					df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_{}/{}/excl_level_2/pkl_{}/{}".format(polarity, directory, suffix, runnum))
				except:
					continue
				for i in range(qadb_range):
					filenum    = 5*i
					qadb       = qaTree.loc[qaTree.index == runnum, filenum]
					if not isinstance(qadb.values[0], dict):
						break    
					chargedb   = chargeTree.loc[qaTree.index == runnum, filenum]
					evnumMin   = qadb.values[0]['evnumMin']
					if (filenum>0) and (evnumMin == 0):
						evnumMin = evnumMax+1
					evnumMax   = qadb.values[0]['evnumMax']
					defect     = qadb.values[0]['defect']
					chargemin  = chargedb.values[0]['fcChargeMin']
					chargemax  = chargedb.values[0]['fcChargeMax']
					df.loc[(df.EventNum>=evnumMin) & (df.EventNum<=evnumMax), "defect"]  = defect
					df.loc[(df.EventNum>=evnumMin) & (df.EventNum<=evnumMax), "chargemin"] = chargemin
					df.loc[(df.EventNum>=evnumMin) & (df.EventNum<=evnumMax), "chargemax"] = chargemax
				df = df.loc[df.defect == 0, :]
				df = df.drop(columns = ["defect"])
				df_merged.append(df)
			df_merged = pd.concat(df_merged)
			df_merged = df_merged.reset_index()
			df_merged = df_merged.loc[:, df_merged.columns[1:]]
			df_merged = assign_efficiency(df_merged, mc = False, pol = pol)
			df_merged.to_pickle('/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_{}/{}/excl_level_2/pkl_{}/fall2018_{}.pkl'.format(polarity, directory, suffix, polarity))
		
			os.makedirs('/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_{}/{}/excl_level_2/restructured_{}'.format(polarity, directory, suffix), exist_ok=True)
			for integrated_binnum in range(1, 147+1):
				df_merged.loc[df_merged.integrated_binnum == integrated_binnum].to_pickle('/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_{}/{}/excl_level_2/restructured_{}/{}.pkl'.format(polarity, directory, suffix, integrated_binnum))
	return


if __name__ == "__main__":

	#for mode in range(6, 13):#range(1, 14):
	for mode in range(6, 13):#range(1, 14):
		main(mode)


