import numpy as np
import pandas as pd
from utils.const import *
from utils.physics import *

import argparse
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
warnings.simplefilter(action='ignore', category=pd.errors.SettingWithCopyWarning)


def main(mode):
	df_summary = pd.DataFrame()
	suffix = schema_suffices[mode]

	if not os.path.exists("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/pi0_1gamma/excl_level_1/restructured_{}".format(suffix)):
		return

	print(suffix)

	for chunk in [1, 2, 3]:
		for integrated_binnum in range(1, 147+1):
			print('Reading /volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/pi0_1gamma/excl_level_1/restructured_{}/{}/{}.pkl'.format(suffix, chunk, integrated_binnum))
			try:
				df = pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/pi0_1gamma/excl_level_1/restructured_{}/{}/{}.pkl'.format(suffix, chunk, integrated_binnum))
			except:
				print('/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/pi0_1gamma/excl_level_1/restructured_{}/{}/{}.pkl does not exist'.format(suffix, chunk, integrated_binnum))
				for phi_binnum in range(24):
					directory    = "sim_rad_rec_fall2018_inb/pi0_1gamma/{}".format(chunk)
					n_entry = 0
					xB_avg  = 0
					Q2_avg  = 0
					t_avg   = 0
					phi_avg = 0
					this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "variation": suffix, "n_entry": n_entry, 
					"xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
					df_summary   = pd.concat([df_summary, this_row])
				continue
			for phi_binnum in range(24):
				df_this_bin  = df.loc[(df.phi_binnum == phi_binnum)]
				directory    = "sim_rad_rec_fall2018_inb/pi0_1gamma/{}".format(chunk)
				if not len(df_this_bin):
					n_entry = 0
					xB_avg  = 0
					Q2_avg  = 0
					t_avg   = 0
					phi_avg = 0
				else:
					n_entry      = len(df_this_bin)
					xB_avg      = np.mean(df_this_bin.xB)
					Q2_avg      = np.mean(df_this_bin.Q2)
					t_avg       = np.mean(df_this_bin.t1)
					phi_avg     = np.mean(df_this_bin.phi1)
				this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "variation": suffix, "n_entry": n_entry, 
					"xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
				df_summary   = pd.concat([df_summary, this_row])

		for integrated_binnum in range(1, 147+1):
			print('Reading /volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/pi0_2gamma/excl_level_1/restructured_{}/{}/{}.pkl'.format(suffix, chunk, integrated_binnum))
			try:
				df = pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/pi0_2gamma/excl_level_1/restructured_{}/{}/{}.pkl'.format(suffix, chunk, integrated_binnum))
			except:
				print('/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/pi0_2gamma/excl_level_1/restructured_{}/{}/{}.pkl does not exist'.format(suffix, chunk, integrated_binnum))
				for phi_binnum in range(24):
					directory    = "sim_rad_rec_fall2018_inb/pi0_2gamma/{}".format(chunk)
					n_entry = 0
					xB_avg  = 0
					Q2_avg  = 0
					t_avg   = 0
					phi_avg = 0
					this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "variation": suffix, "n_entry": n_entry, 
					"xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
					df_summary   = pd.concat([df_summary, this_row])
				continue
			for phi_binnum in range(24):
				df_this_bin  = df.loc[(df.phi_binnum == phi_binnum)]
				directory    = "sim_rad_rec_fall2018_inb/pi0_2gamma/{}".format(chunk)
				if not len(df_this_bin):
					n_entry = 0
					xB_avg  = 0
					Q2_avg  = 0
					t_avg   = 0
					phi_avg = 0
				else:
					n_entry      = len(df_this_bin)
					xB_avg      = np.mean(df_this_bin.xB)
					Q2_avg      = np.mean(df_this_bin.Q2)
					t_avg       = np.mean(df_this_bin.t1)
					phi_avg     = np.mean(df_this_bin.phi1)
				this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "variation": suffix, "n_entry": n_entry, 
					"xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
				df_summary   = pd.concat([df_summary, this_row])


		for integrated_binnum in range(1, 147+1):
			print('Reading /volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/pi0_1gamma/excl_level_1/restructured_{}/{}/{}.pkl'.format(suffix, chunk, integrated_binnum))
			try:
				df = pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/pi0_1gamma/excl_level_1/restructured_{}/{}/{}.pkl'.format(suffix, chunk, integrated_binnum))
			except:
				print('/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/pi0_1gamma/excl_level_1/restructured_{}/{}/{}.pkl does not exist'.format(suffix, chunk, integrated_binnum))
				for phi_binnum in range(24):
					directory    = "sim_rad_rec_fall2018_outb/pi0_1gamma/{}".format(chunk)
					n_entry = 0
					xB_avg  = 0
					Q2_avg  = 0
					t_avg   = 0
					phi_avg = 0
					this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "variation": suffix, "n_entry": n_entry, 
					"xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
					df_summary   = pd.concat([df_summary, this_row])
				continue
			for phi_binnum in range(24):
				df_this_bin  = df.loc[(df.phi_binnum == phi_binnum)]
				directory    = "sim_rad_rec_fall2018_outb/pi0_1gamma/{}".format(chunk)
				if not len(df_this_bin):
					n_entry = 0
					xB_avg  = 0
					Q2_avg  = 0
					t_avg   = 0
					phi_avg = 0
				else:
					n_entry      = len(df_this_bin)
					xB_avg      = np.mean(df_this_bin.xB)
					Q2_avg      = np.mean(df_this_bin.Q2)
					t_avg       = np.mean(df_this_bin.t1)
					phi_avg     = np.mean(df_this_bin.phi1)
				this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "variation": suffix, "n_entry": n_entry, 
					"xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
				df_summary   = pd.concat([df_summary, this_row])

		for integrated_binnum in range(1, 147+1):
			print('Reading /volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/pi0_2gamma/excl_level_1/restructured_{}/{}/{}.pkl'.format(suffix, chunk, integrated_binnum))
			try:
				df = pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/pi0_2gamma/excl_level_1/restructured_{}/{}/{}.pkl'.format(suffix, chunk, integrated_binnum))
			except:
				print('/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/pi0_2gamma/excl_level_1/restructured_{}/{}/{}.pkl does not exist'.format(suffix, chunk, integrated_binnum))
				for phi_binnum in range(24):
					directory    = "sim_rad_rec_fall2018_outb/pi0_2gamma/{}".format(chunk)
					n_entry = 0
					xB_avg  = 0
					Q2_avg  = 0
					t_avg   = 0
					phi_avg = 0
					this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "variation": suffix, "n_entry": n_entry, 
					"xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
					df_summary   = pd.concat([df_summary, this_row])
				continue
			for phi_binnum in range(24):
				df_this_bin  = df.loc[(df.phi_binnum == phi_binnum)]
				directory    = "sim_rad_rec_fall2018_outb/pi0_2gamma/{}".format(chunk)
				if not len(df_this_bin):
					n_entry = 0
					xB_avg  = 0
					Q2_avg  = 0
					t_avg   = 0
					phi_avg = 0
				else:
					n_entry      = len(df_this_bin)
					xB_avg      = np.mean(df_this_bin.xB)
					Q2_avg      = np.mean(df_this_bin.Q2)
					t_avg       = np.mean(df_this_bin.t1)
					phi_avg     = np.mean(df_this_bin.phi1)
				this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "variation": suffix, "n_entry": n_entry, 
					"xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
				df_summary   = pd.concat([df_summary, this_row])

	df_summary = df_summary.reset_index()
	df_summary = df_summary.loc[:, df_summary.columns[1:]]
	return df_summary

if __name__ == "__main__":

	dfs = []
	for mode in range(14):
		result =  main(mode)
		if result is not None:
			dfs.append(result)
	
	dfs = pd.concat(dfs)
	dfs.to_pickle("summary_table.bkg.pkl")


