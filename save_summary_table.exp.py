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

	if not os.path.exists("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_inb/dvcs/excl_level_2/restructured_{}".format(suffix)):
		return

	print(suffix)

	for polarity in ["inb", "outb"]:
		for file_directory in ["dvcs", "pi0"]:
			for integrated_binnum in range(1, 147+1):
				print('Reading /volatile/clas12/sangbaek/dvcs_related/exp_fall2018_{}/{}/excl_level_2/restructured_{}/{}.pkl'.format(polarity, file_directory, suffix, integrated_binnum))
				try:
					df = pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_{}/{}/excl_level_2/restructured_{}/{}.pkl'.format(polarity, file_directory, suffix, integrated_binnum))
				except:
					print('/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_{}/{}/excl_level_2/restructured_{}/{}.pkl does not exist.'.format(polarity, file_directory, suffix, integrated_binnum))
					for phi_binnum in range(24):
						directory    = "exp_fall2018_{}/{}".format(polarity, file_directory)
						n_entry = 0
						n_entry_eff_corrected = 0
						n_entry_eff_corrected_err = 0
						n_entry_CDFT = 0
						n_entry_eff_corrected_CDFT = 0
						n_entry_CD = 0
						n_entry_eff_corrected_CD = 0
						n_entry_FD = 0
						n_entry_eff_corrected_FD = 0
						xB_avg  = 0
						Q2_avg  = 0
						t_avg   = 0
						phi_avg = 0
						this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "variation": suffix, "n_entry": n_entry, "n_entry_eff_corrected": n_entry_eff_corrected, "n_entry_eff_corrected_err": n_entry_eff_corrected_err,
						"n_entry_CDFT": n_entry_CDFT, "n_entry_eff_corrected_CDFT": n_entry_eff_corrected_CDFT,
						"n_entry_CD": n_entry_CD, "n_entry_eff_corrected_CD": n_entry_eff_corrected_CD,
						"n_entry_FD": n_entry_FD, "n_entry_eff_corrected_FD": n_entry_eff_corrected_FD,
							"xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
						df_summary   = pd.concat([df_summary, this_row])
					continue

				for phi_binnum in range(24):
					df_this_bin  = df.loc[(df.phi_binnum == phi_binnum)]
					directory    = "exp_fall2018_{}/{}".format(polarity, file_directory)
					if not len(df_this_bin):
						n_entry = 0
						n_entry_eff_corrected = 0
						n_entry_eff_corrected_err = 0
						n_entry_CDFT = 0
						n_entry_eff_corrected_CDFT = 0
						n_entry_CD = 0
						n_entry_eff_corrected_CD = 0
						n_entry_FD = 0
						n_entry_eff_corrected_FD = 0
						xB_avg  = 0
						Q2_avg  = 0
						t_avg   = 0
						phi_avg = 0
					else:
						efficiency       = inverseHist(df_this_bin.efficiency).to_numpy()
						efficiency[efficiency==0] = 1
						n_entry      = len(df_this_bin)
						n_entry_eff_corrected = np.sum(efficiency)
						n_entry_eff_corrected_err = n_entry_eff_corrected/(n_entry) * np.sqrt(n_entry)
						efficiency_CDFT       = inverseHist(df_this_bin.loc[df_this_bin.config==3].efficiency).to_numpy()
						efficiency_CDFT[efficiency_CDFT==0] = 1
						n_entry_CDFT      = len(df_this_bin.loc[df_this_bin.config==3])
						n_entry_eff_corrected_CDFT = np.sum(efficiency_CDFT)
						efficiency_CD       = inverseHist(df_this_bin.loc[df_this_bin.config==2].efficiency).to_numpy()
						efficiency_CD[efficiency_CD==0] = 1
						n_entry_CD      = len(df_this_bin.loc[df_this_bin.config==2])
						n_entry_eff_corrected_CD = np.sum(efficiency_CD)
						efficiency_FD       = inverseHist(df_this_bin.loc[df_this_bin.config==1].efficiency).to_numpy()
						efficiency_FD[efficiency_FD==0] = 1
						n_entry_FD      = len(df_this_bin.loc[df_this_bin.config==1])
						n_entry_eff_corrected_FD = np.sum(efficiency_FD)
						xB_avg      = np.sum(efficiency * df_this_bin.xB)/n_entry_eff_corrected
						Q2_avg      = np.sum(efficiency * df_this_bin.Q2)/n_entry_eff_corrected
						t_avg       = np.sum(efficiency * df_this_bin.t1)/n_entry_eff_corrected
						phi_avg     = np.sum(efficiency * df_this_bin.phi1)/n_entry_eff_corrected

					this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "variation": suffix, "n_entry": n_entry, "n_entry_eff_corrected": n_entry_eff_corrected, "n_entry_eff_corrected_err": n_entry_eff_corrected_err,
						"n_entry_CDFT": n_entry_CDFT, "n_entry_eff_corrected_CDFT": n_entry_eff_corrected_CDFT,
						"n_entry_CD": n_entry_CD, "n_entry_eff_corrected_CD": n_entry_eff_corrected_CD,
						"n_entry_FD": n_entry_FD, "n_entry_eff_corrected_FD": n_entry_eff_corrected_FD,
						"xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
					df_summary   = pd.concat([df_summary, this_row])

	df_summary = df_summary.reset_index()
	df_summary = df_summary.loc[:, df_summary.columns[1:]]
	return df_summary

if __name__ == "__main__":
	
	dfs = []
	for mode in range(6, 14):
		result =  main(mode)
		if result is not None:
			dfs.append(result)
	
	dfs = pd.concat(dfs)
	dfs.to_pickle("summary_table.exp.pkl")
