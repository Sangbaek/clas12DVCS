import numpy as np
import pandas as pd
from utils.const import *
from utils.physics import *

import argparse
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
warnings.simplefilter(action='ignore', category=pd.errors.SettingWithCopyWarning)


def main(mode, model, polarity):
	df_summary = pd.DataFrame()
	suffix = schema_suffices[mode]

	directory = "{}/fall2018_{}".format(model, polarity)
	location = ""
	if model == "pureBH" and polarity == "inb":
		directory = "{}/fall2018_{}3".format(model, polarity)
		location = "_3"
	if model == "dvcs_km15" and mode <13:
		directory = "{}/fall2018_{}3".format(model, polarity)
	if model == "dvcs_km15" and mode == 13 and polarity == "inb":
		directory = "{}/fall2018_{}3_45nA".format(model, polarity)
	if model == "dvcs_km15" and mode == 13 and polarity == "outb":
		directory = "{}/fall2018_{}3_50nA".format(model, polarity)


	if not os.path.exists("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_{}/{}/excl_level_1/restructured{}_{}".format(polarity, model, location, suffix)):
		print("{} {} invalid.".format(model, suffix))
		return

	print(suffix)

	for integrated_binnum in range(1, 147+1):
		print('Reading /volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_{}/{}/excl_level_1/restructured{}_{}/{}.pkl'.format(polarity, model, location, suffix, integrated_binnum))

		try:
			df = pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_{}/{}/excl_level_1/restructured{}_{}/{}.pkl'.format(polarity, model, location, suffix, integrated_binnum))
		except:
			print('/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_{}/{}/excl_level_1/restructured{}_{}/{}.pkl does not exist'.format(polarity, model, location, suffix, integrated_binnum))
			for phi_binnum in range(24):
				directory  = directory
				n_entry = 0
				n_entry_corrected = 0
				n_entry_corrected_err = 0
				n_entry_CDFT = 0
				n_entry_corrected_CDFT = 0
				n_entry_corrected_err_CDFT = 0
				n_entry_CD = 0
				n_entry_corrected_CD = 0
				n_entry_corrected_err_CD = 0
				n_entry_FD = 0
				n_entry_corrected_FD = 0
				n_entry_corrected_err_FD = 0
				xB_avg = 0
				Q2_avg = 0
				t_avg  = 0
				phi_avg = 0
				this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "variation": suffix, "n_entry": n_entry,
					"n_entry_corrected": n_entry_corrected, "n_entry_corrected_err": n_entry_corrected_err,
					"n_entry_CDFT": n_entry_CDFT, "n_entry_corrected_CDFT": n_entry_corrected_CDFT, "n_entry_corrected_err_CDFT": n_entry_corrected_err_CDFT,
					"n_entry_CD": n_entry_CD, "n_entry_corrected_CD": n_entry_corrected_CD, "n_entry_corrected_err_CD": n_entry_corrected_err_CD,
					"n_entry_FD": n_entry_FD, "n_entry_corrected_FD": n_entry_corrected_FD, "n_entry_corrected_err_FD": n_entry_corrected_err_FD,
					 "xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
				df_summary   = pd.concat([df_summary, this_row])
			continue
		for phi_binnum in range(24):
			df_this_bin  = df.loc[(df.phi_binnum == phi_binnum)]
			directory    = directory
			if not len(df_this_bin):
				n_entry = 0
				n_entry_corrected = 0
				n_entry_corrected_err = 0
				n_entry_CDFT = 0
				n_entry_corrected_CDFT = 0
				n_entry_corrected_err_CDFT = 0
				n_entry_CD = 0
				n_entry_corrected_CD = 0
				n_entry_corrected_err_CD = 0
				n_entry_FD = 0
				n_entry_corrected_FD = 0
				n_entry_corrected_err_FD = 0
				xB_avg  = 0
				Q2_avg  = 0
				t_avg   = 0
				phi_avg = 0
			else:
				n_entry      = len(df_this_bin)
				n_entry_corrected = np.sum(df_this_bin.weights)
				n_entry_corrected_err = np.sqrt(np.sum(df_this_bin.weights**2))
				n_entry_CDFT               = len(df_this_bin.loc[df_this_bin.config==3])
				n_entry_corrected_CDFT     = np.sum(df_this_bin.loc[df_this_bin.config==3].weights)
				n_entry_corrected_err_CDFT = np.sqrt(np.sum(df_this_bin.loc[df_this_bin.config==3].weights**2))
				n_entry_CD                 = len(df_this_bin.loc[df_this_bin.config==2])
				n_entry_corrected_CD       = np.sum(df_this_bin.loc[df_this_bin.config==2].weights)
				n_entry_corrected_err_CD   = np.sqrt(np.sum(df_this_bin.loc[df_this_bin.config==2].weights**2))
				n_entry_FD                 = len(df_this_bin.loc[df_this_bin.config==1])
				n_entry_corrected_FD       = np.sum(df_this_bin.loc[df_this_bin.config==1].weights)
				n_entry_corrected_err_FD   = np.sqrt(np.sum(df_this_bin.loc[df_this_bin.config==1].weights**2))
				xB_avg      = np.sum(df_this_bin.weights * df_this_bin.xB)/n_entry_corrected
				Q2_avg      = np.sum(df_this_bin.weights * df_this_bin.Q2)/n_entry_corrected
				t_avg       = np.sum(df_this_bin.weights * df_this_bin.t1)/n_entry_corrected
				phi_avg     = np.sum(df_this_bin.weights * df_this_bin.phi1)/n_entry_corrected
			this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "variation": suffix, "n_entry": n_entry,
				"n_entry_corrected": n_entry_corrected, "n_entry_corrected_err": n_entry_corrected_err,
				"n_entry_CDFT": n_entry_CDFT, "n_entry_corrected_CDFT": n_entry_corrected_CDFT, "n_entry_corrected_err_CDFT": n_entry_corrected_err_CDFT,
				"n_entry_CD": n_entry_CD, "n_entry_corrected_CD": n_entry_corrected_CD, "n_entry_corrected_err_CD": n_entry_corrected_err_CD,
				"n_entry_FD": n_entry_FD, "n_entry_corrected_FD": n_entry_corrected_FD, "n_entry_corrected_err_FD": n_entry_corrected_err_FD,
				 "xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])

			df_summary   = pd.concat([df_summary, this_row])


	df_summary = df_summary.reset_index()
	df_summary = df_summary.loc[:, df_summary.columns[1:]]
	return df_summary

if __name__ == "__main__":
	dfs = []
	for model in ["dvcs_km15", "dvcs_vgg", "pureBH"]:
		for polarity in ["inb", "outb"]:
			for mode in range(14):
				result =  main(mode, model, polarity)
				if result is not None:
					dfs.append(result)
	
	dfs = pd.concat(dfs)
	dfs.to_pickle("summary_table.sig.pkl")

