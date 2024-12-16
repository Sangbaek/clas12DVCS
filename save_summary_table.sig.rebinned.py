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
	df_summary_table_rebinned = pd.read_csv("df_summary_table_rebinned.csv")
	df_summary_table          = pd.DataFrame()
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
			for phi_binnum, phi_width in zip(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"], df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"]):
				directory  = directory
				n_entry = 0
				n_entry_corrected = 0
				n_entry_corrected_err = 0
				n_entry_eff_corrected = 0
				n_entry_eff_corrected_err = 0
				n_entry_eff_bh_corrected = 0
				n_entry_eff_bh_corrected_err = 0
				n_entry_eff_vgg_corrected = 0
				n_entry_eff_vgg_corrected_err = 0
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
				this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "phi_width": phi_width, "directory": directory, "variation": suffix, "n_entry": n_entry,
					"n_entry_corrected": n_entry_corrected, "n_entry_corrected_err": n_entry_corrected_err,
					"n_entry_eff_corrected": n_entry_eff_corrected, "n_entry_eff_corrected_err": n_entry_eff_corrected_err,
					"n_entry_eff_bh_corrected": n_entry_eff_bh_corrected, "n_entry_eff_bh_corrected_err": n_entry_eff_bh_corrected_err,
					"n_entry_eff_vgg_corrected": n_entry_eff_vgg_corrected, "n_entry_eff_vgg_corrected_err": n_entry_eff_vgg_corrected_err,
					"n_entry_CDFT": n_entry_CDFT, "n_entry_corrected_CDFT": n_entry_corrected_CDFT, "n_entry_corrected_err_CDFT": n_entry_corrected_err_CDFT,
					"n_entry_CD": n_entry_CD, "n_entry_corrected_CD": n_entry_corrected_CD, "n_entry_corrected_err_CD": n_entry_corrected_err_CD,
					"n_entry_FD": n_entry_FD, "n_entry_corrected_FD": n_entry_corrected_FD, "n_entry_corrected_err_FD": n_entry_corrected_err_FD,
					 "xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
				df_summary_table   = pd.concat([df_summary_table, this_row])
			continue
		df = df.loc[df.weights<1000, :]
		for phi_binnum, phi_width in zip(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"], df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"]):
			df_this_bin  = df.loc[(df.phi_binnum >= phi_binnum) & (df.phi_binnum < phi_binnum + phi_width), :]
			directory    = directory
			if not len(df_this_bin):
				n_entry = 0
				n_entry_corrected = 0
				n_entry_corrected_err = 0
				n_entry_eff_corrected = 0
				n_entry_eff_corrected_err = 0
				n_entry_eff_bh_corrected = 0
				n_entry_eff_bh_corrected_err = 0
				n_entry_eff_vgg_corrected = 0
				n_entry_eff_vgg_corrected_err = 0
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
				n_entry_eff_corrected = np.sum(df_this_bin.weights*df_this_bin.efficiency)
				n_entry_eff_corrected_err = np.sqrt(np.sum(df_this_bin.efficiency**2 * df_this_bin.weights**2))
				n_entry_eff_bh_corrected = np.sum(df_this_bin.weights*df_this_bin.efficiency_bh)
				n_entry_eff_bh_corrected_err = np.sqrt(np.sum(df_this_bin.efficiency_bh**2 * df_this_bin.weights**2))
				n_entry_eff_vgg_corrected = np.sum(df_this_bin.weights*df_this_bin.efficiency_vgg)
				n_entry_eff_vgg_corrected_err = np.sqrt(np.sum(df_this_bin.efficiency_vgg**2 * df_this_bin.weights**2))
				n_entry_CDFT               = len(df_this_bin.loc[df_this_bin.config==3])
				n_entry_corrected_CDFT     = np.sum(df_this_bin.loc[df_this_bin.config==3].weights)
				n_entry_corrected_err_CDFT = np.sqrt(np.sum(df_this_bin.loc[df_this_bin.config==3].weights**2))
				n_entry_CD                 = len(df_this_bin.loc[df_this_bin.config==2])
				n_entry_corrected_CD       = np.sum(df_this_bin.loc[df_this_bin.config==2].weights)
				n_entry_corrected_err_CD   = np.sqrt(np.sum(df_this_bin.loc[df_this_bin.config==2].weights**2))
				n_entry_FD                 = len(df_this_bin.loc[df_this_bin.config==1])
				n_entry_corrected_FD       = np.sum(df_this_bin.loc[df_this_bin.config==1].weights)
				n_entry_corrected_err_FD   = np.sqrt(np.sum(df_this_bin.loc[df_this_bin.config==1].weights**2))
				xB_avg      = np.sum(df_this_bin.efficiency * df_this_bin.weights * df_this_bin.xB)/n_entry_corrected
				Q2_avg      = np.sum(df_this_bin.efficiency * df_this_bin.weights * df_this_bin.Q2)/n_entry_corrected
				t_avg       = np.sum(df_this_bin.efficiency * df_this_bin.weights * df_this_bin.t1)/n_entry_corrected
				phi_avg     = np.sum(df_this_bin.efficiency * df_this_bin.weights * df_this_bin.phi1)/n_entry_corrected
			this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "phi_width": phi_width, "directory": directory, "variation": suffix, "n_entry": n_entry,
				"n_entry_corrected": n_entry_corrected, "n_entry_corrected_err": n_entry_corrected_err,
				"n_entry_eff_corrected": n_entry_eff_corrected, "n_entry_eff_corrected_err": n_entry_eff_corrected_err,
				"n_entry_eff_bh_corrected": n_entry_eff_bh_corrected, "n_entry_eff_bh_corrected_err": n_entry_eff_bh_corrected_err,
				"n_entry_eff_vgg_corrected": n_entry_eff_vgg_corrected, "n_entry_eff_vgg_corrected_err": n_entry_eff_vgg_corrected_err,
				"n_entry_CDFT": n_entry_CDFT, "n_entry_corrected_CDFT": n_entry_corrected_CDFT, "n_entry_corrected_err_CDFT": n_entry_corrected_err_CDFT,
				"n_entry_CD": n_entry_CD, "n_entry_corrected_CD": n_entry_corrected_CD, "n_entry_corrected_err_CD": n_entry_corrected_err_CD,
				"n_entry_FD": n_entry_FD, "n_entry_corrected_FD": n_entry_corrected_FD, "n_entry_corrected_err_FD": n_entry_corrected_err_FD,
				 "xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])

			df_summary_table   = pd.concat([df_summary_table, this_row])


	df_summary_table = df_summary_table.reset_index()
	df_summary_table = df_summary_table.loc[:, df_summary_table.columns[1:]]
	return df_summary_table

if __name__ == "__main__":
	dfs = []
	for model in ["dvcs_km15", "dvcs_vgg", "pureBH"]:
		for polarity in ["inb", "outb"]:
			for mode in range(6, 14):
				result =  main(mode, model, polarity)
				if result is not None:
					dfs.append(result)
	
	dfs = pd.concat(dfs)
	dfs.to_pickle("summary_table.sig.rebinned.pkl")

