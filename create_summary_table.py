import numpy as np
import pandas as pd

import argparse

df_summary = pd.DataFrame()

#Gen - bulk
for bin in range(1, 147+1):
	df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_data/sim_rad_gen/dvcs_km15/pkl/fall2018_inb3/{}.pkl".format(binnum))
	df.loc[:, "integrated_binnum_gen"] = bin
	df = df.astype({"integrated_binnum_gen": int, "phi_binnum_gen": int})
	for phi_bin in range(24):
		df_this_bin  = df.loc[df.phi_binnum_gen == phi_bin]
		generator    = "km15gen"
		polarity     = "inbending"
		n_entry      = len(df_this_bin)
		weight_sum   = np.sum(df_this_bin.GenWeight)
		weight_mean  = np.mean(df_this_bin.GenWeight)
		xB_mean      = np.sum(df_this_bin.GenWeight * df_this_bin.GenxB)/weight_sum
		Q2_mean      = np.sum(df_this_bin.GenWeight * df_this_bin.GenQ2)/weight_sum
		t_mean       = np.sum(df_this_bin.GenWeight * df_this_bin.Gent)/weight_sum
		phi_mean     = np.sum(df_this_bin.GenWeight * df_this_bin.Genphi)/weight_sum
		this_row     = pd.DataFrame([{"generator": generator, "polarity": polarity, "n_entry": n_entry,
			"weight_sum": weight_sum, "weight_mean": weight_mean, "xB_mean": xB_mean,
			"Q2_mean": Q2_mean, "t_mean": t_mean, "phi_mean": phi_mean}])
		df_summary   = pd.concat([df_summary, this_row])

for bin in range(1, 159+1):
	df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_data/sim_rad_gen/dvcs_km15_fringe/pkl/fall2018_inb3/{}.pkl".format(binnum))
	df.loc[:, "integrated_binnum_gen"] = bin + 147
	df = df.astype({"integrated_binnum_gen": int, "phi_binnum_gen": int})
	for phi_bin in range(24):
		df_this_bin  = df.loc[df.phi_binnum_gen == phi_bin]
		generator    = "km15gen"
		polarity     = "inbending"
		n_entry      = len(df_this_bin)
		weight_sum   = np.sum(df_this_bin.GenWeight)
		weight_mean  = np.mean(df_this_bin.GenWeight)
		xB_mean      = np.sum(df_this_bin.GenWeight * df_this_bin.GenxB)/weight_sum
		Q2_mean      = np.sum(df_this_bin.GenWeight * df_this_bin.GenQ2)/weight_sum
		t_mean       = np.sum(df_this_bin.GenWeight * df_this_bin.Gent)/weight_sum
		phi_mean     = np.sum(df_this_bin.GenWeight * df_this_bin.Genphi)/weight_sum
		this_row     = pd.DataFrame([{"generator": generator, "polarity": polarity, "n_entry": n_entry,
			"weight_sum": weight_sum, "weight_mean": weight_mean, "xB_mean": xB_mean,
			"Q2_mean": Q2_mean, "t_mean": t_mean, "phi_mean": phi_mean}])
		df_summary   = pd.concat([df_summary, this_row])


for bin in range(1, 147+1):
	df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_data/sim_rad_gen/dvcs_km15/pkl/fall2018_outb3/{}.pkl".format(binnum))
	df.loc[:, "integrated_binnum_gen"] = bin
	df = df.astype({"integrated_binnum_gen": int, "phi_binnum_gen": int})
	for phi_bin in range(24):
		df_this_bin  = df.loc[df.phi_binnum_gen == phi_bin]
		generator    = "km15gen"
		polarity     = "outbending"
		n_entry      = len(df_this_bin)
		weight_sum   = np.sum(df_this_bin.GenWeight)
		weight_mean  = np.mean(df_this_bin.GenWeight)
		xB_mean      = np.sum(df_this_bin.GenWeight * df_this_bin.GenxB)/weight_sum
		Q2_mean      = np.sum(df_this_bin.GenWeight * df_this_bin.GenQ2)/weight_sum
		t_mean       = np.sum(df_this_bin.GenWeight * df_this_bin.Gent)/weight_sum
		phi_mean     = np.sum(df_this_bin.GenWeight * df_this_bin.Genphi)/weight_sum
		this_row     = pd.DataFrame([{"generator": generator, "polarity": polarity, "n_entry": n_entry,
			"weight_sum": weight_sum, "weight_mean": weight_mean, "xB_mean": xB_mean,
			"Q2_mean": Q2_mean, "t_mean": t_mean, "phi_mean": phi_mean}])
		df_summary   = pd.concat([df_summary, this_row])

for bin in range(1, 159+1):
	df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_data/sim_rad_gen/dvcs_km15_fringe/pkl/fall2018_outb3/{}.pkl".format(binnum))
	df.loc[:, "integrated_binnum_gen"] = bin + 147
	df = df.astype({"integrated_binnum_gen": int, "phi_binnum_gen": int})
	for phi_bin in range(24):
		df_this_bin  = df.loc[df.phi_binnum_gen == phi_bin]
		generator    = "km15gen"
		polarity     = "outbending"
		n_entry      = len(df_this_bin)
		weight_sum   = np.sum(df_this_bin.GenWeight)
		weight_mean  = np.mean(df_this_bin.GenWeight)
		xB_mean      = np.sum(df_this_bin.GenWeight * df_this_bin.GenxB)/weight_sum
		Q2_mean      = np.sum(df_this_bin.GenWeight * df_this_bin.GenQ2)/weight_sum
		t_mean       = np.sum(df_this_bin.GenWeight * df_this_bin.Gent)/weight_sum
		phi_mean     = np.sum(df_this_bin.GenWeight * df_this_bin.Genphi)/weight_sum
		this_row     = pd.DataFrame([{"generator": generator, "polarity": polarity, "n_entry": n_entry,
			"weight_sum": weight_sum, "weight_mean": weight_mean, "xB_mean": xB_mean,
			"Q2_mean": Q2_mean, "t_mean": t_mean, "phi_mean": phi_mean}])
		df_summary   = pd.concat([df_summary, this_row])
