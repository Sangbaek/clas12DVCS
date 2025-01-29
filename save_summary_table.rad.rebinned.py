import numpy as np
import pandas as pd
from utils.const import *
from utils.physics import *

import argparse
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
warnings.simplefilter(action='ignore', category=pd.errors.SettingWithCopyWarning)


binnum, bin_volume = np.loadtxt('volume_list.csv', skiprows = 1, delimiter = ',').T # bin_volume is deprecated
bin_volume = {int(binnum[i]): bin_volume[i] for i in range(len(binnum))}
phibins = [-1] + list(np.linspace(0, 360, 24+1)[1:-1]) + [361]
df_summary_table_rebinned = pd.read_csv("df_summary_table_rebinned.csv")

def pureBH_dvcsgen_rad():
	df_summary_table      = pd.DataFrame()
	for integrated_binnum in range(1, 147+1):
		print('Reading /volatile/clas12/sangbaek/regular_backup/dvcs_data_root_only/rad_correction/pkl/pureBH/rad/pureBH_{}.pkl'.format(integrated_binnum))

		df = pd.read_pickle('/volatile/clas12/sangbaek/regular_backup/dvcs_data_root_only/rad_correction/pkl/pureBH/rad/pureBH_{}.pkl'.format(integrated_binnum))
		df.loc[:, "singularity"] = 0
		for phi_binnum in range(24):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+1]
			df_this_bin  = df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), :]
			df_this_bin.loc[df_this_bin.GenWeight.isin(df_this_bin.GenWeight.sort_values().to_numpy()[-100:]), "singularity"] = 1
			df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "singularity"] = df_this_bin.singularity
		for phi_binnum, phi_width in zip(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"], df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"]):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+phi_width]
			if len(df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), :]):
				df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "phi_binnum"] = phi_binnum

		xB_avg_this_integrated_bin = np.sum(df.GenWeight * df.xB)/np.sum(df.GenWeight)
		Q2_avg_this_integrated_bin = np.sum(df.GenWeight * df.Q2)/np.sum(df.GenWeight)
		t_avg_this_integrated_bin = np.sum(df.GenWeight * df.t1)/np.sum(df.GenWeight)

		for phi_binnum, phi_width in zip(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"], df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"]):
			df_this_bin  = df.loc[(df.phi_binnum == phi_binnum)]
			directory    = "dvcsgen/pureBH/rad"
			if len(df_this_bin) <= 100:
				n_entry    = 0
				weight_avg = 0
				xB_avg     = 0
				Q2_avg     = 0
				t_avg      = 0
				phi_avg    = 0				
			else:
				df_this_bin  = df_this_bin.loc[df_this_bin.singularity == 0, :]
				n_entry      = len(df_this_bin)
				weight_mean  = np.mean(df_this_bin.GenWeight)
				weight_mean_err = np.sqrt(np.sum(df_this_bin.GenWeight**2))/n_entry
				xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.xB)/np.sum(df_this_bin.GenWeight)
				Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.Q2)/np.sum(df_this_bin.GenWeight)
				t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.t1)/np.sum(df_this_bin.GenWeight)
				phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.phi1)/np.sum(df_this_bin.GenWeight)
			xB_avg_this_point  = xB_avg_this_integrated_bin
			Q2_avg_this_point  = Q2_avg_this_integrated_bin
			t_avg_this_point   = t_avg_this_integrated_bin
			phi_avg_this_point = phi_avg
			cross_section_this_point = printBHonly(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			tmin_this_point          = tmin(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			tcol_this_point          = tcol(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P1_this_point            = P1(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P2_this_point            = P2(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P1P2                     = P1(df_this_bin.xB, df_this_bin.Q2, df_this_bin.t1, np.radians(df_this_bin.phi1)) * P2(df_this_bin.xB, df_this_bin.Q2, df_this_bin.t1, np.radians(df_this_bin.phi1))
			maxP1P2                  = np.max(P1P2)
			minP1P2                  = np.min(P1P2)
			bin_volume_this_bin      = bin_volume[integrated_binnum]/24.*phi_width

			this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "phi_width": phi_width, "directory": directory, "n_entry": n_entry,
				"weight_mean": weight_mean, "weight_mean_err": weight_mean_err, "xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg,
				"xB_avg_this_point": xB_avg_this_point, "Q2_avg_this_point": Q2_avg_this_point, "t_avg_this_point": t_avg_this_point, "phi_avg_this_point": phi_avg_this_point, "cross_section_this_point": cross_section_this_point,
				"tmin_this_point": tmin_this_point, "tcol_this_point": tcol_this_point, "P1_this_point": P1_this_point, "P2_this_point": P2_this_point, "maxP1P2": maxP1P2, "minP1P2": minP1P2, "bin_volume_this_bin": bin_volume_this_bin}])

			df_summary_table   = pd.concat([df_summary_table, this_row])

	df_summary_table = df_summary_table.reset_index()
	df_summary_table = df_summary_table.loc[:, df_summary_table.columns[1:]]
	return df_summary_table


def pureBH_dvcsgen_norad():
	df_summary_table      = pd.DataFrame()

	for integrated_binnum in range(1, 147+1):
		print('Reading /volatile/clas12/sangbaek/regular_backup/dvcs_data_root_only/rad_correction/pkl/pureBH/norad/pureBH_{}.pkl'.format(integrated_binnum))

		df = pd.read_pickle('/volatile/clas12/sangbaek/regular_backup/dvcs_data_root_only/rad_correction/pkl/pureBH/norad/pureBH_{}.pkl'.format(integrated_binnum))
		df.loc[:, "singularity"] = 0
		for phi_binnum in range(24):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+1]
			df_this_bin  = df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), :]
			df_this_bin.loc[df_this_bin.GenWeight.isin(df_this_bin.GenWeight.sort_values().to_numpy()[-100:]), "singularity"] = 1
			df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "singularity"] = df_this_bin.singularity
		for phi_binnum, phi_width in zip(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"], df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"]):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+phi_width]
			if len(df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), :]):
				df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "phi_binnum"] = phi_binnum

		xB_avg_this_integrated_bin = np.sum(df.GenWeight * df.xB)/np.sum(df.GenWeight)
		Q2_avg_this_integrated_bin = np.sum(df.GenWeight * df.Q2)/np.sum(df.GenWeight)
		t_avg_this_integrated_bin = np.sum(df.GenWeight * df.t1)/np.sum(df.GenWeight)

		for phi_binnum, phi_width in zip(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"], df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"]):
			df_this_bin  = df.loc[(df.phi_binnum == phi_binnum)]
			directory    = "dvcsgen/pureBH/norad"
			if len(df_this_bin) <= 100:
				n_entry    = 0
				weight_avg = 0
				xB_avg     = 0
				Q2_avg     = 0
				t_avg      = 0
				phi_avg    = 0
			else:
				df_this_bin  = df_this_bin.loc[df_this_bin.singularity == 0, :]
				n_entry      = len(df_this_bin)
				weight_mean  = np.mean(df_this_bin.GenWeight)
				weight_mean_err = np.sqrt(np.sum(df_this_bin.GenWeight**2))/n_entry
				xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.xB)/np.sum(df_this_bin.GenWeight)
				Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.Q2)/np.sum(df_this_bin.GenWeight)
				t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.t1)/np.sum(df_this_bin.GenWeight)
				phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.phi1)/np.sum(df_this_bin.GenWeight)
			xB_avg_this_point  = xB_avg_this_integrated_bin
			Q2_avg_this_point  = Q2_avg_this_integrated_bin
			t_avg_this_point   = t_avg_this_integrated_bin
			phi_avg_this_point = phi_avg
			cross_section_this_point = printBHonly(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			tmin_this_point          = tmin(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			tcol_this_point          = tcol(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P1_this_point            = P1(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P2_this_point            = P2(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P1P2                     = P1(df_this_bin.xB, df_this_bin.Q2, df_this_bin.t1, np.radians(df_this_bin.phi1)) * P2(df_this_bin.xB, df_this_bin.Q2, df_this_bin.t1, np.radians(df_this_bin.phi1))
			maxP1P2                  = np.max(P1P2)
			minP1P2                  = np.min(P1P2)
			bin_volume_this_bin      = bin_volume[integrated_binnum]/24.*phi_width

			this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "phi_width": phi_width, "directory": directory, "n_entry": n_entry,
				"weight_mean": weight_mean, "weight_mean_err": weight_mean_err, "xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg,
				"xB_avg_this_point": xB_avg_this_point, "Q2_avg_this_point": Q2_avg_this_point, "t_avg_this_point": t_avg_this_point, "phi_avg_this_point": phi_avg_this_point, "cross_section_this_point": cross_section_this_point,
				"tmin_this_point": tmin_this_point, "tcol_this_point": tcol_this_point, "P1_this_point": P1_this_point, "P2_this_point": P2_this_point, "maxP1P2": maxP1P2, "minP1P2": minP1P2, "bin_volume_this_bin": bin_volume_this_bin}])

			df_summary_table   = pd.concat([df_summary_table, this_row])

	df_summary_table = df_summary_table.reset_index()
	df_summary_table = df_summary_table.loc[:, df_summary_table.columns[1:]]
	return df_summary_table

def dvcs_vgg_dvcsgen_rad():
	df_summary_table      = pd.DataFrame()

	for integrated_binnum in range(1, 147+1):
		print('Reading /volatile/clas12/sangbaek/regular_backup/dvcs_data_root_only/rad_correction/pkl/dvcs_vgg/rad/dvcs_vgg_{}.pkl'.format(integrated_binnum))

		df = pd.read_pickle('/volatile/clas12/sangbaek/regular_backup/dvcs_data_root_only/rad_correction/pkl/dvcs_vgg/rad/dvcs_vgg_{}.pkl'.format(integrated_binnum))
		df.loc[:, "singularity"] = 0
		for phi_binnum in range(24):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+1]
			df_this_bin  = df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), :]
			df_this_bin.loc[df_this_bin.GenWeight.isin(df_this_bin.GenWeight.sort_values().to_numpy()[-100:]), "singularity"] = 1
			df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "singularity"] = df_this_bin.singularity
		for phi_binnum, phi_width in zip(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"], df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"]):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+phi_width]
			if len(df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), :]):
				df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "phi_binnum"] = phi_binnum

		xB_avg_this_integrated_bin = np.sum(df.GenWeight * df.xB)/np.sum(df.GenWeight)
		Q2_avg_this_integrated_bin = np.sum(df.GenWeight * df.Q2)/np.sum(df.GenWeight)
		t_avg_this_integrated_bin = np.sum(df.GenWeight * df.t1)/np.sum(df.GenWeight)

		for phi_binnum, phi_width in zip(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"], df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"]):
			df_this_bin  = df.loc[(df.phi_binnum == phi_binnum)]
			directory    = "dvcsgen/dvcs_vgg/rad"
			if len(df_this_bin) <= 100:
				n_entry    = 0
				weight_avg = 0
				xB_avg     = 0
				Q2_avg     = 0
				t_avg      = 0
				phi_avg    = 0
			else:
				df_this_bin  = df_this_bin.loc[df_this_bin.singularity == 0, :]
				n_entry      = len(df_this_bin)
				weight_mean  = np.mean(df_this_bin.GenWeight)
				weight_mean_err = np.sqrt(np.sum(df_this_bin.GenWeight**2))/n_entry
				xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.xB)/np.sum(df_this_bin.GenWeight)
				Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.Q2)/np.sum(df_this_bin.GenWeight)
				t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.t1)/np.sum(df_this_bin.GenWeight)
				phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.phi1)/np.sum(df_this_bin.GenWeight)
			xB_avg_this_point  = xB_avg_this_integrated_bin
			Q2_avg_this_point  = Q2_avg_this_integrated_bin
			t_avg_this_point   = t_avg_this_integrated_bin
			phi_avg_this_point = phi_avg
			cross_section_this_point = printVGG(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			tmin_this_point          = tmin(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			tcol_this_point          = tcol(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P1_this_point            = P1(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P2_this_point            = P2(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P1P2                     = P1(df_this_bin.xB, df_this_bin.Q2, df_this_bin.t1, np.radians(df_this_bin.phi1)) * P2(df_this_bin.xB, df_this_bin.Q2, df_this_bin.t1, np.radians(df_this_bin.phi1))
			maxP1P2                  = np.max(P1P2)
			minP1P2                  = np.min(P1P2)
			bin_volume_this_bin      = bin_volume[integrated_binnum]/24.*phi_width

			this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "phi_width": phi_width, "directory": directory, "n_entry": n_entry,
				"weight_mean": weight_mean, "weight_mean_err": weight_mean_err, "xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg,
				"xB_avg_this_point": xB_avg_this_point, "Q2_avg_this_point": Q2_avg_this_point, "t_avg_this_point": t_avg_this_point, "phi_avg_this_point": phi_avg_this_point, "cross_section_this_point": cross_section_this_point,
				"tmin_this_point": tmin_this_point, "tcol_this_point": tcol_this_point, "P1_this_point": P1_this_point, "P2_this_point": P2_this_point, "maxP1P2": maxP1P2, "minP1P2": minP1P2, "bin_volume_this_bin": bin_volume_this_bin}])

			df_summary_table   = pd.concat([df_summary_table, this_row])

	df_summary_table = df_summary_table.reset_index()
	df_summary_table = df_summary_table.loc[:, df_summary_table.columns[1:]]
	return df_summary_table


def dvcs_vgg_dvcsgen_norad():
	df_summary_table      = pd.DataFrame()

	for integrated_binnum in range(1, 147+1):
		print('Reading /volatile/clas12/sangbaek/regular_backup/dvcs_data_root_only/rad_correction/pkl/dvcs_vgg/norad/dvcs_vgg_{}.pkl'.format(integrated_binnum))

		df = pd.read_pickle('/volatile/clas12/sangbaek/regular_backup/dvcs_data_root_only/rad_correction/pkl/dvcs_vgg/norad/dvcs_vgg_{}.pkl'.format(integrated_binnum))
		df.loc[:, "singularity"] = 0
		for phi_binnum in range(24):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+1]
			df_this_bin  = df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), :]
			df_this_bin.loc[df_this_bin.GenWeight.isin(df_this_bin.GenWeight.sort_values().to_numpy()[-100:]), "singularity"] = 1
			df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "singularity"] = df_this_bin.singularity
		for phi_binnum, phi_width in zip(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"], df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"]):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+phi_width]
			if len(df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), :]):
				df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "phi_binnum"] = phi_binnum

		xB_avg_this_integrated_bin = np.sum(df.GenWeight * df.xB)/np.sum(df.GenWeight)
		Q2_avg_this_integrated_bin = np.sum(df.GenWeight * df.Q2)/np.sum(df.GenWeight)
		t_avg_this_integrated_bin = np.sum(df.GenWeight * df.t1)/np.sum(df.GenWeight)

		for phi_binnum, phi_width in zip(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"], df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"]):
			df_this_bin  = df.loc[(df.phi_binnum == phi_binnum)]
			directory    = "dvcsgen/dvcs_vgg/norad"
			if len(df_this_bin) <= 100:
				n_entry    = 0
				weight_avg = 0
				xB_avg     = 0
				Q2_avg     = 0
				t_avg      = 0
				phi_avg    = 0
			else:
				df_this_bin  = df_this_bin.loc[df_this_bin.singularity == 0, :]
				n_entry      = len(df_this_bin)
				weight_mean  = np.mean(df_this_bin.GenWeight)
				weight_mean_err = np.sqrt(np.sum(df_this_bin.GenWeight**2))/n_entry
				xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.xB)/np.sum(df_this_bin.GenWeight)
				Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.Q2)/np.sum(df_this_bin.GenWeight)
				t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.t1)/np.sum(df_this_bin.GenWeight)
				phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.phi1)/np.sum(df_this_bin.GenWeight)
			xB_avg_this_point  = xB_avg_this_integrated_bin
			Q2_avg_this_point  = Q2_avg_this_integrated_bin
			t_avg_this_point   = t_avg_this_integrated_bin
			phi_avg_this_point = phi_avg
			cross_section_this_point = printVGG(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			tmin_this_point          = tmin(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			tcol_this_point          = tcol(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P1_this_point            = P1(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P2_this_point            = P2(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P1P2                     = P1(df_this_bin.xB, df_this_bin.Q2, df_this_bin.t1, np.radians(df_this_bin.phi1)) * P2(df_this_bin.xB, df_this_bin.Q2, df_this_bin.t1, np.radians(df_this_bin.phi1))
			maxP1P2                  = np.max(P1P2)
			minP1P2                  = np.min(P1P2)
			bin_volume_this_bin      = bin_volume[integrated_binnum]/24.*phi_width

			this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "phi_width": phi_width, "directory": directory, "n_entry": n_entry,
				"weight_mean": weight_mean, "weight_mean_err": weight_mean_err, "xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg,
				"xB_avg_this_point": xB_avg_this_point, "Q2_avg_this_point": Q2_avg_this_point, "t_avg_this_point": t_avg_this_point, "phi_avg_this_point": phi_avg_this_point, "cross_section_this_point": cross_section_this_point,
				"tmin_this_point": tmin_this_point, "tcol_this_point": tcol_this_point, "P1_this_point": P1_this_point, "P2_this_point": P2_this_point, "maxP1P2": maxP1P2, "minP1P2": minP1P2, "bin_volume_this_bin": bin_volume_this_bin}])

			df_summary_table   = pd.concat([df_summary_table, this_row])

	df_summary_table = df_summary_table.reset_index()
	df_summary_table = df_summary_table.loc[:, df_summary_table.columns[1:]]
	return df_summary_table


def dvcs_km15_rad():
	df_summary_table      = pd.DataFrame()

	for integrated_binnum in range(1, 147+1):

		df = []
		for i in range(1, 51):
			print('Reading /volatile/clas12/sangbaek/regular_backup/dvcs_data_root_only/rad_correction/pkl/dvcs_km15/rad/dvcs_km15_{0}/dvcs_km15_{0}_{1}.pkl'.format(integrated_binnum, i))
			df.append(pd.read_pickle('/volatile/clas12/sangbaek/regular_backup/dvcs_data_root_only/rad_correction/pkl/dvcs_km15/rad/dvcs_km15_{0}/dvcs_km15_{0}_{1}.pkl'.format(integrated_binnum, i)))
		df = pd.concat(df)
		df = df.reset_index()
		df = df.loc[:, df.columns[1:]]
		df.loc[:, "singularity"] = 0
		for phi_binnum in range(24):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+1]
			df_this_bin  = df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), :]
			df_this_bin.loc[df_this_bin.GenWeight.isin(df_this_bin.GenWeight.sort_values().to_numpy()[-100:]), "singularity"] = 1
			df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "singularity"] = df_this_bin.singularity
		for phi_binnum, phi_width in zip(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"], df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"]):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+phi_width]
			if len(df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), :]):
				df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "phi_binnum"] = phi_binnum

		xB_avg_this_integrated_bin = np.sum(df.GenWeight * df.xB)/np.sum(df.GenWeight)
		Q2_avg_this_integrated_bin = np.sum(df.GenWeight * df.Q2)/np.sum(df.GenWeight)
		t_avg_this_integrated_bin = np.sum(df.GenWeight * df.t1)/np.sum(df.GenWeight)

		for phi_binnum, phi_width in zip(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"], df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"]):
			df_this_bin  = df.loc[(df.phi_binnum == phi_binnum)]
			directory    = "km15gen/dvcs_km15/rad"
			if len(df_this_bin) <= 100:
				n_entry    = 0
				weight_avg = 0
				xB_avg     = 0
				Q2_avg     = 0
				t_avg      = 0
				phi_avg    = 0
			else:
				df_this_bin  = df_this_bin.loc[df_this_bin.singularity == 0, :]
				n_entry      = len(df_this_bin)
				weight_mean  = np.mean(df_this_bin.GenWeight)
				weight_mean_err = np.sqrt(np.sum(df_this_bin.GenWeight**2))/n_entry
				xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.xB)/np.sum(df_this_bin.GenWeight)
				Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.Q2)/np.sum(df_this_bin.GenWeight)
				t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.t1)/np.sum(df_this_bin.GenWeight)
				phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.phi1)/np.sum(df_this_bin.GenWeight)
			xB_avg_this_point  = xB_avg_this_integrated_bin
			Q2_avg_this_point  = Q2_avg_this_integrated_bin
			t_avg_this_point   = t_avg_this_integrated_bin
			phi_avg_this_point = phi_avg
			cross_section_this_point = printKM(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			tmin_this_point          = tmin(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			tcol_this_point          = tcol(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P1_this_point            = P1(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P2_this_point            = P2(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P1P2                     = P1(df_this_bin.xB, df_this_bin.Q2, df_this_bin.t1, np.radians(df_this_bin.phi1)) * P2(df_this_bin.xB, df_this_bin.Q2, df_this_bin.t1, np.radians(df_this_bin.phi1))
			maxP1P2                  = np.max(P1P2)
			minP1P2                  = np.min(P1P2)
			bin_volume_this_bin      = bin_volume[integrated_binnum]/24.*phi_width

			this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "phi_width": phi_width, "directory": directory, "n_entry": n_entry,
				"weight_mean": weight_mean, "weight_mean_err": weight_mean_err, "xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg,
				"xB_avg_this_point": xB_avg_this_point, "Q2_avg_this_point": Q2_avg_this_point, "t_avg_this_point": t_avg_this_point, "phi_avg_this_point": phi_avg_this_point, "cross_section_this_point": cross_section_this_point,
				"tmin_this_point": tmin_this_point, "tcol_this_point": tcol_this_point, "P1_this_point": P1_this_point, "P2_this_point": P2_this_point, "maxP1P2": maxP1P2, "minP1P2": minP1P2, "bin_volume_this_bin": bin_volume_this_bin}])

			df_summary_table   = pd.concat([df_summary_table, this_row])

	df_summary_table = df_summary_table.reset_index()
	df_summary_table = df_summary_table.loc[:, df_summary_table.columns[1:]]
	return df_summary_table


def dvcs_km15_norad():
	df_summary_table      = pd.DataFrame()

	for integrated_binnum in range(1, 147+1):

		df = []
		for i in range(1, 51):
			print('Reading /volatile/clas12/sangbaek/regular_backup/dvcs_data_root_only/rad_correction/pkl/dvcs_km15/norad/dvcs_km15_{0}/dvcs_km15_{0}_{1}.pkl'.format(integrated_binnum, i))
			df.append(pd.read_pickle('/volatile/clas12/sangbaek/regular_backup/dvcs_data_root_only/rad_correction/pkl/dvcs_km15/norad/dvcs_km15_{0}/dvcs_km15_{0}_{1}.pkl'.format(integrated_binnum, i)))
		df = pd.concat(df)
		df = df.reset_index()
		df = df.loc[:, df.columns[1:]]
		df.loc[:, "singularity"] = 0
		for phi_binnum in range(24):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+1]
			df_this_bin  = df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), :]
			df_this_bin.loc[df_this_bin.GenWeight.isin(df_this_bin.GenWeight.sort_values().to_numpy()[-100:]), "singularity"] = 1
			df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "singularity"] = df_this_bin.singularity
		for phi_binnum, phi_width in zip(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"], df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"]):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+phi_width]
			if len(df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), :]):
				df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "phi_binnum"] = phi_binnum

		xB_avg_this_integrated_bin = np.sum(df.GenWeight * df.xB)/np.sum(df.GenWeight)
		Q2_avg_this_integrated_bin = np.sum(df.GenWeight * df.Q2)/np.sum(df.GenWeight)
		t_avg_this_integrated_bin = np.sum(df.GenWeight * df.t1)/np.sum(df.GenWeight)

		for phi_binnum, phi_width in zip(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"], df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"]):
			df_this_bin  = df.loc[(df.phi_binnum == phi_binnum)]
			directory    = "km15gen/dvcs_km15/norad"
			if len(df_this_bin) <= 100:
				n_entry    = 0
				weight_avg = 0
				xB_avg     = 0
				Q2_avg     = 0
				t_avg      = 0
				phi_avg    = 0
			else:
				df_this_bin  = df_this_bin.loc[df_this_bin.singularity == 0, :]
				n_entry      = len(df_this_bin)
				weight_mean  = np.mean(df_this_bin.GenWeight)
				weight_mean_err = np.sqrt(np.sum(df_this_bin.GenWeight**2))/n_entry
				xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.xB)/np.sum(df_this_bin.GenWeight)
				Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.Q2)/np.sum(df_this_bin.GenWeight)
				t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.t1)/np.sum(df_this_bin.GenWeight)
				phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.phi1)/np.sum(df_this_bin.GenWeight)
			xB_avg_this_point  = xB_avg_this_integrated_bin
			Q2_avg_this_point  = Q2_avg_this_integrated_bin
			t_avg_this_point   = t_avg_this_integrated_bin
			phi_avg_this_point = phi_avg
			cross_section_this_point = printKM(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			tmin_this_point          = tmin(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			tcol_this_point          = tcol(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P1_this_point            = P1(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P2_this_point            = P2(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P1P2                     = P1(df_this_bin.xB, df_this_bin.Q2, df_this_bin.t1, np.radians(df_this_bin.phi1)) * P2(df_this_bin.xB, df_this_bin.Q2, df_this_bin.t1, np.radians(df_this_bin.phi1))
			maxP1P2                  = np.max(P1P2)
			minP1P2                  = np.min(P1P2)
			bin_volume_this_bin      = bin_volume[integrated_binnum]/24.*phi_width

			this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "phi_width": phi_width, "directory": directory, "n_entry": n_entry,
				"weight_mean": weight_mean, "weight_mean_err": weight_mean_err, "xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg,
				"xB_avg_this_point": xB_avg_this_point, "Q2_avg_this_point": Q2_avg_this_point, "t_avg_this_point": t_avg_this_point, "phi_avg_this_point": phi_avg_this_point, "cross_section_this_point": cross_section_this_point,
				"tmin_this_point": tmin_this_point, "tcol_this_point": tcol_this_point, "P1_this_point": P1_this_point, "P2_this_point": P2_this_point, "maxP1P2": maxP1P2, "minP1P2": minP1P2, "bin_volume_this_bin": bin_volume_this_bin}])

			df_summary_table   = pd.concat([df_summary_table, this_row])

	df_summary_table = df_summary_table.reset_index()
	df_summary_table = df_summary_table.loc[:, df_summary_table.columns[1:]]
	return df_summary_table

def pureBH_km15_rad():
	df_summary_table      = pd.DataFrame()

	for integrated_binnum in range(1, 147+1):

		df = []
		for i in range(1, 51):
			print('Reading /volatile/clas12/sangbaek/regular_backup/dvcs_data_root_only/rad_correction/pkl/pureBH_km15/rad/pureBH_km15_{0}/pureBH_km15_{0}_{1}.pkl'.format(integrated_binnum, i))
			df.append(pd.read_pickle('/volatile/clas12/sangbaek/regular_backup/dvcs_data_root_only/rad_correction/pkl/pureBH_km15/rad/pureBH_km15_{0}/pureBH_km15_{0}_{1}.pkl'.format(integrated_binnum, i)))
		df = pd.concat(df)
		df = df.reset_index()
		df = df.loc[:, df.columns[1:]]
		df.loc[:, "singularity"] = 0
		for phi_binnum in range(24):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+1]
			df_this_bin  = df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), :]
			df_this_bin.loc[df_this_bin.GenWeight.isin(df_this_bin.GenWeight.sort_values().to_numpy()[-100:]), "singularity"] = 1
			df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "singularity"] = df_this_bin.singularity
		for phi_binnum, phi_width in zip(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"], df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"]):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+phi_width]
			if len(df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), :]):
				df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "phi_binnum"] = phi_binnum

		xB_avg_this_integrated_bin = np.sum(df.GenWeight * df.xB)/np.sum(df.GenWeight)
		Q2_avg_this_integrated_bin = np.sum(df.GenWeight * df.Q2)/np.sum(df.GenWeight)
		t_avg_this_integrated_bin = np.sum(df.GenWeight * df.t1)/np.sum(df.GenWeight)

		for phi_binnum, phi_width in zip(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"], df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"]):
			df_this_bin  = df.loc[(df.phi_binnum == phi_binnum)]
			directory    = "km15gen/pureBH_km15/rad"
			if len(df_this_bin) <= 100:
				n_entry    = 0
				weight_avg = 0
				xB_avg     = 0
				Q2_avg     = 0
				t_avg      = 0
				phi_avg    = 0
			else:
				df_this_bin  = df_this_bin.loc[df_this_bin.singularity == 0, :]
				n_entry      = len(df_this_bin)
				weight_mean  = np.mean(df_this_bin.GenWeight)
				weight_mean_err = np.sqrt(np.sum(df_this_bin.GenWeight**2))/n_entry
				xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.xB)/np.sum(df_this_bin.GenWeight)
				Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.Q2)/np.sum(df_this_bin.GenWeight)
				t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.t1)/np.sum(df_this_bin.GenWeight)
				phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.phi1)/np.sum(df_this_bin.GenWeight)
			xB_avg_this_point  = xB_avg_this_integrated_bin
			Q2_avg_this_point  = Q2_avg_this_integrated_bin
			t_avg_this_point   = t_avg_this_integrated_bin
			phi_avg_this_point = phi_avg
			cross_section_this_point = printKM(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point), mode = 1)
			tmin_this_point          = tmin(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			tcol_this_point          = tcol(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P1_this_point            = P1(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P2_this_point            = P2(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P1P2                     = P1(df_this_bin.xB, df_this_bin.Q2, df_this_bin.t1, np.radians(df_this_bin.phi1)) * P2(df_this_bin.xB, df_this_bin.Q2, df_this_bin.t1, np.radians(df_this_bin.phi1))
			maxP1P2                  = np.max(P1P2)
			minP1P2                  = np.min(P1P2)
			bin_volume_this_bin      = bin_volume[integrated_binnum]/24.*phi_width

			this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "phi_width": phi_width, "directory": directory, "n_entry": n_entry,
				"weight_mean": weight_mean, "weight_mean_err": weight_mean_err, "xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg,
				"xB_avg_this_point": xB_avg_this_point, "Q2_avg_this_point": Q2_avg_this_point, "t_avg_this_point": t_avg_this_point, "phi_avg_this_point": phi_avg_this_point, "cross_section_this_point": cross_section_this_point,
				"tmin_this_point": tmin_this_point, "tcol_this_point": tcol_this_point, "P1_this_point": P1_this_point, "P2_this_point": P2_this_point, "maxP1P2": maxP1P2, "minP1P2": minP1P2, "bin_volume_this_bin": bin_volume_this_bin}])

			df_summary_table   = pd.concat([df_summary_table, this_row])

	df_summary_table = df_summary_table.reset_index()
	df_summary_table = df_summary_table.loc[:, df_summary_table.columns[1:]]
	return df_summary_table


def pureBH_km15_norad():
	df_summary_table      = pd.DataFrame()

	for integrated_binnum in range(1, 147+1):

		df = []
		for i in range(1, 51):
			print('Reading /volatile/clas12/sangbaek/regular_backup/dvcs_data_root_only/rad_correction/pkl/pureBH_km15/norad/pureBH_km15_{0}/pureBH_km15_{0}_{1}.pkl'.format(integrated_binnum, i))
			df.append(pd.read_pickle('/volatile/clas12/sangbaek/regular_backup/dvcs_data_root_only/rad_correction/pkl/pureBH_km15/norad/pureBH_km15_{0}/pureBH_km15_{0}_{1}.pkl'.format(integrated_binnum, i)))
		df = pd.concat(df)
		df = df.reset_index()
		df = df.loc[:, df.columns[1:]]
		df.loc[:, "singularity"] = 0
		for phi_binnum in range(24):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+1]
			df_this_bin  = df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), :]
			df_this_bin.loc[df_this_bin.GenWeight.isin(df_this_bin.GenWeight.sort_values().to_numpy()[-100:]), "singularity"] = 1
			df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "singularity"] = df_this_bin.singularity

		for phi_binnum, phi_width in zip(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"], df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"]):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+phi_width]
			if len(df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), :]):
				df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "phi_binnum"] = phi_binnum

		xB_avg_this_integrated_bin = np.sum(df.GenWeight * df.xB)/np.sum(df.GenWeight)
		Q2_avg_this_integrated_bin = np.sum(df.GenWeight * df.Q2)/np.sum(df.GenWeight)
		t_avg_this_integrated_bin = np.sum(df.GenWeight * df.t1)/np.sum(df.GenWeight)

		for phi_binnum, phi_width in zip(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"], df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"]):
			df_this_bin  = df.loc[(df.phi_binnum == phi_binnum)]
			directory    = "km15gen/pureBH_km15/norad"
			if len(df_this_bin) <= 100:
				n_entry    = 0
				weight_avg = 0
				xB_avg     = 0
				Q2_avg     = 0
				t_avg      = 0
				phi_avg    = 0
			else:
				df_this_bin  = df_this_bin.loc[df_this_bin.singularity == 0, :]
				n_entry      = len(df_this_bin)
				weight_mean  = np.mean(df_this_bin.GenWeight)
				weight_mean_err = np.sqrt(np.sum(df_this_bin.GenWeight**2))/n_entry
				xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.xB)/np.sum(df_this_bin.GenWeight)
				Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.Q2)/np.sum(df_this_bin.GenWeight)
				t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.t1)/np.sum(df_this_bin.GenWeight)
				phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.phi1)/np.sum(df_this_bin.GenWeight)
			xB_avg_this_point  = xB_avg_this_integrated_bin
			Q2_avg_this_point  = Q2_avg_this_integrated_bin
			t_avg_this_point   = t_avg_this_integrated_bin
			phi_avg_this_point = phi_avg
			cross_section_this_point = printKM(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point), mode = 1)
			tmin_this_point          = tmin(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			tcol_this_point          = tcol(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P1_this_point            = P1(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P2_this_point            = P2(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))
			P1P2                     = P1(df_this_bin.xB, df_this_bin.Q2, df_this_bin.t1, np.radians(df_this_bin.phi1)) * P2(df_this_bin.xB, df_this_bin.Q2, df_this_bin.t1, np.radians(df_this_bin.phi1))
			maxP1P2                  = np.max(P1P2)
			minP1P2                  = np.min(P1P2)
			bin_volume_this_bin      = bin_volume[integrated_binnum]/24.*phi_width

			this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "phi_width": phi_width, "directory": directory, "n_entry": n_entry,
				"weight_mean": weight_mean, "weight_mean_err": weight_mean_err, "xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg,
				"xB_avg_this_point": xB_avg_this_point, "Q2_avg_this_point": Q2_avg_this_point, "t_avg_this_point": t_avg_this_point, "phi_avg_this_point": phi_avg_this_point, "cross_section_this_point": cross_section_this_point,
				"tmin_this_point": tmin_this_point, "tcol_this_point": tcol_this_point, "P1_this_point": P1_this_point, "P2_this_point": P2_this_point, "maxP1P2": maxP1P2, "minP1P2": minP1P2, "bin_volume_this_bin": bin_volume_this_bin}])

			df_summary_table   = pd.concat([df_summary_table, this_row])

	df_summary_table = df_summary_table.reset_index()
	df_summary_table = df_summary_table.loc[:, df_summary_table.columns[1:]]
	return df_summary_table


if __name__ == "__main__":
	dfs = []
	result =  pureBH_dvcsgen_rad()
	dfs.append(result)
	df = pd.concat(dfs)
	df = df.reset_index()
	df = df.loc[:, df.columns[1:]]
	df.to_pickle("summary_table.rad.suppressed.rebinned.pkl")
	result =  pureBH_dvcsgen_norad()
	dfs.append(result)
	df = pd.concat(dfs)
	df = df.reset_index()
	df = df.loc[:, df.columns[1:]]
	df.to_pickle("summary_table.rad.suppressed.rebinned.pkl")
	result =  dvcs_vgg_dvcsgen_rad()
	dfs.append(result)
	df = pd.concat(dfs)
	df = df.reset_index()
	df = df.loc[:, df.columns[1:]]
	df.to_pickle("summary_table.rad.suppressed.rebinned.pkl")
	result =  dvcs_vgg_dvcsgen_norad()
	dfs.append(result)
	df = pd.concat(dfs)
	df = df.reset_index()
	df = df.loc[:, df.columns[1:]]
	df.to_pickle("summary_table.rad.suppressed.rebinned.pkl")
	# # df_existing = pd.read_pickle("summary_table.rad.pkl")
	# # dfs.append(df_existing)
	result =  dvcs_km15_rad()
	dfs.append(result)
	df = pd.concat(dfs)
	df = df.reset_index()
	df = df.loc[:, df.columns[1:]]
	df.to_pickle("summary_table.rad.suppressed.rebinned.pkl")
	result =  dvcs_km15_norad()
	dfs.append(result)
	df = pd.concat(dfs)
	df = df.reset_index()
	df = df.loc[:, df.columns[1:]]
	df.to_pickle("summary_table.rad.suppressed.rebinned.pkl")
	result =  pureBH_km15_rad()
	dfs.append(result)
	df = pd.concat(dfs)
	df = df.reset_index()
	df = df.loc[:, df.columns[1:]]
	df.to_pickle("summary_table.rad.suppressed.rebinned.pkl")
	result =  pureBH_km15_norad()
	dfs.append(result)
	df = pd.concat(dfs)
	df = df.reset_index()
	df = df.loc[:, df.columns[1:]]
	df.to_pickle("summary_table.rad.suppressed.rebinned.pkl")
