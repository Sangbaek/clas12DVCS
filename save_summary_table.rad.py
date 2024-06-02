import numpy as np
import pandas as pd
from utils.const import *
from utils.physics import *

import argparse
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
warnings.simplefilter(action='ignore', category=pd.errors.SettingWithCopyWarning)


def pureBH_dvcsgen_rad():
	df_summary = pd.DataFrame()

	for integrated_binnum in range(1, 147+1):
		print('Reading /volatile/clas12/sangbaek/dvcs_related/rad_correction/pkl/pureBH/rad/pureBH_{}.pkl'.format(integrated_binnum))

		df = pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/rad_correction/pkl/pureBH/rad/pureBH_{}.pkl'.format(integrated_binnum))
		phibins = [-1] + list(np.linspace(0, 360, 24+1)[1:-1]) + [361]
		for phi_binnum in range(24):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+1]
			df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "phi_binnum"] = phi_binnum

		xB_avg_this_integrated_bin = np.sum(df.GenWeight * df.xB)/np.sum(df.GenWeight)
		Q2_avg_this_integrated_bin = np.sum(df.GenWeight * df.Q2)/np.sum(df.GenWeight)
		t_avg_this_integrated_bin = np.sum(df.GenWeight * df.t1)/np.sum(df.GenWeight)

		for phi_binnum in range(24):
			df_this_bin  = df.loc[(df.phi_binnum == phi_binnum)]
			directory    = "dvcsgen/pureBH/rad"
			if not len(df_this_bin):
				n_entry    = 0
				weight_avg = 0
				xB_avg     = 0
				Q2_avg     = 0
				t_avg      = 0
				phi_avg    = 0				
			else:
				n_entry      = len(df_this_bin)
				weight_mean  = np.mean(df_this_bin.GenWeight)
				weight_mean_err = np.sqrt(np.sum(df_this_bin.GenWeight**2))
				xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.xB)/np.sum(df_this_bin.GenWeight)
				Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.Q2)/np.sum(df_this_bin.GenWeight)
				t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.t1)/np.sum(df_this_bin.GenWeight)
				phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.phi1)/np.sum(df_this_bin.GenWeight)
			xB_avg_this_point  = xB_avg_this_integrated_bin
			Q2_avg_this_point  = Q2_avg_this_integrated_bin
			t_avg_this_point   = t_avg_this_integrated_bin
			phi_avg_this_point = phi_avg
			cross_section_this_point = printBHonly(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))

			this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "n_entry": n_entry,
				"weight_mean": weight_mean, "weight_mean_err": weight_mean_err, "xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg,
				"xB_avg_this_point": xB_avg_this_point, "Q2_avg_this_point": Q2_avg_this_point, "t_avg_this_point": t_avg_this_point, "phi_avg_this_point": phi_avg_this_point, "cross_section_this_point": cross_section_this_point}])

			df_summary   = pd.concat([df_summary, this_row])


	df_summary = df_summary.reset_index()
	df_summary = df_summary.loc[:, df_summary.columns[1:]]
	return df_summary


def pureBH_dvcsgen_norad():
	df_summary = pd.DataFrame()

	for integrated_binnum in range(1, 147+1):
		print('Reading /volatile/clas12/sangbaek/dvcs_related/rad_correction/pkl/pureBH/norad/pureBH_{}.pkl'.format(integrated_binnum))

		df = pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/rad_correction/pkl/pureBH/norad/pureBH_{}.pkl'.format(integrated_binnum))
		phibins = [-1] + list(np.linspace(0, 360, 24+1)[1:-1]) + [361]
		for phi_binnum in range(24):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+1]
			df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "phi_binnum"] = phi_binnum

		xB_avg_this_integrated_bin = np.sum(df.GenWeight * df.xB)/np.sum(df.GenWeight)
		Q2_avg_this_integrated_bin = np.sum(df.GenWeight * df.Q2)/np.sum(df.GenWeight)
		t_avg_this_integrated_bin = np.sum(df.GenWeight * df.t1)/np.sum(df.GenWeight)

		for phi_binnum in range(24):
			df_this_bin  = df.loc[(df.phi_binnum == phi_binnum)]
			directory    = "dvcsgen/pureBH/norad"
			if not len(df_this_bin):
				n_entry    = 0
				weight_avg = 0
				xB_avg     = 0
				Q2_avg     = 0
				t_avg      = 0
				phi_avg    = 0
			else:
				n_entry      = len(df_this_bin)
				weight_mean  = np.mean(df_this_bin.GenWeight)
				weight_mean_err = np.sqrt(np.sum(df_this_bin.GenWeight**2))
				xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.xB)/np.sum(df_this_bin.GenWeight)
				Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.Q2)/np.sum(df_this_bin.GenWeight)
				t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.t1)/np.sum(df_this_bin.GenWeight)
				phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.phi1)/np.sum(df_this_bin.GenWeight)
			xB_avg_this_point  = xB_avg_this_integrated_bin
			Q2_avg_this_point  = Q2_avg_this_integrated_bin
			t_avg_this_point   = t_avg_this_integrated_bin
			phi_avg_this_point = phi_avg
			cross_section_this_point = printBHonly(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))

			this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "n_entry": n_entry,
				"weight_mean": weight_mean, "weight_mean_err": weight_mean_err, "xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg,
				"xB_avg_this_point": xB_avg_this_point, "Q2_avg_this_point": Q2_avg_this_point, "t_avg_this_point": t_avg_this_point, "phi_avg_this_point": phi_avg_this_point, "cross_section_this_point": cross_section_this_point}])

			df_summary   = pd.concat([df_summary, this_row])


	df_summary = df_summary.reset_index()
	df_summary = df_summary.loc[:, df_summary.columns[1:]]
	return df_summary

def dvcs_vgg_dvcsgen_rad():
	df_summary = pd.DataFrame()

	for integrated_binnum in range(1, 147+1):
		print('Reading /volatile/clas12/sangbaek/dvcs_related/rad_correction/pkl/dvcs_vgg/rad/dvcs_vgg_{}.pkl'.format(integrated_binnum))

		df = pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/rad_correction/pkl/dvcs_vgg/rad/dvcs_vgg_{}.pkl'.format(integrated_binnum))
		phibins = [-1] + list(np.linspace(0, 360, 24+1)[1:-1]) + [361]
		for phi_binnum in range(24):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+1]
			df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "phi_binnum"] = phi_binnum

		xB_avg_this_integrated_bin = np.sum(df.GenWeight * df.xB)/np.sum(df.GenWeight)
		Q2_avg_this_integrated_bin = np.sum(df.GenWeight * df.Q2)/np.sum(df.GenWeight)
		t_avg_this_integrated_bin = np.sum(df.GenWeight * df.t1)/np.sum(df.GenWeight)

		for phi_binnum in range(24):
			df_this_bin  = df.loc[(df.phi_binnum == phi_binnum)]
			directory    = "dvcsgen/dvcs_vgg/rad"
			if not len(df_this_bin):
				n_entry    = 0
				weight_avg = 0
				xB_avg     = 0
				Q2_avg     = 0
				t_avg      = 0
				phi_avg    = 0
			else:
				n_entry      = len(df_this_bin)
				weight_mean  = np.mean(df_this_bin.GenWeight)
				weight_mean_err = np.sqrt(np.sum(df_this_bin.GenWeight**2))
				xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.xB)/np.sum(df_this_bin.GenWeight)
				Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.Q2)/np.sum(df_this_bin.GenWeight)
				t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.t1)/np.sum(df_this_bin.GenWeight)
				phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.phi1)/np.sum(df_this_bin.GenWeight)
			xB_avg_this_point  = xB_avg_this_integrated_bin
			Q2_avg_this_point  = Q2_avg_this_integrated_bin
			t_avg_this_point   = t_avg_this_integrated_bin
			phi_avg_this_point = phi_avg
			cross_section_this_point = printVGG(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))

			this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "n_entry": n_entry,
				"weight_mean": weight_mean, "weight_mean_err": weight_mean_err, "xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg,
				"xB_avg_this_point": xB_avg_this_point, "Q2_avg_this_point": Q2_avg_this_point, "t_avg_this_point": t_avg_this_point, "phi_avg_this_point": phi_avg_this_point, "cross_section_this_point": cross_section_this_point}])

			df_summary   = pd.concat([df_summary, this_row])


	df_summary = df_summary.reset_index()
	df_summary = df_summary.loc[:, df_summary.columns[1:]]
	return df_summary


def dvcs_vgg_dvcsgen_norad():
	df_summary = pd.DataFrame()

	for integrated_binnum in range(1, 147+1):
		print('Reading /volatile/clas12/sangbaek/dvcs_related/rad_correction/pkl/dvcs_vgg/norad/dvcs_vgg_{}.pkl'.format(integrated_binnum))

		df = pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/rad_correction/pkl/dvcs_vgg/norad/dvcs_vgg_{}.pkl'.format(integrated_binnum))
		phibins = [-1] + list(np.linspace(0, 360, 24+1)[1:-1]) + [361]
		for phi_binnum in range(24):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+1]
			df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "phi_binnum"] = phi_binnum

		xB_avg_this_integrated_bin = np.sum(df.GenWeight * df.xB)/np.sum(df.GenWeight)
		Q2_avg_this_integrated_bin = np.sum(df.GenWeight * df.Q2)/np.sum(df.GenWeight)
		t_avg_this_integrated_bin = np.sum(df.GenWeight * df.t1)/np.sum(df.GenWeight)

		for phi_binnum in range(24):
			df_this_bin  = df.loc[(df.phi_binnum == phi_binnum)]
			directory    = "dvcsgen/dvcs_vgg/norad"
			if not len(df_this_bin):
				n_entry    = 0
				weight_avg = 0
				xB_avg     = 0
				Q2_avg     = 0
				t_avg      = 0
				phi_avg    = 0
			else:
				n_entry      = len(df_this_bin)
				weight_mean  = np.mean(df_this_bin.GenWeight)
				weight_mean_err = np.sqrt(np.sum(df_this_bin.GenWeight**2))
				xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.xB)/np.sum(df_this_bin.GenWeight)
				Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.Q2)/np.sum(df_this_bin.GenWeight)
				t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.t1)/np.sum(df_this_bin.GenWeight)
				phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.phi1)/np.sum(df_this_bin.GenWeight)
			xB_avg_this_point  = xB_avg_this_integrated_bin
			Q2_avg_this_point  = Q2_avg_this_integrated_bin
			t_avg_this_point   = t_avg_this_integrated_bin
			phi_avg_this_point = phi_avg
			cross_section_this_point = printVGG(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))

			this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "n_entry": n_entry,
				"weight_mean": weight_mean, "weight_mean_err": weight_mean_err, "xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg,
				"xB_avg_this_point": xB_avg_this_point, "Q2_avg_this_point": Q2_avg_this_point, "t_avg_this_point": t_avg_this_point, "phi_avg_this_point": phi_avg_this_point, "cross_section_this_point": cross_section_this_point}])

			df_summary   = pd.concat([df_summary, this_row])


	df_summary = df_summary.reset_index()
	df_summary = df_summary.loc[:, df_summary.columns[1:]]
	return df_summary


def dvcs_km15_rad():
	df_summary = pd.DataFrame()

	for integrated_binnum in range(1, 147+1):

		df = []
		for i in range(1, 51):
			print('Reading /volatile/clas12/sangbaek/dvcs_related/rad_correction/pkl/dvcs_km15/rad/dvcs_km15_{0}/dvcs_km15_{0}_{1}.pkl'.format(integrated_binnum, i))
			df.append(pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/rad_correction/pkl/dvcs_km15/rad/dvcs_km15_{0}/dvcs_km15_{0}_{1}.pkl'.format(integrated_binnum, i)))
		df = pd.concat(df)
		df = df.reset_index()
		df = df.loc[:, df.columns[1:]]
		phibins = [-1] + list(np.linspace(0, 360, 24+1)[1:-1]) + [361]
		for phi_binnum in range(24):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+1]
			df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "phi_binnum"] = phi_binnum

		xB_avg_this_integrated_bin = np.sum(df.GenWeight * df.xB)/np.sum(df.GenWeight)
		Q2_avg_this_integrated_bin = np.sum(df.GenWeight * df.Q2)/np.sum(df.GenWeight)
		t_avg_this_integrated_bin = np.sum(df.GenWeight * df.t1)/np.sum(df.GenWeight)

		for phi_binnum in range(24):
			df_this_bin  = df.loc[(df.phi_binnum == phi_binnum)]
			directory    = "km15gen/dvcs_km15/rad"
			if not len(df_this_bin):
				n_entry    = 0
				weight_avg = 0
				xB_avg     = 0
				Q2_avg     = 0
				t_avg      = 0
				phi_avg    = 0
			else:
				n_entry      = len(df_this_bin)
				weight_mean  = np.mean(df_this_bin.GenWeight)
				weight_mean_err = np.sqrt(np.sum(df_this_bin.GenWeight**2))
				xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.xB)/np.sum(df_this_bin.GenWeight)
				Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.Q2)/np.sum(df_this_bin.GenWeight)
				t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.t1)/np.sum(df_this_bin.GenWeight)
				phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.phi1)/np.sum(df_this_bin.GenWeight)
			xB_avg_this_point  = xB_avg_this_integrated_bin
			Q2_avg_this_point  = Q2_avg_this_integrated_bin
			t_avg_this_point   = t_avg_this_integrated_bin
			phi_avg_this_point = phi_avg
			cross_section_this_point = printKM(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))

			this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "n_entry": n_entry,
				"weight_mean": weight_mean, "weight_mean_err": weight_mean_err, "xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg,
				"xB_avg_this_point": xB_avg_this_point, "Q2_avg_this_point": Q2_avg_this_point, "t_avg_this_point": t_avg_this_point, "phi_avg_this_point": phi_avg_this_point, "cross_section_this_point": cross_section_this_point}])

			df_summary   = pd.concat([df_summary, this_row])


	df_summary = df_summary.reset_index()
	df_summary = df_summary.loc[:, df_summary.columns[1:]]
	return df_summary


def dvcs_km15_norad():
	df_summary = pd.DataFrame()

	for integrated_binnum in range(1, 147+1):

		df = []
		for i in range(1, 51):
			print('Reading /volatile/clas12/sangbaek/dvcs_related/rad_correction/pkl/dvcs_km15/norad/dvcs_km15_{0}/dvcs_km15_{0}_{1}.pkl'.format(integrated_binnum, i))
			df.append(pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/rad_correction/pkl/dvcs_km15/norad/dvcs_km15_{0}/dvcs_km15_{0}_{1}.pkl'.format(integrated_binnum, i)))
		df = pd.concat(df)
		df = df.reset_index()
		df = df.loc[:, df.columns[1:]]
		phibins = [-1] + list(np.linspace(0, 360, 24+1)[1:-1]) + [361]
		for phi_binnum in range(24):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+1]
			df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "phi_binnum"] = phi_binnum

		xB_avg_this_integrated_bin = np.sum(df.GenWeight * df.xB)/np.sum(df.GenWeight)
		Q2_avg_this_integrated_bin = np.sum(df.GenWeight * df.Q2)/np.sum(df.GenWeight)
		t_avg_this_integrated_bin = np.sum(df.GenWeight * df.t1)/np.sum(df.GenWeight)

		for phi_binnum in range(24):
			df_this_bin  = df.loc[(df.phi_binnum == phi_binnum)]
			directory    = "km15gen/dvcs_km15/norad"
			if not len(df_this_bin):
				n_entry    = 0
				weight_avg = 0
				xB_avg     = 0
				Q2_avg     = 0
				t_avg      = 0
				phi_avg    = 0
			else:
				n_entry      = len(df_this_bin)
				weight_mean  = np.mean(df_this_bin.GenWeight)
				weight_mean_err = np.sqrt(np.sum(df_this_bin.GenWeight**2))
				xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.xB)/np.sum(df_this_bin.GenWeight)
				Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.Q2)/np.sum(df_this_bin.GenWeight)
				t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.t1)/np.sum(df_this_bin.GenWeight)
				phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.phi1)/np.sum(df_this_bin.GenWeight)
			xB_avg_this_point  = xB_avg_this_integrated_bin
			Q2_avg_this_point  = Q2_avg_this_integrated_bin
			t_avg_this_point   = t_avg_this_integrated_bin
			phi_avg_this_point = phi_avg
			cross_section_this_point = printKM(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point))

			this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "n_entry": n_entry,
				"weight_mean": weight_mean, "weight_mean_err": weight_mean_err, "xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg,
				"xB_avg_this_point": xB_avg_this_point, "Q2_avg_this_point": Q2_avg_this_point, "t_avg_this_point": t_avg_this_point, "phi_avg_this_point": phi_avg_this_point, "cross_section_this_point": cross_section_this_point}])

			df_summary   = pd.concat([df_summary, this_row])


	df_summary = df_summary.reset_index()
	df_summary = df_summary.loc[:, df_summary.columns[1:]]
	return df_summary

def pureBH_km15_rad():
	df_summary = pd.DataFrame()

	for integrated_binnum in range(1, 147+1):

		df = []
		for i in range(1, 51):
			print('Reading /volatile/clas12/sangbaek/dvcs_related/rad_correction/pkl/pureBH_km15/rad/pureBH_km15_{0}/pureBH_km15_{0}_{1}.pkl'.format(integrated_binnum, i))
			df.append(pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/rad_correction/pkl/pureBH_km15/rad/pureBH_km15_{0}/pureBH_km15_{0}_{1}.pkl'.format(integrated_binnum, i)))
		df = pd.concat(df)
		df = df.reset_index()
		df = df.loc[:, df.columns[1:]]
		phibins = [-1] + list(np.linspace(0, 360, 24+1)[1:-1]) + [361]
		for phi_binnum in range(24):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+1]
			df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "phi_binnum"] = phi_binnum

		xB_avg_this_integrated_bin = np.sum(df.GenWeight * df.xB)/np.sum(df.GenWeight)
		Q2_avg_this_integrated_bin = np.sum(df.GenWeight * df.Q2)/np.sum(df.GenWeight)
		t_avg_this_integrated_bin = np.sum(df.GenWeight * df.t1)/np.sum(df.GenWeight)

		for phi_binnum in range(24):
			df_this_bin  = df.loc[(df.phi_binnum == phi_binnum)]
			directory    = "km15gen/pureBH_km15/rad"
			if not len(df_this_bin):
				n_entry    = 0
				weight_avg = 0
				xB_avg     = 0
				Q2_avg     = 0
				t_avg      = 0
				phi_avg    = 0
			else:
				n_entry      = len(df_this_bin)
				weight_mean  = np.mean(df_this_bin.GenWeight)
				weight_mean_err = np.sqrt(np.sum(df_this_bin.GenWeight**2))
				xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.xB)/np.sum(df_this_bin.GenWeight)
				Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.Q2)/np.sum(df_this_bin.GenWeight)
				t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.t1)/np.sum(df_this_bin.GenWeight)
				phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.phi1)/np.sum(df_this_bin.GenWeight)
			xB_avg_this_point  = xB_avg_this_integrated_bin
			Q2_avg_this_point  = Q2_avg_this_integrated_bin
			t_avg_this_point   = t_avg_this_integrated_bin
			phi_avg_this_point = phi_avg
			cross_section_this_point = printKM(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point), mode = 1)

			this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "n_entry": n_entry,
				"weight_mean": weight_mean, "weight_mean_err": weight_mean_err, "xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg,
				"xB_avg_this_point": xB_avg_this_point, "Q2_avg_this_point": Q2_avg_this_point, "t_avg_this_point": t_avg_this_point, "phi_avg_this_point": phi_avg_this_point, "cross_section_this_point": cross_section_this_point}])

			df_summary   = pd.concat([df_summary, this_row])


	df_summary = df_summary.reset_index()
	df_summary = df_summary.loc[:, df_summary.columns[1:]]
	return df_summary


def pureBH_km15_norad():
	df_summary = pd.DataFrame()

	for integrated_binnum in range(1, 147+1):

		df = []
		for i in range(1, 51):
			print('Reading /volatile/clas12/sangbaek/dvcs_related/rad_correction/pkl/pureBH_km15/norad/pureBH_km15_{0}/pureBH_km15_{0}_{1}.pkl'.format(integrated_binnum, i))
			df.append(pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/rad_correction/pkl/pureBH_km15/norad/pureBH_km15_{0}/pureBH_km15_{0}_{1}.pkl'.format(integrated_binnum, i)))
		df = pd.concat(df)
		df = df.reset_index()
		df = df.loc[:, df.columns[1:]]
		phibins = [-1] + list(np.linspace(0, 360, 24+1)[1:-1]) + [361]
		for phi_binnum in range(24):
			phimin = phibins[phi_binnum]
			phimax = phibins[phi_binnum+1]
			df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), "phi_binnum"] = phi_binnum

		xB_avg_this_integrated_bin = np.sum(df.GenWeight * df.xB)/np.sum(df.GenWeight)
		Q2_avg_this_integrated_bin = np.sum(df.GenWeight * df.Q2)/np.sum(df.GenWeight)
		t_avg_this_integrated_bin = np.sum(df.GenWeight * df.t1)/np.sum(df.GenWeight)

		for phi_binnum in range(24):
			df_this_bin  = df.loc[(df.phi_binnum == phi_binnum)]
			directory    = "km15gen/pureBH_km15/norad"
			if not len(df_this_bin):
				n_entry    = 0
				weight_avg = 0
				xB_avg     = 0
				Q2_avg     = 0
				t_avg      = 0
				phi_avg    = 0
			else:
				n_entry      = len(df_this_bin)
				weight_mean  = np.mean(df_this_bin.GenWeight)
				weight_mean_err = np.sqrt(np.sum(df_this_bin.GenWeight**2))
				xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.xB)/np.sum(df_this_bin.GenWeight)
				Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.Q2)/np.sum(df_this_bin.GenWeight)
				t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.t1)/np.sum(df_this_bin.GenWeight)
				phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.phi1)/np.sum(df_this_bin.GenWeight)
			xB_avg_this_point  = xB_avg_this_integrated_bin
			Q2_avg_this_point  = Q2_avg_this_integrated_bin
			t_avg_this_point   = t_avg_this_integrated_bin
			phi_avg_this_point = phi_avg
			cross_section_this_point = printKM(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point), mode = 1)

			this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "n_entry": n_entry,
				"weight_mean": weight_mean, "weight_mean_err": weight_mean_err, "xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg,
				"xB_avg_this_point": xB_avg_this_point, "Q2_avg_this_point": Q2_avg_this_point, "t_avg_this_point": t_avg_this_point, "phi_avg_this_point": phi_avg_this_point, "cross_section_this_point": cross_section_this_point}])

			df_summary   = pd.concat([df_summary, this_row])


	df_summary = df_summary.reset_index()
	df_summary = df_summary.loc[:, df_summary.columns[1:]]
	return df_summary


if __name__ == "__main__":
	dfs = []
	# result =  pureBH_dvcsgen_rad()
	# dfs.append(result)
	# df = pd.concat(dfs)
	# df = df.reset_index()
	# df = df.loc[:, df.columns[1:]]
	# df.to_pickle("summary_table.rad.pkl")
	# result =  pureBH_dvcsgen_norad()
	# dfs.append(result)
	# df = pd.concat(dfs)
	# df = df.reset_index()
	# df = df.loc[:, df.columns[1:]]
	# df.to_pickle("summary_table.rad.pkl")
	# result =  dvcs_vgg_dvcsgen_rad()
	# dfs.append(result)
	# df = pd.concat(dfs)
	# df = df.reset_index()
	# df = df.loc[:, df.columns[1:]]
	# df.to_pickle("summary_table.rad.pkl")
	# result =  dvcs_vgg_dvcsgen_norad()
	# dfs.append(result)
	# df = pd.concat(dfs)
	# df = df.reset_index()
	# df = df.loc[:, df.columns[1:]]
	# df.to_pickle("summary_table.rad.pkl")
	# dfs.append(df_existing)
	# result =  dvcs_km15_rad()
	# dfs.append(result)
	# df = pd.concat(dfs)
	# df = df.reset_index()
	# df = df.loc[:, df.columns[1:]]
	# df.to_pickle("summary_table.rad.pkl")
	# result =  dvcs_km15_norad()
	# dfs.append(result)
	# df = pd.concat(dfs)
	# df = df.reset_index()
	# df = df.loc[:, df.columns[1:]]
	# df.to_pickle("summary_table.rad.pkl")
	df_existing = pd.read_pickle("summary_table.rad.pkl")
	df_existing.loc[df_existing.directory == "dvcsgen/dvcs_km15/rad", "directory"]   = "km15gen/dvcs_km15/rad"
	df_existing.loc[df_existing.directory == "dvcsgen/dvcs_km15/norad", "directory"] = "km15gen/dvcs_km15/rad"
	result =  pureBH_km15_rad()
	dfs.append(result)
	df = pd.concat(dfs)
	df = df.reset_index()
	df = df.loc[:, df.columns[1:]]
	df.to_pickle("summary_table.rad.pkl")
	result =  pureBH_km15_norad()
	dfs.append(result)
	df = pd.concat(dfs)
	df = df.reset_index()
	df = df.loc[:, df.columns[1:]]
	df.to_pickle("summary_table.rad.pkl")
