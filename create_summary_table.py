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
	binnum, bin_volume = np.loadtxt('volume_list.csv', skiprows = 1, delimiter = ',').T
	bin_volume = {int(binnum[i]): bin_volume[i] for i in range(len(binnum))}
	models = ["dvcs_vgg", "pureBH", "dvcs_km15"]
	if mode == "Rec":
		df_gen = pd.read_pickle("summary_table.gen.pkl")
		df_merged = {}

		for integrated_binnum in range(1, 147+1):
			df_merged[integrated_binnum] = []
		for integrated_binnum_gen in range(1, 147+1):
			print(integrated_binnum_gen)
			try:
				df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/dvcs_km15/excl_level_1/pkl/{}".format(integrated_binnum_gen))
			except:
				continue
			df.loc[:, "integrated_binnum_gen"] = integrated_binnum_gen
			df = df.astype({"integrated_binnum_gen": int, "phi_binnum_gen": int})
			df.loc[:, "n_gen"] = 0
			df.loc[:, "bin_volume_gen"] = 0
			for phi_binnum_gen in range(24):
				df.loc[df.phi_binnum_gen == phi_binnum_gen, "n_gen"] = df_gen.loc[(df_gen.directory == "dvcs_km15/fall2018_inb3") & (df_gen.integrated_binnum_gen == integrated_binnum_gen) & (df_gen.phi_binnum_gen == phi_binnum_gen)].n_entry.to_numpy()[0]
				df.loc[df.phi_binnum_gen == phi_binnum_gen, "bin_volume_gen"] = bin_volume[integrated_binnum_gen]/24.
				df.loc[df.phi_binnum_gen == phi_binnum_gen, "weights"] = df.loc[df.phi_binnum_gen == phi_binnum_gen].GenWeight \
				 * df.loc[df.phi_binnum_gen == phi_binnum_gen, "bin_volume_gen"] * luminosity_inb / df.loc[df.phi_binnum_gen == phi_binnum_gen].n_gen
			for integrated_binnum in range(1,147+1):
				df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum), :]
				df_merged[integrated_binnum].append(df_this_bin)
		for integrated_binnum_gen in range(148, 147+159+1):
			print(integrated_binnum_gen)
			try:
				df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/dvcs_km15_fringe/excl_level_1/pkl/{}".format(integrated_binnum_gen - 147))
			except:
				continue
			df.loc[:, "integrated_binnum_gen"] = integrated_binnum_gen
			df = df.astype({"integrated_binnum_gen": int, "phi_binnum_gen": int})
			df.loc[:, "n_gen"] = 0
			df.loc[:, "bin_volume_gen"] = 0
			for phi_binnum_gen in range(24):
				df.loc[df.phi_binnum_gen == phi_binnum_gen, "n_gen"] = df_gen.loc[(df_gen.directory == "dvcs_km15/fall2018_inb3") & (df_gen.integrated_binnum_gen == integrated_binnum_gen) & (df_gen.phi_binnum_gen == phi_binnum_gen)].n_entry.to_numpy()[0]
				df.loc[df.phi_binnum_gen == phi_binnum_gen, "bin_volume_gen"] = bin_volume[integrated_binnum_gen]/24.
				df.loc[df.phi_binnum_gen == phi_binnum_gen, "weights"] = df.loc[df.phi_binnum_gen == phi_binnum_gen].GenWeight \
				 * df.loc[df.phi_binnum_gen == phi_binnum_gen, "bin_volume_gen"] * luminosity_inb / df.loc[df.phi_binnum_gen == phi_binnum_gen].n_gen
			for integrated_binnum in range(1,147+1):
				df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum), :]
				df_merged[integrated_binnum].append(df_this_bin)

		for integrated_binnum in range(1, 147+1):
			df_merged[integrated_binnum] = pd.concat(df_merged[integrated_binnum])
			df_merged[integrated_binnum].to_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/dvcs_km15/excl_level_1/restructured/{}.pkl".format(integrated_binnum))


		df_merged = {}

		for integrated_binnum in range(1, 147+1):
			df_merged[integrated_binnum] = []
		for integrated_binnum_gen in range(1, 147+1):
			print(integrated_binnum_gen)
			try:
				df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/dvcs_km15/excl_level_1/pkl/{}".format(integrated_binnum_gen))
			except:
				continue
			df.loc[:, "integrated_binnum_gen"] = integrated_binnum_gen
			df = df.astype({"integrated_binnum_gen": int, "phi_binnum_gen": int})
			df.loc[:, "n_gen"] = 0
			df.loc[:, "bin_volume_gen"] = 0
			for phi_binnum_gen in range(24):
				df.loc[df.phi_binnum_gen == phi_binnum_gen, "n_gen"] = df_gen.loc[(df_gen.directory == "dvcs_km15/fall2018_outb3") & (df_gen.integrated_binnum_gen == integrated_binnum_gen) & (df_gen.phi_binnum_gen == phi_binnum_gen)].n_entry.to_numpy()[0]
				df.loc[df.phi_binnum_gen == phi_binnum_gen, "bin_volume_gen"] = bin_volume[integrated_binnum_gen]/24.
				df.loc[df.phi_binnum_gen == phi_binnum_gen, "weights"] = df.loc[df.phi_binnum_gen == phi_binnum_gen].GenWeight \
				 * df.loc[df.phi_binnum_gen == phi_binnum_gen, "bin_volume_gen"] * luminosity_outb / df.loc[df.phi_binnum_gen == phi_binnum_gen].n_gen
			for integrated_binnum in range(1,147+1):
				df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum), :]
				df_merged[integrated_binnum].append(df_this_bin)
		for integrated_binnum_gen in range(148, 147+159+1):
			print(integrated_binnum_gen)
			try:
				df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/dvcs_km15_fringe/excl_level_1/pkl/{}".format(integrated_binnum_gen - 147))
			except:
				continue
			df.loc[:, "integrated_binnum_gen"] = integrated_binnum_gen
			df = df.astype({"integrated_binnum_gen": int, "phi_binnum_gen": int})
			df.loc[:, "n_gen"] = 0
			df.loc[:, "bin_volume_gen"] = 0
			for phi_binnum_gen in range(24):
				df.loc[df.phi_binnum_gen == phi_binnum_gen, "n_gen"] = df_gen.loc[(df_gen.directory == "dvcs_km15/fall2018_outb3") & (df_gen.integrated_binnum_gen == integrated_binnum_gen) & (df_gen.phi_binnum_gen == phi_binnum_gen)].n_entry.to_numpy()[0]
				df.loc[df.phi_binnum_gen == phi_binnum_gen, "bin_volume_gen"] = bin_volume[integrated_binnum_gen]/24.
				df.loc[df.phi_binnum_gen == phi_binnum_gen, "weights"] = df.loc[df.phi_binnum_gen == phi_binnum_gen].GenWeight \
				 * df.loc[df.phi_binnum_gen == phi_binnum_gen, "bin_volume_gen"] * luminosity_outb / df.loc[df.phi_binnum_gen == phi_binnum_gen].n_gen
			for integrated_binnum in range(1,147+1):
				df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum), :]
				df_merged[integrated_binnum].append(df_this_bin)

		for integrated_binnum in range(1, 147+1):
			df_merged[integrated_binnum] = pd.concat(df_merged[integrated_binnum])
			df_merged[integrated_binnum].to_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/dvcs_km15/excl_level_1/restructured/{}.pkl".format(integrated_binnum))


		for integrated_binnum in range(1, 147+1):
			print(integrated_binnum)
			df = pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/dvcs_km15/excl_level_1/restructured/{}.pkl'.format(integrated_binnum))
			for phi_binnum in range(24):
				df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum) & (df.phi_binnum == phi_binnum)]
				# if not len(df_this_bin):
				# 	continue
				directory    = "sim_rad_rec_fall2018_inb/dvcs_km15"
				n_entry                      = len(df_this_bin)
				# n_entry_with_efficiency      = np.sum(df_this_bin.weight)
				n_entry_with_weight          = np.sum(df_this_bin.weights)
				weight_sum  = np.sum(df_this_bin.GenWeight)
				weight_avg  = np.mean(df_this_bin.GenWeight)
				if weight_sum:
					xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.GenxB)/weight_sum
					Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.GenQ2)/weight_sum
					t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.Gent)/weight_sum
					phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.Genphi)/weight_sum
				else:
					xB_avg      = 0
					Q2_avg      = 0
					t_avg       = 0
					phi_avg     = 0
				this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "n_entry": n_entry,
					"n_entry_with_weight": n_entry_with_weight, "xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
				df_summary   = pd.concat([df_summary, this_row])

		for integrated_binnum in range(1, 147+1):
			print(integrated_binnum)
			df = pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/dvcs_km15/excl_level_1/restructured/{}.pkl'.format(integrated_binnum))
			for phi_binnum in range(24):
				df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum) & (df.phi_binnum == phi_binnum)]
				# if not len(df_this_bin):
				# 	continue
				directory    = "sim_rad_rec_fall2018_outb/dvcs_km15"
				n_entry      = len(df_this_bin)
				# n_entry_with_efficiency      = np.sum(df_this_bin.weight)
				n_entry_with_weight          = np.sum(df_this_bin.weights)
				weight_sum   = np.sum(df_this_bin.GenWeight)
				weight_avg  = np.mean(df_this_bin.GenWeight)
				if weight_sum:
					xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.GenxB)/weight_sum
					Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.GenQ2)/weight_sum
					t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.Gent)/weight_sum
					phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.Genphi)/weight_sum
				else:
					xB_avg      = 0
					Q2_avg      = 0
					t_avg       = 0
					phi_avg     = 0
				this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "n_entry": n_entry,
					"n_entry_with_weight": n_entry_with_weight, "xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
				df_summary   = pd.concat([df_summary, this_row])

		df_summary.to_pickle("/volatile/clas12/sangbaek/clas12DVCS/summary_table.dvcs_km15.pkl")
	if mode == "Pi0_2Gamma":
		# #Rec - fall 2018 inbending ep->epg
		# df_merged = {}
		# for integrated_binnum in range(1, 147+1):
		# 	df_merged[integrated_binnum] = []
		
		# for directory in [1, 2, 3]:
		# 	for chunk in np.linspace(1, 20, 20).astype(int):
		# 		df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/pi0_2gamma/excl_level_1/pkl/{}/{}".format(directory, chunk))
		# 		for integrated_binnum in range(1, 147+1):
		# 			df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum), :]
		# 			df_merged[integrated_binnum].append(df_this_bin)
		# 			print(directory, chunk, integrated_binnum)
		# for integrated_binnum in range(1, 147+1):
		# 	df_merged[integrated_binnum] = pd.concat(df_merged[integrated_binnum])
		# 	df_merged[integrated_binnum].to_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/pi0_2gamma/excl_level_1/pkl/{}.pkl".format(integrated_binnum))

		# #Rec - fall 2018 outbending ep->epg
		# df_merged = {}
		# for integrated_binnum in range(1, 147+1):
		# 	df_merged[integrated_binnum] = []
		
		# for directory in [1, 2, 3]:
		# 	for chunk in np.linspace(1, 20, 20).astype(int):
		# 		df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/pi0_2gamma/excl_level_1/pkl/{}/{}".format(directory, chunk))
		# 		for integrated_binnum in range(1, 147+1):
		# 			df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum), :]
		# 			if directory == 1:
		# 				if df_this_bin.Q2.min()<1.4:
		# 					print(integrated_binnum)
		# 					continue
		# 			df_merged[integrated_binnum].append(df_this_bin)
		# 			# print(directory, chunk, integrated_binnum)
		# for integrated_binnum in range(1, 147+1):
		# 	df_merged[integrated_binnum] = pd.concat(df_merged[integrated_binnum])
		# 	df_merged[integrated_binnum].to_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/pi0_2gamma/excl_level_1/pkl/{}.pkl".format(integrated_binnum))


		for integrated_binnum in range(1, 147+1):
			print(integrated_binnum)
			df = pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/pi0_2gamma/excl_level_1/pkl/{}.pkl'.format(integrated_binnum))
			# df.loc[df.weight==0, "weight"] = 0.98 * df.loc[df.weight==0, "EFtof1bEfficiency"] * df.loc[df.weight==0, "PFtof1bEfficiency"]
			for phi_binnum in range(24):
				df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum) & (df.phi_binnum == phi_binnum)]
				# if not len(df_this_bin):
				# 	continue
				directory    = "sim_rad_rec_fall2018_inb/pi0_2gamma"
				n_entry                      = len(df_this_bin)
				n_entry_with_efficiency      = np.sum(df_this_bin.weight)
				xB_avg      = np.mean(df_this_bin.xB)
				Q2_avg      = np.mean(df_this_bin.Q2)
				t_avg       = np.mean(df_this_bin.t1)
				phi_avg     = np.mean(df_this_bin.phi1)
				this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "n_entry": n_entry,
					"xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
				df_summary   = pd.concat([df_summary, this_row])

		for integrated_binnum in range(1, 147+1):
			df = pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/pi0_2gamma/excl_level_1/pkl/{}.pkl'.format(integrated_binnum))
			# df.loc[df.weight==0, "weight"] = 0.98 * df.loc[df.weight==0, "EFtof1bEfficiency"] * df.loc[df.weight==0, "PFtof1bEfficiency"]
			for phi_binnum in range(24):
				df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum) & (df.phi_binnum == phi_binnum)]
				# if not len(df_this_bin):
				# 	continue
				directory    = "sim_rad_rec_fall2018_outb/pi0_2gamma"
				n_entry      = len(df_this_bin)
				n_entry_with_efficiency      = np.sum(df_this_bin.weight)
				xB_avg      = np.mean(df_this_bin.xB)
				Q2_avg      = np.mean(df_this_bin.Q2)
				t_avg       = np.mean(df_this_bin.t1)
				phi_avg     = np.mean(df_this_bin.phi1)
				this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "n_entry": n_entry,
					"xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
				df_summary   = pd.concat([df_summary, this_row])

		df_summary.to_pickle("/volatile/clas12/sangbaek/clas12DVCS/summary_table.pi0_2gamma.pkl")
	if mode == "Pi0_1Gamma":
		# #Rec - fall 2018 inbending ep->epg
		# df_merged = {}
		# for integrated_binnum in range(1, 147+1):
		# 	df_merged[integrated_binnum] = []
		
		# for directory in [1, 2, 3]:
		# 	for chunk in np.linspace(1, 20, 20).astype(int):
		# 		df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/pi0_1gamma/excl_level_1/pkl/{}/{}".format(directory, chunk))
		# 		for integrated_binnum in range(1, 147+1):
		# 			df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum), :]
		# 			df_merged[integrated_binnum].append(df_this_bin)
		# 			print(directory, chunk, integrated_binnum)
		# for integrated_binnum in range(1, 147+1):
		# 	df_merged[integrated_binnum] = pd.concat(df_merged[integrated_binnum])
		# 	df_merged[integrated_binnum].to_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/pi0_1gamma/excl_level_1/pkl/{}.pkl".format(integrated_binnum))

		# #Rec - fall 2018 outbending ep->epg
		# df_merged = {}
		# for integrated_binnum in range(1, 147+1):
		# 	df_merged[integrated_binnum] = []
		
		# for directory in [1, 2, 3]:
		# 	for chunk in np.linspace(1, 20, 20).astype(int):
		# 		df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/pi0_1gamma/excl_level_1/pkl/{}/{}".format(directory, chunk))
		# 		for integrated_binnum in range(1, 147+1):
		# 			df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum), :]
		# 			if directory == 1:
		# 				if df_this_bin.Q2.min()<1.4:
		# 					print(integrated_binnum)
		# 					continue
		# 			df_merged[integrated_binnum].append(df_this_bin)
		# 			# print(directory, chunk, integrated_binnum)
		# for integrated_binnum in range(1, 147+1):
		# 	df_merged[integrated_binnum] = pd.concat(df_merged[integrated_binnum])
		# 	df_merged[integrated_binnum].to_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/pi0_1gamma/excl_level_1/pkl/{}.pkl".format(integrated_binnum))


		for integrated_binnum in range(1, 147+1):
			print(integrated_binnum)
			df = pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/pi0_1gamma/excl_level_1/pkl/{}.pkl'.format(integrated_binnum))
			# df.loc[df.weight==0, "weight"] = 0.98 * df.loc[df.weight==0, "EFtof1bEfficiency"] * df.loc[df.weight==0, "PFtof1bEfficiency"]
			for phi_binnum in range(24):
				df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum) & (df.phi_binnum == phi_binnum)]
				# if not len(df_this_bin):
				# 	continue
				directory    = "sim_rad_rec_fall2018_inb/pi0_1gamma"
				n_entry                      = len(df_this_bin)
				# n_entry_with_efficiency      = np.sum(df_this_bin.weight)
				xB_avg      = np.mean(df_this_bin.xB)
				Q2_avg      = np.mean(df_this_bin.Q2)
				t_avg       = np.mean(df_this_bin.t1)
				phi_avg     = np.mean(df_this_bin.phi1)
				this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "n_entry": n_entry,
					"xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
				df_summary   = pd.concat([df_summary, this_row])

		for integrated_binnum in range(1, 147+1):
			df = pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/pi0_1gamma/excl_level_1/pkl/{}.pkl'.format(integrated_binnum))
			# df.loc[df.weight==0, "weight"] = 0.98 * df.loc[df.weight==0, "EFtof1bEfficiency"] * df.loc[df.weight==0, "PFtof1bEfficiency"]
			for phi_binnum in range(24):
				df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum) & (df.phi_binnum == phi_binnum)]
				# if not len(df_this_bin):
				# 	continue
				directory    = "sim_rad_rec_fall2018_outb/pi0_1gamma"
				n_entry      = len(df_this_bin)
				# n_entry_with_efficiency      = np.sum(df_this_bin.weight)
				xB_avg      = np.mean(df_this_bin.xB)
				Q2_avg      = np.mean(df_this_bin.Q2)
				t_avg       = np.mean(df_this_bin.t1)
				phi_avg     = np.mean(df_this_bin.phi1)
				this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "n_entry": n_entry,
					"xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
				df_summary   = pd.concat([df_summary, this_row])

		df_summary.to_pickle("/volatile/clas12/sangbaek/clas12DVCS/summary_table.pi0_1gamma.pkl")
	if mode == "Exp":
		# #Exp - fall 2018 inbending ep->epg
		# polarity = "inb"
		# qaTree     = pd.read_json("/volatile/clas12/sangbaek/clasqaDB/qadb/qa.rga_{}ending/qaTree.json".format(polarity)).T
		# chargeTree = pd.read_json("/volatile/clas12/sangbaek/clasqaDB/qadb/qa.rga_{}ending/chargeTree.json".format(polarity)).T
		# df_merged = []
		# for runnum in runlist_inb:
		# 	print(runnum)
		# 	df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_inb/dvcs/excl_level_2/pkl/{}".format(runnum))
		# 	for i in range(339):
		# 		filenum    = 5*i
		# 		qadb       = qaTree.loc[qaTree.index == runnum, filenum]
		# 		if not isinstance(qadb.values[0], dict):
		# 			break    
		# 		chargedb   = chargeTree.loc[qaTree.index == runnum, filenum]
		# 		evnumMin   = qadb.values[0]['evnumMin']
		# 		if (filenum>0) and (evnumMin == 0):
		# 			evnumMin = evnumMax+1
		# 		evnumMax   = qadb.values[0]['evnumMax']
		# 		defect     = qadb.values[0]['defect']
		# 		chargemin  = chargedb.values[0]['fcChargeMin']
		# 		chargemax  = chargedb.values[0]['fcChargeMax']
		# 		df.loc[(df.EventNum>=evnumMin) & (df.EventNum<=evnumMax), "defect"]  = defect
		# 		df.loc[(df.EventNum>=evnumMin) & (df.EventNum<=evnumMax), "chargemin"] = chargemin
		# 		df.loc[(df.EventNum>=evnumMin) & (df.EventNum<=evnumMax), "chargemax"] = chargemax
		# 	df = df.loc[df.defect == 0, :]
		# 	df = df.drop(columns = ["defect"])
		# 	df_merged.append(df)
		# df_merged = pd.concat(df_merged)
		# df_merged = df_merged.reset_index()
		# df_merged = df_merged.loc[:, df_merged.columns[1:]]
		# df_merged.to_pickle('/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_inb/dvcs/excl_level_2/pkl/fall2018_inb.pkl')

		# #Exp - fall 2018 inbending ep->epgg
		# df_merged = []
		# for runnum in runlist_inb:
		# 	print(runnum)
		# 	df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_inb/pi0/excl_level_2/pkl/{}".format(runnum))
		# 	for i in range(339):
		# 		filenum    = 5*i
		# 		qadb       = qaTree.loc[qaTree.index == runnum, filenum]
		# 		if not isinstance(qadb.values[0], dict):
		# 			break    
		# 		chargedb   = chargeTree.loc[qaTree.index == runnum, filenum]
		# 		evnumMin   = qadb.values[0]['evnumMin']
		# 		if (filenum>0) and (evnumMin == 0):
		# 			evnumMin = evnumMax+1
		# 		evnumMax   = qadb.values[0]['evnumMax']
		# 		defect     = qadb.values[0]['defect']
		# 		chargemin  = chargedb.values[0]['fcChargeMin']
		# 		chargemax  = chargedb.values[0]['fcChargeMax']
		# 		df.loc[(df.EventNum>=evnumMin) & (df.EventNum<=evnumMax), "defect"]  = defect
		# 		df.loc[(df.EventNum>=evnumMin) & (df.EventNum<=evnumMax), "chargemin"] = chargemin
		# 		df.loc[(df.EventNum>=evnumMin) & (df.EventNum<=evnumMax), "chargemax"] = chargemax
		# 	df = df.loc[df.defect == 0, :]
		# 	df = df.drop(columns = ["defect"])
		# 	df_merged.append(df)
		# df_merged = pd.concat(df_merged)
		# df_merged = df_merged.reset_index()
		# df_merged = df_merged.loc[:, df_merged.columns[1:]]
		# df_merged.to_pickle('/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_inb/pi0/excl_level_2/pkl/fall2018_inb.pkl')

		# #Exp - fall 2018 outbending ep->epg
		# polarity = "outb"
		# qaTree     = pd.read_json("/volatile/clas12/sangbaek/clasqaDB/qadb/qa.rga_{}ending/qaTree.json".format(polarity)).T
		# chargeTree = pd.read_json("/volatile/clas12/sangbaek/clasqaDB/qadb/qa.rga_{}ending/chargeTree.json".format(polarity)).T
		# df_merged = []
		# for runnum in runlist_outb:
		# 	print(runnum)
		# 	try:
		# 		df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_outb/dvcs/excl_level_2/pkl/{}".format(runnum))
		# 	except:
		# 		continue
		# 	for i in range(339):
		# 		filenum    = 5*i
		# 		qadb       = qaTree.loc[qaTree.index == runnum, filenum]
		# 		if not isinstance(qadb.values[0], dict):
		# 			break    
		# 		chargedb   = chargeTree.loc[qaTree.index == runnum, filenum]
		# 		evnumMin   = qadb.values[0]['evnumMin']
		# 		if (filenum>0) and (evnumMin == 0):
		# 			evnumMin = evnumMax+1
		# 		evnumMax   = qadb.values[0]['evnumMax']
		# 		defect     = qadb.values[0]['defect']
		# 		chargemin  = chargedb.values[0]['fcChargeMin']
		# 		chargemax  = chargedb.values[0]['fcChargeMax']
		# 		df.loc[(df.EventNum>=evnumMin) & (df.EventNum<=evnumMax), "defect"]  = defect
		# 		df.loc[(df.EventNum>=evnumMin) & (df.EventNum<=evnumMax), "chargemin"] = chargemin
		# 		df.loc[(df.EventNum>=evnumMin) & (df.EventNum<=evnumMax), "chargemax"] = chargemax
		# 	df = df.loc[df.defect == 0, :]
		# 	df = df.drop(columns = ["defect"])
		# 	df_merged.append(df)
		# df_merged = pd.concat(df_merged)
		# df_merged = df_merged.reset_index()
		# df_merged = df_merged.loc[:, df_merged.columns[1:]]
		# df_merged.to_pickle('/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_outb/dvcs/excl_level_2/pkl/fall2018_outb.pkl')

		# #Exp - fall 2018 outbending ep->epgg
		# df_merged = []
		# for runnum in runlist_outb:
		# 	print(runnum)
		# 	try:
		# 		df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_outb/pi0/excl_level_2/pkl/{}".format(runnum))
		# 	except:
		# 		continue
		# 	if not len(df):
		# 		continue
		# 	for i in range(339):
		# 		filenum    = 5*i
		# 		qadb       = qaTree.loc[qaTree.index == runnum, filenum]
		# 		if not isinstance(qadb.values[0], dict):
		# 			break    
		# 		chargedb   = chargeTree.loc[qaTree.index == runnum, filenum]
		# 		evnumMin   = qadb.values[0]['evnumMin']
		# 		if (filenum>0) and (evnumMin == 0):
		# 			evnumMin = evnumMax+1
		# 		evnumMax   = qadb.values[0]['evnumMax']
		# 		defect     = qadb.values[0]['defect']
		# 		chargemin  = chargedb.values[0]['fcChargeMin']
		# 		chargemax  = chargedb.values[0]['fcChargeMax']
		# 		df.loc[(df.EventNum>=evnumMin) & (df.EventNum<=evnumMax), "defect"]  = defect
		# 		df.loc[(df.EventNum>=evnumMin) & (df.EventNum<=evnumMax), "chargemin"] = chargemin
		# 		df.loc[(df.EventNum>=evnumMin) & (df.EventNum<=evnumMax), "chargemax"] = chargemax
		# 	df = df.loc[df.defect == 0, :]
		# 	df = df.drop(columns = ["defect"])
		# 	df_merged.append(df)
		# df_merged = pd.concat(df_merged)
		# df_merged = df_merged.reset_index()
		# df_merged = df_merged.loc[:, df_merged.columns[1:]]
		# df_merged.to_pickle('/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_outb/pi0/excl_level_2/pkl/fall2018_outb.pkl')

		df = pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_inb/dvcs/excl_level_2/pkl/fall2018_inb.pkl')
		for integrated_binnum in range(1, 147+1):
			for phi_binnum in range(24):
				df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum) & (df.phi_binnum == phi_binnum) & (df.weight != 0)]
				# if not len(df_this_bin):
				# 	continue
				directory    = "exp_fall2018_inb/dvcs"
				n_entry      = len(df_this_bin)
				n_entry_corrected = np.sum(1/df_this_bin.weight)
				xB_avg      = np.mean(df_this_bin.xB)
				Q2_avg      = np.mean(df_this_bin.Q2)
				t_avg       = np.mean(df_this_bin.t1)
				phi_avg     = np.mean(df_this_bin.phi1)
				this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "n_entry": n_entry, "n_entry_corrected": n_entry_corrected,
					"xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
				df_summary   = pd.concat([df_summary, this_row])

		df = pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_inb/pi0/excl_level_2/pkl/fall2018_inb.pkl')
		for integrated_binnum in range(1, 147+1):
			for phi_binnum in range(24):
				df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum) & (df.phi_binnum == phi_binnum) & (df.weight != 0)]
				# if not len(df_this_bin):
				# 	continue
				directory    = "exp_fall2018_inb/pi0"
				n_entry      = len(df_this_bin)
				n_entry_corrected = np.sum(1/df_this_bin.weight)
				xB_avg      = np.mean(df_this_bin.xB)
				Q2_avg      = np.mean(df_this_bin.Q2)
				t_avg       = np.mean(df_this_bin.t1)
				phi_avg     = np.mean(df_this_bin.phi1)
				this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "n_entry": n_entry, "n_entry_corrected": n_entry_corrected,
					"xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
				df_summary   = pd.concat([df_summary, this_row])


		df = pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_outb/dvcs/excl_level_2/pkl/fall2018_outb.pkl')
		for integrated_binnum in range(1, 147+1):
			for phi_binnum in range(24):
				df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum) & (df.phi_binnum == phi_binnum) & (df.weight != 0)]
				# if not len(df_this_bin):
				# 	continue
				directory    = "exp_fall2018_outb/dvcs"
				n_entry      = len(df_this_bin)
				n_entry_corrected = np.sum(1/df_this_bin.weight)
				xB_avg      = np.mean(df_this_bin.xB)
				Q2_avg      = np.mean(df_this_bin.Q2)
				t_avg       = np.mean(df_this_bin.t1)
				phi_avg     = np.mean(df_this_bin.phi1)
				this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "n_entry": n_entry, "n_entry_corrected": n_entry_corrected,
					"xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
				df_summary   = pd.concat([df_summary, this_row])

		df = pd.read_pickle('/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_outb/pi0/excl_level_2/pkl/fall2018_outb.pkl')
		for integrated_binnum in range(1, 147+1):
			for phi_binnum in range(24):
				df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum) & (df.phi_binnum == phi_binnum) & (df.weight != 0)]
				# if not len(df_this_bin):
				# 	continue
				directory    = "exp_fall2018_outb/pi0"
				n_entry      = len(df_this_bin)
				n_entry_corrected = np.sum(1/df_this_bin.weight)
				xB_avg      = np.mean(df_this_bin.xB)
				Q2_avg      = np.mean(df_this_bin.Q2)
				t_avg       = np.mean(df_this_bin.t1)
				phi_avg     = np.mean(df_this_bin.phi1)
				this_row     = pd.DataFrame([{"integrated_binnum": integrated_binnum, "phi_binnum": phi_binnum, "directory": directory, "n_entry": n_entry, "n_entry_corrected": n_entry_corrected,
					"xB_avg": xB_avg, "Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
				df_summary   = pd.concat([df_summary, this_row])

		df_summary.to_pickle("/volatile/clas12/sangbaek/clas12DVCS/summary_table.exp.pkl")
	if mode == "Gen_norad":
		#Gen - bulk
		phibins = [-1] + list(np.linspace(0, 360, 24+1)[1:-1]) + [361]

		for model in models:
			print(model, "inbending")
			if model == "dvcs_km15":
				string_suffix = "3"
				root_suffix   = "_3"
			else:
				string_suffix = ""
				root_suffix   = ""
			for integrated_binnum_gen in range(1, 147+1):
				df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_data/sim_norad_gen/{}/root/fall2018_inb{}/{}.pkl".format(model, string_suffix, integrated_binnum_gen))
				df.loc[:, "integrated_binnum_gen"] = integrated_binnum_gen
				for phi_binnum_gen in range(24):
					phimin = phibins[phi_binnum_gen]
					phimax = phibins[phi_binnum_gen+1]
					df_this_bin = df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), :]
					directory    = "{}/fall2018_inb{}".format(model, string_suffix)
					n_entry      = len(df_this_bin)
					weight_sum   = np.sum(df_this_bin.GenWeight)
					weight_avg  = np.mean(df_this_bin.GenWeight)
					xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.xB)/weight_sum
					Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.Q2)/weight_sum
					t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.t1)/weight_sum
					phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.phi1)/weight_sum
					df_this_bin.loc[:, "bin_volume_gen"] = bin_volume[integrated_binnum_gen]/24.
					df_this_bin.loc[:, "weights"] = df_this_bin.GenWeight * df_this_bin.loc[:, "bin_volume_gen"] * luminosity_inb / len(df_this_bin)
					n_generated = np.sum(df_this_bin.weights)
					km15_this_point  = printKM(xB_avg, Q2_avg, t_avg, np.radians(phi_avg))
					bh_this_point1   = printBHonly(xB_avg, Q2_avg, t_avg, np.radians(phi_avg))
					bh_this_point2   = printKM(xB_avg, Q2_avg, t_avg, np.radians(phi_avg), mode = 1)
					vgg_this_point  = printVGG(xB_avg, Q2_avg, t_avg, np.radians(phi_avg))
					this_row     = pd.DataFrame([{"integrated_binnum_gen": integrated_binnum_gen, "phi_binnum_gen": phi_binnum_gen, "directory": directory, "n_entry": n_entry,
						"bin_volume_gen": bin_volume[integrated_binnum_gen]/24., "n_generated": n_generated, "weight_sum": weight_sum, "weight_avg": weight_avg, "xB_avg": xB_avg,
						"Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg, "xsec_km15": km15_this_point, "xsec_bh": bh_this_point1,  "xsec_bh_2": bh_this_point2, "xsec_vgg": vgg_this_point}])
					print(this_row)
					df_summary   = pd.concat([df_summary, this_row])

			for integrated_binnum_gen in range(1, 147+1):
				print(model, "outbending")
				df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_data/sim_norad_gen/{}/root/fall2018_outb{}/{}.pkl".format(model, string_suffix, integrated_binnum_gen))
				df.loc[:, "integrated_binnum_gen"] = integrated_binnum_gen
				for phi_binnum_gen in range(24):
					phimin = phibins[phi_binnum_gen]
					phimax = phibins[phi_binnum_gen+1]
					df_this_bin = df.loc[ (df.phi1>=phimin) & (df.phi1<phimax), :]
					directory    = "{}/fall2018_out{}".format(model, string_suffix)
					n_entry      = len(df_this_bin)
					weight_sum   = np.sum(df_this_bin.GenWeight)
					weight_avg  = np.mean(df_this_bin.GenWeight)
					df_this_bin.loc[:, "bin_volume_gen"] = bin_volume[integrated_binnum_gen]/24.
					df_this_bin.loc[:, "weights"] =df_this_bin.GenWeight * df_this_bin.loc[:, "bin_volume_gen"] * luminosity_inb / len(df_this_bin)
					n_generated = np.sum(df_this_bin.weights)
					xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.xB)/weight_sum
					Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.Q2)/weight_sum
					t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.t1)/weight_sum
					phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.phi1)/weight_sum
					df_this_bin.loc[:, "bin_volume_gen"] = bin_volume[integrated_binnum_gen]/24.
					df_this_bin.loc[:, "weights"] = df_this_bin.GenWeight * df_this_bin.loc[:, "bin_volume_gen"] * luminosity_inb / len(df_this_bin)
					n_generated = np.sum(df_this_bin.weights)
					km15_this_point  = printKM(xB_avg, Q2_avg, t_avg, np.radians(phi_avg))
					bh_this_point1   = printBHonly(xB_avg, Q2_avg, t_avg, np.radians(phi_avg))
					bh_this_point2   = printKM(xB_avg, Q2_avg, t_avg, np.radians(phi_avg), mode = 1)
					vgg_this_point  = printVGG(xB_avg, Q2_avg, t_avg, np.radians(phi_avg))
					this_row     = pd.DataFrame([{"integrated_binnum_gen": integrated_binnum_gen, "phi_binnum_gen": phi_binnum_gen, "directory": directory, "n_entry": n_entry,
						"bin_volume_gen": bin_volume[integrated_binnum_gen]/24., "n_generated": n_generated, "weight_sum": weight_sum, "weight_avg": weight_avg, "xB_avg": xB_avg,
						"Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg, "xsec_km15": km15_this_point, "xsec_bh": bh_this_point1,  "xsec_bh_2": bh_this_point2, "xsec_vgg": vgg_this_point}])
					print(this_row)
					df_summary   = pd.concat([df_summary, this_row])

		df_summary.to_pickle("/volatile/clas12/sangbaek/clas12DVCS/summary_table.gen.norad.pkl")
	if mode == "Gen":
		#Gen - bulk
		for model in models:
			print(model, "inbending")
			if model == "dvcs_km15":
				string_suffix = "3"
				root_suffix   = "_3"
			else:
				string_suffix = ""
				root_suffix   = ""
			for integrated_binnum_gen in range(1, 147+1):
				df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_data/sim_rad_gen/{}/pkl/fall2018_inb{}/{}.pkl".format(model, string_suffix, integrated_binnum_gen))
				df.loc[:, "integrated_binnum_gen"] = integrated_binnum_gen
				df = df.astype({"integrated_binnum_gen": int, "phi_binnum_gen": int})
				for phi_binnum_gen in range(24):
					df_this_bin  = df.loc[df.phi_binnum_gen == phi_binnum_gen]
					directory    = "{}/fall2018_inb{}".format(model, string_suffix)
					n_entry      = len(df_this_bin)
					if n_entry:
						weight_sum   = np.sum(df_this_bin.GenWeight)
						weight_avg  = np.mean(df_this_bin.GenWeight)
						xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.GenxB)/weight_sum
						Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.GenQ2)/weight_sum
						t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.Gent)/weight_sum
						phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.Genphi)/weight_sum
						df_this_bin.loc[:, "bin_volume_gen"] = bin_volume[integrated_binnum_gen]/24.
						df_this_bin.loc[:, "weights"] = df_this_bin.GenWeight * df_this_bin.loc[:, "bin_volume_gen"] * luminosity_inb / len(df_this_bin)
						n_generated = np.sum(df_this_bin.weights)
					else:
						weight_sum  = 0#np.sum(df_this_bin.GenWeight)
						weight_avg  = 0#np.mean(df_this_bin.GenWeight)
						xB_avg      = 0#np.sum(df_this_bin.GenWeight * df_this_bin.GenxB)/weight_sum
						Q2_avg      = 0#np.sum(df_this_bin.GenWeight * df_this_bin.GenQ2)/weight_sum
						t_avg       = 0#np.sum(df_this_bin.GenWeight * df_this_bin.Gent)/weight_sum
						phi_avg     = 0#np.sum(df_this_bin.GenWeight * df_this_bin.Genphi)/weight_sum
						n_generated = 0#np.sum(df_this_bin.weights)
					this_row     = pd.DataFrame([{"integrated_binnum_gen": integrated_binnum_gen, "phi_binnum_gen": phi_binnum_gen, "directory": directory, "n_entry": n_entry,
						"bin_volume_gen": bin_volume[integrated_binnum_gen]/24., "n_generated": n_generated, "weight_sum": weight_sum, "weight_avg": weight_avg, "xB_avg": xB_avg,
						"Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
					print(this_row)
					df_summary   = pd.concat([df_summary, this_row])

			for integrated_binnum_gen in range(1, 159+1):
				df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_data/sim_rad_gen/{}_fringe/pkl/fall2018_inb{}/{}.pkl".format(model, string_suffix, integrated_binnum_gen))
				df.loc[:, "integrated_binnum_gen"] = integrated_binnum_gen + 147
				df = df.astype({"integrated_binnum_gen": int, "phi_binnum_gen": int})
				for phi_binnum_gen in range(24):
					df_this_bin  = df.loc[df.phi_binnum_gen == phi_binnum_gen]
					directory    = "{}/fall2018_inb{}".format(model, string_suffix)
					n_entry      = len(df_this_bin)
					if n_entry:
						weight_sum   = np.sum(df_this_bin.GenWeight)
						weight_avg  = np.mean(df_this_bin.GenWeight)
						xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.GenxB)/weight_sum
						Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.GenQ2)/weight_sum
						t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.Gent)/weight_sum
						phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.Genphi)/weight_sum
						df_this_bin.loc[:, "bin_volume_gen"] = bin_volume[integrated_binnum_gen]/24.
						df_this_bin.loc[:, "weights"] = df_this_bin.GenWeight * df_this_bin.loc[:, "bin_volume_gen"] * luminosity_inb / len(df_this_bin)
						n_generated = np.sum(df_this_bin.weights)
					else:
						weight_sum  = 0#np.sum(df_this_bin.GenWeight)
						weight_avg  = 0#np.mean(df_this_bin.GenWeight)
						xB_avg      = 0#np.sum(df_this_bin.GenWeight * df_this_bin.GenxB)/weight_sum
						Q2_avg      = 0#np.sum(df_this_bin.GenWeight * df_this_bin.GenQ2)/weight_sum
						t_avg       = 0#np.sum(df_this_bin.GenWeight * df_this_bin.Gent)/weight_sum
						phi_avg     = 0#np.sum(df_this_bin.GenWeight * df_this_bin.Genphi)/weight_sum
						n_generated = 0#np.sum(df_this_bin.weights)
					this_row     = pd.DataFrame([{"integrated_binnum_gen": integrated_binnum_gen + 147, "phi_binnum_gen": phi_binnum_gen, "directory": directory, "n_entry": n_entry,
						"bin_volume_gen": bin_volume[integrated_binnum_gen]/24., "n_generated": n_generated, "weight_sum": weight_sum, "weight_avg": weight_avg, "xB_avg": xB_avg,
						"Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
					print(this_row)
					df_summary   = pd.concat([df_summary, this_row])


			for integrated_binnum_gen in range(1, 147+1):
				df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_data/sim_rad_gen/{}/pkl/fall2018_outb{}/{}.pkl".format(model, string_suffix, integrated_binnum_gen))
				df.loc[:, "integrated_binnum_gen"] = integrated_binnum_gen
				df = df.astype({"integrated_binnum_gen": int, "phi_binnum_gen": int})
				for phi_binnum_gen in range(24):
					df_this_bin  = df.loc[df.phi_binnum_gen == phi_binnum_gen]
					directory    = "{}/fall2018_outb{}".format(model, string_suffix)
					n_entry      = len(df_this_bin)
					if n_entry:
						weight_sum   = np.sum(df_this_bin.GenWeight)
						weight_avg  = np.mean(df_this_bin.GenWeight)
						xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.GenxB)/weight_sum
						Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.GenQ2)/weight_sum
						t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.Gent)/weight_sum
						phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.Genphi)/weight_sum
						df_this_bin.loc[:, "bin_volume_gen"] = bin_volume[integrated_binnum_gen]/24.
						df_this_bin.loc[:, "weights"] = df_this_bin.GenWeight * df_this_bin.loc[:, "bin_volume_gen"] * luminosity_inb / len(df_this_bin)
						n_generated = np.sum(df_this_bin.weights)
					else:
						weight_sum  = 0#np.sum(df_this_bin.GenWeight)
						weight_avg  = 0#np.mean(df_this_bin.GenWeight)
						xB_avg      = 0#np.sum(df_this_bin.GenWeight * df_this_bin.GenxB)/weight_sum
						Q2_avg      = 0#np.sum(df_this_bin.GenWeight * df_this_bin.GenQ2)/weight_sum
						t_avg       = 0#np.sum(df_this_bin.GenWeight * df_this_bin.Gent)/weight_sum
						phi_avg     = 0#np.sum(df_this_bin.GenWeight * df_this_bin.Genphi)/weight_sum
						n_generated = 0#np.sum(df_this_bin.weights)
					this_row     = pd.DataFrame([{"integrated_binnum_gen": integrated_binnum_gen, "phi_binnum_gen": phi_binnum_gen, "directory": directory, "n_entry": n_entry,
						"bin_volume_gen": bin_volume[integrated_binnum_gen]/24., "n_generated": n_generated, "weight_sum": weight_sum, "weight_avg": weight_avg, "xB_avg": xB_avg,
						"Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
					print(this_row)
					df_summary   = pd.concat([df_summary, this_row])

			for integrated_binnum_gen in range(1, 159+1):
				df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_data/sim_rad_gen/{}_fringe/pkl/fall2018_outb{}/{}.pkl".format(model, string_suffix, integrated_binnum_gen))
				df.loc[:, "integrated_binnum_gen"] = integrated_binnum_gen + 147
				df = df.astype({"integrated_binnum_gen": int, "phi_binnum_gen": int})
				for phi_binnum_gen in range(24):
					df_this_bin  = df.loc[df.phi_binnum_gen == phi_binnum_gen]
					directory    = "{}/fall2018_outb{}".format(model, string_suffix)
					n_entry      = len(df_this_bin)
					if n_entry:
						weight_sum   = np.sum(df_this_bin.GenWeight)
						weight_avg  = np.mean(df_this_bin.GenWeight)
						xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.GenxB)/weight_sum
						Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.GenQ2)/weight_sum
						t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.Gent)/weight_sum
						phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.Genphi)/weight_sum
						df_this_bin.loc[:, "bin_volume_gen"] = bin_volume[integrated_binnum_gen]/24.
						df_this_bin.loc[:, "weights"] = df_this_bin.GenWeight * df_this_bin.loc[:, "bin_volume_gen"] * luminosity_inb / len(df_this_bin)
						n_generated = np.sum(df_this_bin.weights)
					else:
						weight_sum  = 0#np.sum(df_this_bin.GenWeight)
						weight_avg  = 0#np.mean(df_this_bin.GenWeight)
						xB_avg      = 0#np.sum(df_this_bin.GenWeight * df_this_bin.GenxB)/weight_sum
						Q2_avg      = 0#np.sum(df_this_bin.GenWeight * df_this_bin.GenQ2)/weight_sum
						t_avg       = 0#np.sum(df_this_bin.GenWeight * df_this_bin.Gent)/weight_sum
						phi_avg     = 0#np.sum(df_this_bin.GenWeight * df_this_bin.Genphi)/weight_sum
						n_generated = 0#np.sum(df_this_bin.weights)
					this_row     = pd.DataFrame([{"integrated_binnum_gen": integrated_binnum_gen + 147, "phi_binnum_gen": phi_binnum_gen, "directory": directory, "n_entry": n_entry,
						"bin_volume_gen": bin_volume[integrated_binnum_gen]/24., "n_generated": n_generated, "weight_sum": weight_sum, "weight_avg": weight_avg, "xB_avg": xB_avg,
						"Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
					print(this_row)
					df_summary   = pd.concat([df_summary, this_row])
			if model =="pureBH":
				for integrated_binnum_gen in range(1, 147+1):
					df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_data/sim_rad_gen/{}/pkl/fall2018_inb{}/{}.pkl".format(model, "3", integrated_binnum_gen))
					df.loc[:, "integrated_binnum_gen"] = integrated_binnum_gen
					df = df.astype({"integrated_binnum_gen": int, "phi_binnum_gen": int})
					for phi_binnum_gen in range(24):
						df_this_bin  = df.loc[df.phi_binnum_gen == phi_binnum_gen]
						directory    = "{}/fall2018_inb{}".format(model, "3")
						n_entry      = len(df_this_bin)
					if n_entry:
						weight_sum   = np.sum(df_this_bin.GenWeight)
						weight_avg  = np.mean(df_this_bin.GenWeight)
						xB_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.GenxB)/weight_sum
						Q2_avg      = np.sum(df_this_bin.GenWeight * df_this_bin.GenQ2)/weight_sum
						t_avg       = np.sum(df_this_bin.GenWeight * df_this_bin.Gent)/weight_sum
						phi_avg     = np.sum(df_this_bin.GenWeight * df_this_bin.Genphi)/weight_sum
						df_this_bin.loc[:, "bin_volume_gen"] = bin_volume[integrated_binnum_gen]/24.
						df_this_bin.loc[:, "weights"] = df_this_bin.GenWeight * df_this_bin.loc[:, "bin_volume_gen"] * luminosity_inb / len(df_this_bin)
						n_generated = np.sum(df_this_bin.weights)
					else:
						weight_sum  = 0#np.sum(df_this_bin.GenWeight)
						weight_avg  = 0#np.mean(df_this_bin.GenWeight)
						xB_avg      = 0#np.sum(df_this_bin.GenWeight * df_this_bin.GenxB)/weight_sum
						Q2_avg      = 0#np.sum(df_this_bin.GenWeight * df_this_bin.GenQ2)/weight_sum
						t_avg       = 0#np.sum(df_this_bin.GenWeight * df_this_bin.Gent)/weight_sum
						phi_avg     = 0#np.sum(df_this_bin.GenWeight * df_this_bin.Genphi)/weight_sum
						n_generated = 0#np.sum(df_this_bin.weights)
					this_row     = pd.DataFrame([{"integrated_binnum_gen": integrated_binnum_gen, "phi_binnum_gen": phi_binnum_gen, "directory": directory, "n_entry": n_entry,
						"bin_volume_gen": bin_volume[integrated_binnum_gen]/24., "n_generated": n_generated, "weight_sum": weight_sum, "weight_avg": weight_avg, "xB_avg": xB_avg,
						"Q2_avg": Q2_avg, "t_avg": t_avg, "phi_avg": phi_avg}])
					print(this_row)
					df_summary   = pd.concat([df_summary, this_row])

		df_summary.to_pickle("/volatile/clas12/sangbaek/clas12DVCS/summary_table.gen.pkl")

if __name__ == "__main__":

	for mode in ["Gen"]:
		main(mode)


