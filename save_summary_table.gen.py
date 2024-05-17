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
	if mode == "Gen":
		#Gen - bulk
		for model in models:
			print(model, "inbending")
			if model == "dvcs_km15":
				string_suffix = "3"
			else:
				string_suffix = ""
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

			if model == "dvcs_km15":
				string_suffix = "3_45nA"
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

				string_suffix = "3_50nA"

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

	df_summary = df_summary.reset_index()
	df_summary = df_summary.loc[:, df_summary.columns[1:]]
	df_summary.to_pickle("/volatile/clas12/sangbaek/clas12DVCS/summary_table.gen.pkl")

if __name__ == "__main__":

	main("Gen")


