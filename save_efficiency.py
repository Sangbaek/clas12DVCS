import numpy as np
import pandas as pd
from utils.const import *
from utils.physics import *
import os

import argparse
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
warnings.simplefilter(action='ignore', category=pd.errors.SettingWithCopyWarning)

def rebin(epg_exp, new_bins):
    epg_exp_rebinned = np.zeros(len(new_bins)-1)
    assert new_bins[0] == 0
    assert new_bins[-1] == 24
    for i in range(len(new_bins)-1):
        bin_start = new_bins[i]
        bin_end   = new_bins[i+1]
        epg_exp_rebinned[i] = np.sum(epg_exp[bin_start: bin_end])
    return epg_exp_rebinned, new_bins

def assign_Pthetaefficiency(df, ptheta_efficiency_sector_1, ptheta_efficiency_sector_2, ptheta_efficiency_sector_3, ptheta_efficiency_sector_4, ptheta_efficiency_sector_5, ptheta_efficiency_sector_6, ptheta_efficiency_cd, ptheta_bins_fd, ptheta_bins_cd):
    df.loc[:, "Pthetaefficiency"] = 1
    for i in range(len(ptheta_bins_fd)-1):
        df.loc[(df.Psector == 1) & (df.Ptheta>ptheta_bins_fd[i]) & (df.Ptheta<ptheta_bins_fd[i+1]), "Pthetaefficiency"] = ptheta_efficiency_sector_1[i]
        df.loc[(df.Psector == 2) & (df.Ptheta>ptheta_bins_fd[i]) & (df.Ptheta<ptheta_bins_fd[i+1]), "Pthetaefficiency"] = ptheta_efficiency_sector_2[i]
        df.loc[(df.Psector == 3) & (df.Ptheta>ptheta_bins_fd[i]) & (df.Ptheta<ptheta_bins_fd[i+1]), "Pthetaefficiency"] = ptheta_efficiency_sector_3[i]
        df.loc[(df.Psector == 4) & (df.Ptheta>ptheta_bins_fd[i]) & (df.Ptheta<ptheta_bins_fd[i+1]), "Pthetaefficiency"] = ptheta_efficiency_sector_4[i]
        df.loc[(df.Psector == 5) & (df.Ptheta>ptheta_bins_fd[i]) & (df.Ptheta<ptheta_bins_fd[i+1]), "Pthetaefficiency"] = ptheta_efficiency_sector_5[i]
        df.loc[(df.Psector == 6) & (df.Ptheta>ptheta_bins_fd[i]) & (df.Ptheta<ptheta_bins_fd[i+1]), "Pthetaefficiency"] = ptheta_efficiency_sector_6[i]
    for i in range(len(ptheta_bins_cd)-1):
        df.loc[(df.Psector >  7) & (df.Ptheta>ptheta_bins_cd[i]) & (df.Ptheta<ptheta_bins_cd[i+1]), "Pthetaefficiency"] = ptheta_efficiency_cd[i]
    return df

def assign_Ppefficiency(df, pp_efficiency_sector_1, pp_efficiency_sector_2, pp_efficiency_sector_3, pp_efficiency_sector_4, pp_efficiency_sector_5, pp_efficiency_sector_6, pp_efficiency_cd, pp_bins_fd, pp_bins_cd):
    df.loc[:, "Ppefficiency"] = 1
    for i in range(len(pp_bins_fd)-1):
        df.loc[(df.Psector == 1) & (df.Pp>pp_bins_fd[i]) & (df.Pp<pp_bins_fd[i+1]), "Ppefficiency"] = pp_efficiency_sector_1[i]
        df.loc[(df.Psector == 2) & (df.Pp>pp_bins_fd[i]) & (df.Pp<pp_bins_fd[i+1]), "Ppefficiency"] = pp_efficiency_sector_2[i]
        df.loc[(df.Psector == 3) & (df.Pp>pp_bins_fd[i]) & (df.Pp<pp_bins_fd[i+1]), "Ppefficiency"] = pp_efficiency_sector_3[i]
        df.loc[(df.Psector == 4) & (df.Pp>pp_bins_fd[i]) & (df.Pp<pp_bins_fd[i+1]), "Ppefficiency"] = pp_efficiency_sector_4[i]
        df.loc[(df.Psector == 5) & (df.Pp>pp_bins_fd[i]) & (df.Pp<pp_bins_fd[i+1]), "Ppefficiency"] = pp_efficiency_sector_5[i]
        df.loc[(df.Psector == 6) & (df.Pp>pp_bins_fd[i]) & (df.Pp<pp_bins_fd[i+1]), "Ppefficiency"] = pp_efficiency_sector_6[i]
    for i in range(len(pp_bins_cd)-1):
        df.loc[(df.Psector >  7) & (df.Pp>pp_bins_cd[i]) & (df.Pp<pp_bins_cd[i+1]), "Ppefficiency"] = pp_efficiency_cd[i]
    return df

def assign_Pphiefficiency(df, pphi_efficiency_sector_1, pphi_efficiency_sector_2, pphi_efficiency_sector_3, pphi_efficiency_sector_4, pphi_efficiency_sector_5, pphi_efficiency_sector_6, pphi_efficiency_cd, pphi_bins_fd, pphi_bins_cd):
    df.loc[:, "Pphiefficiency"] = 1
    for i in range(len(pphi_bins_fd)-1):
        df.loc[(df.Psector == 1) & (df.Pphi>pphi_bins_fd[i]) & (df.Pphi<pphi_bins_fd[i+1]), "Pphiefficiency"] = pphi_efficiency_sector_1[i]
        df.loc[(df.Psector == 2) & (df.Pphi>pphi_bins_fd[i]) & (df.Pphi<pphi_bins_fd[i+1]), "Pphiefficiency"] = pphi_efficiency_sector_2[i]
        df.loc[(df.Psector == 3) & (df.Pphi>pphi_bins_fd[i]) & (df.Pphi<pphi_bins_fd[i+1]), "Pphiefficiency"] = pphi_efficiency_sector_3[i]
        df.loc[(df.Psector == 4) & (df.Pphi>pphi_bins_fd[i]) & (df.Pphi<pphi_bins_fd[i+1]), "Pphiefficiency"] = pphi_efficiency_sector_4[i]
        df.loc[(df.Psector == 5) & (df.Pphi>pphi_bins_fd[i]) & (df.Pphi<pphi_bins_fd[i+1]), "Pphiefficiency"] = pphi_efficiency_sector_5[i]
        df.loc[(df.Psector == 6) & (df.Pphi>pphi_bins_fd[i]) & (df.Pphi<pphi_bins_fd[i+1]), "Pphiefficiency"] = pphi_efficiency_sector_6[i]
    for i in range(len(pphi_bins_cd)-1):
        df.loc[(df.Psector >  7) & (df.Pphi>pphi_bins_cd[i]) & (df.Pphi<pphi_bins_cd[i+1]), "Pphiefficiency"] = pphi_efficiency_cd[i]
    return df


# df_summary_table_rebinned = pd.read_csv("df_summary_table_rebinned.csv")
df_summary_table_rebinned = pd.read_csv("summary_table.rebinned.08212024.csv")

ptheta_bins_inb_fd             = np.linspace(7,  40, 8)
ptheta_bins_inb_cd             = np.linspace(40, 65, 11)

pp_bins_inb_fd             = np.linspace(0.442,  1.133, 8)
pp_bins_inb_cd             = np.linspace(0.336,  1.133, 11)

pphi_bins_inb_fd             = np.linspace(-180, 180, 37)
pphi_bins_inb_cd             = np.linspace(-180, 180, 37)

#initialize the efficiency map
ptheta_efficiency_inb_sector_1 = np.ones(7)
ptheta_efficiency_inb_sector_2 = np.ones(7)
ptheta_efficiency_inb_sector_3 = np.ones(7)
ptheta_efficiency_inb_sector_4 = np.ones(7)
ptheta_efficiency_inb_sector_5 = np.ones(7)
ptheta_efficiency_inb_sector_6 = np.ones(7)
ptheta_efficiency_inb_cd       = np.ones(10)

pp_efficiency_inb_sector_1 = np.ones(7)
pp_efficiency_inb_sector_2 = np.ones(7)
pp_efficiency_inb_sector_3 = np.ones(7)
pp_efficiency_inb_sector_4 = np.ones(7)
pp_efficiency_inb_sector_5 = np.ones(7)
pp_efficiency_inb_sector_6 = np.ones(7)
pp_efficiency_inb_cd       = np.ones(10)


pphi_efficiency_inb_sector_1 = np.ones(36)
pphi_efficiency_inb_sector_2 = np.ones(36)
pphi_efficiency_inb_sector_3 = np.ones(36)
pphi_efficiency_inb_sector_4 = np.ones(36)
pphi_efficiency_inb_sector_5 = np.ones(36)
pphi_efficiency_inb_sector_6 = np.ones(36)
pphi_efficiency_inb_cd       = np.ones(36)

effective_current_inb = (40*charge_inb_40nA + 45* charge_inb_45nA + 50*charge_inb_50nA + 55*charge_inb_55nA)/charge_inb
chunk_inb             = [1, 2, 3]

#inbending
for inbending_trial in range(10):
    for ptheta_trial in range(3):

        df_exp_epg_inbs        = []
        df_sim_dvcs_inbs       = []
        print(ptheta_trial, "ptheta_trial", "inb sector 1", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_1], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_1], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_1])
        print(ptheta_trial, "ptheta_trial", "inb sector 2", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_2], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_2], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_2])
        print(ptheta_trial, "ptheta_trial", "inb sector 3", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_3], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_3], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_3])
        print(ptheta_trial, "ptheta_trial", "inb sector 4", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_4], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_4], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_4])
        print(ptheta_trial, "ptheta_trial", "inb sector 5", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_5], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_5], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_5])
        print(ptheta_trial, "ptheta_trial", "inb sector 6", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_6], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_6], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_6])
        print(ptheta_trial, "ptheta_trial", "inb CD"      , ["{:.3f}".format(i) for i in ptheta_efficiency_inb_cd      ], ["{:.3f}".format(i) for i in pp_efficiency_inb_cd      ], ["{:.3f}".format(i) for i in pphi_efficiency_inb_cd      ])

        for integrated_binnum in range(1, 147+1):
            this_rebinned_phi_binnums = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"].to_numpy()
            this_rebinned_phi_widths  = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"].to_numpy()

            try:
                df_exp_epg_inb           =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_inb/dvcs/excl_level_2/restructured_7_nominal/{}.pkl".format(integrated_binnum))
                df_exp_epg_inb.loc[df_exp_epg_inb.efficiency == 0, "efficiency"] =  1
                df_exp_epg_inb.loc[:, "contamination"] = 1
            except:
                # print("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_inb/dvcs/excl_level_2/restructured_7_nominal/{}.pkl does not exist.".format(integrated_binnum))
                continue

            # if len(df_exp_epg_inb)<50:
            if df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum), "active_bin_inb"].sum() == 0:
                continue

            try:
                df_exp_pi0_inb            =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_inb/pi0/excl_level_2/restructured_7_nominal/{}.pkl".format(integrated_binnum))
            except:
                df_exp_pi0_inb            =  pd.DataFrame({var: [] for var in ['Pp', 'Ptheta', 'Pphi', 'efficiency']})
            # df_exp_pi0_inb.loc[df_exp_pi0_inb.efficiency == 0, "efficiency"]  = 1
            try:
                df_sim_pi0_inb            =  pd.concat([pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/pi0_2gamma/excl_level_1/restructured_7_nominal/{}/{}.pkl".format(chunk, integrated_binnum)) for chunk in chunk_inb])
                df_sim_pi0_inb            = assign_Pthetaefficiency(df_sim_pi0_inb, ptheta_efficiency_inb_sector_1, ptheta_efficiency_inb_sector_2, ptheta_efficiency_inb_sector_3, ptheta_efficiency_inb_sector_4, ptheta_efficiency_inb_sector_5, ptheta_efficiency_inb_sector_6, ptheta_efficiency_inb_cd, ptheta_bins_inb_fd, ptheta_bins_inb_cd)
                df_sim_pi0_inb            = assign_Ppefficiency(df_sim_pi0_inb, pp_efficiency_inb_sector_1, pp_efficiency_inb_sector_2, pp_efficiency_inb_sector_3, pp_efficiency_inb_sector_4, pp_efficiency_inb_sector_5, pp_efficiency_inb_sector_6, pp_efficiency_inb_cd, pp_bins_inb_fd, pp_bins_inb_cd)
                df_sim_pi0_inb            = assign_Pphiefficiency(df_sim_pi0_inb, pphi_efficiency_inb_sector_1, pphi_efficiency_inb_sector_2, pphi_efficiency_inb_sector_3, pphi_efficiency_inb_sector_4, pphi_efficiency_inb_sector_5, pphi_efficiency_inb_sector_6, pphi_efficiency_inb_cd, pphi_bins_inb_fd, pphi_bins_inb_cd)
            except:
                df_sim_pi0_inb            = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})

            try:
                df_sim_pi0_1gamma_inb     =  pd.concat([pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/pi0_1gamma/excl_level_1/restructured_7_nominal/{}/{}.pkl".format(chunk, integrated_binnum)) for chunk in chunk_inb])
                df_sim_pi0_1gamma_inb     = assign_Pthetaefficiency(df_sim_pi0_1gamma_inb, ptheta_efficiency_inb_sector_1, ptheta_efficiency_inb_sector_2, ptheta_efficiency_inb_sector_3, ptheta_efficiency_inb_sector_4, ptheta_efficiency_inb_sector_5, ptheta_efficiency_inb_sector_6, ptheta_efficiency_inb_cd, ptheta_bins_inb_fd, ptheta_bins_inb_cd)
                df_sim_pi0_1gamma_inb     = assign_Ppefficiency(df_sim_pi0_1gamma_inb, pp_efficiency_inb_sector_1, pp_efficiency_inb_sector_2, pp_efficiency_inb_sector_3, pp_efficiency_inb_sector_4, pp_efficiency_inb_sector_5, pp_efficiency_inb_sector_6, pp_efficiency_inb_cd, pp_bins_inb_fd, pp_bins_inb_cd)
                df_sim_pi0_1gamma_inb     = assign_Pphiefficiency(df_sim_pi0_1gamma_inb, pphi_efficiency_inb_sector_1, pphi_efficiency_inb_sector_2, pphi_efficiency_inb_sector_3, pphi_efficiency_inb_sector_4, pphi_efficiency_inb_sector_5, pphi_efficiency_inb_sector_6, pphi_efficiency_inb_cd, pphi_bins_inb_fd, pphi_bins_inb_cd)
            except:
                df_sim_pi0_1gamma_inb      = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})
            try:
                df_sim_nobkgmerging_inb   =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/dvcs_km15/excl_level_1/restructured_7_nominal/{}.pkl".format(integrated_binnum))
                df_sim_nobkgmerging_inb   = assign_Pthetaefficiency(df_sim_nobkgmerging_inb, ptheta_efficiency_inb_sector_1, ptheta_efficiency_inb_sector_2, ptheta_efficiency_inb_sector_3, ptheta_efficiency_inb_sector_4, ptheta_efficiency_inb_sector_5, ptheta_efficiency_inb_sector_6, ptheta_efficiency_inb_cd, ptheta_bins_inb_fd, ptheta_bins_inb_cd)
                df_sim_nobkgmerging_inb   = assign_Ppefficiency(df_sim_nobkgmerging_inb, pp_efficiency_inb_sector_1, pp_efficiency_inb_sector_2, pp_efficiency_inb_sector_3, pp_efficiency_inb_sector_4, pp_efficiency_inb_sector_5, pp_efficiency_inb_sector_6, pp_efficiency_inb_cd, pp_bins_inb_fd, pp_bins_inb_cd)
                df_sim_nobkgmerging_inb   = assign_Pphiefficiency(df_sim_nobkgmerging_inb, pphi_efficiency_inb_sector_1, pphi_efficiency_inb_sector_2, pphi_efficiency_inb_sector_3, pphi_efficiency_inb_sector_4, pphi_efficiency_inb_sector_5, pphi_efficiency_inb_sector_6, pphi_efficiency_inb_cd, pphi_bins_inb_fd, pphi_bins_inb_cd)
            except:
                df_sim_nobkgmerging_inb      = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})
            try:
                df_sim_bkgmerging_inb     =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/dvcs_km15/excl_level_1/restructured_14_bkgmerging/{}.pkl".format(integrated_binnum))
                df_sim_bkgmerging_inb     = assign_Pthetaefficiency(df_sim_bkgmerging_inb, ptheta_efficiency_inb_sector_1, ptheta_efficiency_inb_sector_2, ptheta_efficiency_inb_sector_3, ptheta_efficiency_inb_sector_4, ptheta_efficiency_inb_sector_5, ptheta_efficiency_inb_sector_6, ptheta_efficiency_inb_cd, ptheta_bins_inb_fd, ptheta_bins_inb_cd)
                df_sim_bkgmerging_inb     = assign_Ppefficiency(df_sim_bkgmerging_inb, pp_efficiency_inb_sector_1, pp_efficiency_inb_sector_2, pp_efficiency_inb_sector_3, pp_efficiency_inb_sector_4, pp_efficiency_inb_sector_5, pp_efficiency_inb_sector_6, pp_efficiency_inb_cd, pp_bins_inb_fd, pp_bins_inb_cd)
                df_sim_bkgmerging_inb     = assign_Pphiefficiency(df_sim_bkgmerging_inb, pphi_efficiency_inb_sector_1, pphi_efficiency_inb_sector_2, pphi_efficiency_inb_sector_3, pphi_efficiency_inb_sector_4, pphi_efficiency_inb_sector_5, pphi_efficiency_inb_sector_6, pphi_efficiency_inb_cd, pphi_bins_inb_fd, pphi_bins_inb_cd)
            except:
                df_sim_bkgmerging_inb      = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})
            # df_sim_nobkgmerging_inb.loc[:, "efficiency"] = 0

            # bkg_to_nobkg    = np.sum(df_sim_bkgmerging_inb.weights * df_sim_bkgmerging_inb.Pthetaefficiency * df_sim_bkgmerging_inb.Ppefficiency * df_sim_bkgmerging_inb.Pphiefficiency)/np.sum(df_sim_nobkgmerging_inb.weights * df_sim_nobkgmerging_inb.Pthetaefficiency * df_sim_nobkgmerging_inb.Ppefficiency * df_sim_nobkgmerging_inb.Pphiefficiency)
            # eff_bkg_merging = ( 1 + effective_current_inb/45 * ( -1 + bkg_to_nobkg))

            # df_sim_nobkgmerging_inb.loc[:, "efficiency"] = eff_bkg_merging
            # epg_exp_inb_this_integrated_bin           = np.sum(1/df_exp_epg_inb.efficiency)
            # pi0_exp_inb_this_integrated_bin           = np.sum(1/df_exp_pi0_inb.efficiency)
            # pi0_sim_inb_this_integrated_bin           = np.sum(df_sim_pi0_inb.Pthetaefficiency * df_sim_pi0_inb.Ppefficiency * df_sim_pi0_inb.Pphiefficiency)
            # bkg_sim_inb_this_integrated_bin           = np.sum(df_sim_pi0_1gamma_inb.Pthetaefficiency * df_sim_pi0_1gamma_inb.Ppefficiency * df_sim_pi0_1gamma_inb.Pphiefficiency)
            # bkg_exp_inb_this_integrated_bin           = bkg_sim_inb_this_integrated_bin  * pi0_exp_inb_this_integrated_bin/pi0_sim_inb_this_integrated_bin
            # contamination_ratio                       = np.minimum(bkg_exp_inb_this_integrated_bin , epg_exp_inb_this_integrated_bin )/epg_exp_inb_this_integrated_bin 

            # df_exp_epg_inb.loc[:, "contamination"]   = contamination_ratio

            for i, this_rebinned_phi_binnum in enumerate(this_rebinned_phi_binnums):
                df_exp_epg_inb_this_bin = df_exp_epg_inb.loc[(df_exp_epg_inb.integrated_binnum == integrated_binnum) & (df_exp_epg_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                # if len(df_exp_epg_inb_this_bin) < 10:
                if df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == this_rebinned_phi_binnum), "active_bin_inb"].sum() == 0:
                    # print("{} inactive bin between {} and {}".format(integrated_binnum, this_rebinned_phi_binnum, this_rebinned_phi_binnum+this_rebinned_phi_widths[i]))
                    continue
                else:
                    df_exp_pi0_inb_this_bin = df_exp_pi0_inb.loc[(df_exp_pi0_inb.integrated_binnum == integrated_binnum) & (df_exp_pi0_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_pi0_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_pi0_inb_this_bin = df_sim_pi0_inb.loc[(df_sim_pi0_inb.integrated_binnum == integrated_binnum) & (df_sim_pi0_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_pi0_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_pi0_1gamma_inb_this_bin = df_sim_pi0_1gamma_inb.loc[(df_sim_pi0_1gamma_inb.integrated_binnum == integrated_binnum) & (df_sim_pi0_1gamma_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_pi0_1gamma_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_nobkgmerging_inb_this_bin = df_sim_nobkgmerging_inb.loc[(df_sim_nobkgmerging_inb.integrated_binnum == integrated_binnum) & (df_sim_nobkgmerging_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_nobkgmerging_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_bkgmerging_inb_this_bin = df_sim_bkgmerging_inb.loc[(df_sim_bkgmerging_inb.integrated_binnum == integrated_binnum) & (df_sim_bkgmerging_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_bkgmerging_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]

                    bkg_to_nobkg    = np.sum(df_sim_bkgmerging_inb_this_bin.weights * df_sim_bkgmerging_inb_this_bin.Pthetaefficiency * df_sim_bkgmerging_inb_this_bin.Ppefficiency * df_sim_bkgmerging_inb_this_bin.Pphiefficiency)/np.sum(df_sim_nobkgmerging_inb_this_bin.weights * df_sim_nobkgmerging_inb_this_bin.Pthetaefficiency * df_sim_nobkgmerging_inb_this_bin.Ppefficiency * df_sim_nobkgmerging_inb_this_bin.Pphiefficiency)
                    eff_bkg_merging = ( 1 + effective_current_inb/45 * ( -1 + bkg_to_nobkg))
                    
                    epg_exp_inb_this_bin               = np.sum(1/df_exp_epg_inb_this_bin.efficiency)
                    try:
                        pi0_exp_inb_this_integrated_bin    = np.sum(1/df_exp_pi0_inb.efficiency)
                    except:
                        pi0_exp_inb_this_integrated_bin    = 0
                    try:
                        pi0_sim_inb_this_integrated_bin    = np.sum(df_sim_pi0_inb.Pthetaefficiency)
                    except:
                        pi0_sim_inb_this_integrated_bin    = 0
                    try:
                        bkg_sim_inb_this_bin               = np.sum(df_sim_pi0_1gamma_inb_this_bin.Pthetaefficiency)
                    except:
                        bkg_sim_inb_this_bin               = 0
                    try:
                        bkg_exp_inb_this_bin               = bkg_sim_inb_this_bin * pi0_exp_inb_this_integrated_bin/pi0_sim_inb_this_integrated_bin
                    except:
                        bkg_exp_inb_this_bin               = 0
                    contamination_ratio                = np.minimum(bkg_exp_inb_this_bin, epg_exp_inb_this_bin)/epg_exp_inb_this_bin

                    df_exp_epg_inb.loc[(df_exp_epg_inb.integrated_binnum == integrated_binnum) & (df_exp_epg_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "contamination"]   = contamination_ratio
                    df_exp_epg_inb.loc[(df_exp_epg_inb.integrated_binnum == integrated_binnum) & (df_exp_epg_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "phi_rebinnum"]    = this_rebinned_phi_binnum
                    df_exp_epg_inb.loc[(df_exp_epg_inb.integrated_binnum == integrated_binnum) & (df_exp_epg_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "phi_rebinwidth"]  = this_rebinned_phi_widths[i]
                    df_sim_nobkgmerging_inb.loc[(df_sim_nobkgmerging_inb.integrated_binnum == integrated_binnum) & (df_sim_nobkgmerging_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_nobkgmerging_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "efficiency"] = eff_bkg_merging

            df_exp_epg_inb.loc[:, "signal"]          = (1 - df_exp_epg_inb.contamination)/df_exp_epg_inb.efficiency
            df_sim_nobkgmerging_inb.loc[:, "signal"] = df_sim_nobkgmerging_inb.weights * df_sim_nobkgmerging_inb.efficiency * df_sim_nobkgmerging_inb.Pthetaefficiency  * df_sim_nobkgmerging_inb.Ppefficiency * df_sim_nobkgmerging_inb.Pphiefficiency
            
            df_exp_epg_inbs.append(df_exp_epg_inb)
            df_sim_dvcs_inbs.append(df_sim_nobkgmerging_inb)

        df_exp_epg_inbs  = pd.concat(df_exp_epg_inbs)
        df_sim_dvcs_inbs = pd.concat(df_sim_dvcs_inbs)

        # FD sector 1
        ptheta_exp_sector_1, _ = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 1, "Ptheta"], weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 1, "signal"], bins = ptheta_bins_inb_fd)
        ptheta_sim_sector_1, _ = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 1, "Ptheta"], weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 1, "signal"], bins = ptheta_bins_inb_fd)
        ptheta_efficiency_inb_sector_1 = ptheta_efficiency_inb_sector_1*divideHist(ptheta_exp_sector_1, ptheta_sim_sector_1)

        # FD sector 2
        ptheta_exp_sector_2, _ = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 2, "Ptheta"], weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 2, "signal"], bins = ptheta_bins_inb_fd)
        ptheta_sim_sector_2, _ = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 2, "Ptheta"], weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 2, "signal"], bins = ptheta_bins_inb_fd)
        ptheta_efficiency_inb_sector_2 = ptheta_efficiency_inb_sector_2*divideHist(ptheta_exp_sector_2, ptheta_sim_sector_2)

        # FD sector 3
        ptheta_exp_sector_3, _ = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 3, "Ptheta"], weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 3, "signal"], bins = ptheta_bins_inb_fd)
        ptheta_sim_sector_3, _ = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 3, "Ptheta"], weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 3, "signal"], bins = ptheta_bins_inb_fd)
        ptheta_efficiency_inb_sector_3 = ptheta_efficiency_inb_sector_3*divideHist(ptheta_exp_sector_3, ptheta_sim_sector_3)

        # FD sector 4
        ptheta_exp_sector_4, _ = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 4, "Ptheta"], weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 4, "signal"], bins = ptheta_bins_inb_fd)
        ptheta_sim_sector_4, _ = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 4, "Ptheta"], weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 4, "signal"], bins = ptheta_bins_inb_fd)
        ptheta_efficiency_inb_sector_4 = ptheta_efficiency_inb_sector_4*divideHist(ptheta_exp_sector_4, ptheta_sim_sector_4)

        # FD sector 5
        ptheta_exp_sector_5, _ = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 5, "Ptheta"], weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 5, "signal"], bins = ptheta_bins_inb_fd)
        ptheta_sim_sector_5, _ = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 5, "Ptheta"], weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 5, "signal"], bins = ptheta_bins_inb_fd)
        ptheta_efficiency_inb_sector_5 = ptheta_efficiency_inb_sector_5*divideHist(ptheta_exp_sector_5, ptheta_sim_sector_5)

        # FD sector 6
        ptheta_exp_sector_6, _ = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 6, "Ptheta"], weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 6, "signal"], bins = ptheta_bins_inb_fd)
        ptheta_sim_sector_6, _ = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 6, "Ptheta"], weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 6, "signal"], bins = ptheta_bins_inb_fd)
        ptheta_efficiency_inb_sector_6 = ptheta_efficiency_inb_sector_6*divideHist(ptheta_exp_sector_6, ptheta_sim_sector_6)

        # CD
        ptheta_exp_cd, _ = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector > 7, "Ptheta"], weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector > 7, "signal"], bins = ptheta_bins_inb_cd)
        ptheta_sim_cd, _ = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector > 7, "Ptheta"], weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector > 7, "signal"], bins = ptheta_bins_inb_cd)
        ptheta_efficiency_inb_cd = ptheta_efficiency_inb_cd*divideHist(ptheta_exp_cd, ptheta_sim_cd)

        ptheta_efficiency_inb_sector_1[ptheta_efficiency_inb_sector_1<0.1] = 0.1
        ptheta_efficiency_inb_sector_2[ptheta_efficiency_inb_sector_2<0.1] = 0.1
        ptheta_efficiency_inb_sector_3[ptheta_efficiency_inb_sector_3<0.1] = 0.1
        ptheta_efficiency_inb_sector_4[ptheta_efficiency_inb_sector_4<0.1] = 0.1
        ptheta_efficiency_inb_sector_5[ptheta_efficiency_inb_sector_5<0.1] = 0.1
        ptheta_efficiency_inb_sector_6[ptheta_efficiency_inb_sector_6<0.1] = 0.1
        ptheta_efficiency_inb_cd      [ptheta_efficiency_inb_cd      <0.1] = 0.1

        ptheta_efficiency_inb_sector_1[ptheta_efficiency_inb_sector_1>1.5] = 1.5
        ptheta_efficiency_inb_sector_2[ptheta_efficiency_inb_sector_2>1.5] = 1.5
        ptheta_efficiency_inb_sector_3[ptheta_efficiency_inb_sector_3>1.5] = 1.5
        ptheta_efficiency_inb_sector_4[ptheta_efficiency_inb_sector_4>1.5] = 1.5
        ptheta_efficiency_inb_sector_5[ptheta_efficiency_inb_sector_5>1.5] = 1.5
        ptheta_efficiency_inb_sector_6[ptheta_efficiency_inb_sector_6>1.5] = 1.5
        ptheta_efficiency_inb_cd      [ptheta_efficiency_inb_cd      >1.5] = 1.5

    for pp_trial in range(3):
        df_exp_epg_inbs        = []
        df_sim_dvcs_inbs       = []
        print(pp_trial, "pp_trial", "inb sector 1", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_1], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_1], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_1])
        print(pp_trial, "pp_trial", "inb sector 2", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_2], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_2], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_2])
        print(pp_trial, "pp_trial", "inb sector 3", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_3], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_3], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_3])
        print(pp_trial, "pp_trial", "inb sector 4", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_4], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_4], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_4])
        print(pp_trial, "pp_trial", "inb sector 5", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_5], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_5], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_5])
        print(pp_trial, "pp_trial", "inb sector 6", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_6], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_6], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_6])
        print(pp_trial, "pp_trial", "inb CD"      , ["{:.3f}".format(i) for i in ptheta_efficiency_inb_cd      ], ["{:.3f}".format(i) for i in pp_efficiency_inb_cd      ], ["{:.3f}".format(i) for i in pphi_efficiency_inb_cd      ])

        for integrated_binnum in range(1, 147+1):
            this_rebinned_phi_binnums = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"].to_numpy()
            this_rebinned_phi_widths  = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"].to_numpy()

            try:
                df_exp_epg_inb           =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_inb/dvcs/excl_level_2/restructured_7_nominal/{}.pkl".format(integrated_binnum))
                df_exp_epg_inb.loc[df_exp_epg_inb.efficiency == 0, "efficiency"] =  1
                df_exp_epg_inb.loc[:, "contamination"] = 1
            except:
                # print("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_inb/dvcs/excl_level_2/restructured_7_nominal/{}.pkl does not exist.".format(integrated_binnum))
                continue

            # if len(df_exp_epg_inb)<50:
            if df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum), "active_bin_inb"].sum() == 0:
                continue

            try:
                df_exp_pi0_inb            =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_inb/pi0/excl_level_2/restructured_7_nominal/{}.pkl".format(integrated_binnum))
            except:
                df_exp_pi0_inb            =  pd.DataFrame({var: [] for var in ['Pp', 'Ptheta', 'Pphi', 'efficiency']})
            # df_exp_pi0_inb.loc[df_exp_pi0_inb.efficiency == 0, "efficiency"]  = 1
            try:
                df_sim_pi0_inb            =  pd.concat([pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/pi0_2gamma/excl_level_1/restructured_7_nominal/{}/{}.pkl".format(chunk, integrated_binnum)) for chunk in chunk_inb])
                df_sim_pi0_inb            = assign_Pthetaefficiency(df_sim_pi0_inb, ptheta_efficiency_inb_sector_1, ptheta_efficiency_inb_sector_2, ptheta_efficiency_inb_sector_3, ptheta_efficiency_inb_sector_4, ptheta_efficiency_inb_sector_5, ptheta_efficiency_inb_sector_6, ptheta_efficiency_inb_cd, ptheta_bins_inb_fd, ptheta_bins_inb_cd)
                df_sim_pi0_inb            = assign_Ppefficiency(df_sim_pi0_inb, pp_efficiency_inb_sector_1, pp_efficiency_inb_sector_2, pp_efficiency_inb_sector_3, pp_efficiency_inb_sector_4, pp_efficiency_inb_sector_5, pp_efficiency_inb_sector_6, pp_efficiency_inb_cd, pp_bins_inb_fd, pp_bins_inb_cd)
                df_sim_pi0_inb            = assign_Pphiefficiency(df_sim_pi0_inb, pphi_efficiency_inb_sector_1, pphi_efficiency_inb_sector_2, pphi_efficiency_inb_sector_3, pphi_efficiency_inb_sector_4, pphi_efficiency_inb_sector_5, pphi_efficiency_inb_sector_6, pphi_efficiency_inb_cd, pphi_bins_inb_fd, pphi_bins_inb_cd)
            except:
                df_sim_pi0_inb            = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})

            try:
                df_sim_pi0_1gamma_inb     =  pd.concat([pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/pi0_1gamma/excl_level_1/restructured_7_nominal/{}/{}.pkl".format(chunk, integrated_binnum)) for chunk in chunk_inb])
                df_sim_pi0_1gamma_inb     = assign_Pthetaefficiency(df_sim_pi0_1gamma_inb, ptheta_efficiency_inb_sector_1, ptheta_efficiency_inb_sector_2, ptheta_efficiency_inb_sector_3, ptheta_efficiency_inb_sector_4, ptheta_efficiency_inb_sector_5, ptheta_efficiency_inb_sector_6, ptheta_efficiency_inb_cd, ptheta_bins_inb_fd, ptheta_bins_inb_cd)
                df_sim_pi0_1gamma_inb     = assign_Ppefficiency(df_sim_pi0_1gamma_inb, pp_efficiency_inb_sector_1, pp_efficiency_inb_sector_2, pp_efficiency_inb_sector_3, pp_efficiency_inb_sector_4, pp_efficiency_inb_sector_5, pp_efficiency_inb_sector_6, pp_efficiency_inb_cd, pp_bins_inb_fd, pp_bins_inb_cd)
                df_sim_pi0_1gamma_inb     = assign_Pphiefficiency(df_sim_pi0_1gamma_inb, pphi_efficiency_inb_sector_1, pphi_efficiency_inb_sector_2, pphi_efficiency_inb_sector_3, pphi_efficiency_inb_sector_4, pphi_efficiency_inb_sector_5, pphi_efficiency_inb_sector_6, pphi_efficiency_inb_cd, pphi_bins_inb_fd, pphi_bins_inb_cd)
            except:
                df_sim_pi0_1gamma_inb      = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})
            try:
                df_sim_nobkgmerging_inb   =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/dvcs_km15/excl_level_1/restructured_7_nominal/{}.pkl".format(integrated_binnum))
                df_sim_nobkgmerging_inb   = assign_Pthetaefficiency(df_sim_nobkgmerging_inb, ptheta_efficiency_inb_sector_1, ptheta_efficiency_inb_sector_2, ptheta_efficiency_inb_sector_3, ptheta_efficiency_inb_sector_4, ptheta_efficiency_inb_sector_5, ptheta_efficiency_inb_sector_6, ptheta_efficiency_inb_cd, ptheta_bins_inb_fd, ptheta_bins_inb_cd)
                df_sim_nobkgmerging_inb   = assign_Ppefficiency(df_sim_nobkgmerging_inb, pp_efficiency_inb_sector_1, pp_efficiency_inb_sector_2, pp_efficiency_inb_sector_3, pp_efficiency_inb_sector_4, pp_efficiency_inb_sector_5, pp_efficiency_inb_sector_6, pp_efficiency_inb_cd, pp_bins_inb_fd, pp_bins_inb_cd)
                df_sim_nobkgmerging_inb   = assign_Pphiefficiency(df_sim_nobkgmerging_inb, pphi_efficiency_inb_sector_1, pphi_efficiency_inb_sector_2, pphi_efficiency_inb_sector_3, pphi_efficiency_inb_sector_4, pphi_efficiency_inb_sector_5, pphi_efficiency_inb_sector_6, pphi_efficiency_inb_cd, pphi_bins_inb_fd, pphi_bins_inb_cd)
            except:
                df_sim_nobkgmerging_inb      = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})
            try:
                df_sim_bkgmerging_inb     =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/dvcs_km15/excl_level_1/restructured_14_bkgmerging/{}.pkl".format(integrated_binnum))
                df_sim_bkgmerging_inb     = assign_Pthetaefficiency(df_sim_bkgmerging_inb, ptheta_efficiency_inb_sector_1, ptheta_efficiency_inb_sector_2, ptheta_efficiency_inb_sector_3, ptheta_efficiency_inb_sector_4, ptheta_efficiency_inb_sector_5, ptheta_efficiency_inb_sector_6, ptheta_efficiency_inb_cd, ptheta_bins_inb_fd, ptheta_bins_inb_cd)
                df_sim_bkgmerging_inb     = assign_Ppefficiency(df_sim_bkgmerging_inb, pp_efficiency_inb_sector_1, pp_efficiency_inb_sector_2, pp_efficiency_inb_sector_3, pp_efficiency_inb_sector_4, pp_efficiency_inb_sector_5, pp_efficiency_inb_sector_6, pp_efficiency_inb_cd, pp_bins_inb_fd, pp_bins_inb_cd)
                df_sim_bkgmerging_inb     = assign_Pphiefficiency(df_sim_bkgmerging_inb, pphi_efficiency_inb_sector_1, pphi_efficiency_inb_sector_2, pphi_efficiency_inb_sector_3, pphi_efficiency_inb_sector_4, pphi_efficiency_inb_sector_5, pphi_efficiency_inb_sector_6, pphi_efficiency_inb_cd, pphi_bins_inb_fd, pphi_bins_inb_cd)
            except:
                df_sim_bkgmerging_inb      = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})
            # df_sim_nobkgmerging_inb.loc[:, "efficiency"] = 0

            # bkg_to_nobkg    = np.sum(df_sim_bkgmerging_inb.weights * df_sim_bkgmerging_inb.Pthetaefficiency * df_sim_bkgmerging_inb.Ppefficiency * df_sim_bkgmerging_inb.Pphiefficiency)/np.sum(df_sim_nobkgmerging_inb.weights * df_sim_nobkgmerging_inb.Pthetaefficiency * df_sim_nobkgmerging_inb.Ppefficiency * df_sim_nobkgmerging_inb.Pphiefficiency)
            # eff_bkg_merging = ( 1 + effective_current_inb/45 * ( -1 + bkg_to_nobkg))

            # df_sim_nobkgmerging_inb.loc[:, "efficiency"] = eff_bkg_merging
            # epg_exp_inb_this_integrated_bin           = np.sum(1/df_exp_epg_inb.efficiency)
            # pi0_exp_inb_this_integrated_bin           = np.sum(1/df_exp_pi0_inb.efficiency)
            # pi0_sim_inb_this_integrated_bin           = np.sum(df_sim_pi0_inb.Pthetaefficiency * df_sim_pi0_inb.Ppefficiency * df_sim_pi0_inb.Pphiefficiency)
            # bkg_sim_inb_this_integrated_bin           = np.sum(df_sim_pi0_1gamma_inb.Pthetaefficiency * df_sim_pi0_1gamma_inb.Ppefficiency * df_sim_pi0_1gamma_inb.Pphiefficiency)
            # bkg_exp_inb_this_integrated_bin           = bkg_sim_inb_this_integrated_bin  * pi0_exp_inb_this_integrated_bin/pi0_sim_inb_this_integrated_bin
            # contamination_ratio                       = np.minimum(bkg_exp_inb_this_integrated_bin , epg_exp_inb_this_integrated_bin )/epg_exp_inb_this_integrated_bin 

            # df_exp_epg_inb.loc[:, "contamination"]   = contamination_ratio

            for i, this_rebinned_phi_binnum in enumerate(this_rebinned_phi_binnums):
                df_exp_epg_inb_this_bin = df_exp_epg_inb.loc[(df_exp_epg_inb.integrated_binnum == integrated_binnum) & (df_exp_epg_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                # if len(df_exp_epg_inb_this_bin) < 10:
                if df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == this_rebinned_phi_binnum), "active_bin_inb"].sum() == 0:
                    # print("{} inactive bin between {} and {}".format(integrated_binnum, this_rebinned_phi_binnum, this_rebinned_phi_binnum+this_rebinned_phi_widths[i]))
                    continue
                else:
                    df_exp_pi0_inb_this_bin = df_exp_pi0_inb.loc[(df_exp_pi0_inb.integrated_binnum == integrated_binnum) & (df_exp_pi0_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_pi0_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_pi0_inb_this_bin = df_sim_pi0_inb.loc[(df_sim_pi0_inb.integrated_binnum == integrated_binnum) & (df_sim_pi0_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_pi0_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_pi0_1gamma_inb_this_bin = df_sim_pi0_1gamma_inb.loc[(df_sim_pi0_1gamma_inb.integrated_binnum == integrated_binnum) & (df_sim_pi0_1gamma_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_pi0_1gamma_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_nobkgmerging_inb_this_bin = df_sim_nobkgmerging_inb.loc[(df_sim_nobkgmerging_inb.integrated_binnum == integrated_binnum) & (df_sim_nobkgmerging_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_nobkgmerging_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_bkgmerging_inb_this_bin = df_sim_bkgmerging_inb.loc[(df_sim_bkgmerging_inb.integrated_binnum == integrated_binnum) & (df_sim_bkgmerging_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_bkgmerging_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]

                    bkg_to_nobkg    = np.sum(df_sim_bkgmerging_inb_this_bin.weights * df_sim_bkgmerging_inb_this_bin.Pthetaefficiency * df_sim_bkgmerging_inb_this_bin.Ppefficiency * df_sim_bkgmerging_inb_this_bin.Pphiefficiency)/np.sum(df_sim_nobkgmerging_inb_this_bin.weights * df_sim_nobkgmerging_inb_this_bin.Pthetaefficiency * df_sim_nobkgmerging_inb_this_bin.Ppefficiency * df_sim_nobkgmerging_inb_this_bin.Pphiefficiency)
                    eff_bkg_merging = ( 1 + effective_current_inb/45 * ( -1 + bkg_to_nobkg))
                    
                    epg_exp_inb_this_bin               = np.sum(1/df_exp_epg_inb_this_bin.efficiency)
                    try:
                        pi0_exp_inb_this_integrated_bin    = np.sum(1/df_exp_pi0_inb.efficiency)
                    except:
                        pi0_exp_inb_this_integrated_bin    = 0
                    try:
                        pi0_sim_inb_this_integrated_bin    = np.sum(df_sim_pi0_inb.Pthetaefficiency)
                    except:
                        pi0_sim_inb_this_integrated_bin    = 0
                    try:
                        bkg_sim_inb_this_bin               = np.sum(df_sim_pi0_1gamma_inb_this_bin.Pthetaefficiency)
                    except:
                        bkg_sim_inb_this_bin               = 0
                    try:
                        bkg_exp_inb_this_bin               = bkg_sim_inb_this_bin * pi0_exp_inb_this_integrated_bin/pi0_sim_inb_this_integrated_bin
                    except:
                        bkg_exp_inb_this_bin               = 0
                    contamination_ratio                = np.minimum(bkg_exp_inb_this_bin, epg_exp_inb_this_bin)/epg_exp_inb_this_bin

                    df_exp_epg_inb.loc[(df_exp_epg_inb.integrated_binnum == integrated_binnum) & (df_exp_epg_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "contamination"]   = contamination_ratio
                    df_exp_epg_inb.loc[(df_exp_epg_inb.integrated_binnum == integrated_binnum) & (df_exp_epg_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "phi_rebinnum"]    = this_rebinned_phi_binnum
                    df_exp_epg_inb.loc[(df_exp_epg_inb.integrated_binnum == integrated_binnum) & (df_exp_epg_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "phi_rebinwidth"]  = this_rebinned_phi_widths[i]
                    df_sim_nobkgmerging_inb.loc[(df_sim_nobkgmerging_inb.integrated_binnum == integrated_binnum) & (df_sim_nobkgmerging_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_nobkgmerging_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "efficiency"] = eff_bkg_merging

            df_exp_epg_inb.loc[:, "signal"]          = (1 - df_exp_epg_inb.contamination)/df_exp_epg_inb.efficiency
            df_sim_nobkgmerging_inb.loc[:, "signal"] = df_sim_nobkgmerging_inb.weights * df_sim_nobkgmerging_inb.efficiency * df_sim_nobkgmerging_inb.Pthetaefficiency  * df_sim_nobkgmerging_inb.Ppefficiency * df_sim_nobkgmerging_inb.Pphiefficiency
            
            df_exp_epg_inbs.append(df_exp_epg_inb)
            df_sim_dvcs_inbs.append(df_sim_nobkgmerging_inb)

        df_exp_epg_inbs  = pd.concat(df_exp_epg_inbs)
        df_sim_dvcs_inbs = pd.concat(df_sim_dvcs_inbs)

        # FD sector 1
        pp_exp_sector_1, _ = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 1, "Pp"], weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 1, "signal"], bins = pp_bins_inb_fd)
        pp_sim_sector_1, _ = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 1, "Pp"], weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 1, "signal"], bins = pp_bins_inb_fd)
        pp_efficiency_inb_sector_1 = pp_efficiency_inb_sector_1*divideHist(pp_exp_sector_1, pp_sim_sector_1)

        # FD sector 2
        pp_exp_sector_2, _ = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 2, "Pp"], weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 2, "signal"], bins = pp_bins_inb_fd)
        pp_sim_sector_2, _ = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 2, "Pp"], weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 2, "signal"], bins = pp_bins_inb_fd)
        pp_efficiency_inb_sector_2 = pp_efficiency_inb_sector_2*divideHist(pp_exp_sector_2, pp_sim_sector_2)

        # FD sector 3
        pp_exp_sector_3, _ = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 3, "Pp"], weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 3, "signal"], bins = pp_bins_inb_fd)
        pp_sim_sector_3, _ = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 3, "Pp"], weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 3, "signal"], bins = pp_bins_inb_fd)
        pp_efficiency_inb_sector_3 = pp_efficiency_inb_sector_3*divideHist(pp_exp_sector_3, pp_sim_sector_3)

        # FD sector 4
        pp_exp_sector_4, _ = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 4, "Pp"], weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 4, "signal"], bins = pp_bins_inb_fd)
        pp_sim_sector_4, _ = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 4, "Pp"], weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 4, "signal"], bins = pp_bins_inb_fd)
        pp_efficiency_inb_sector_4 = pp_efficiency_inb_sector_4*divideHist(pp_exp_sector_4, pp_sim_sector_4)

        # FD sector 5
        pp_exp_sector_5, _ = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 5, "Pp"], weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 5, "signal"], bins = pp_bins_inb_fd)
        pp_sim_sector_5, _ = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 5, "Pp"], weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 5, "signal"], bins = pp_bins_inb_fd)
        pp_efficiency_inb_sector_5 = pp_efficiency_inb_sector_5*divideHist(pp_exp_sector_5, pp_sim_sector_5)

        # FD sector 6
        pp_exp_sector_6, _ = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 6, "Pp"], weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 6, "signal"], bins = pp_bins_inb_fd)
        pp_sim_sector_6, _ = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 6, "Pp"], weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 6, "signal"], bins = pp_bins_inb_fd)
        pp_efficiency_inb_sector_6 = pp_efficiency_inb_sector_6*divideHist(pp_exp_sector_6, pp_sim_sector_6)

        # CD
        pp_exp_cd, _ = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector > 7, "Pp"], weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector > 7, "signal"], bins = pp_bins_inb_cd)
        pp_sim_cd, _ = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector > 7, "Pp"], weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector > 7, "signal"], bins = pp_bins_inb_cd)
        pp_efficiency_inb_cd = pp_efficiency_inb_cd*divideHist(pp_exp_cd, pp_sim_cd)

        pp_efficiency_inb_sector_1[pp_efficiency_inb_sector_1<0.5] = 0.5
        pp_efficiency_inb_sector_2[pp_efficiency_inb_sector_2<0.5] = 0.5
        pp_efficiency_inb_sector_3[pp_efficiency_inb_sector_3<0.5] = 0.5
        pp_efficiency_inb_sector_4[pp_efficiency_inb_sector_4<0.5] = 0.5
        pp_efficiency_inb_sector_5[pp_efficiency_inb_sector_5<0.5] = 0.5
        pp_efficiency_inb_sector_6[pp_efficiency_inb_sector_6<0.5] = 0.5
        pp_efficiency_inb_cd      [pp_efficiency_inb_cd      <0.5] = 0.5

        pp_efficiency_inb_sector_1[pp_efficiency_inb_sector_1>1.5] = 1.5
        pp_efficiency_inb_sector_2[pp_efficiency_inb_sector_2>1.5] = 1.5
        pp_efficiency_inb_sector_3[pp_efficiency_inb_sector_3>1.5] = 1.5
        pp_efficiency_inb_sector_4[pp_efficiency_inb_sector_4>1.5] = 1.5
        pp_efficiency_inb_sector_5[pp_efficiency_inb_sector_5>1.5] = 1.5
        pp_efficiency_inb_sector_6[pp_efficiency_inb_sector_6>1.5] = 1.5
        pp_efficiency_inb_cd      [pp_efficiency_inb_cd      >1.5] = 1.5

    for pphi_trial in range(3):
        df_exp_epg_inbs        = []
        df_sim_dvcs_inbs       = []
        print(pphi_trial, "pphi_trial", "inb sector 1", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_1], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_1], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_1])
        print(pphi_trial, "pphi_trial", "inb sector 2", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_2], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_2], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_2])
        print(pphi_trial, "pphi_trial", "inb sector 3", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_3], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_3], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_3])
        print(pphi_trial, "pphi_trial", "inb sector 4", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_4], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_4], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_4])
        print(pphi_trial, "pphi_trial", "inb sector 5", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_5], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_5], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_5])
        print(pphi_trial, "pphi_trial", "inb sector 6", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_6], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_6], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_6])
        print(pphi_trial, "pphi_trial", "inb CD"      , ["{:.3f}".format(i) for i in ptheta_efficiency_inb_cd      ], ["{:.3f}".format(i) for i in pp_efficiency_inb_cd      ], ["{:.3f}".format(i) for i in pphi_efficiency_inb_cd      ])

        for integrated_binnum in range(1, 147+1):
            this_rebinned_phi_binnums = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"].to_numpy()
            this_rebinned_phi_widths  = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"].to_numpy()

            try:
                df_exp_epg_inb           =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_inb/dvcs/excl_level_2/restructured_7_nominal/{}.pkl".format(integrated_binnum))
                df_exp_epg_inb.loc[df_exp_epg_inb.efficiency == 0, "efficiency"] =  1
                df_exp_epg_inb.loc[:, "contamination"] = 1
            except:
                # print("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_inb/dvcs/excl_level_2/restructured_7_nominal/{}.pkl does not exist.".format(integrated_binnum))
                continue

            # if len(df_exp_epg_inb)<50:
            if df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum), "active_bin_inb"].sum() == 0:
                continue

            try:
                df_exp_pi0_inb            =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_inb/pi0/excl_level_2/restructured_7_nominal/{}.pkl".format(integrated_binnum))
            except:
                df_exp_pi0_inb            =  pd.DataFrame({var: [] for var in ['Pp', 'Ptheta', 'Pphi', 'efficiency']})
            # df_exp_pi0_inb.loc[df_exp_pi0_inb.efficiency == 0, "efficiency"]  = 1
            try:
                df_sim_pi0_inb            =  pd.concat([pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/pi0_2gamma/excl_level_1/restructured_7_nominal/{}/{}.pkl".format(chunk, integrated_binnum)) for chunk in chunk_inb])
                df_sim_pi0_inb            = assign_Pthetaefficiency(df_sim_pi0_inb, ptheta_efficiency_inb_sector_1, ptheta_efficiency_inb_sector_2, ptheta_efficiency_inb_sector_3, ptheta_efficiency_inb_sector_4, ptheta_efficiency_inb_sector_5, ptheta_efficiency_inb_sector_6, ptheta_efficiency_inb_cd, ptheta_bins_inb_fd, ptheta_bins_inb_cd)
                df_sim_pi0_inb            = assign_Ppefficiency(df_sim_pi0_inb, pp_efficiency_inb_sector_1, pp_efficiency_inb_sector_2, pp_efficiency_inb_sector_3, pp_efficiency_inb_sector_4, pp_efficiency_inb_sector_5, pp_efficiency_inb_sector_6, pp_efficiency_inb_cd, pp_bins_inb_fd, pp_bins_inb_cd)
                df_sim_pi0_inb            = assign_Pphiefficiency(df_sim_pi0_inb, pphi_efficiency_inb_sector_1, pphi_efficiency_inb_sector_2, pphi_efficiency_inb_sector_3, pphi_efficiency_inb_sector_4, pphi_efficiency_inb_sector_5, pphi_efficiency_inb_sector_6, pphi_efficiency_inb_cd, pphi_bins_inb_fd, pphi_bins_inb_cd)
            except:
                df_sim_pi0_inb            = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})

            try:
                df_sim_pi0_1gamma_inb     =  pd.concat([pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/pi0_1gamma/excl_level_1/restructured_7_nominal/{}/{}.pkl".format(chunk, integrated_binnum)) for chunk in chunk_inb])
                df_sim_pi0_1gamma_inb     = assign_Pthetaefficiency(df_sim_pi0_1gamma_inb, ptheta_efficiency_inb_sector_1, ptheta_efficiency_inb_sector_2, ptheta_efficiency_inb_sector_3, ptheta_efficiency_inb_sector_4, ptheta_efficiency_inb_sector_5, ptheta_efficiency_inb_sector_6, ptheta_efficiency_inb_cd, ptheta_bins_inb_fd, ptheta_bins_inb_cd)
                df_sim_pi0_1gamma_inb     = assign_Ppefficiency(df_sim_pi0_1gamma_inb, pp_efficiency_inb_sector_1, pp_efficiency_inb_sector_2, pp_efficiency_inb_sector_3, pp_efficiency_inb_sector_4, pp_efficiency_inb_sector_5, pp_efficiency_inb_sector_6, pp_efficiency_inb_cd, pp_bins_inb_fd, pp_bins_inb_cd)
                df_sim_pi0_1gamma_inb     = assign_Pphiefficiency(df_sim_pi0_1gamma_inb, pphi_efficiency_inb_sector_1, pphi_efficiency_inb_sector_2, pphi_efficiency_inb_sector_3, pphi_efficiency_inb_sector_4, pphi_efficiency_inb_sector_5, pphi_efficiency_inb_sector_6, pphi_efficiency_inb_cd, pphi_bins_inb_fd, pphi_bins_inb_cd)
            except:
                df_sim_pi0_1gamma_inb      = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})
            try:
                df_sim_nobkgmerging_inb   =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/dvcs_km15/excl_level_1/restructured_7_nominal/{}.pkl".format(integrated_binnum))
                df_sim_nobkgmerging_inb   = assign_Pthetaefficiency(df_sim_nobkgmerging_inb, ptheta_efficiency_inb_sector_1, ptheta_efficiency_inb_sector_2, ptheta_efficiency_inb_sector_3, ptheta_efficiency_inb_sector_4, ptheta_efficiency_inb_sector_5, ptheta_efficiency_inb_sector_6, ptheta_efficiency_inb_cd, ptheta_bins_inb_fd, ptheta_bins_inb_cd)
                df_sim_nobkgmerging_inb   = assign_Ppefficiency(df_sim_nobkgmerging_inb, pp_efficiency_inb_sector_1, pp_efficiency_inb_sector_2, pp_efficiency_inb_sector_3, pp_efficiency_inb_sector_4, pp_efficiency_inb_sector_5, pp_efficiency_inb_sector_6, pp_efficiency_inb_cd, pp_bins_inb_fd, pp_bins_inb_cd)
                df_sim_nobkgmerging_inb   = assign_Pphiefficiency(df_sim_nobkgmerging_inb, pphi_efficiency_inb_sector_1, pphi_efficiency_inb_sector_2, pphi_efficiency_inb_sector_3, pphi_efficiency_inb_sector_4, pphi_efficiency_inb_sector_5, pphi_efficiency_inb_sector_6, pphi_efficiency_inb_cd, pphi_bins_inb_fd, pphi_bins_inb_cd)
            except:
                df_sim_nobkgmerging_inb      = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})
            try:
                df_sim_bkgmerging_inb     =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/dvcs_km15/excl_level_1/restructured_14_bkgmerging/{}.pkl".format(integrated_binnum))
                df_sim_bkgmerging_inb     = assign_Pthetaefficiency(df_sim_bkgmerging_inb, ptheta_efficiency_inb_sector_1, ptheta_efficiency_inb_sector_2, ptheta_efficiency_inb_sector_3, ptheta_efficiency_inb_sector_4, ptheta_efficiency_inb_sector_5, ptheta_efficiency_inb_sector_6, ptheta_efficiency_inb_cd, ptheta_bins_inb_fd, ptheta_bins_inb_cd)
                df_sim_bkgmerging_inb     = assign_Ppefficiency(df_sim_bkgmerging_inb, pp_efficiency_inb_sector_1, pp_efficiency_inb_sector_2, pp_efficiency_inb_sector_3, pp_efficiency_inb_sector_4, pp_efficiency_inb_sector_5, pp_efficiency_inb_sector_6, pp_efficiency_inb_cd, pp_bins_inb_fd, pp_bins_inb_cd)
                df_sim_bkgmerging_inb     = assign_Pphiefficiency(df_sim_bkgmerging_inb, pphi_efficiency_inb_sector_1, pphi_efficiency_inb_sector_2, pphi_efficiency_inb_sector_3, pphi_efficiency_inb_sector_4, pphi_efficiency_inb_sector_5, pphi_efficiency_inb_sector_6, pphi_efficiency_inb_cd, pphi_bins_inb_fd, pphi_bins_inb_cd)
            except:
                df_sim_bkgmerging_inb      = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})
            # df_sim_nobkgmerging_inb.loc[:, "efficiency"] = 0

            # bkg_to_nobkg    = np.sum(df_sim_bkgmerging_inb.weights * df_sim_bkgmerging_inb.Pthetaefficiency * df_sim_bkgmerging_inb.Ppefficiency * df_sim_bkgmerging_inb.Pphiefficiency)/np.sum(df_sim_nobkgmerging_inb.weights * df_sim_nobkgmerging_inb.Pthetaefficiency * df_sim_nobkgmerging_inb.Ppefficiency * df_sim_nobkgmerging_inb.Pphiefficiency)
            # eff_bkg_merging = ( 1 + effective_current_inb/45 * ( -1 + bkg_to_nobkg))

            # df_sim_nobkgmerging_inb.loc[:, "efficiency"] = eff_bkg_merging
            # epg_exp_inb_this_integrated_bin           = np.sum(1/df_exp_epg_inb.efficiency)
            # pi0_exp_inb_this_integrated_bin           = np.sum(1/df_exp_pi0_inb.efficiency)
            # pi0_sim_inb_this_integrated_bin           = np.sum(df_sim_pi0_inb.Pthetaefficiency * df_sim_pi0_inb.Ppefficiency * df_sim_pi0_inb.Pphiefficiency)
            # bkg_sim_inb_this_integrated_bin           = np.sum(df_sim_pi0_1gamma_inb.Pthetaefficiency * df_sim_pi0_1gamma_inb.Ppefficiency * df_sim_pi0_1gamma_inb.Pphiefficiency)
            # bkg_exp_inb_this_integrated_bin           = bkg_sim_inb_this_integrated_bin  * pi0_exp_inb_this_integrated_bin/pi0_sim_inb_this_integrated_bin
            # contamination_ratio                       = np.minimum(bkg_exp_inb_this_integrated_bin , epg_exp_inb_this_integrated_bin )/epg_exp_inb_this_integrated_bin 

            # df_exp_epg_inb.loc[:, "contamination"]   = contamination_ratio

            for i, this_rebinned_phi_binnum in enumerate(this_rebinned_phi_binnums):
                df_exp_epg_inb_this_bin = df_exp_epg_inb.loc[(df_exp_epg_inb.integrated_binnum == integrated_binnum) & (df_exp_epg_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                # if len(df_exp_epg_inb_this_bin) < 10:
                if df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == this_rebinned_phi_binnum), "active_bin_inb"].sum() == 0:
                    # print("{} inactive bin between {} and {}".format(integrated_binnum, this_rebinned_phi_binnum, this_rebinned_phi_binnum+this_rebinned_phi_widths[i]))
                    continue
                else:
                    df_exp_pi0_inb_this_bin = df_exp_pi0_inb.loc[(df_exp_pi0_inb.integrated_binnum == integrated_binnum) & (df_exp_pi0_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_pi0_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_pi0_inb_this_bin = df_sim_pi0_inb.loc[(df_sim_pi0_inb.integrated_binnum == integrated_binnum) & (df_sim_pi0_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_pi0_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_pi0_1gamma_inb_this_bin = df_sim_pi0_1gamma_inb.loc[(df_sim_pi0_1gamma_inb.integrated_binnum == integrated_binnum) & (df_sim_pi0_1gamma_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_pi0_1gamma_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_nobkgmerging_inb_this_bin = df_sim_nobkgmerging_inb.loc[(df_sim_nobkgmerging_inb.integrated_binnum == integrated_binnum) & (df_sim_nobkgmerging_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_nobkgmerging_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_bkgmerging_inb_this_bin = df_sim_bkgmerging_inb.loc[(df_sim_bkgmerging_inb.integrated_binnum == integrated_binnum) & (df_sim_bkgmerging_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_bkgmerging_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]

                    bkg_to_nobkg    = np.sum(df_sim_bkgmerging_inb_this_bin.weights * df_sim_bkgmerging_inb_this_bin.Pthetaefficiency * df_sim_bkgmerging_inb_this_bin.Ppefficiency * df_sim_bkgmerging_inb_this_bin.Pphiefficiency)/np.sum(df_sim_nobkgmerging_inb_this_bin.weights * df_sim_nobkgmerging_inb_this_bin.Pthetaefficiency * df_sim_nobkgmerging_inb_this_bin.Ppefficiency * df_sim_nobkgmerging_inb_this_bin.Pphiefficiency)
                    eff_bkg_merging = ( 1 + effective_current_inb/45 * ( -1 + bkg_to_nobkg))
                    
                    epg_exp_inb_this_bin               = np.sum(1/df_exp_epg_inb_this_bin.efficiency)
                    try:
                        pi0_exp_inb_this_integrated_bin    = np.sum(1/df_exp_pi0_inb.efficiency)
                    except:
                        pi0_exp_inb_this_integrated_bin    = 0
                    try:
                        pi0_sim_inb_this_integrated_bin    = np.sum(df_sim_pi0_inb.Pthetaefficiency)
                    except:
                        pi0_sim_inb_this_integrated_bin    = 0
                    try:
                        bkg_sim_inb_this_bin               = np.sum(df_sim_pi0_1gamma_inb_this_bin.Pthetaefficiency)
                    except:
                        bkg_sim_inb_this_bin               = 0
                    try:
                        bkg_exp_inb_this_bin               = bkg_sim_inb_this_bin * pi0_exp_inb_this_integrated_bin/pi0_sim_inb_this_integrated_bin
                    except:
                        bkg_exp_inb_this_bin               = 0
                    contamination_ratio                = np.minimum(bkg_exp_inb_this_bin, epg_exp_inb_this_bin)/epg_exp_inb_this_bin

                    df_exp_epg_inb.loc[(df_exp_epg_inb.integrated_binnum == integrated_binnum) & (df_exp_epg_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "contamination"]   = contamination_ratio
                    df_exp_epg_inb.loc[(df_exp_epg_inb.integrated_binnum == integrated_binnum) & (df_exp_epg_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "phi_rebinnum"]    = this_rebinned_phi_binnum
                    df_exp_epg_inb.loc[(df_exp_epg_inb.integrated_binnum == integrated_binnum) & (df_exp_epg_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "phi_rebinwidth"]  = this_rebinned_phi_widths[i]
                    df_sim_nobkgmerging_inb.loc[(df_sim_nobkgmerging_inb.integrated_binnum == integrated_binnum) & (df_sim_nobkgmerging_inb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_nobkgmerging_inb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "efficiency"] = eff_bkg_merging

            df_exp_epg_inb.loc[:, "signal"]          = (1 - df_exp_epg_inb.contamination)/df_exp_epg_inb.efficiency
            df_sim_nobkgmerging_inb.loc[:, "signal"] = df_sim_nobkgmerging_inb.weights * df_sim_nobkgmerging_inb.efficiency * df_sim_nobkgmerging_inb.Pthetaefficiency  * df_sim_nobkgmerging_inb.Ppefficiency * df_sim_nobkgmerging_inb.Pphiefficiency
            
            df_exp_epg_inbs.append(df_exp_epg_inb)
            df_sim_dvcs_inbs.append(df_sim_nobkgmerging_inb)

        df_exp_epg_inbs  = pd.concat(df_exp_epg_inbs)
        df_sim_dvcs_inbs = pd.concat(df_sim_dvcs_inbs)

        # FD sector 1
        pphi_exp_sector_1, _ = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 1, "Pphi"], weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 1, "signal"], bins = pphi_bins_inb_fd)
        pphi_sim_sector_1, _ = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 1, "Pphi"], weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 1, "signal"], bins = pphi_bins_inb_fd)
        pphi_efficiency_inb_sector_1 = pphi_efficiency_inb_sector_1*divideHist(pphi_exp_sector_1, pphi_sim_sector_1)

        # FD sector 2
        pphi_exp_sector_2, _ = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 2, "Pphi"], weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 2, "signal"], bins = pphi_bins_inb_fd)
        pphi_sim_sector_2, _ = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 2, "Pphi"], weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 2, "signal"], bins = pphi_bins_inb_fd)
        pphi_efficiency_inb_sector_2 = pphi_efficiency_inb_sector_2*divideHist(pphi_exp_sector_2, pphi_sim_sector_2)

        # FD sector 3
        pphi_exp_sector_3, _ = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 3, "Pphi"], weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 3, "signal"], bins = pphi_bins_inb_fd)
        pphi_sim_sector_3, _ = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 3, "Pphi"], weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 3, "signal"], bins = pphi_bins_inb_fd)
        pphi_efficiency_inb_sector_3 = pphi_efficiency_inb_sector_3*divideHist(pphi_exp_sector_3, pphi_sim_sector_3)

        # FD sector 4
        pphi_exp_sector_4, _ = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 4, "Pphi"], weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 4, "signal"], bins = pphi_bins_inb_fd)
        pphi_sim_sector_4, _ = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 4, "Pphi"], weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 4, "signal"], bins = pphi_bins_inb_fd)
        pphi_efficiency_inb_sector_4 = pphi_efficiency_inb_sector_4*divideHist(pphi_exp_sector_4, pphi_sim_sector_4)

        # FD sector 5
        pphi_exp_sector_5, _ = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 5, "Pphi"], weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 5, "signal"], bins = pphi_bins_inb_fd)
        pphi_sim_sector_5, _ = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 5, "Pphi"], weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 5, "signal"], bins = pphi_bins_inb_fd)
        pphi_efficiency_inb_sector_5 = pphi_efficiency_inb_sector_5*divideHist(pphi_exp_sector_5, pphi_sim_sector_5)

        # FD sector 6
        pphi_exp_sector_6, _ = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 6, "Pphi"], weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector == 6, "signal"], bins = pphi_bins_inb_fd)
        pphi_sim_sector_6, _ = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 6, "Pphi"], weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector == 6, "signal"], bins = pphi_bins_inb_fd)
        pphi_efficiency_inb_sector_6 = pphi_efficiency_inb_sector_6*divideHist(pphi_exp_sector_6, pphi_sim_sector_6)

        # CD
        pphi_exp_cd, _ = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector > 7, "Pphi"], weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector > 7, "signal"], bins = pphi_bins_inb_cd)
        pphi_sim_cd, _ = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector > 7, "Pphi"], weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector > 7, "signal"], bins = pphi_bins_inb_cd)
        pphi_efficiency_inb_cd = pphi_efficiency_inb_cd*divideHist(pphi_exp_cd, pphi_sim_cd)

        pphi_efficiency_inb_sector_1[pphi_efficiency_inb_sector_1<0.5] = 0.5
        pphi_efficiency_inb_sector_2[pphi_efficiency_inb_sector_2<0.5] = 0.5
        pphi_efficiency_inb_sector_3[pphi_efficiency_inb_sector_3<0.5] = 0.5
        pphi_efficiency_inb_sector_4[pphi_efficiency_inb_sector_4<0.5] = 0.5
        pphi_efficiency_inb_sector_5[pphi_efficiency_inb_sector_5<0.5] = 0.5
        pphi_efficiency_inb_sector_6[pphi_efficiency_inb_sector_6<0.5] = 0.5
        pphi_efficiency_inb_cd      [pphi_efficiency_inb_cd      <0.5] = 0.5

        pphi_efficiency_inb_sector_1[pphi_efficiency_inb_sector_1>1.5] = 1.5
        pphi_efficiency_inb_sector_2[pphi_efficiency_inb_sector_2>1.5] = 1.5
        pphi_efficiency_inb_sector_3[pphi_efficiency_inb_sector_3>1.5] = 1.5
        pphi_efficiency_inb_sector_4[pphi_efficiency_inb_sector_4>1.5] = 1.5
        pphi_efficiency_inb_sector_5[pphi_efficiency_inb_sector_5>1.5] = 1.5
        pphi_efficiency_inb_sector_6[pphi_efficiency_inb_sector_6>1.5] = 1.5
        pphi_efficiency_inb_cd      [pphi_efficiency_inb_cd      >1.5] = 1.5


    print("inbending", inbending_trial, "inb sector 1", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_1], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_1], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_1])
    print("inbending", inbending_trial, "inb sector 2", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_2], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_2], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_2])
    print("inbending", inbending_trial, "inb sector 3", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_3], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_3], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_3])
    print("inbending", inbending_trial, "inb sector 4", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_4], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_4], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_4])
    print("inbending", inbending_trial, "inb sector 5", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_5], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_5], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_5])
    print("inbending", inbending_trial, "inb sector 6", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_6], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_6], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_6])
    print("inbending", inbending_trial, "inb CD"      , ["{:.3f}".format(i) for i in ptheta_efficiency_inb_cd      ], ["{:.3f}".format(i) for i in pp_efficiency_inb_cd      ], ["{:.3f}".format(i) for i in pphi_efficiency_inb_cd      ])
#outbending
ptheta_bins_outb_fd             = np.linspace(17, 40, 11)
ptheta_bins_outb_cd             = np.linspace(40, 65, 11)

pp_bins_outb_fd             = np.linspace(0.520,  1.133, 11)
pp_bins_outb_cd             = np.linspace(0.336,  1.133, 11)

pphi_bins_outb_fd             = np.linspace(-180, 180, 37)
pphi_bins_outb_cd             = np.linspace(-180, 180, 37)


#initialize the efficiency map
ptheta_efficiency_outb_sector_1 = np.ones(10)
ptheta_efficiency_outb_sector_2 = np.ones(10)
ptheta_efficiency_outb_sector_3 = np.ones(10)
ptheta_efficiency_outb_sector_4 = np.ones(10)
ptheta_efficiency_outb_sector_5 = np.ones(10)
ptheta_efficiency_outb_sector_6 = np.ones(10)
ptheta_efficiency_outb_cd       = np.ones(10)

pp_efficiency_outb_sector_1 = np.ones(10)
pp_efficiency_outb_sector_2 = np.ones(10)
pp_efficiency_outb_sector_3 = np.ones(10)
pp_efficiency_outb_sector_4 = np.ones(10)
pp_efficiency_outb_sector_5 = np.ones(10)
pp_efficiency_outb_sector_6 = np.ones(10)
pp_efficiency_outb_cd       = np.ones(10)


pphi_efficiency_outb_sector_1 = np.ones(36)
pphi_efficiency_outb_sector_2 = np.ones(36)
pphi_efficiency_outb_sector_3 = np.ones(36)
pphi_efficiency_outb_sector_4 = np.ones(36)
pphi_efficiency_outb_sector_5 = np.ones(36)
pphi_efficiency_outb_sector_6 = np.ones(36)
pphi_efficiency_outb_cd       = np.ones(36)



effective_current_outb = (5*charge_outb_5nA + 40* charge_outb_40nA + 50*charge_outb_50nA)/charge_outb
chunk_outb             = [2, 3]

for outbending_trial in range(10):
    for ptheta_trial in range(3):

        df_exp_epg_outbs        = []
        df_sim_dvcs_outbs       = []
        print(ptheta_trial, "ptheta_trial", "outb sector 1", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_1], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_1], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_1])
        print(ptheta_trial, "ptheta_trial", "outb sector 2", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_2], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_2], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_2])
        print(ptheta_trial, "ptheta_trial", "outb sector 3", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_3], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_3], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_3])
        print(ptheta_trial, "ptheta_trial", "outb sector 4", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_4], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_4], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_4])
        print(ptheta_trial, "ptheta_trial", "outb sector 5", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_5], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_5], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_5])
        print(ptheta_trial, "ptheta_trial", "outb sector 6", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_6], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_6], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_6])
        print(ptheta_trial, "ptheta_trial", "outb CD"      , ["{:.3f}".format(i) for i in ptheta_efficiency_outb_cd      ], ["{:.3f}".format(i) for i in pp_efficiency_outb_cd      ], ["{:.3f}".format(i) for i in pphi_efficiency_outb_cd      ])

        for integrated_binnum in range(1, 147+1):
            this_rebinned_phi_binnums = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"].to_numpy()
            this_rebinned_phi_widths  = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"].to_numpy()

            try:
                df_exp_epg_outb           =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_outb/dvcs/excl_level_2/restructured_7_nominal/{}.pkl".format(integrated_binnum))
                df_exp_epg_outb.loc[df_exp_epg_outb.efficiency == 0, "efficiency"] =  1
                df_exp_epg_outb.loc[:, "contamination"] = 1
            except:
                # print("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_outb/dvcs/excl_level_2/restructured_7_nominal/{}.pkl does not exist.".format(integrated_binnum))
                continue

            # if len(df_exp_epg_outb)<50:
            if df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum), "active_bin_outb"].sum() == 0:
                continue

            try:
                df_exp_pi0_outb            =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_outb/pi0/excl_level_2/restructured_7_nominal/{}.pkl".format(integrated_binnum))
            except:
                df_exp_pi0_outb            =  pd.DataFrame({var: [] for var in ['Pp', 'Ptheta', 'Pphi', 'efficiency']})
            # df_exp_pi0_outb.loc[df_exp_pi0_outb.efficiency == 0, "efficiency"]  = 1
            try:
                df_sim_pi0_outb            =  pd.concat([pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/pi0_2gamma/excl_level_1/restructured_7_nominal/{}/{}.pkl".format(chunk, integrated_binnum)) for chunk in chunk_outb])
                df_sim_pi0_outb            = assign_Pthetaefficiency(df_sim_pi0_outb, ptheta_efficiency_outb_sector_1, ptheta_efficiency_outb_sector_2, ptheta_efficiency_outb_sector_3, ptheta_efficiency_outb_sector_4, ptheta_efficiency_outb_sector_5, ptheta_efficiency_outb_sector_6, ptheta_efficiency_outb_cd, ptheta_bins_outb_fd, ptheta_bins_outb_cd)
                df_sim_pi0_outb            = assign_Ppefficiency(df_sim_pi0_outb, pp_efficiency_outb_sector_1, pp_efficiency_outb_sector_2, pp_efficiency_outb_sector_3, pp_efficiency_outb_sector_4, pp_efficiency_outb_sector_5, pp_efficiency_outb_sector_6, pp_efficiency_outb_cd, pp_bins_outb_fd, pp_bins_outb_cd)
                df_sim_pi0_outb            = assign_Pphiefficiency(df_sim_pi0_outb, pphi_efficiency_outb_sector_1, pphi_efficiency_outb_sector_2, pphi_efficiency_outb_sector_3, pphi_efficiency_outb_sector_4, pphi_efficiency_outb_sector_5, pphi_efficiency_outb_sector_6, pphi_efficiency_outb_cd, pphi_bins_outb_fd, pphi_bins_outb_cd)
            except:
                df_sim_pi0_outb            = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})

            try:
                df_sim_pi0_1gamma_outb     =  pd.concat([pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/pi0_1gamma/excl_level_1/restructured_7_nominal/{}/{}.pkl".format(chunk, integrated_binnum)) for chunk in chunk_outb])
                df_sim_pi0_1gamma_outb     = assign_Pthetaefficiency(df_sim_pi0_1gamma_outb, ptheta_efficiency_outb_sector_1, ptheta_efficiency_outb_sector_2, ptheta_efficiency_outb_sector_3, ptheta_efficiency_outb_sector_4, ptheta_efficiency_outb_sector_5, ptheta_efficiency_outb_sector_6, ptheta_efficiency_outb_cd, ptheta_bins_outb_fd, ptheta_bins_outb_cd)
                df_sim_pi0_1gamma_outb     = assign_Ppefficiency(df_sim_pi0_1gamma_outb, pp_efficiency_outb_sector_1, pp_efficiency_outb_sector_2, pp_efficiency_outb_sector_3, pp_efficiency_outb_sector_4, pp_efficiency_outb_sector_5, pp_efficiency_outb_sector_6, pp_efficiency_outb_cd, pp_bins_outb_fd, pp_bins_outb_cd)
                df_sim_pi0_1gamma_outb     = assign_Pphiefficiency(df_sim_pi0_1gamma_outb, pphi_efficiency_outb_sector_1, pphi_efficiency_outb_sector_2, pphi_efficiency_outb_sector_3, pphi_efficiency_outb_sector_4, pphi_efficiency_outb_sector_5, pphi_efficiency_outb_sector_6, pphi_efficiency_outb_cd, pphi_bins_outb_fd, pphi_bins_outb_cd)
            except:
                df_sim_pi0_1gamma_outb      = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})
            try:
                df_sim_nobkgmerging_outb   =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/dvcs_km15/excl_level_1/restructured_7_nominal/{}.pkl".format(integrated_binnum))
                df_sim_nobkgmerging_outb   = assign_Pthetaefficiency(df_sim_nobkgmerging_outb, ptheta_efficiency_outb_sector_1, ptheta_efficiency_outb_sector_2, ptheta_efficiency_outb_sector_3, ptheta_efficiency_outb_sector_4, ptheta_efficiency_outb_sector_5, ptheta_efficiency_outb_sector_6, ptheta_efficiency_outb_cd, ptheta_bins_outb_fd, ptheta_bins_outb_cd)
                df_sim_nobkgmerging_outb   = assign_Ppefficiency(df_sim_nobkgmerging_outb, pp_efficiency_outb_sector_1, pp_efficiency_outb_sector_2, pp_efficiency_outb_sector_3, pp_efficiency_outb_sector_4, pp_efficiency_outb_sector_5, pp_efficiency_outb_sector_6, pp_efficiency_outb_cd, pp_bins_outb_fd, pp_bins_outb_cd)
                df_sim_nobkgmerging_outb   = assign_Pphiefficiency(df_sim_nobkgmerging_outb, pphi_efficiency_outb_sector_1, pphi_efficiency_outb_sector_2, pphi_efficiency_outb_sector_3, pphi_efficiency_outb_sector_4, pphi_efficiency_outb_sector_5, pphi_efficiency_outb_sector_6, pphi_efficiency_outb_cd, pphi_bins_outb_fd, pphi_bins_outb_cd)
            except:
                df_sim_nobkgmerging_outb      = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})
            try:
                df_sim_bkgmerging_outb     =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/dvcs_km15/excl_level_1/restructured_14_bkgmerging/{}.pkl".format(integrated_binnum))
                df_sim_bkgmerging_outb     = assign_Pthetaefficiency(df_sim_bkgmerging_outb, ptheta_efficiency_outb_sector_1, ptheta_efficiency_outb_sector_2, ptheta_efficiency_outb_sector_3, ptheta_efficiency_outb_sector_4, ptheta_efficiency_outb_sector_5, ptheta_efficiency_outb_sector_6, ptheta_efficiency_outb_cd, ptheta_bins_outb_fd, ptheta_bins_outb_cd)
                df_sim_bkgmerging_outb     = assign_Ppefficiency(df_sim_bkgmerging_outb, pp_efficiency_outb_sector_1, pp_efficiency_outb_sector_2, pp_efficiency_outb_sector_3, pp_efficiency_outb_sector_4, pp_efficiency_outb_sector_5, pp_efficiency_outb_sector_6, pp_efficiency_outb_cd, pp_bins_outb_fd, pp_bins_outb_cd)
                df_sim_bkgmerging_outb     = assign_Pphiefficiency(df_sim_bkgmerging_outb, pphi_efficiency_outb_sector_1, pphi_efficiency_outb_sector_2, pphi_efficiency_outb_sector_3, pphi_efficiency_outb_sector_4, pphi_efficiency_outb_sector_5, pphi_efficiency_outb_sector_6, pphi_efficiency_outb_cd, pphi_bins_outb_fd, pphi_bins_outb_cd)
            except:
                df_sim_bkgmerging_outb      = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})
            # df_sim_nobkgmerging_outb.loc[:, "efficiency"] = 0

            # bkg_to_nobkg    = np.sum(df_sim_bkgmerging_outb.weights * df_sim_bkgmerging_outb.Pthetaefficiency * df_sim_bkgmerging_outb.Ppefficiency * df_sim_bkgmerging_outb.Pphiefficiency)/np.sum(df_sim_nobkgmerging_outb.weights * df_sim_nobkgmerging_outb.Pthetaefficiency * df_sim_nobkgmerging_outb.Ppefficiency * df_sim_nobkgmerging_outb.Pphiefficiency)
            # eff_bkg_merging = ( 1 + effective_current_outb/50 * ( -1 + bkg_to_nobkg))

            # df_sim_nobkgmerging_outb.loc[:, "efficiency"] = eff_bkg_merging
            # epg_exp_outb_this_integrated_bin           = np.sum(1/df_exp_epg_outb.efficiency)
            # pi0_exp_outb_this_integrated_bin           = np.sum(1/df_exp_pi0_outb.efficiency)
            # pi0_sim_outb_this_integrated_bin           = np.sum(df_sim_pi0_outb.Pthetaefficiency * df_sim_pi0_outb.Ppefficiency * df_sim_pi0_outb.Pphiefficiency)
            # bkg_sim_outb_this_integrated_bin           = np.sum(df_sim_pi0_1gamma_outb.Pthetaefficiency * df_sim_pi0_1gamma_outb.Ppefficiency * df_sim_pi0_1gamma_outb.Pphiefficiency)
            # bkg_exp_outb_this_integrated_bin           = bkg_sim_outb_this_integrated_bin  * pi0_exp_outb_this_integrated_bin/pi0_sim_outb_this_integrated_bin
            # contamination_ratio                       = np.minimum(bkg_exp_outb_this_integrated_bin , epg_exp_outb_this_integrated_bin )/epg_exp_outb_this_integrated_bin 

            # df_exp_epg_outb.loc[:, "contamination"]   = contamination_ratio

            for i, this_rebinned_phi_binnum in enumerate(this_rebinned_phi_binnums):
                df_exp_epg_outb_this_bin = df_exp_epg_outb.loc[(df_exp_epg_outb.integrated_binnum == integrated_binnum) & (df_exp_epg_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                # if len(df_exp_epg_outb_this_bin) < 10:
                if df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == this_rebinned_phi_binnum), "active_bin_outb"].sum() == 0:
                    # print("{} inactive bin between {} and {}".format(integrated_binnum, this_rebinned_phi_binnum, this_rebinned_phi_binnum+this_rebinned_phi_widths[i]))
                    continue
                else:
                    df_exp_pi0_outb_this_bin = df_exp_pi0_outb.loc[(df_exp_pi0_outb.integrated_binnum == integrated_binnum) & (df_exp_pi0_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_pi0_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_pi0_outb_this_bin = df_sim_pi0_outb.loc[(df_sim_pi0_outb.integrated_binnum == integrated_binnum) & (df_sim_pi0_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_pi0_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_pi0_1gamma_outb_this_bin = df_sim_pi0_1gamma_outb.loc[(df_sim_pi0_1gamma_outb.integrated_binnum == integrated_binnum) & (df_sim_pi0_1gamma_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_pi0_1gamma_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_nobkgmerging_outb_this_bin = df_sim_nobkgmerging_outb.loc[(df_sim_nobkgmerging_outb.integrated_binnum == integrated_binnum) & (df_sim_nobkgmerging_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_nobkgmerging_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_bkgmerging_outb_this_bin = df_sim_bkgmerging_outb.loc[(df_sim_bkgmerging_outb.integrated_binnum == integrated_binnum) & (df_sim_bkgmerging_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_bkgmerging_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]

                    bkg_to_nobkg    = np.sum(df_sim_bkgmerging_outb_this_bin.weights * df_sim_bkgmerging_outb_this_bin.Pthetaefficiency * df_sim_bkgmerging_outb_this_bin.Ppefficiency * df_sim_bkgmerging_outb_this_bin.Pphiefficiency)/np.sum(df_sim_nobkgmerging_outb_this_bin.weights * df_sim_nobkgmerging_outb_this_bin.Pthetaefficiency * df_sim_nobkgmerging_outb_this_bin.Ppefficiency * df_sim_nobkgmerging_outb_this_bin.Pphiefficiency)
                    eff_bkg_merging = ( 1 + effective_current_outb/50 * ( -1 + bkg_to_nobkg))
                    
                    epg_exp_outb_this_bin               = np.sum(1/df_exp_epg_outb_this_bin.efficiency)
                    try:
                        pi0_exp_outb_this_integrated_bin    = np.sum(1/df_exp_pi0_outb.efficiency)
                    except:
                        pi0_exp_outb_this_integrated_bin    = 0
                    try:
                        pi0_sim_outb_this_integrated_bin    = np.sum(df_sim_pi0_outb.Pthetaefficiency)
                    except:
                        pi0_sim_outb_this_integrated_bin    = 0
                    try:
                        bkg_sim_outb_this_bin               = np.sum(df_sim_pi0_1gamma_outb_this_bin.Pthetaefficiency)
                    except:
                        bkg_sim_outb_this_bin               = 0
                    try:
                        bkg_exp_outb_this_bin               = bkg_sim_outb_this_bin * pi0_exp_outb_this_integrated_bin/pi0_sim_outb_this_integrated_bin
                    except:
                        bkg_exp_outb_this_bin               = 0
                    contamination_ratio                = np.minimum(bkg_exp_outb_this_bin, epg_exp_outb_this_bin)/epg_exp_outb_this_bin

                    df_exp_epg_outb.loc[(df_exp_epg_outb.integrated_binnum == integrated_binnum) & (df_exp_epg_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "contamination"]   = contamination_ratio
                    df_exp_epg_outb.loc[(df_exp_epg_outb.integrated_binnum == integrated_binnum) & (df_exp_epg_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "phi_rebinnum"]    = this_rebinned_phi_binnum
                    df_exp_epg_outb.loc[(df_exp_epg_outb.integrated_binnum == integrated_binnum) & (df_exp_epg_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "phi_rebinwidth"]  = this_rebinned_phi_widths[i]
                    df_sim_nobkgmerging_outb.loc[(df_sim_nobkgmerging_outb.integrated_binnum == integrated_binnum) & (df_sim_nobkgmerging_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_nobkgmerging_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "efficiency"] = eff_bkg_merging

            df_exp_epg_outb.loc[:, "signal"]          = (1 - df_exp_epg_outb.contamination)/df_exp_epg_outb.efficiency
            df_sim_nobkgmerging_outb.loc[:, "signal"] = df_sim_nobkgmerging_outb.weights * df_sim_nobkgmerging_outb.efficiency * df_sim_nobkgmerging_outb.Pthetaefficiency  * df_sim_nobkgmerging_outb.Ppefficiency * df_sim_nobkgmerging_outb.Pphiefficiency
            
            df_exp_epg_outbs.append(df_exp_epg_outb)
            df_sim_dvcs_outbs.append(df_sim_nobkgmerging_outb)

        df_exp_epg_outbs  = pd.concat(df_exp_epg_outbs)
        df_sim_dvcs_outbs = pd.concat(df_sim_dvcs_outbs)

        # FD sector 1
        ptheta_exp_sector_1, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 1, "Ptheta"], weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 1, "signal"], bins = ptheta_bins_outb_fd)
        ptheta_sim_sector_1, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 1, "Ptheta"], weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 1, "signal"], bins = ptheta_bins_outb_fd)
        ptheta_efficiency_outb_sector_1 = ptheta_efficiency_outb_sector_1*divideHist(ptheta_exp_sector_1, ptheta_sim_sector_1)

        # FD sector 2
        ptheta_exp_sector_2, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 2, "Ptheta"], weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 2, "signal"], bins = ptheta_bins_outb_fd)
        ptheta_sim_sector_2, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 2, "Ptheta"], weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 2, "signal"], bins = ptheta_bins_outb_fd)
        ptheta_efficiency_outb_sector_2 = ptheta_efficiency_outb_sector_2*divideHist(ptheta_exp_sector_2, ptheta_sim_sector_2)

        # FD sector 3
        ptheta_exp_sector_3, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 3, "Ptheta"], weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 3, "signal"], bins = ptheta_bins_outb_fd)
        ptheta_sim_sector_3, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 3, "Ptheta"], weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 3, "signal"], bins = ptheta_bins_outb_fd)
        ptheta_efficiency_outb_sector_3 = ptheta_efficiency_outb_sector_3*divideHist(ptheta_exp_sector_3, ptheta_sim_sector_3)

        # FD sector 4
        ptheta_exp_sector_4, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 4, "Ptheta"], weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 4, "signal"], bins = ptheta_bins_outb_fd)
        ptheta_sim_sector_4, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 4, "Ptheta"], weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 4, "signal"], bins = ptheta_bins_outb_fd)
        ptheta_efficiency_outb_sector_4 = ptheta_efficiency_outb_sector_4*divideHist(ptheta_exp_sector_4, ptheta_sim_sector_4)

        # FD sector 5
        ptheta_exp_sector_5, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 5, "Ptheta"], weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 5, "signal"], bins = ptheta_bins_outb_fd)
        ptheta_sim_sector_5, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 5, "Ptheta"], weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 5, "signal"], bins = ptheta_bins_outb_fd)
        ptheta_efficiency_outb_sector_5 = ptheta_efficiency_outb_sector_5*divideHist(ptheta_exp_sector_5, ptheta_sim_sector_5)

        # FD sector 6
        ptheta_exp_sector_6, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 6, "Ptheta"], weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 6, "signal"], bins = ptheta_bins_outb_fd)
        ptheta_sim_sector_6, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 6, "Ptheta"], weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 6, "signal"], bins = ptheta_bins_outb_fd)
        ptheta_efficiency_outb_sector_6 = ptheta_efficiency_outb_sector_6*divideHist(ptheta_exp_sector_6, ptheta_sim_sector_6)

        # CD
        ptheta_exp_cd, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector > 7, "Ptheta"], weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector > 7, "signal"], bins = ptheta_bins_outb_cd)
        ptheta_sim_cd, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector > 7, "Ptheta"], weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector > 7, "signal"], bins = ptheta_bins_outb_cd)
        ptheta_efficiency_outb_cd = ptheta_efficiency_outb_cd*divideHist(ptheta_exp_cd, ptheta_sim_cd)

        ptheta_efficiency_outb_sector_1[ptheta_efficiency_outb_sector_1<0.1] = 0.1
        ptheta_efficiency_outb_sector_2[ptheta_efficiency_outb_sector_2<0.1] = 0.1
        ptheta_efficiency_outb_sector_3[ptheta_efficiency_outb_sector_3<0.1] = 0.1
        ptheta_efficiency_outb_sector_4[ptheta_efficiency_outb_sector_4<0.1] = 0.1
        ptheta_efficiency_outb_sector_5[ptheta_efficiency_outb_sector_5<0.1] = 0.1
        ptheta_efficiency_outb_sector_6[ptheta_efficiency_outb_sector_6<0.1] = 0.1
        ptheta_efficiency_outb_cd      [ptheta_efficiency_outb_cd      <0.1] = 0.1

        ptheta_efficiency_outb_sector_1[ptheta_efficiency_outb_sector_1>1.5] = 1.5
        ptheta_efficiency_outb_sector_2[ptheta_efficiency_outb_sector_2>1.5] = 1.5
        ptheta_efficiency_outb_sector_3[ptheta_efficiency_outb_sector_3>1.5] = 1.5
        ptheta_efficiency_outb_sector_4[ptheta_efficiency_outb_sector_4>1.5] = 1.5
        ptheta_efficiency_outb_sector_5[ptheta_efficiency_outb_sector_5>1.5] = 1.5
        ptheta_efficiency_outb_sector_6[ptheta_efficiency_outb_sector_6>1.5] = 1.5
        ptheta_efficiency_outb_cd      [ptheta_efficiency_outb_cd      >1.5] = 1.5

    for pp_trial in range(3):

        df_exp_epg_outbs        = []
        df_sim_dvcs_outbs       = []
        print(pp_trial, "pp_trial", "outb sector 1", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_1], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_1], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_1])
        print(pp_trial, "pp_trial", "outb sector 2", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_2], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_2], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_2])
        print(pp_trial, "pp_trial", "outb sector 3", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_3], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_3], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_3])
        print(pp_trial, "pp_trial", "outb sector 4", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_4], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_4], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_4])
        print(pp_trial, "pp_trial", "outb sector 5", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_5], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_5], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_5])
        print(pp_trial, "pp_trial", "outb sector 6", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_6], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_6], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_6])
        print(pp_trial, "pp_trial", "outb CD"      , ["{:.3f}".format(i) for i in ptheta_efficiency_outb_cd      ], ["{:.3f}".format(i) for i in pp_efficiency_outb_cd      ], ["{:.3f}".format(i) for i in pphi_efficiency_outb_cd      ])

        for integrated_binnum in range(1, 147+1):
            this_rebinned_phi_binnums = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"].to_numpy()
            this_rebinned_phi_widths  = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"].to_numpy()

            try:
                df_exp_epg_outb           =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_outb/dvcs/excl_level_2/restructured_7_nominal/{}.pkl".format(integrated_binnum))
                df_exp_epg_outb.loc[df_exp_epg_outb.efficiency == 0, "efficiency"] =  1
                df_exp_epg_outb.loc[:, "contamination"] = 1
            except:
                # print("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_outb/dvcs/excl_level_2/restructured_7_nominal/{}.pkl does not exist.".format(integrated_binnum))
                continue

            # if len(df_exp_epg_outb)<50:
            if df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum), "active_bin_outb"].sum() == 0:
                continue

            try:
                df_exp_pi0_outb            =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_outb/pi0/excl_level_2/restructured_7_nominal/{}.pkl".format(integrated_binnum))
            except:
                df_exp_pi0_outb            =  pd.DataFrame({var: [] for var in ['Pp', 'Ptheta', 'Pphi', 'efficiency']})
            # df_exp_pi0_outb.loc[df_exp_pi0_outb.efficiency == 0, "efficiency"]  = 1
            try:
                df_sim_pi0_outb            =  pd.concat([pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/pi0_2gamma/excl_level_1/restructured_7_nominal/{}/{}.pkl".format(chunk, integrated_binnum)) for chunk in chunk_outb])
                df_sim_pi0_outb            = assign_Pthetaefficiency(df_sim_pi0_outb, ptheta_efficiency_outb_sector_1, ptheta_efficiency_outb_sector_2, ptheta_efficiency_outb_sector_3, ptheta_efficiency_outb_sector_4, ptheta_efficiency_outb_sector_5, ptheta_efficiency_outb_sector_6, ptheta_efficiency_outb_cd, ptheta_bins_outb_fd, ptheta_bins_outb_cd)
                df_sim_pi0_outb            = assign_Ppefficiency(df_sim_pi0_outb, pp_efficiency_outb_sector_1, pp_efficiency_outb_sector_2, pp_efficiency_outb_sector_3, pp_efficiency_outb_sector_4, pp_efficiency_outb_sector_5, pp_efficiency_outb_sector_6, pp_efficiency_outb_cd, pp_bins_outb_fd, pp_bins_outb_cd)
                df_sim_pi0_outb            = assign_Pphiefficiency(df_sim_pi0_outb, pphi_efficiency_outb_sector_1, pphi_efficiency_outb_sector_2, pphi_efficiency_outb_sector_3, pphi_efficiency_outb_sector_4, pphi_efficiency_outb_sector_5, pphi_efficiency_outb_sector_6, pphi_efficiency_outb_cd, pphi_bins_outb_fd, pphi_bins_outb_cd)
            except:
                df_sim_pi0_outb            = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})

            try:
                df_sim_pi0_1gamma_outb     =  pd.concat([pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/pi0_1gamma/excl_level_1/restructured_7_nominal/{}/{}.pkl".format(chunk, integrated_binnum)) for chunk in chunk_outb])
                df_sim_pi0_1gamma_outb     = assign_Pthetaefficiency(df_sim_pi0_1gamma_outb, ptheta_efficiency_outb_sector_1, ptheta_efficiency_outb_sector_2, ptheta_efficiency_outb_sector_3, ptheta_efficiency_outb_sector_4, ptheta_efficiency_outb_sector_5, ptheta_efficiency_outb_sector_6, ptheta_efficiency_outb_cd, ptheta_bins_outb_fd, ptheta_bins_outb_cd)
                df_sim_pi0_1gamma_outb     = assign_Ppefficiency(df_sim_pi0_1gamma_outb, pp_efficiency_outb_sector_1, pp_efficiency_outb_sector_2, pp_efficiency_outb_sector_3, pp_efficiency_outb_sector_4, pp_efficiency_outb_sector_5, pp_efficiency_outb_sector_6, pp_efficiency_outb_cd, pp_bins_outb_fd, pp_bins_outb_cd)
                df_sim_pi0_1gamma_outb     = assign_Pphiefficiency(df_sim_pi0_1gamma_outb, pphi_efficiency_outb_sector_1, pphi_efficiency_outb_sector_2, pphi_efficiency_outb_sector_3, pphi_efficiency_outb_sector_4, pphi_efficiency_outb_sector_5, pphi_efficiency_outb_sector_6, pphi_efficiency_outb_cd, pphi_bins_outb_fd, pphi_bins_outb_cd)
            except:
                df_sim_pi0_1gamma_outb      = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})
            try:
                df_sim_nobkgmerging_outb   =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/dvcs_km15/excl_level_1/restructured_7_nominal/{}.pkl".format(integrated_binnum))
                df_sim_nobkgmerging_outb   = assign_Pthetaefficiency(df_sim_nobkgmerging_outb, ptheta_efficiency_outb_sector_1, ptheta_efficiency_outb_sector_2, ptheta_efficiency_outb_sector_3, ptheta_efficiency_outb_sector_4, ptheta_efficiency_outb_sector_5, ptheta_efficiency_outb_sector_6, ptheta_efficiency_outb_cd, ptheta_bins_outb_fd, ptheta_bins_outb_cd)
                df_sim_nobkgmerging_outb   = assign_Ppefficiency(df_sim_nobkgmerging_outb, pp_efficiency_outb_sector_1, pp_efficiency_outb_sector_2, pp_efficiency_outb_sector_3, pp_efficiency_outb_sector_4, pp_efficiency_outb_sector_5, pp_efficiency_outb_sector_6, pp_efficiency_outb_cd, pp_bins_outb_fd, pp_bins_outb_cd)
                df_sim_nobkgmerging_outb   = assign_Pphiefficiency(df_sim_nobkgmerging_outb, pphi_efficiency_outb_sector_1, pphi_efficiency_outb_sector_2, pphi_efficiency_outb_sector_3, pphi_efficiency_outb_sector_4, pphi_efficiency_outb_sector_5, pphi_efficiency_outb_sector_6, pphi_efficiency_outb_cd, pphi_bins_outb_fd, pphi_bins_outb_cd)
            except:
                df_sim_nobkgmerging_outb      = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})
            try:
                df_sim_bkgmerging_outb     =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/dvcs_km15/excl_level_1/restructured_14_bkgmerging/{}.pkl".format(integrated_binnum))
                df_sim_bkgmerging_outb     = assign_Pthetaefficiency(df_sim_bkgmerging_outb, ptheta_efficiency_outb_sector_1, ptheta_efficiency_outb_sector_2, ptheta_efficiency_outb_sector_3, ptheta_efficiency_outb_sector_4, ptheta_efficiency_outb_sector_5, ptheta_efficiency_outb_sector_6, ptheta_efficiency_outb_cd, ptheta_bins_outb_fd, ptheta_bins_outb_cd)
                df_sim_bkgmerging_outb     = assign_Ppefficiency(df_sim_bkgmerging_outb, pp_efficiency_outb_sector_1, pp_efficiency_outb_sector_2, pp_efficiency_outb_sector_3, pp_efficiency_outb_sector_4, pp_efficiency_outb_sector_5, pp_efficiency_outb_sector_6, pp_efficiency_outb_cd, pp_bins_outb_fd, pp_bins_outb_cd)
                df_sim_bkgmerging_outb     = assign_Pphiefficiency(df_sim_bkgmerging_outb, pphi_efficiency_outb_sector_1, pphi_efficiency_outb_sector_2, pphi_efficiency_outb_sector_3, pphi_efficiency_outb_sector_4, pphi_efficiency_outb_sector_5, pphi_efficiency_outb_sector_6, pphi_efficiency_outb_cd, pphi_bins_outb_fd, pphi_bins_outb_cd)
            except:
                df_sim_bkgmerging_outb      = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})
            # df_sim_nobkgmerging_outb.loc[:, "efficiency"] = 0

            # bkg_to_nobkg    = np.sum(df_sim_bkgmerging_outb.weights * df_sim_bkgmerging_outb.Pthetaefficiency * df_sim_bkgmerging_outb.Ppefficiency * df_sim_bkgmerging_outb.Pphiefficiency)/np.sum(df_sim_nobkgmerging_outb.weights * df_sim_nobkgmerging_outb.Pthetaefficiency * df_sim_nobkgmerging_outb.Ppefficiency * df_sim_nobkgmerging_outb.Pphiefficiency)
            # eff_bkg_merging = ( 1 + effective_current_outb/50 * ( -1 + bkg_to_nobkg))

            # df_sim_nobkgmerging_outb.loc[:, "efficiency"] = eff_bkg_merging
            # epg_exp_outb_this_integrated_bin           = np.sum(1/df_exp_epg_outb.efficiency)
            # pi0_exp_outb_this_integrated_bin           = np.sum(1/df_exp_pi0_outb.efficiency)
            # pi0_sim_outb_this_integrated_bin           = np.sum(df_sim_pi0_outb.Pthetaefficiency * df_sim_pi0_outb.Ppefficiency * df_sim_pi0_outb.Pphiefficiency)
            # bkg_sim_outb_this_integrated_bin           = np.sum(df_sim_pi0_1gamma_outb.Pthetaefficiency * df_sim_pi0_1gamma_outb.Ppefficiency * df_sim_pi0_1gamma_outb.Pphiefficiency)
            # bkg_exp_outb_this_integrated_bin           = bkg_sim_outb_this_integrated_bin  * pi0_exp_outb_this_integrated_bin/pi0_sim_outb_this_integrated_bin
            # contamination_ratio                       = np.minimum(bkg_exp_outb_this_integrated_bin , epg_exp_outb_this_integrated_bin )/epg_exp_outb_this_integrated_bin 

            # df_exp_epg_outb.loc[:, "contamination"]   = contamination_ratio

            for i, this_rebinned_phi_binnum in enumerate(this_rebinned_phi_binnums):
                df_exp_epg_outb_this_bin = df_exp_epg_outb.loc[(df_exp_epg_outb.integrated_binnum == integrated_binnum) & (df_exp_epg_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                # if len(df_exp_epg_outb_this_bin) < 10:
                if df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == this_rebinned_phi_binnum), "active_bin_outb"].sum() == 0:
                    # print("{} inactive bin between {} and {}".format(integrated_binnum, this_rebinned_phi_binnum, this_rebinned_phi_binnum+this_rebinned_phi_widths[i]))
                    continue
                else:
                    df_exp_pi0_outb_this_bin = df_exp_pi0_outb.loc[(df_exp_pi0_outb.integrated_binnum == integrated_binnum) & (df_exp_pi0_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_pi0_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_pi0_outb_this_bin = df_sim_pi0_outb.loc[(df_sim_pi0_outb.integrated_binnum == integrated_binnum) & (df_sim_pi0_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_pi0_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_pi0_1gamma_outb_this_bin = df_sim_pi0_1gamma_outb.loc[(df_sim_pi0_1gamma_outb.integrated_binnum == integrated_binnum) & (df_sim_pi0_1gamma_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_pi0_1gamma_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_nobkgmerging_outb_this_bin = df_sim_nobkgmerging_outb.loc[(df_sim_nobkgmerging_outb.integrated_binnum == integrated_binnum) & (df_sim_nobkgmerging_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_nobkgmerging_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_bkgmerging_outb_this_bin = df_sim_bkgmerging_outb.loc[(df_sim_bkgmerging_outb.integrated_binnum == integrated_binnum) & (df_sim_bkgmerging_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_bkgmerging_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]

                    bkg_to_nobkg    = np.sum(df_sim_bkgmerging_outb_this_bin.weights * df_sim_bkgmerging_outb_this_bin.Pthetaefficiency * df_sim_bkgmerging_outb_this_bin.Ppefficiency * df_sim_bkgmerging_outb_this_bin.Pphiefficiency)/np.sum(df_sim_nobkgmerging_outb_this_bin.weights * df_sim_nobkgmerging_outb_this_bin.Pthetaefficiency * df_sim_nobkgmerging_outb_this_bin.Ppefficiency * df_sim_nobkgmerging_outb_this_bin.Pphiefficiency)
                    eff_bkg_merging = ( 1 + effective_current_outb/50 * ( -1 + bkg_to_nobkg))
                    
                    epg_exp_outb_this_bin               = np.sum(1/df_exp_epg_outb_this_bin.efficiency)
                    try:
                        pi0_exp_outb_this_integrated_bin    = np.sum(1/df_exp_pi0_outb.efficiency)
                    except:
                        pi0_exp_outb_this_integrated_bin    = 0
                    try:
                        pi0_sim_outb_this_integrated_bin    = np.sum(df_sim_pi0_outb.Pthetaefficiency)
                    except:
                        pi0_sim_outb_this_integrated_bin    = 0
                    try:
                        bkg_sim_outb_this_bin               = np.sum(df_sim_pi0_1gamma_outb_this_bin.Pthetaefficiency)
                    except:
                        bkg_sim_outb_this_bin               = 0
                    try:
                        bkg_exp_outb_this_bin               = bkg_sim_outb_this_bin * pi0_exp_outb_this_integrated_bin/pi0_sim_outb_this_integrated_bin
                    except:
                        bkg_exp_outb_this_bin               = 0
                    contamination_ratio                = np.minimum(bkg_exp_outb_this_bin, epg_exp_outb_this_bin)/epg_exp_outb_this_bin

                    df_exp_epg_outb.loc[(df_exp_epg_outb.integrated_binnum == integrated_binnum) & (df_exp_epg_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "contamination"]   = contamination_ratio
                    df_exp_epg_outb.loc[(df_exp_epg_outb.integrated_binnum == integrated_binnum) & (df_exp_epg_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "phi_rebinnum"]    = this_rebinned_phi_binnum
                    df_exp_epg_outb.loc[(df_exp_epg_outb.integrated_binnum == integrated_binnum) & (df_exp_epg_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "phi_rebinwidth"]  = this_rebinned_phi_widths[i]
                    df_sim_nobkgmerging_outb.loc[(df_sim_nobkgmerging_outb.integrated_binnum == integrated_binnum) & (df_sim_nobkgmerging_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_nobkgmerging_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "efficiency"] = eff_bkg_merging

            df_exp_epg_outb.loc[:, "signal"]          = (1 - df_exp_epg_outb.contamination)/df_exp_epg_outb.efficiency
            df_sim_nobkgmerging_outb.loc[:, "signal"] = df_sim_nobkgmerging_outb.weights * df_sim_nobkgmerging_outb.efficiency * df_sim_nobkgmerging_outb.Pthetaefficiency  * df_sim_nobkgmerging_outb.Ppefficiency * df_sim_nobkgmerging_outb.Pphiefficiency
            
            df_exp_epg_outbs.append(df_exp_epg_outb)
            df_sim_dvcs_outbs.append(df_sim_nobkgmerging_outb)

        df_exp_epg_outbs  = pd.concat(df_exp_epg_outbs)
        df_sim_dvcs_outbs = pd.concat(df_sim_dvcs_outbs)

        # FD sector 1
        pp_exp_sector_1, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 1, "Pp"], weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 1, "signal"], bins = pp_bins_outb_fd)
        pp_sim_sector_1, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 1, "Pp"], weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 1, "signal"], bins = pp_bins_outb_fd)
        pp_efficiency_outb_sector_1 = pp_efficiency_outb_sector_1*divideHist(pp_exp_sector_1, pp_sim_sector_1)

        # FD sector 2
        pp_exp_sector_2, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 2, "Pp"], weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 2, "signal"], bins = pp_bins_outb_fd)
        pp_sim_sector_2, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 2, "Pp"], weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 2, "signal"], bins = pp_bins_outb_fd)
        pp_efficiency_outb_sector_2 = pp_efficiency_outb_sector_2*divideHist(pp_exp_sector_2, pp_sim_sector_2)

        # FD sector 3
        pp_exp_sector_3, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 3, "Pp"], weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 3, "signal"], bins = pp_bins_outb_fd)
        pp_sim_sector_3, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 3, "Pp"], weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 3, "signal"], bins = pp_bins_outb_fd)
        pp_efficiency_outb_sector_3 = pp_efficiency_outb_sector_3*divideHist(pp_exp_sector_3, pp_sim_sector_3)

        # FD sector 4
        pp_exp_sector_4, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 4, "Pp"], weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 4, "signal"], bins = pp_bins_outb_fd)
        pp_sim_sector_4, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 4, "Pp"], weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 4, "signal"], bins = pp_bins_outb_fd)
        pp_efficiency_outb_sector_4 = pp_efficiency_outb_sector_4*divideHist(pp_exp_sector_4, pp_sim_sector_4)

        # FD sector 5
        pp_exp_sector_5, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 5, "Pp"], weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 5, "signal"], bins = pp_bins_outb_fd)
        pp_sim_sector_5, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 5, "Pp"], weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 5, "signal"], bins = pp_bins_outb_fd)
        pp_efficiency_outb_sector_5 = pp_efficiency_outb_sector_5*divideHist(pp_exp_sector_5, pp_sim_sector_5)

        # FD sector 6
        pp_exp_sector_6, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 6, "Pp"], weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 6, "signal"], bins = pp_bins_outb_fd)
        pp_sim_sector_6, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 6, "Pp"], weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 6, "signal"], bins = pp_bins_outb_fd)
        pp_efficiency_outb_sector_6 = pp_efficiency_outb_sector_6*divideHist(pp_exp_sector_6, pp_sim_sector_6)

        # CD
        pp_exp_cd, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector > 7, "Pp"], weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector > 7, "signal"], bins = pp_bins_outb_cd)
        pp_sim_cd, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector > 7, "Pp"], weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector > 7, "signal"], bins = pp_bins_outb_cd)
        pp_efficiency_outb_cd = pp_efficiency_outb_cd*divideHist(pp_exp_cd, pp_sim_cd)

        pp_efficiency_outb_sector_1[pp_efficiency_outb_sector_1<0.5] = 0.5
        pp_efficiency_outb_sector_2[pp_efficiency_outb_sector_2<0.5] = 0.5
        pp_efficiency_outb_sector_3[pp_efficiency_outb_sector_3<0.5] = 0.5
        pp_efficiency_outb_sector_4[pp_efficiency_outb_sector_4<0.5] = 0.5
        pp_efficiency_outb_sector_5[pp_efficiency_outb_sector_5<0.5] = 0.5
        pp_efficiency_outb_sector_6[pp_efficiency_outb_sector_6<0.5] = 0.5
        pp_efficiency_outb_cd      [pp_efficiency_outb_cd      <0.5] = 0.5

        pp_efficiency_outb_sector_1[pp_efficiency_outb_sector_1>1.5] = 1.5
        pp_efficiency_outb_sector_2[pp_efficiency_outb_sector_2>1.5] = 1.5
        pp_efficiency_outb_sector_3[pp_efficiency_outb_sector_3>1.5] = 1.5
        pp_efficiency_outb_sector_4[pp_efficiency_outb_sector_4>1.5] = 1.5
        pp_efficiency_outb_sector_5[pp_efficiency_outb_sector_5>1.5] = 1.5
        pp_efficiency_outb_sector_6[pp_efficiency_outb_sector_6>1.5] = 1.5
        pp_efficiency_outb_cd      [pp_efficiency_outb_cd      >1.5] = 1.5

    for pphi_trial in range(3):

        df_exp_epg_outbs        = []
        df_sim_dvcs_outbs       = []
        print(pphi_trial, "pphi_trial", "outb sector 1", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_1], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_1], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_1])
        print(pphi_trial, "pphi_trial", "outb sector 2", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_2], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_2], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_2])
        print(pphi_trial, "pphi_trial", "outb sector 3", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_3], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_3], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_3])
        print(pphi_trial, "pphi_trial", "outb sector 4", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_4], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_4], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_4])
        print(pphi_trial, "pphi_trial", "outb sector 5", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_5], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_5], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_5])
        print(pphi_trial, "pphi_trial", "outb sector 6", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_6], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_6], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_6])
        print(pphi_trial, "pphi_trial", "outb CD"      , ["{:.3f}".format(i) for i in ptheta_efficiency_outb_cd      ], ["{:.3f}".format(i) for i in pp_efficiency_outb_cd      ], ["{:.3f}".format(i) for i in pphi_efficiency_outb_cd      ])

        for integrated_binnum in range(1, 147+1):
            this_rebinned_phi_binnums = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_binnum"].to_numpy()
            this_rebinned_phi_widths  = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "phi_width"].to_numpy()

            try:
                df_exp_epg_outb           =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_outb/dvcs/excl_level_2/restructured_7_nominal/{}.pkl".format(integrated_binnum))
                df_exp_epg_outb.loc[df_exp_epg_outb.efficiency == 0, "efficiency"] =  1
                df_exp_epg_outb.loc[:, "contamination"] = 1
            except:
                # print("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_outb/dvcs/excl_level_2/restructured_7_nominal/{}.pkl does not exist.".format(integrated_binnum))
                continue

            # if len(df_exp_epg_outb)<50:
            if df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum), "active_bin_outb"].sum() == 0:
                continue

            try:
                df_exp_pi0_outb            =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_outb/pi0/excl_level_2/restructured_7_nominal/{}.pkl".format(integrated_binnum))
            except:
                df_exp_pi0_outb            =  pd.DataFrame({var: [] for var in ['Pp', 'Ptheta', 'Pphi', 'efficiency']})
            # df_exp_pi0_outb.loc[df_exp_pi0_outb.efficiency == 0, "efficiency"]  = 1
            try:
                df_sim_pi0_outb            =  pd.concat([pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/pi0_2gamma/excl_level_1/restructured_7_nominal/{}/{}.pkl".format(chunk, integrated_binnum)) for chunk in chunk_outb])
                df_sim_pi0_outb            = assign_Pthetaefficiency(df_sim_pi0_outb, ptheta_efficiency_outb_sector_1, ptheta_efficiency_outb_sector_2, ptheta_efficiency_outb_sector_3, ptheta_efficiency_outb_sector_4, ptheta_efficiency_outb_sector_5, ptheta_efficiency_outb_sector_6, ptheta_efficiency_outb_cd, ptheta_bins_outb_fd, ptheta_bins_outb_cd)
                df_sim_pi0_outb            = assign_Ppefficiency(df_sim_pi0_outb, pp_efficiency_outb_sector_1, pp_efficiency_outb_sector_2, pp_efficiency_outb_sector_3, pp_efficiency_outb_sector_4, pp_efficiency_outb_sector_5, pp_efficiency_outb_sector_6, pp_efficiency_outb_cd, pp_bins_outb_fd, pp_bins_outb_cd)
                df_sim_pi0_outb            = assign_Pphiefficiency(df_sim_pi0_outb, pphi_efficiency_outb_sector_1, pphi_efficiency_outb_sector_2, pphi_efficiency_outb_sector_3, pphi_efficiency_outb_sector_4, pphi_efficiency_outb_sector_5, pphi_efficiency_outb_sector_6, pphi_efficiency_outb_cd, pphi_bins_outb_fd, pphi_bins_outb_cd)
            except:
                df_sim_pi0_outb            = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})

            try:
                df_sim_pi0_1gamma_outb     =  pd.concat([pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/pi0_1gamma/excl_level_1/restructured_7_nominal/{}/{}.pkl".format(chunk, integrated_binnum)) for chunk in chunk_outb])
                df_sim_pi0_1gamma_outb     = assign_Pthetaefficiency(df_sim_pi0_1gamma_outb, ptheta_efficiency_outb_sector_1, ptheta_efficiency_outb_sector_2, ptheta_efficiency_outb_sector_3, ptheta_efficiency_outb_sector_4, ptheta_efficiency_outb_sector_5, ptheta_efficiency_outb_sector_6, ptheta_efficiency_outb_cd, ptheta_bins_outb_fd, ptheta_bins_outb_cd)
                df_sim_pi0_1gamma_outb     = assign_Ppefficiency(df_sim_pi0_1gamma_outb, pp_efficiency_outb_sector_1, pp_efficiency_outb_sector_2, pp_efficiency_outb_sector_3, pp_efficiency_outb_sector_4, pp_efficiency_outb_sector_5, pp_efficiency_outb_sector_6, pp_efficiency_outb_cd, pp_bins_outb_fd, pp_bins_outb_cd)
                df_sim_pi0_1gamma_outb     = assign_Pphiefficiency(df_sim_pi0_1gamma_outb, pphi_efficiency_outb_sector_1, pphi_efficiency_outb_sector_2, pphi_efficiency_outb_sector_3, pphi_efficiency_outb_sector_4, pphi_efficiency_outb_sector_5, pphi_efficiency_outb_sector_6, pphi_efficiency_outb_cd, pphi_bins_outb_fd, pphi_bins_outb_cd)
            except:
                df_sim_pi0_1gamma_outb      = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})
            try:
                df_sim_nobkgmerging_outb   =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/dvcs_km15/excl_level_1/restructured_7_nominal/{}.pkl".format(integrated_binnum))
                df_sim_nobkgmerging_outb   = assign_Pthetaefficiency(df_sim_nobkgmerging_outb, ptheta_efficiency_outb_sector_1, ptheta_efficiency_outb_sector_2, ptheta_efficiency_outb_sector_3, ptheta_efficiency_outb_sector_4, ptheta_efficiency_outb_sector_5, ptheta_efficiency_outb_sector_6, ptheta_efficiency_outb_cd, ptheta_bins_outb_fd, ptheta_bins_outb_cd)
                df_sim_nobkgmerging_outb   = assign_Ppefficiency(df_sim_nobkgmerging_outb, pp_efficiency_outb_sector_1, pp_efficiency_outb_sector_2, pp_efficiency_outb_sector_3, pp_efficiency_outb_sector_4, pp_efficiency_outb_sector_5, pp_efficiency_outb_sector_6, pp_efficiency_outb_cd, pp_bins_outb_fd, pp_bins_outb_cd)
                df_sim_nobkgmerging_outb   = assign_Pphiefficiency(df_sim_nobkgmerging_outb, pphi_efficiency_outb_sector_1, pphi_efficiency_outb_sector_2, pphi_efficiency_outb_sector_3, pphi_efficiency_outb_sector_4, pphi_efficiency_outb_sector_5, pphi_efficiency_outb_sector_6, pphi_efficiency_outb_cd, pphi_bins_outb_fd, pphi_bins_outb_cd)
            except:
                df_sim_nobkgmerging_outb      = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})
            try:
                df_sim_bkgmerging_outb     =  pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/dvcs_km15/excl_level_1/restructured_14_bkgmerging/{}.pkl".format(integrated_binnum))
                df_sim_bkgmerging_outb     = assign_Pthetaefficiency(df_sim_bkgmerging_outb, ptheta_efficiency_outb_sector_1, ptheta_efficiency_outb_sector_2, ptheta_efficiency_outb_sector_3, ptheta_efficiency_outb_sector_4, ptheta_efficiency_outb_sector_5, ptheta_efficiency_outb_sector_6, ptheta_efficiency_outb_cd, ptheta_bins_outb_fd, ptheta_bins_outb_cd)
                df_sim_bkgmerging_outb     = assign_Ppefficiency(df_sim_bkgmerging_outb, pp_efficiency_outb_sector_1, pp_efficiency_outb_sector_2, pp_efficiency_outb_sector_3, pp_efficiency_outb_sector_4, pp_efficiency_outb_sector_5, pp_efficiency_outb_sector_6, pp_efficiency_outb_cd, pp_bins_outb_fd, pp_bins_outb_cd)
                df_sim_bkgmerging_outb     = assign_Pphiefficiency(df_sim_bkgmerging_outb, pphi_efficiency_outb_sector_1, pphi_efficiency_outb_sector_2, pphi_efficiency_outb_sector_3, pphi_efficiency_outb_sector_4, pphi_efficiency_outb_sector_5, pphi_efficiency_outb_sector_6, pphi_efficiency_outb_cd, pphi_bins_outb_fd, pphi_bins_outb_cd)
            except:
                df_sim_bkgmerging_outb      = pd.DataFrame({var: [] for var in ['integrated_binnum', 'phi_binnum', 'Pp', 'Ptheta', 'Pphi', 'efficiency']})
            # df_sim_nobkgmerging_outb.loc[:, "efficiency"] = 0

            # bkg_to_nobkg    = np.sum(df_sim_bkgmerging_outb.weights * df_sim_bkgmerging_outb.Pthetaefficiency * df_sim_bkgmerging_outb.Ppefficiency * df_sim_bkgmerging_outb.Pphiefficiency)/np.sum(df_sim_nobkgmerging_outb.weights * df_sim_nobkgmerging_outb.Pthetaefficiency * df_sim_nobkgmerging_outb.Ppefficiency * df_sim_nobkgmerging_outb.Pphiefficiency)
            # eff_bkg_merging = ( 1 + effective_current_outb/50 * ( -1 + bkg_to_nobkg))

            # df_sim_nobkgmerging_outb.loc[:, "efficiency"] = eff_bkg_merging
            # epg_exp_outb_this_integrated_bin           = np.sum(1/df_exp_epg_outb.efficiency)
            # pi0_exp_outb_this_integrated_bin           = np.sum(1/df_exp_pi0_outb.efficiency)
            # pi0_sim_outb_this_integrated_bin           = np.sum(df_sim_pi0_outb.Pthetaefficiency * df_sim_pi0_outb.Ppefficiency * df_sim_pi0_outb.Pphiefficiency)
            # bkg_sim_outb_this_integrated_bin           = np.sum(df_sim_pi0_1gamma_outb.Pthetaefficiency * df_sim_pi0_1gamma_outb.Ppefficiency * df_sim_pi0_1gamma_outb.Pphiefficiency)
            # bkg_exp_outb_this_integrated_bin           = bkg_sim_outb_this_integrated_bin  * pi0_exp_outb_this_integrated_bin/pi0_sim_outb_this_integrated_bin
            # contamination_ratio                       = np.minimum(bkg_exp_outb_this_integrated_bin , epg_exp_outb_this_integrated_bin )/epg_exp_outb_this_integrated_bin 

            # df_exp_epg_outb.loc[:, "contamination"]   = contamination_ratio

            for i, this_rebinned_phi_binnum in enumerate(this_rebinned_phi_binnums):
                df_exp_epg_outb_this_bin = df_exp_epg_outb.loc[(df_exp_epg_outb.integrated_binnum == integrated_binnum) & (df_exp_epg_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                # if len(df_exp_epg_outb_this_bin) < 10:
                if df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == this_rebinned_phi_binnum), "active_bin_outb"].sum() == 0:
                    # print("{} inactive bin between {} and {}".format(integrated_binnum, this_rebinned_phi_binnum, this_rebinned_phi_binnum+this_rebinned_phi_widths[i]))
                    continue
                else:
                    df_exp_pi0_outb_this_bin = df_exp_pi0_outb.loc[(df_exp_pi0_outb.integrated_binnum == integrated_binnum) & (df_exp_pi0_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_pi0_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_pi0_outb_this_bin = df_sim_pi0_outb.loc[(df_sim_pi0_outb.integrated_binnum == integrated_binnum) & (df_sim_pi0_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_pi0_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_pi0_1gamma_outb_this_bin = df_sim_pi0_1gamma_outb.loc[(df_sim_pi0_1gamma_outb.integrated_binnum == integrated_binnum) & (df_sim_pi0_1gamma_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_pi0_1gamma_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_nobkgmerging_outb_this_bin = df_sim_nobkgmerging_outb.loc[(df_sim_nobkgmerging_outb.integrated_binnum == integrated_binnum) & (df_sim_nobkgmerging_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_nobkgmerging_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]
                    df_sim_bkgmerging_outb_this_bin = df_sim_bkgmerging_outb.loc[(df_sim_bkgmerging_outb.integrated_binnum == integrated_binnum) & (df_sim_bkgmerging_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_bkgmerging_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), :]

                    bkg_to_nobkg    = np.sum(df_sim_bkgmerging_outb_this_bin.weights * df_sim_bkgmerging_outb_this_bin.Pthetaefficiency * df_sim_bkgmerging_outb_this_bin.Ppefficiency * df_sim_bkgmerging_outb_this_bin.Pphiefficiency)/np.sum(df_sim_nobkgmerging_outb_this_bin.weights * df_sim_nobkgmerging_outb_this_bin.Pthetaefficiency * df_sim_nobkgmerging_outb_this_bin.Ppefficiency * df_sim_nobkgmerging_outb_this_bin.Pphiefficiency)
                    eff_bkg_merging = ( 1 + effective_current_outb/50 * ( -1 + bkg_to_nobkg))
                    
                    epg_exp_outb_this_bin               = np.sum(1/df_exp_epg_outb_this_bin.efficiency)
                    try:
                        pi0_exp_outb_this_integrated_bin    = np.sum(1/df_exp_pi0_outb.efficiency)
                    except:
                        pi0_exp_outb_this_integrated_bin    = 0
                    try:
                        pi0_sim_outb_this_integrated_bin    = np.sum(df_sim_pi0_outb.Pthetaefficiency)
                    except:
                        pi0_sim_outb_this_integrated_bin    = 0
                    try:
                        bkg_sim_outb_this_bin               = np.sum(df_sim_pi0_1gamma_outb_this_bin.Pthetaefficiency)
                    except:
                        bkg_sim_outb_this_bin               = 0
                    try:
                        bkg_exp_outb_this_bin               = bkg_sim_outb_this_bin * pi0_exp_outb_this_integrated_bin/pi0_sim_outb_this_integrated_bin
                    except:
                        bkg_exp_outb_this_bin               = 0
                    contamination_ratio                = np.minimum(bkg_exp_outb_this_bin, epg_exp_outb_this_bin)/epg_exp_outb_this_bin

                    df_exp_epg_outb.loc[(df_exp_epg_outb.integrated_binnum == integrated_binnum) & (df_exp_epg_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "contamination"]   = contamination_ratio
                    df_exp_epg_outb.loc[(df_exp_epg_outb.integrated_binnum == integrated_binnum) & (df_exp_epg_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "phi_rebinnum"]    = this_rebinned_phi_binnum
                    df_exp_epg_outb.loc[(df_exp_epg_outb.integrated_binnum == integrated_binnum) & (df_exp_epg_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_exp_epg_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "phi_rebinwidth"]  = this_rebinned_phi_widths[i]
                    df_sim_nobkgmerging_outb.loc[(df_sim_nobkgmerging_outb.integrated_binnum == integrated_binnum) & (df_sim_nobkgmerging_outb.phi_binnum >= this_rebinned_phi_binnum) & (df_sim_nobkgmerging_outb.phi_binnum < this_rebinned_phi_binnum + this_rebinned_phi_widths[i]), "efficiency"] = eff_bkg_merging

            df_exp_epg_outb.loc[:, "signal"]          = (1 - df_exp_epg_outb.contamination)/df_exp_epg_outb.efficiency
            df_sim_nobkgmerging_outb.loc[:, "signal"] = df_sim_nobkgmerging_outb.weights * df_sim_nobkgmerging_outb.efficiency * df_sim_nobkgmerging_outb.Pthetaefficiency  * df_sim_nobkgmerging_outb.Ppefficiency * df_sim_nobkgmerging_outb.Pphiefficiency
            
            df_exp_epg_outbs.append(df_exp_epg_outb)
            df_sim_dvcs_outbs.append(df_sim_nobkgmerging_outb)

        df_exp_epg_outbs  = pd.concat(df_exp_epg_outbs)
        df_sim_dvcs_outbs = pd.concat(df_sim_dvcs_outbs)

        # FD sector 1
        pphi_exp_sector_1, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 1, "Pphi"], weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 1, "signal"], bins = pphi_bins_outb_fd)
        pphi_sim_sector_1, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 1, "Pphi"], weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 1, "signal"], bins = pphi_bins_outb_fd)
        pphi_efficiency_outb_sector_1 = pphi_efficiency_outb_sector_1*divideHist(pphi_exp_sector_1, pphi_sim_sector_1)

        # FD sector 2
        pphi_exp_sector_2, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 2, "Pphi"], weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 2, "signal"], bins = pphi_bins_outb_fd)
        pphi_sim_sector_2, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 2, "Pphi"], weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 2, "signal"], bins = pphi_bins_outb_fd)
        pphi_efficiency_outb_sector_2 = pphi_efficiency_outb_sector_2*divideHist(pphi_exp_sector_2, pphi_sim_sector_2)

        # FD sector 3
        pphi_exp_sector_3, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 3, "Pphi"], weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 3, "signal"], bins = pphi_bins_outb_fd)
        pphi_sim_sector_3, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 3, "Pphi"], weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 3, "signal"], bins = pphi_bins_outb_fd)
        pphi_efficiency_outb_sector_3 = pphi_efficiency_outb_sector_3*divideHist(pphi_exp_sector_3, pphi_sim_sector_3)

        # FD sector 4
        pphi_exp_sector_4, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 4, "Pphi"], weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 4, "signal"], bins = pphi_bins_outb_fd)
        pphi_sim_sector_4, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 4, "Pphi"], weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 4, "signal"], bins = pphi_bins_outb_fd)
        pphi_efficiency_outb_sector_4 = pphi_efficiency_outb_sector_4*divideHist(pphi_exp_sector_4, pphi_sim_sector_4)

        # FD sector 5
        pphi_exp_sector_5, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 5, "Pphi"], weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 5, "signal"], bins = pphi_bins_outb_fd)
        pphi_sim_sector_5, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 5, "Pphi"], weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 5, "signal"], bins = pphi_bins_outb_fd)
        pphi_efficiency_outb_sector_5 = pphi_efficiency_outb_sector_5*divideHist(pphi_exp_sector_5, pphi_sim_sector_5)

        # FD sector 6
        pphi_exp_sector_6, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 6, "Pphi"], weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector == 6, "signal"], bins = pphi_bins_outb_fd)
        pphi_sim_sector_6, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 6, "Pphi"], weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector == 6, "signal"], bins = pphi_bins_outb_fd)
        pphi_efficiency_outb_sector_6 = pphi_efficiency_outb_sector_6*divideHist(pphi_exp_sector_6, pphi_sim_sector_6)

        # CD
        pphi_exp_cd, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector > 7, "Pphi"], weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector > 7, "signal"], bins = pphi_bins_outb_cd)
        pphi_sim_cd, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector > 7, "Pphi"], weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector > 7, "signal"], bins = pphi_bins_outb_cd)
        pphi_efficiency_outb_cd = pphi_efficiency_outb_cd*divideHist(pphi_exp_cd, pphi_sim_cd)

        pphi_efficiency_outb_sector_1[pphi_efficiency_outb_sector_1<0.5] = 0.5
        pphi_efficiency_outb_sector_2[pphi_efficiency_outb_sector_2<0.5] = 0.5
        pphi_efficiency_outb_sector_3[pphi_efficiency_outb_sector_3<0.5] = 0.5
        pphi_efficiency_outb_sector_4[pphi_efficiency_outb_sector_4<0.5] = 0.5
        pphi_efficiency_outb_sector_5[pphi_efficiency_outb_sector_5<0.5] = 0.5
        pphi_efficiency_outb_sector_6[pphi_efficiency_outb_sector_6<0.5] = 0.5
        pphi_efficiency_outb_cd      [pphi_efficiency_outb_cd      <0.5] = 0.5

        pphi_efficiency_outb_sector_1[pphi_efficiency_outb_sector_1>1.5] = 1.5
        pphi_efficiency_outb_sector_2[pphi_efficiency_outb_sector_2>1.5] = 1.5
        pphi_efficiency_outb_sector_3[pphi_efficiency_outb_sector_3>1.5] = 1.5
        pphi_efficiency_outb_sector_4[pphi_efficiency_outb_sector_4>1.5] = 1.5
        pphi_efficiency_outb_sector_5[pphi_efficiency_outb_sector_5>1.5] = 1.5
        pphi_efficiency_outb_sector_6[pphi_efficiency_outb_sector_6>1.5] = 1.5
        pphi_efficiency_outb_cd      [pphi_efficiency_outb_cd      >1.5] = 1.5

    print("outbending", outbending_trial, "outb sector 1", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_1], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_1], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_1])
    print("outbending", outbending_trial, "outb sector 2", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_2], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_2], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_2])
    print("outbending", outbending_trial, "outb sector 3", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_3], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_3], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_3])
    print("outbending", outbending_trial, "outb sector 4", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_4], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_4], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_4])
    print("outbending", outbending_trial, "outb sector 5", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_5], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_5], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_5])
    print("outbending", outbending_trial, "outb sector 6", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_6], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_6], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_6])
    print("outbending", outbending_trial, "outb CD"      , ["{:.3f}".format(i) for i in ptheta_efficiency_outb_cd      ], ["{:.3f}".format(i) for i in pp_efficiency_outb_cd      ], ["{:.3f}".format(i) for i in pphi_efficiency_outb_cd      ])

print("inb sector 1", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_1], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_1], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_1])
print("inb sector 2", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_2], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_2], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_2])
print("inb sector 3", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_3], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_3], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_3])
print("inb sector 4", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_4], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_4], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_4])
print("inb sector 5", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_5], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_5], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_5])
print("inb sector 6", ["{:.3f}".format(i) for i in ptheta_efficiency_inb_sector_6], ["{:.3f}".format(i) for i in pp_efficiency_inb_sector_6], ["{:.3f}".format(i) for i in pphi_efficiency_inb_sector_6])
print("inb CD"      , ["{:.3f}".format(i) for i in ptheta_efficiency_inb_cd      ], ["{:.3f}".format(i) for i in pp_efficiency_inb_cd      ], ["{:.3f}".format(i) for i in pphi_efficiency_inb_cd      ])

print("outb sector 1", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_1], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_1], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_1])
print("outb sector 2", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_2], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_2], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_2])
print("outb sector 3", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_3], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_3], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_3])
print("outb sector 4", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_4], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_4], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_4])
print("outb sector 5", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_5], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_5], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_5])
print("outb sector 6", ["{:.3f}".format(i) for i in ptheta_efficiency_outb_sector_6], ["{:.3f}".format(i) for i in pp_efficiency_outb_sector_6], ["{:.3f}".format(i) for i in pphi_efficiency_outb_sector_6])
print("outb CD"      , ["{:.3f}".format(i) for i in ptheta_efficiency_outb_cd      ], ["{:.3f}".format(i) for i in pp_efficiency_outb_cd      ], ["{:.3f}".format(i) for i in pphi_efficiency_outb_cd      ])


