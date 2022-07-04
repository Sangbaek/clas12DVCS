from utils.const import *
from utils.physics import *

def electronFiducial(df_electronRec, pol = "inbending", mc = False):
	df_electronRec.loc[:, "EFid"] = 1

	# #PCAL dead wires
	# exclusion1_1 = (df_electronRec.EcalW1 > 74) & (df_electronRec.EcalW1 < 79.8)
	# exclusion1_2 = (df_electronRec.EcalW1 > 83.6) & (df_electronRec.EcalW1 < 92.2)
	# exclusion1_3 = (df_electronRec.EcalW1 > 212.5) & (df_electronRec.EcalW1 < 230)
	# exclusion1 = exclusion1_1 | exclusion1_2 | exclusion1_3
	# df_electronRec.loc[(df_electronRec.Esector == 1) & exclusion1, "EFid"] = 0
	# exclusion2_1 = (df_electronRec.EcalW1 < 14)
	# exclusion2_2 = (df_electronRec.EcalU1 > 111.2) & (df_electronRec.EcalU1 < 119.3)
	# exclusion2_3 = (df_electronRec.EcalV1 > 113) & (df_electronRec.EcalV1 < 118.7)
	# exclusion2 = exclusion2_1 | exclusion2_2 | exclusion2_3
	# df_electronRec.loc[(df_electronRec.Esector == 2) & exclusion2, "EFid"] = 0
	# exclusion3 = df_electronRec.EcalW1 < 14
	# df_electronRec.loc[(df_electronRec.Esector == 3) & exclusion3, "EFid"] = 0
	# exclusion4_1 = (df_electronRec.EcalV1 < 14)
	# exclusion4_2 = (df_electronRec.EcalV1 > 229.4) & (df_electronRec.EcalV1 < 240.7)
	# exclusion4_3 = (df_electronRec.EcalW1 > 135) & (df_electronRec.EcalW1 < 150)
	# exclusion4 = exclusion4_1 | exclusion4_2 | exclusion4_3
	# df_electronRec.loc[(df_electronRec.Esector == 4) & exclusion4, "EFid"] = 0
	# exclusion6 = (df_electronRec.EcalW1 > 170) & (df_electronRec.EcalW1 < 192)
	# df_electronRec.loc[(df_electronRec.Esector == 6) & exclusion6, "EFid"] = 0

	# passElectronTrackQualityCut (pass)
	sector_cond = [df_electronRec.Esector ==1, df_electronRec.Esector ==2, df_electronRec.Esector ==3, df_electronRec.Esector ==4, df_electronRec.Esector ==5, df_electronRec.Esector ==6]

	# passElectronSamplingFractionCut
	ecal_e_sampl_mu_0 = np.select(sector_cond, ecal_e_sampl_mu[0])
	ecal_e_sampl_mu_1 = np.select(sector_cond, ecal_e_sampl_mu[1])
	ecal_e_sampl_mu_2 = np.select(sector_cond, ecal_e_sampl_mu[2])
	ecal_e_sampl_sigm_0 = np.select(sector_cond, ecal_e_sampl_sigm[0])
	ecal_e_sampl_sigm_1 = np.select(sector_cond, ecal_e_sampl_sigm[1])
	ecal_e_sampl_sigm_2 = np.select(sector_cond, ecal_e_sampl_sigm[2])

	if mc:
		ecal_e_sampl_mu_0 = np.select(sector_cond, ecal_e_sampl_mu_mc[0])
		ecal_e_sampl_mu_1 = np.select(sector_cond, ecal_e_sampl_mu_mc[1])
		ecal_e_sampl_mu_2 = np.select(sector_cond, ecal_e_sampl_mu_mc[2])
		ecal_e_sampl_sigm_0 = np.select(sector_cond, ecal_e_sampl_sigm_mc[0])
		ecal_e_sampl_sigm_1 = np.select(sector_cond, ecal_e_sampl_sigm_mc[1])
		ecal_e_sampl_sigm_2 = np.select(sector_cond, ecal_e_sampl_sigm_mc[2])
	mean = ecal_e_sampl_mu_0 + ecal_e_sampl_mu_1/1000*pow(df_electronRec.Ep-ecal_e_sampl_mu_2,2)
	sigma = ecal_e_sampl_sigm_0 + ecal_e_sampl_sigm_1/(10*(df_electronRec.Ep-ecal_e_sampl_sigm_2))
	df_electronRec.loc[df_electronRec.ESamplFrac < mean - e_sampl_sigma_range*sigma, "EFid"]  = 0
	df_electronRec.loc[df_electronRec.ESamplFrac > mean + e_sampl_sigma_range*sigma, "EFid"]  = 0


	#passElectronNpheCut
	df_electronRec.loc[df_electronRec.Enphe < min_nphe, "EFid"] = 0

	#passElectronVertexCut
	if pol == 'inbending':
		min_vz = vz_min_inb
		max_vz = vz_max_inb
	if pol == 'outbending':
		min_vz = vz_min_outb
		max_vz = vz_max_outb
	df_electronRec.loc[df_electronRec.Evz < min_vz, "EFid"] = 0
	df_electronRec.loc[df_electronRec.Evz > max_vz, "EFid"] = 0

	# passElectronPCALFiducialCut
	df_electronRec.loc[df_electronRec.EcalV1 < min_v, "EFid"] = 0
	df_electronRec.loc[df_electronRec.EcalW1 < min_w, "EFid"] = 0

	#passElectronPCALEdepCut
	df_electronRec.loc[df_electronRec.Eedep1 < min_pcal_dep, "EFid"] = 0

	#passElectronDCR1
	if pol == 'inbending':
		minparams = e_dc_minparams_in
		maxparams = e_dc_maxparams_in
	if pol == 'outbending':
		minparams = e_dc_minparams_out
		maxparams = e_dc_maxparams_out

	dcsec = determineSector(df_electronRec.EDc1Hitx, df_electronRec.EDc1Hity)
	x_rot, y_rot = rotateDCHitPosition(df_electronRec.EDc1Hitx, df_electronRec.EDc1Hity, dcsec)
	calc_min, calc_max = e_DC_fiducial_cut_XY(x_rot, dcsec, 0, minparams, maxparams)
	df_electronRec.loc[y_rot<calc_min, "EFid"] = 0
	df_electronRec.loc[y_rot>calc_max, "EFid"] = 0
	#passElectronDCR2
	dcsec = determineSector(df_electronRec.EDc2Hitx, df_electronRec.EDc2Hity)
	x_rot, y_rot = rotateDCHitPosition(df_electronRec.EDc2Hitx, df_electronRec.EDc2Hity, dcsec)
	calc_min, calc_max = e_DC_fiducial_cut_XY(x_rot, dcsec, 1, minparams, maxparams)
	df_electronRec.loc[y_rot<calc_min, "EFid"] = 0
	df_electronRec.loc[y_rot>calc_max, "EFid"] = 0

	#passElectronDCR3
	dcsec = determineSector(df_electronRec.EDc3Hitx, df_electronRec.EDc3Hity)
	x_rot, y_rot = rotateDCHitPosition(df_electronRec.EDc3Hitx, df_electronRec.EDc3Hity, dcsec)
	calc_min, calc_max = e_DC_fiducial_cut_XY(x_rot, dcsec, 2, minparams, maxparams)
	df_electronRec.loc[y_rot<calc_min, "EFid"] = 0
	df_electronRec.loc[y_rot>calc_max, "EFid"] = 0

	# #passElectronAntiPionCut
	# df_electronRec.loc[-df_electronRec.Edep1/df_electronRec.Ep + anti_pion_threshold > df_electronRec.Edep2/event.p[index], "EFid"] = 0
	return df_electronRec

def gammaFiducial(df_gammaRec):
	#passGammaPCALFiducialCut
	df_gammaRec.loc[(df_gammaRec.GcalV1 < min_v) & (df_gammaRec.Gsector<7), "GFid"] = 0
	df_gammaRec.loc[(df_gammaRec.GcalW1 < min_w) & (df_gammaRec.Gsector<7), "GFid"] = 0
	#passGammaBetaCut
	df_gammaRec.loc[df_gammaRec.Gbeta < min_Gbeta, "GFid"] = 0
	df_gammaRec.loc[df_gammaRec.Gbeta > max_Gbeta, "GFid"] = 0

	return df_gammaRec

def protonFiducial(df_protonRec, pol = 'inbending'):

	dcsec = determineSector(df_protonRec.PDc1Hitx, df_protonRec.PDc1Hity)
	if pol == 'inbending':
		minparams = p_dc_minparams_in
		maxparams = p_dc_maxparams_in

		theta_DC, phi_DC = thetaphifromhit(df_protonRec.PDc1Hitx, df_protonRec.PDc1Hity, df_protonRec.PDc1Hitz)
		phi_DC_min, phi_DC_max = p_DC_fiducial_cut_thetaphi(theta_DC, dcsec, 0, minparams, maxparams)
		df_protonRec.loc[(phi_DC<phi_DC_min) & (df_protonRec.Psector<7), "PFid"] = 0
		df_protonRec.loc[(phi_DC>phi_DC_max) & (df_protonRec.Psector<7), "PFid"] = 0

		theta_DC, phi_DC = thetaphifromhit(df_protonRec.PDc2Hitx, df_protonRec.PDc2Hity, df_protonRec.PDc2Hitz)
		phi_DC_min, phi_DC_max = p_DC_fiducial_cut_thetaphi(theta_DC, dcsec, 1, minparams, maxparams)
		df_protonRec.loc[(phi_DC<phi_DC_min) & (df_protonRec.Psector<7), "PFid"] = 0
		df_protonRec.loc[(phi_DC>phi_DC_max) & (df_protonRec.Psector<7), "PFid"] = 0

		theta_DC, phi_DC = thetaphifromhit(df_protonRec.PDc3Hitx, df_protonRec.PDc3Hity, df_protonRec.PDc3Hitz)
		phi_DC_min, phi_DC_max = p_DC_fiducial_cut_thetaphi(theta_DC, dcsec, 2, minparams, maxparams)
		df_protonRec.loc[(phi_DC<phi_DC_min) & (df_protonRec.Psector<7), "PFid"] = 0
		df_protonRec.loc[(phi_DC>phi_DC_max) & (df_protonRec.Psector<7), "PFid"] = 0

	if pol == 'outbending':
		minparams = p_dc_minparams_out
		maxparams = p_dc_maxparams_out

		dcsec = determineSector(df_protonRec.PDc1Hitx, df_protonRec.PDc1Hity)
		x_rot, y_rot = rotateDCHitPosition(df_protonRec.PDc1Hitx, df_protonRec.PDc1Hity, dcsec)
		calc_min, calc_max = p_DC_fiducial_cut_XY(x_rot, dcsec, 0, minparams, maxparams)
		df_protonRec.loc[(y_rot<calc_min) & (df_protonRec.Psector<7), "PFid"] = 0
		df_protonRec.loc[(y_rot>calc_max) & (df_protonRec.Psector<7), "PFid"] = 0
		#passElectronDCR2
		dcsec = determineSector(df_protonRec.PDc2Hitx, df_protonRec.PDc2Hity)
		x_rot, y_rot = rotateDCHitPosition(df_protonRec.PDc2Hitx, df_protonRec.PDc2Hity, dcsec)
		calc_min, calc_max = p_DC_fiducial_cut_XY(x_rot, dcsec, 1, minparams, maxparams)
		df_protonRec.loc[(y_rot<calc_min) & (df_protonRec.Psector<7), "PFid"] = 0
		df_protonRec.loc[(y_rot>calc_max) & (df_protonRec.Psector<7), "PFid"] = 0

		#passElectronDCR3
		dcsec = determineSector(df_protonRec.PDc3Hitx, df_protonRec.PDc3Hity)
		x_rot, y_rot = rotateDCHitPosition(df_protonRec.PDc3Hitx, df_protonRec.PDc3Hity, dcsec)
		calc_min, calc_max = p_DC_fiducial_cut_XY(x_rot, dcsec, 2, minparams, maxparams)
		df_protonRec.loc[(y_rot<calc_min) & (df_protonRec.Psector<7), "PFid"] = 0
		df_protonRec.loc[(y_rot>calc_max) & (df_protonRec.Psector<7), "PFid"] = 0

	return df_protonRec