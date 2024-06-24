from utils.const import *
from utils.physics import *
from copy import copy

def assign_efficiency(df_Rec, mc = False):
	if not mc:
		# HTCC efficiency map
		htcc_eff_map = np.loadtxt("/work/clas12/sangbaek/Inclusive/HTCCEfficiencyData.dat").reshape(250,250)
		EhtccXBin = (df_Rec.loc[:, "EhtcctrajX"] + 125 ).astype(int).to_numpy()
		EhtccYBin = (df_Rec.loc[:, "EhtcctrajY"] + 125 ).astype(int).to_numpy()
		EhtccEfficiency = []
		for i in range(len(EhtccXBin)):
			EhtccEfficiency.append(htcc_eff_map[EhtccXBin[i], EhtccYBin[i]])
		df_Rec.loc[:, "EhtccEfficiency"]  = EhtccEfficiency
		# df_Rec.loc[df_Rec.EhtccEfficiency < 0.95, "EFid"] = 0

		df_Rec = copy(df_Rec)
		df_Rec.loc[:, "EFtof1bEfficiency"] = 1
		df_Rec.loc[(df_Rec.Esector==6) & (df_Rec.EFtof1bComponent>=33) & (df_Rec.EFtof1bComponent<=48), "EFtof1bEfficiency"]= 1/1.013
		df_Rec.loc[:, "PFtof1bEfficiency"] = 1
		df_Rec.loc[(df_Rec.PFtof1bSector==6) & (df_Rec.PFtof1bComponent>=33) & (df_Rec.PFtof1bComponent<=48), "PFtof1bEfficiency"]= 1/1.013
		if "weight" in df_Rec.columns:
			df_Rec.loc[:, "weight"] = df_Rec.weight * df_Rec.EhtccEfficiency * df_Rec.EFtof1bEfficiency * df_Rec.PFtof1bEfficiency
		else:
			df_Rec.loc[:, "weight"] = df_Rec.EhtccEfficiency * df_Rec.EFtof1bEfficiency * df_Rec.PFtof1bEfficiency
	return df_Rec

def electronFiducial(df_electronRec, mc = False, fidlevel = 'mid'):
	df_electronRec = copy(df_electronRec)
	# following inclusive analysis note
	df_electronRec.loc[:, "EFid"] = 1
	# D. nphe cut
	df_electronRec.loc[df_electronRec.Enphe <= min_nphe, "EFid"] = 0

	# E. vz cut
	if fidlevel == 'mid':
		df_electronRec.loc[df_electronRec.Evz < -8, "EFid"] = 0
		df_electronRec.loc[df_electronRec.Evz >  2, "EFid"] = 0
	elif fidlevel == 'loose':
		df_electronRec.loc[df_electronRec.Evz < -8.5, "EFid"] = 0
		df_electronRec.loc[df_electronRec.Evz >  2.5, "EFid"] = 0
	elif fidlevel == 'tight':
		df_electronRec.loc[df_electronRec.Evz < -7.5, "EFid"] = 0
		df_electronRec.loc[df_electronRec.Evz >  1.5, "EFid"] = 0
	else:
		print("check fidlevel {}".format(fidlevel))
	# # F. Minimum PCAL energy Threshold cut # removed at the latest inclusive analysis
	# if fidlevel == 'mid':
	# 	df_electronRec.loc[df_electronRec.Eedep1 < 0.07, "EFid"] = 0
	# elif fidlevel == 'loose':
	# 	df_electronRec.loc[df_electronRec.Eedep1 < 0.06, "EFid"] = 0
	# elif fidlevel == 'tight':
	# 	df_electronRec.loc[df_electronRec.Eedep1 < 0.08, "EFid"] = 0
	# else:
	# 	print("check fidlevel {}".format(fidlevel))
	# F. DC Fiducial Cuts
	if fidlevel == 'mid':
		adjustment_layer1 = 0
		adjustment_layer2 = 0
		adjustment_layer3 = 0
	elif fidlevel == 'loose':
		adjustment_layer1 = 0.6*1
		adjustment_layer2 = 0.6*2
		adjustment_layer3 = 0.6*3
	elif fidlevel == 'tight':
		adjustment_layer1 = -0.6*1
		adjustment_layer2 = -0.6*2
		adjustment_layer3 = -0.6*3
	else:
		print("check fidlevel {}".format(fidlevel))

	dcsec_l1 = determineSector(df_electronRec.EDc1Hitx, df_electronRec.EDc1Hity)
	x_rot_l1, y_rot_l1 = rotateDCHitPosition(df_electronRec.EDc1Hitx, df_electronRec.EDc1Hity, dcsec_l1)
	x_rot_l1, _ = rotateDCHitPosition_alongY(x_rot_l1, df_electronRec.EDc1Hitz)
	calc_min_l1 = -0.50 * (x_rot_l1 + 72 + adjustment_layer1)
	calc_max_l1 =  0.50 * (x_rot_l1 + 72 + adjustment_layer1)
	df_electronRec.loc[y_rot_l1 < calc_min_l1, "EFid"] = 0
	df_electronRec.loc[y_rot_l1 > calc_max_l1, "EFid"] = 0

	dcsec_l2 = determineSector(df_electronRec.EDc2Hitx, df_electronRec.EDc2Hity)
	x_rot_l2, y_rot_l2 = rotateDCHitPosition(df_electronRec.EDc2Hitx, df_electronRec.EDc2Hity, dcsec_l2)
	x_rot_l2, _ = rotateDCHitPosition_alongY(x_rot_l2, df_electronRec.EDc2Hitz)
	calc_min_l2 = -0.505 * (x_rot_l2 + 114 + adjustment_layer2)
	calc_max_l2 =  0.505 * (x_rot_l2 + 114 + adjustment_layer2)
	df_electronRec.loc[y_rot_l2 < calc_min_l2, "EFid"] = 0
	df_electronRec.loc[y_rot_l2 > calc_max_l2, "EFid"] = 0

	dcsec_l3 = determineSector(df_electronRec.EDc3Hitx, df_electronRec.EDc3Hity)
	x_rot_l3, y_rot_l3 = rotateDCHitPosition(df_electronRec.EDc3Hitx, df_electronRec.EDc3Hity, dcsec_l3)
	x_rot_l3, _ = rotateDCHitPosition_alongY(x_rot_l3, df_electronRec.EDc3Hitz)
	calc_min_l3 = -0.495 * (x_rot_l3 + 180 + adjustment_layer3)
	calc_max_l3 =  0.495 * (x_rot_l3 + 180 + adjustment_layer3)
	df_electronRec.loc[y_rot_l3 < calc_min_l3, "EFid"] = 0
	df_electronRec.loc[y_rot_l3 > calc_max_l3, "EFid"] = 0
	# G. PCAL Fid Cuts
	if fidlevel == 'mid':
		df_electronRec.loc[df_electronRec.EcalV1<19, "EFid"] = 0
		df_electronRec.loc[df_electronRec.EcalW1<19, "EFid"] = 0
		df_electronRec.loc[df_electronRec.EcalU1>395, "EFid"] = 0
	elif fidlevel == 'loose':
		df_electronRec.loc[df_electronRec.EcalV1<19-2.5, "EFid"] = 0
		df_electronRec.loc[df_electronRec.EcalW1<19-2.5, "EFid"] = 0
		df_electronRec.loc[df_electronRec.EcalU1>395+2.5, "EFid"] = 0
	elif fidlevel == 'tight':
		df_electronRec.loc[df_electronRec.EcalV1<19+2.5, "EFid"] = 0
		df_electronRec.loc[df_electronRec.EcalW1<19+2.5, "EFid"] = 0
		df_electronRec.loc[df_electronRec.EcalU1>395-2.5, "EFid"] = 0
	else:
		print("check fidlevel {}".format(fidlevel))
	# H. ECAL SF Cut
	A = [0.286, 0.280, 0.275, 0.273, 0.271, 0.276]
	B = [-0.040, -0.038, -0.034, -0.033, -0.032, -0.034]
	C = [-0.0030, -0.0012, -0.0014, -0.0007, 0.0005, -0.0014]
	A_sim = [0.29]*6
	B_sim = [-0.040]*6
	C_sim = [-0.0029]*6
	D = [0.017, 0.019, 0.017, 0.0157, 0.016, 0.017]
	E = [-0.0012, -0.003, -0.002, 0.0003, -0.00135, -0.002]
	F = [-0.0012, -0.00135, -0.00129, -0.0013, -0.001, -0.001]
	D_sim = [0.015]*6
	E_sim = [-0.00053]*6
	F_sim = [-0.0014]*6

	pcal_sf_mu    = [A, B, C]
	pcal_sf_sigma = [D, E, F]
	if mc:
		pcal_sf_mu    = [A_sim, B_sim, C_sim]
		pcal_sf_sigma = [D_sim, E_sim, F_sim]

	sector_cond = [df_electronRec.Esector ==1, df_electronRec.Esector ==2, df_electronRec.Esector ==3, df_electronRec.Esector ==4, df_electronRec.Esector ==5, df_electronRec.Esector ==6]

	ecal_e_sampl_mu_0 = np.select(sector_cond, pcal_sf_mu[0])
	ecal_e_sampl_mu_1 = np.select(sector_cond, pcal_sf_mu[1])
	ecal_e_sampl_mu_2 = np.select(sector_cond, pcal_sf_mu[2])
	ecal_e_sampl_sigm_0 = np.select(sector_cond, pcal_sf_sigma[0])
	ecal_e_sampl_sigm_1 = np.select(sector_cond, pcal_sf_sigma[1])
	ecal_e_sampl_sigm_2 = np.select(sector_cond, pcal_sf_sigma[2])

	mean =  ecal_e_sampl_mu_0   + ecal_e_sampl_mu_1  /df_electronRec.Eedep + ecal_e_sampl_mu_2  /df_electronRec.Eedep/df_electronRec.Eedep
	sigma = ecal_e_sampl_sigm_0 + ecal_e_sampl_sigm_1/df_electronRec.Eedep + ecal_e_sampl_sigm_2/df_electronRec.Eedep/df_electronRec.Eedep

	if fidlevel == 'mid':
		df_electronRec.loc[df_electronRec.ESamplFrac < mean - 3.5*sigma, "EFid"]  = 0
		# df_electronRec.loc[df_electronRec.ESamplFrac > mean + 3.5*sigma, "EFid"]  = 0
	elif fidlevel == 'loose':
		df_electronRec.loc[df_electronRec.ESamplFrac < mean - (3.5+0.5)*sigma, "EFid"]  = 0
		# df_electronRec.loc[df_electronRec.ESamplFrac > mean + (3.5+0.5)*sigma, "EFid"]  = 0
	elif fidlevel == 'tight':
		df_electronRec.loc[df_electronRec.ESamplFrac < mean - (3.5-0.5)*sigma, "EFid"]  = 0
		# df_electronRec.loc[df_electronRec.ESamplFrac > mean + (3.5-0.5)*sigma, "EFid"]  = 0
	else:
		print("check fidlevel {}".format(fidlevel))
	# I. Pion Separtaion Cut
	eleFidCut = df_electronRec.loc[:, ["Ep", "Esector", "Eedep1", "Eedep2"]]
	eleFidCut.loc[:, "a"] = 0
	eleFidCut.loc[:, "b"] = 0

	if not mc:
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep < 2)       , "a"] =  0.201232 
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 2) &  3 , "a"] =  0.207553 
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 3) &  4 , "a"] =  0.213829 
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 4) &  5 , "a"] =  0.217145 
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 5) &  6 , "a"] =  0.220458 
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 6) &  7 , "a"] =  0.22359  
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 7) &  8 , "a"] =  0.226479 
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 8) &  9 , "a"] =  0.226668 
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 9)      , "a"] =  0.225488 

		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep < 2)       , "b"] = -1.00947
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 2) &  3 , "b"] = -1.07716
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 3) &  4 , "b"] = -1.10165
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 4) &  5 , "b"] = -1.12186
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 5) &  6 , "b"] = -1.16741
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 6) &  7 , "b"] = -1.22825
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 7) &  8 , "b"] = -1.30402
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 8) &  9 , "b"] = -1.37011
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 9)      , "b"] = -1.45688

		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep < 2)       , "a"] =  0.197423
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 2) &  3 , "a"] =  0.20437 
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 3) &  4 , "a"] =  0.215049
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 4) &  5 , "a"] =  0.218486
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 5) &  6 , "a"] =  0.218645
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 6) &  7 , "a"] =  0.219856
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 7) &  8 , "a"] =  0.219286
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 8) &  9 , "a"] =  0.219087
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 9)      , "a"] =  0.221795

		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep < 2)       , "b"] = -0.975077
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 2) &  3 , "b"] = -1.01246 
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 3) &  4 , "b"] = -1.08631 
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 4) &  5 , "b"] = -1.1006  
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 5) &  6 , "b"] = -1.10554 
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 6) &  7 , "b"] = -1.13945 
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 7) &  8 , "b"] = -1.16559 
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 8) &  9 , "b"] = -1.21059 
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 9)      , "b"] = -1.30088 

		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep < 2)       , "a"] =  0.197381
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 2) &  3 , "a"] =  0.209949
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 3) &  4 , "a"] =  0.215857
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 4) &  5 , "a"] =  0.21988 
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 5) &  6 , "a"] =  0.220504
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 6) &  7 , "a"] =  0.226815
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 7) &  8 , "a"] =  0.22881 
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 8) &  9 , "a"] =  0.226924
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 9)      , "a"] =  0.21997 

		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep < 2)       , "b"] = -0.988142
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 2) &  3 , "b"] = -1.10285 
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 3) &  4 , "b"] = -1.1188  
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 4) &  5 , "b"] = -1.14603 
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 5) &  6 , "b"] = -1.15965 
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 6) &  7 , "b"] = -1.26374 
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 7) &  8 , "b"] = -1.33568 
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 8) &  9 , "b"] = -1.38956 
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 9)      , "b"] = -1.4082  

		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep < 2)       , "a"] =  0.18777 
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 2) &  3 , "a"] =  0.198804
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 3) &  4 , "a"] =  0.209816
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 4) &  5 , "a"] =  0.215048
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 5) &  6 , "a"] =  0.218401
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 6) &  7 , "a"] =  0.221764
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 7) &  8 , "a"] =  0.225656
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 8) &  9 , "a"] =  0.228833
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 9)      , "a"] =  0.228161

		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep < 2)       , "b"] = -0.892892
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 2) &  3 , "b"] = -1.00481 
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 3) &  4 , "b"] = -1.0906  
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 4) &  5 , "b"] = -1.11891 
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 5) &  6 , "b"] = -1.13694 
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 6) &  7 , "b"] = -1.16584 
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 7) &  8 , "b"] = -1.2161  
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 8) &  9 , "b"] = -1.28444 
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 9)      , "b"] = -1.34519 

		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep < 2)       , "a"] = 0.197711
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 2) &  3 , "a"] = 0.208762
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 3) &  4 , "a"] = 0.218481
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 4) &  5 , "a"] = 0.221863
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 5) &  6 , "a"] = 0.222802
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 6) &  7 , "a"] = 0.222474
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 7) &  8 , "a"] = 0.22175 
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 8) &  9 , "a"] = 0.21954 
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 9)      , "a"] = 0.216218

		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep < 2)       , "b"] = -0.97406
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 2) &  3 , "b"] = -1.07982
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 3) &  4 , "b"] = -1.15156
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 4) &  5 , "b"] = -1.16233
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 5) &  6 , "b"] = -1.1562 
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 6) &  7 , "b"] = -1.14528
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 7) &  8 , "b"] = -1.13725
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 8) &  9 , "b"] = -1.11712
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 9)      , "b"] = -1.09426

		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep < 2)       , "a"] = 0.198306
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 2) &  3 , "a"] = 0.209544
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 3) &  4 , "a"] = 0.221003
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 4) &  5 , "a"] = 0.22733 
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 5) &  6 , "a"] = 0.230317
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 6) &  7 , "a"] = 0.23842 
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 7) &  8 , "a"] = 0.242428
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 8) &  9 , "a"] = 0.24472 
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 9)      , "a"] = 0.246185

		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep < 2)       , "b"] = -0.981174
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 2) &  3 , "b"] = -1.08578 
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 3) &  4 , "b"] = -1.17655 
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 4) &  5 , "b"] = -1.22691 
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 5) &  6 , "b"] = -1.26216 
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 6) &  7 , "b"] = -1.38061 
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 7) &  8 , "b"] = -1.47588 
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 8) &  9 , "b"] = -1.57429 
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 9)      , "b"] = -1.67625 
	else:
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep < 2)       , "a"] = 0.20233  
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 2) &  3 , "a"] = 0.212437 
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 3) &  4 , "a"] = 0.219554 
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 4) &  5 , "a"] = 0.224078 
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 5) &  6 , "a"] = 0.22785  
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 6) &  7 , "a"] = 0.230326 
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 7) &  8 , "a"] = 0.23258  
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 8) &  9 , "a"] = 0.232341 
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 9)      , "a"] = 0.22484  

		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep < 2)       , "b"] = -0.949695
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 2) &  3 , "b"] = -1.04219
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 3) &  4 , "b"] = -1.08307
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 4) &  5 , "b"] = -1.11223
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 5) &  6 , "b"] = -1.13177
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 6) &  7 , "b"] = -1.14596
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 7) &  8 , "b"] = -1.16058
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 8) &  9 , "b"] = -1.15733
		eleFidCut.loc[ (eleFidCut.Esector == 1) & (eleFidCut.Ep >= 9)      , "b"] = -1.08707

		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep < 2)       , "a"] =  0.202647
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 2) &  3 , "a"] =  0.21288 
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 3) &  4 , "a"] =  0.220162
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 4) &  5 , "a"] =  0.224914
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 5) &  6 , "a"] =  0.228319
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 6) &  7 , "a"] =  0.230837
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 7) &  8 , "a"] =  0.233087
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 8) &  9 , "a"] =  0.233287
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 9)      , "a"] =  0.233669

		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep < 2)       , "b"] = -0.956382
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 2) &  3 , "b"] = -1.04446
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 3) &  4 , "b"] = -1.09327
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 4) &  5 , "b"] = -1.11772
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 5) &  6 , "b"] = -1.13684
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 6) &  7 , "b"] = -1.14755
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 7) &  8 , "b"] = -1.16113
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 8) &  9 , "b"] = -1.16306
		eleFidCut.loc[ (eleFidCut.Esector == 2) & (eleFidCut.Ep >= 9)      , "b"] = -1.16674

		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep < 2)       , "a"] =  0.200986
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 2) &  3 , "a"] =  0.211826
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 3) &  4 , "a"] =  0.219478
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 4) &  5 , "a"] =  0.223814
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 5) &  6 , "a"] =  0.227713
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 6) &  7 , "a"] =  0.230451
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 7) &  8 , "a"] =  0.2319  
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 8) &  9 , "a"] =  0.231747
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 9)      , "a"] =  0.236727

		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep < 2)       , "b"] = -0.934262
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 2) &  3 , "b"] = -1.035
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 3) &  4 , "b"] = -1.08529
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 4) &  5 , "b"] = -1.10889
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 5) &  6 , "b"] = -1.13048
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 6) &  7 , "b"] = -1.14677
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 7) &  8 , "b"] = -1.1519
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 8) &  9 , "b"] = -1.14965
		eleFidCut.loc[ (eleFidCut.Esector == 3) & (eleFidCut.Ep >= 9)      , "b"] = -1.21746

		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep < 2)       , "a"] =  0.201774
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 2) &  3 , "a"] =  0.213158
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 3) &  4 , "a"] =  0.220034
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 4) &  5 , "a"] =  0.224357
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 5) &  6 , "a"] =  0.228105
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 6) &  7 , "a"] =  0.230969
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 7) &  8 , "a"] =  0.233145
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 8) &  9 , "a"] =  0.232848
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 9)      , "a"] =  0.233346

		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep < 2)       , "b"] = -0.940217
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 2) &  3 , "b"] = -1.04832
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 3) &  4 , "b"] = -1.08885
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 4) &  5 , "b"] = -1.11142
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 5) &  6 , "b"] = -1.13521
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 6) &  7 , "b"] = -1.15218
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 7) &  8 , "b"] = -1.1639
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 8) &  9 , "b"] = -1.15888
		eleFidCut.loc[ (eleFidCut.Esector == 4) & (eleFidCut.Ep >= 9)      , "b"] = -1.18771

		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep < 2)       , "a"] = 0.201084 
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 2) &  3 , "a"] = 0.212129 
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 3) &  4 , "a"] = 0.219012 
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 4) &  5 , "a"] = 0.224371 
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 5) &  6 , "a"] = 0.22709  
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 6) &  7 , "a"] = 0.229402 
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 7) &  8 , "a"] = 0.232148 
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 8) &  9 , "a"] = 0.231751 
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 9)      , "a"] = 0.233825 

		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep < 2)       , "b"] = -0.933617
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 2) &  3 , "b"] = -1.03597
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 3) &  4 , "b"] = -1.07962
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 4) &  5 , "b"] = -1.11528
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 5) &  6 , "b"] = -1.12421
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 6) &  7 , "b"] = -1.13488
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 7) &  8 , "b"] = -1.15457
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 8) &  9 , "b"] = -1.15131
		eleFidCut.loc[ (eleFidCut.Esector == 5) & (eleFidCut.Ep >= 9)      , "b"] = -1.17542

		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep < 2)       , "a"] = 0.201648 
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 2) &  3 , "a"] = 0.212346 
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 3) &  4 , "a"] = 0.219218 
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 4) &  5 , "a"] = 0.224192 
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 5) &  6 , "a"] = 0.227346 
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 6) &  7 , "a"] = 0.229665 
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 7) &  8 , "a"] = 0.232072 
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 8) &  9 , "a"] = 0.232098 
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 9)      , "a"] = 0.23524  

		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep < 2)       , "b"] = -0.936363
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 2) &  3 , "b"] = -1.03963
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 3) &  4 , "b"] = -1.08145
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 4) &  5 , "b"] = -1.11228
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 5) &  6 , "b"] = -1.12983
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 6) &  7 , "b"] = -1.14008
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 7) &  8 , "b"] = -1.1564
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 8) &  9 , "b"] = -1.15527
		eleFidCut.loc[ (eleFidCut.Esector == 6) & (eleFidCut.Ep >= 9)      , "b"] = -1.19943

	if fidlevel == 'mid':
		pass
	elif fidlevel == 'loose':
		eleFidCut.loc[:, "a"] = eleFidCut.loc[:, "a"] - 0.005
	elif fidlevel == 'tight':
		eleFidCut.loc[:, "a"] = eleFidCut.loc[:, "a"] + 0.005
	else:
		print("check fidlevel {}".format(fidlevel))

	df_electronRec.loc[eleFidCut.Eedep1/eleFidCut.Ep <= eleFidCut.a + eleFidCut.b * (eleFidCut.Eedep2/eleFidCut.Ep), "EFid"] = 0

	#Table XII.
	if fidlevel == 'mid':
		adjustment = 0
	elif fidlevel == 'loose':
		adjustment = -0.5
	elif fidlevel == 'tight':
		adjustment = +0.5
	else:
		print("check fidlevel {}".format(fidlevel))

	df_electronRec.loc[ (df_electronRec.Esector == 1) & (df_electronRec.EcalHy1 <= 0.56575  * df_electronRec.EcalHx1 -92        + 0.25 + adjustment) & (df_electronRec.EcalHy1 >= 0.56575 * df_electronRec.EcalHx1 -94.4         - 0.25 -adjustment), "EFid"] = 0
	df_electronRec.loc[ (df_electronRec.Esector == 1) & (df_electronRec.EcalHy1 <= 0.56575  * df_electronRec.EcalHx1 -101.1     + 0.25 + adjustment) & (df_electronRec.EcalHy1 >= 0.56575 * df_electronRec.EcalHx1 -103.5        - 0.25 -adjustment), "EFid"] = 0
	df_electronRec.loc[ (df_electronRec.Esector == 1) & (df_electronRec.EcalHy1 <= 0.56575  * df_electronRec.EcalHx1 -219       + 0.25 + adjustment) & (df_electronRec.EcalHy1 >= 0.56575 * df_electronRec.EcalHx1 -221.4        - 0.25 -adjustment), "EFid"] = 0
	df_electronRec.loc[ (df_electronRec.Esector == 1) & (df_electronRec.EcalHy1 <= 0.56575  * df_electronRec.EcalHx1 -227       + 0.25 + adjustment) & (df_electronRec.EcalHy1 >= 0.56575 * df_electronRec.EcalHx1 -229.4        - 0.25 -adjustment), "EFid"] = 0
	df_electronRec.loc[ (df_electronRec.Esector == 2) & (df_electronRec.EcalHy1 <= 0.5897   * df_electronRec.EcalHx1 +120.7937  + 0.25 + adjustment) & (df_electronRec.EcalHy1 >= 0.5913  * df_electronRec.EcalHx1 +114.3872     - 0.25 -adjustment), "EFid"] = 0
	df_electronRec.loc[ (df_electronRec.Esector == 2) & (df_electronRec.EcalHy1 <= 107.2766 * df_electronRec.EcalHx1 -10602.9779+ 0.25 + adjustment) & (df_electronRec.EcalHy1 >= 98.9667 * df_electronRec.EcalHx1 -10262.0167   - 0.25 -adjustment), "EFid"] = 0
	df_electronRec.loc[ (df_electronRec.Esector == 3) & (df_electronRec.EcalHx1 <= -302.38 + adjustment) & (df_electronRec.EcalHx1 >= -313.71 - adjustment), "EFid"] = 0
	df_electronRec.loc[ (df_electronRec.Esector == 4) & (df_electronRec.EcalHx1 <= -122.5  + adjustment) & (df_electronRec.EcalHx1 >= -127.5  - adjustment), "EFid"] = 0
	df_electronRec.loc[ (df_electronRec.Esector == 4) & (df_electronRec.EcalHy1 <= -0.568   * df_electronRec.EcalHx1 -232.8     + 0.25 + adjustment) & (df_electronRec.EcalHy1 >= -0.568  * df_electronRec.EcalHx1 -236.3        - 0.25 -adjustment), "EFid"] = 0
	df_electronRec.loc[ (df_electronRec.Esector == 5) & (df_electronRec.EcalHy1 <= 98.0644  * df_electronRec.EcalHx1 +5825.4023 + 0.25 + adjustment) & (df_electronRec.EcalHy1 >= 99.9337 * df_electronRec.EcalHx1 +5098.3456    - 0.25 -adjustment), "EFid"] = 0
	df_electronRec.loc[ (df_electronRec.Esector == 6) & (df_electronRec.EcalHy1 <= 0.4547   * df_electronRec.EcalHx1 -275.9317  + 0.25 + adjustment) & (df_electronRec.EcalHy1 >= 0.4547  * df_electronRec.EcalHx1 -285.9317     - 0.25 -adjustment), "EFid"] = 0
	df_electronRec.loc[ (df_electronRec.Esector == 6) & (df_electronRec.EcalHy1 <= 0.591377  * df_electronRec.EcalHx1 -185      + 0.25 + adjustment) & (df_electronRec.EcalHy1 >= 0.591377* df_electronRec.EcalHx1 -187          - 0.25 -adjustment), "EFid"] = 0
	df_electronRec.loc[ (df_electronRec.Esector == 6) & (df_electronRec.EcalHy1 <= 0.591377  * df_electronRec.EcalHx1 -193.3    + 0.25 + adjustment) & (df_electronRec.EcalHy1 >= 0.591377* df_electronRec.EcalHx1 -195.5        - 0.25 -adjustment), "EFid"] = 0

	df_electronRec.loc[ (df_electronRec.Esector == 5) & (df_electronRec.EcalHy3 <= -0.5841  * df_electronRec.EcalHx3 -252.11    + 0.25 + adjustment) & (df_electronRec.EcalHy3 >= -0.5775 * df_electronRec.EcalHx3 -263.2072    - 0.25 -adjustment), "EFid"] = 0

	return df_electronRec.loc[df_electronRec.EFid==1, :]

def electronFiducial_legacy(df_electronRec, pol = "inbending", mc = False, fidlevel = 'mid'):
	df_electronRec = copy(df_electronRec)
	df_electronRec.loc[:, "EFid"] = 1

	# #PCAL dead wires
	exclusion1_1 = (df_electronRec.EcalW1 > 74) & (df_electronRec.EcalW1 < 79.8)
	exclusion1_2 = (df_electronRec.EcalW1 > 83.6) & (df_electronRec.EcalW1 < 92.2)
	exclusion1_3 = (df_electronRec.EcalW1 > 212.5) & (df_electronRec.EcalW1 < 230)
	exclusion1 = exclusion1_1 | exclusion1_2 | exclusion1_3
	df_electronRec.loc[(df_electronRec.Esector == 1) & exclusion1, "EFid"] = 0
	exclusion2_1 = (df_electronRec.EcalW1 < 14)
	exclusion2_2 = (df_electronRec.EcalU1 > 111.2) & (df_electronRec.EcalU1 < 119.3)
	exclusion2_3 = (df_electronRec.EcalV1 > 113) & (df_electronRec.EcalV1 < 118.7)
	exclusion2 = exclusion2_1 | exclusion2_2 | exclusion2_3
	df_electronRec.loc[(df_electronRec.Esector == 2) & exclusion2, "EFid"] = 0
	exclusion3 = df_electronRec.EcalW1 < 14
	df_electronRec.loc[(df_electronRec.Esector == 3) & exclusion3, "EFid"] = 0
	exclusion4_1 = (df_electronRec.EcalV1 < 14)
	exclusion4_2 = (df_electronRec.EcalV1 > 229.4) & (df_electronRec.EcalV1 < 240.7)
	exclusion4_3 = (df_electronRec.EcalW1 > 135) & (df_electronRec.EcalW1 < 150)
	exclusion4 = exclusion4_1 | exclusion4_2 | exclusion4_3
	df_electronRec.loc[(df_electronRec.Esector == 4) & exclusion4, "EFid"] = 0
	exclusion6 = (df_electronRec.EcalW1 > 170) & (df_electronRec.EcalW1 < 192)
	df_electronRec.loc[(df_electronRec.Esector == 6) & exclusion6, "EFid"] = 0

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
	if fidlevel == 'mid':
		df_electronRec.loc[df_electronRec.ESamplFrac < mean - e_sampl_sigma_range*sigma, "EFid"]  = 0
		df_electronRec.loc[df_electronRec.ESamplFrac > mean + e_sampl_sigma_range*sigma, "EFid"]  = 0
	elif fidlevel == 'tight':
		df_electronRec.loc[df_electronRec.ESamplFrac < mean - (e_sampl_sigma_range-0.5)*sigma, "EFid"]  = 0
		df_electronRec.loc[df_electronRec.ESamplFrac > mean + (e_sampl_sigma_range-0.5)*sigma, "EFid"]  = 0


	#passElectronNpheCut
	df_electronRec.loc[df_electronRec.Enphe <= min_nphe, "EFid"] = 0

	#passElectronVertexCut
	if pol == 'inbending':
		min_vz = vz_min_inb
		max_vz = vz_max_inb
	if pol == 'outbending':
		min_vz = vz_min_outb
		max_vz = vz_max_outb
	df_electronRec.loc[df_electronRec.Evz <= min_vz, "EFid"] = 0
	df_electronRec.loc[df_electronRec.Evz >= max_vz, "EFid"] = 0

	# passElectronPCALFiducialCut
	if fidlevel == 'mid':
		df_electronRec.loc[df_electronRec.EcalV1 <= min_v, "EFid"] = 0
		df_electronRec.loc[df_electronRec.EcalW1 <= min_w, "EFid"] = 0
	elif fidlevel == 'tight':
		df_electronRec.loc[df_electronRec.EcalV1 <= min_v+5, "EFid"] = 0
		df_electronRec.loc[df_electronRec.EcalW1 <= min_w+5, "EFid"] = 0

	#passElectronPCALEdepCut
	df_electronRec.loc[df_electronRec.Eedep1 <= min_pcal_dep, "EFid"] = 0

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
	df_electronRec.loc[y_rot<=calc_min, "EFid"] = 0
	df_electronRec.loc[y_rot>=calc_max, "EFid"] = 0
	#passElectronDCR2
	dcsec = determineSector(df_electronRec.EDc2Hitx, df_electronRec.EDc2Hity)
	x_rot, y_rot = rotateDCHitPosition(df_electronRec.EDc2Hitx, df_electronRec.EDc2Hity, dcsec)
	calc_min, calc_max = e_DC_fiducial_cut_XY(x_rot, dcsec, 1, minparams, maxparams)
	df_electronRec.loc[y_rot<=calc_min, "EFid"] = 0
	df_electronRec.loc[y_rot>=calc_max, "EFid"] = 0

	#passElectronDCR3
	dcsec = determineSector(df_electronRec.EDc3Hitx, df_electronRec.EDc3Hity)
	x_rot, y_rot = rotateDCHitPosition(df_electronRec.EDc3Hitx, df_electronRec.EDc3Hity, dcsec)
	calc_min, calc_max = e_DC_fiducial_cut_XY(x_rot, dcsec, 2, minparams, maxparams)
	df_electronRec.loc[y_rot<=calc_min, "EFid"] = 0
	df_electronRec.loc[y_rot>=calc_max, "EFid"] = 0

	# #passElectronAntiPionCut
	df_electronRec.loc[(df_electronRec.Ep>4.5)&(-df_electronRec.Eedep1/df_electronRec.Ep + anti_pion_threshold > df_electronRec.Eedep2/df_electronRec.Ep), "EFid"] = 0
	return df_electronRec.loc[df_electronRec.EFid==1, :]

def gammaFiducial(df_gammaRec, fidlevel = 'mid'):
	df_gammaRec = copy(df_gammaRec)
	df_gammaRec.loc[:, "GFid"] = 1
	# H. PCAL Fid Cuts
	if fidlevel == 'mid':
		df_gammaRec.loc[(df_gammaRec.GcalV1<19) & (df_gammaRec.Gsector<7), "GFid"] = 0
		df_gammaRec.loc[(df_gammaRec.GcalW1<19) & (df_gammaRec.Gsector<7), "GFid"] = 0
		df_gammaRec.loc[(df_gammaRec.GcalU1>395) & (df_gammaRec.Gsector<7), "GFid"] = 0
	elif fidlevel == 'loose':
		df_gammaRec.loc[(df_gammaRec.GcalV1<19-2.5) & (df_gammaRec.Gsector<7), "GFid"] = 0
		df_gammaRec.loc[(df_gammaRec.GcalW1<19-2.5) & (df_gammaRec.Gsector<7), "GFid"] = 0
		df_gammaRec.loc[(df_gammaRec.GcalU1>395+2.5) & (df_gammaRec.Gsector<7), "GFid"] = 0
	elif fidlevel == 'tight':
		df_gammaRec.loc[(df_gammaRec.GcalV1<19+2.5) & (df_gammaRec.Gsector<7), "GFid"] = 0
		df_gammaRec.loc[(df_gammaRec.GcalW1<19+2.5) & (df_gammaRec.Gsector<7), "GFid"] = 0
		df_gammaRec.loc[(df_gammaRec.GcalU1>395-2.5) & (df_gammaRec.Gsector<7), "GFid"] = 0
	else:
		print("check fidlevel {}".format(fidlevel))
	#passGammaBetaCut
	df_gammaRec.loc[df_gammaRec.Gbeta <= min_Gbeta, "GFid"] = 0
	df_gammaRec.loc[df_gammaRec.Gbeta >= max_Gbeta, "GFid"] = 0

	df_gammaRec = copy(df_gammaRec.loc[df_gammaRec.GFid==1, :])

	df_gammaRec.loc[df_gammaRec.Gsector<7, "GFid"] = 0

	#apply additional photon fiducial cuts
	sector_cond = [df_gammaRec.Gsector ==1, df_gammaRec.Gsector ==2, df_gammaRec.Gsector ==3, df_gammaRec.Gsector ==4, df_gammaRec.Gsector ==5, df_gammaRec.Gsector ==6]
	psplit = np.select(sector_cond, [87, 82, 85, 77, 78, 82])
	tleft = np.select(sector_cond, [58.7356, 62.8204, 62.2296, 53.7756, 58.2888, 54.5822])
	tright = np.select(sector_cond, [58.7477, 51.2589, 59.2357, 56.2415, 60.8219, 49.8914])
	sleft = np.select(sector_cond, [0.582053, 0.544976, 0.549788, 0.56899, 0.56414, 0.57343])
	sright = np.select(sector_cond, [-0.591876, -0.562926, -0.562246, -0.563726, -0.568902, -0.550729])
	rleft = np.select(sector_cond, [64.9348, 64.7541, 67.832, 55.9324, 55.9225, 60.0997])
	rright = np.select(sector_cond, [65.424, 54.6992, 63.6628, 57.8931, 56.5367, 56.4641])
	qleft = np.select(sector_cond, [0.745578, 0.606081, 0.729202, 0.627239, 0.503674, 0.717899])
	qright = np.select(sector_cond, [-0.775022, -0.633863, -0.678901, -0.612458, -0.455319, -0.692481])
	#first condition
	ang = np.radians((df_gammaRec.loc[df_gammaRec.Gsector<7, "Gsector"]-1) * 60)
	GcX_rot = df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] * np.sin(ang) + df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] * np.cos(ang)
	GcY_rot = df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] * np.cos(ang) - df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] * np.sin(ang)

	df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] = GcX_rot
	df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] = GcY_rot

	cond1_1 = df_gammaRec.GcX >= psplit
	cond1_2 = df_gammaRec.GcY < sleft * (df_gammaRec.GcX - tleft)
	cond1_3 = df_gammaRec.GcY > sright * (df_gammaRec.GcX - tright)
	cond1_4 = df_gammaRec.Gsector < 7
	cond1 = cond1_1 & cond1_2 & cond1_3 & cond1_4
	df_gammaRec.loc[cond1, "GFid"] = 1
	#second condition else if the first
	# cond2_0 = df_gammaRec.GFid == 0 # not necessary, because cond2_1 rules out the first (S. Lee)
	cond2_1 = df_gammaRec.GcX < psplit
	cond2_2 = df_gammaRec.GcY < qleft * (df_gammaRec.GcX - rleft)
	cond2_3 = df_gammaRec.GcY > qright * (df_gammaRec.GcX - rright)
	cond2_4 = df_gammaRec.Gsector < 7
	cond2 = cond2_1 & cond2_2 & cond2_3 & cond2_4
	df_gammaRec.loc[cond2, "GFid"] = 1
    #photon FD fiducial cuts by F.X. Girod

	#FT fiducial cuts
	circleCenterX1 = -8.419
	circleCenterY1 = 9.889
	circleRadius1 = 1.6

	circleCenterX2 = -9.89
	circleCenterY2 = -5.327
	circleRadius2 = 1.6

	circleCenterX3 = -6.15
	circleCenterY3 = -13
	circleRadius3 = 2.3

	circleCenterX4 = 3.7
	circleCenterY4 = -6.5
	circleRadius4 = 2

	circle1 = (df_gammaRec.GcX - circleCenterX1)**2 + (df_gammaRec.GcY - circleCenterY1)**2 < circleRadius1**2
	circle2 = (df_gammaRec.GcX - circleCenterX2)**2 + (df_gammaRec.GcY - circleCenterY2)**2 < circleRadius2**2
	circle3 = (df_gammaRec.GcX - circleCenterX3)**2 + (df_gammaRec.GcY - circleCenterY3)**2 < circleRadius3**2
	circle4 = (df_gammaRec.GcX - circleCenterX4)**2 + (df_gammaRec.GcY - circleCenterY4)**2 < circleRadius4**2

	df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle1, "GFid"] = 0
	df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle2, "GFid"] = 0
	df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle3, "GFid"] = 0
	df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle4, "GFid"] = 0

	#Table XII.
	if fidlevel == 'mid':
		adjustment = 0
	elif fidlevel == 'loose':
		adjustment = -0.5
	elif fidlevel == 'tight':
		adjustment = +0.5
	else:
		print("check fidlevel {}".format(fidlevel))

	df_gammaRec.loc[ (df_gammaRec.Gsector == 1) & (df_gammaRec.GcalY1 <= 0.56575  * df_gammaRec.GcalX1 -92        + 0.25 + adjustment) & (df_gammaRec.GcalY1 >= 0.56575 * df_gammaRec.GcalX1 -94.4         - 0.25 -adjustment), "GFid"] = 0
	df_gammaRec.loc[ (df_gammaRec.Gsector == 1) & (df_gammaRec.GcalY1 <= 0.56575  * df_gammaRec.GcalX1 -101.1     + 0.25 + adjustment) & (df_gammaRec.GcalY1 >= 0.56575 * df_gammaRec.GcalX1 -103.5        - 0.25 -adjustment), "GFid"] = 0
	df_gammaRec.loc[ (df_gammaRec.Gsector == 1) & (df_gammaRec.GcalY1 <= 0.56575  * df_gammaRec.GcalX1 -219       + 0.25 + adjustment) & (df_gammaRec.GcalY1 >= 0.56575 * df_gammaRec.GcalX1 -221.4        - 0.25 -adjustment), "GFid"] = 0
	df_gammaRec.loc[ (df_gammaRec.Gsector == 1) & (df_gammaRec.GcalY1 <= 0.56575  * df_gammaRec.GcalX1 -227       + 0.25 + adjustment) & (df_gammaRec.GcalY1 >= 0.56575 * df_gammaRec.GcalX1 -229.4        - 0.25 -adjustment), "GFid"] = 0
	df_gammaRec.loc[ (df_gammaRec.Gsector == 2) & (df_gammaRec.GcalY1 <= 0.5897   * df_gammaRec.GcalX1 +120.7937  + 0.25 + adjustment) & (df_gammaRec.GcalY1 >= 0.5913  * df_gammaRec.GcalX1 +114.3872     - 0.25 -adjustment), "GFid"] = 0
	df_gammaRec.loc[ (df_gammaRec.Gsector == 2) & (df_gammaRec.GcalY1 <= 107.2766 * df_gammaRec.GcalX1 -10602.9779+ 0.25 + adjustment) & (df_gammaRec.GcalY1 >= 98.9667 * df_gammaRec.GcalX1 -10262.0167   - 0.25 -adjustment), "GFid"] = 0
	df_gammaRec.loc[ (df_gammaRec.Gsector == 3) & (df_gammaRec.GcalX1 <= -302.38 + adjustment) & (df_gammaRec.GcalX1 >= -313.71 - adjustment), "GFid"] = 0
	df_gammaRec.loc[ (df_gammaRec.Gsector == 4) & (df_gammaRec.GcalX1 <= -122.5  + adjustment) & (df_gammaRec.GcalX1 >= -127.5  - adjustment), "GFid"] = 0
	df_gammaRec.loc[ (df_gammaRec.Gsector == 4) & (df_gammaRec.GcalY1 <= -0.568   * df_gammaRec.GcalX1 -232.8     + 0.25 + adjustment) & (df_gammaRec.GcalY1 >= -0.568  * df_gammaRec.GcalX1 -236.3        - 0.25 -adjustment), "GFid"] = 0
	df_gammaRec.loc[ (df_gammaRec.Gsector == 5) & (df_gammaRec.GcalY1 <= 98.0644  * df_gammaRec.GcalX1 +5825.4023 + 0.25 + adjustment) & (df_gammaRec.GcalY1 >= 99.9337 * df_gammaRec.GcalX1 +5098.3456    - 0.25 -adjustment), "GFid"] = 0
	df_gammaRec.loc[ (df_gammaRec.Gsector == 6) & (df_gammaRec.GcalY1 <= 0.4547   * df_gammaRec.GcalX1 -275.9317  + 0.25 + adjustment) & (df_gammaRec.GcalY1 >= 0.4547  * df_gammaRec.GcalX1 -285.9317     - 0.25 -adjustment), "GFid"] = 0
	df_gammaRec.loc[ (df_gammaRec.Gsector == 6) & (df_gammaRec.GcalY1 <= 0.591377  * df_gammaRec.GcalX1 -185      + 0.25 + adjustment) & (df_gammaRec.GcalY1 >= 0.591377* df_gammaRec.GcalX1 -187          - 0.25 -adjustment), "GFid"] = 0
	df_gammaRec.loc[ (df_gammaRec.Gsector == 6) & (df_gammaRec.GcalY1 <= 0.591377  * df_gammaRec.GcalX1 -193.3    + 0.25 + adjustment) & (df_gammaRec.GcalY1 >= 0.591377* df_gammaRec.GcalX1 -195.5        - 0.25 -adjustment), "GFid"] = 0

	df_gammaRec.loc[ (df_gammaRec.Gsector == 5) & (df_gammaRec.GcalY3 <= -0.5841  * df_gammaRec.GcalX3 -252.11    + 0.25 + adjustment) & (df_gammaRec.GcalY3 >= -0.5775 * df_gammaRec.GcalX3 -263.2072    - 0.25 -adjustment), "GFid"] = 0

	return df_gammaRec.loc[df_gammaRec.GFid==1, :]



def gammaFiducialLegacy(df_gammaRec):
	df_gammaRec = copy(df_gammaRec)
	df_gammaRec.loc[:, "GFid"] = 1
	#passGammaPCALFiducialCut
	df_gammaRec.loc[(df_gammaRec.GcalV1 <= g_min_v) & (df_gammaRec.Gsector<7), "GFid"] = 0
	df_gammaRec.loc[(df_gammaRec.GcalW1 <= g_min_w) & (df_gammaRec.Gsector<7), "GFid"] = 0
	#passGammaBetaCut
	df_gammaRec.loc[df_gammaRec.Gbeta <= min_Gbeta, "GFid"] = 0
	df_gammaRec.loc[df_gammaRec.Gbeta >= max_Gbeta, "GFid"] = 0

	df_gammaRec.loc[df_gammaRec.Gsector<7, "GFid"] = 0

	#apply photon fiducial cuts
	sector_cond = [df_gammaRec.Gsector ==1, df_gammaRec.Gsector ==2, df_gammaRec.Gsector ==3, df_gammaRec.Gsector ==4, df_gammaRec.Gsector ==5, df_gammaRec.Gsector ==6]
	psplit = np.select(sector_cond, [87, 82, 85, 77, 78, 82])
	tleft = np.select(sector_cond, [58.7356, 62.8204, 62.2296, 53.7756, 58.2888, 54.5822])
	tright = np.select(sector_cond, [58.7477, 51.2589, 59.2357, 56.2415, 60.8219, 49.8914])
	sleft = np.select(sector_cond, [0.582053, 0.544976, 0.549788, 0.56899, 0.56414, 0.57343])
	sright = np.select(sector_cond, [-0.591876, -0.562926, -0.562246, -0.563726, -0.568902, -0.550729])
	rleft = np.select(sector_cond, [64.9348, 64.7541, 67.832, 55.9324, 55.9225, 60.0997])
	rright = np.select(sector_cond, [65.424, 54.6992, 63.6628, 57.8931, 56.5367, 56.4641])
	qleft = np.select(sector_cond, [0.745578, 0.606081, 0.729202, 0.627239, 0.503674, 0.717899])
	qright = np.select(sector_cond, [-0.775022, -0.633863, -0.678901, -0.612458, -0.455319, -0.692481])
	#first condition
	ang = np.radians((df_gammaRec.loc[df_gammaRec.Gsector<7, "Gsector"]-1) * 60)
	GcX_rot = df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] * np.sin(ang) + df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] * np.cos(ang)
	GcY_rot = df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] * np.cos(ang) - df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] * np.sin(ang)

	df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] = GcX_rot
	df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] = GcY_rot

	cond1_1 = df_gammaRec.GcX >= psplit
	cond1_2 = df_gammaRec.GcY < sleft * (df_gammaRec.GcX - tleft)
	cond1_3 = df_gammaRec.GcY > sright * (df_gammaRec.GcX - tright)
	cond1_4 = df_gammaRec.Gsector < 7
	cond1 = cond1_1 & cond1_2 & cond1_3 & cond1_4
	df_gammaRec.loc[cond1, "GFid"] = 1
	#second condition else if the first
	# cond2_0 = df_gammaRec.GFid == 0 # not necessary, because cond2_1 rules out the first (S. Lee)
	cond2_1 = df_gammaRec.GcX < psplit
	cond2_2 = df_gammaRec.GcY < qleft * (df_gammaRec.GcX - rleft)
	cond2_3 = df_gammaRec.GcY > qright * (df_gammaRec.GcX - rright)
	cond2_4 = df_gammaRec.Gsector < 7
	cond2 = cond2_1 & cond2_2 & cond2_3 & cond2_4
	df_gammaRec.loc[cond2, "GFid"] = 1
    #photon FD fiducial cuts by F.X. Girod

	#FT fiducial cuts
	circleCenterX1 = -8.419
	circleCenterY1 = 9.889
	circleRadius1 = 1.6

	circleCenterX2 = -9.89
	circleCenterY2 = -5.327
	circleRadius2 = 1.6

	circleCenterX3 = -6.15
	circleCenterY3 = -13
	circleRadius3 = 2.3

	circleCenterX4 = 3.7
	circleCenterY4 = -6.5
	circleRadius4 = 2

	circle1 = (df_gammaRec.GcX - circleCenterX1)**2 + (df_gammaRec.GcY - circleCenterY1)**2 < circleRadius1**2
	circle2 = (df_gammaRec.GcX - circleCenterX2)**2 + (df_gammaRec.GcY - circleCenterY2)**2 < circleRadius2**2
	circle3 = (df_gammaRec.GcX - circleCenterX3)**2 + (df_gammaRec.GcY - circleCenterY3)**2 < circleRadius3**2
	circle4 = (df_gammaRec.GcX - circleCenterX4)**2 + (df_gammaRec.GcY - circleCenterY4)**2 < circleRadius4**2

	df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle1, "GFid"] = 0
	df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle2, "GFid"] = 0
	df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle3, "GFid"] = 0
	df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle4, "GFid"] = 0

	exclusion1_1 = (df_gammaRec.GcalW1 > 74) & (df_gammaRec.GcalW1 < 79.8)
	exclusion1_2 = (df_gammaRec.GcalW1 > 83.6) & (df_gammaRec.GcalW1 < 92.2)
	exclusion1_3 = (df_gammaRec.GcalW1 > 212.5) & (df_gammaRec.GcalW1 < 230)
	exclusion1 = exclusion1_1 | exclusion1_2 | exclusion1_3
	df_gammaRec.loc[(df_gammaRec.Gsector == 1) & exclusion1, "GFid"] = 0
	exclusion2_1 = (df_gammaRec.GcalW1 < 14)
	exclusion2_2 = (df_gammaRec.GcalU1 > 111.2) & (df_gammaRec.GcalU1 < 119.3)
	exclusion2_3 = (df_gammaRec.GcalV1 > 113) & (df_gammaRec.GcalV1 < 118.7)
	exclusion2 = exclusion2_1 | exclusion2_2 | exclusion2_3
	df_gammaRec.loc[(df_gammaRec.Gsector == 2) & exclusion2, "GFid"] = 0
	exclusion3 = df_gammaRec.GcalW1 < 14
	df_gammaRec.loc[(df_gammaRec.Gsector == 3) & exclusion3, "GFid"] = 0
	exclusion4_1 = (df_gammaRec.GcalV1 < 14)
	exclusion4_2 = (df_gammaRec.GcalV1 > 229.4) & (df_gammaRec.GcalV1 < 240.7)
	exclusion4_3 = (df_gammaRec.GcalW1 > 135) & (df_gammaRec.GcalW1 < 150)
	exclusion4 = exclusion4_1 | exclusion4_2 | exclusion4_3
	df_gammaRec.loc[(df_gammaRec.Gsector == 4) & exclusion4, "GFid"] = 0
	exclusion6 = (df_gammaRec.GcalW1 > 170) & (df_gammaRec.GcalW1 < 192)
	df_gammaRec.loc[(df_gammaRec.Gsector == 6) & exclusion6, "GFid"] = 0

	return df_gammaRec.loc[df_gammaRec.GFid==1, :]

def protonFiducial(df_protonRec, fidlevel = 'mid'):
	df_protonRec = copy(df_protonRec)
	df_protonRec.loc[:, "PFid"] = 1

	#proton DC fiducial cut
	if fidlevel == 'mid':
		adjustment_layer1 = 0
		adjustment_layer2 = 0
		adjustment_layer3 = 0
	elif fidlevel == 'loose':
		adjustment_layer1 = 0.6*1
		adjustment_layer2 = 0.6*2
		adjustment_layer3 = 0.6*3
	elif fidlevel == 'tight':
		adjustment_layer1 = -0.6*1
		adjustment_layer2 = -0.6*2
		adjustment_layer3 = -0.6*3
	else:
		print("check fidlevel {}".format(fidlevel))

	dcsec_l1 = determineSector(df_protonRec.PDc1Hitx, df_protonRec.PDc1Hity)
	x_rot_l1, y_rot_l1 = rotateDCHitPosition(df_protonRec.PDc1Hitx, df_protonRec.PDc1Hity, dcsec_l1)
	x_rot_l1, _ = rotateDCHitPosition_alongY(x_rot_l1, df_protonRec.PDc1Hitz)
	calc_min_l1 = -0.50 * (x_rot_l1 + 72 + adjustment_layer1)
	calc_max_l1 =  0.50 * (x_rot_l1 + 72 + adjustment_layer1)
	df_protonRec.loc[y_rot_l1 < calc_min_l1, "PFid"] = 0
	df_protonRec.loc[y_rot_l1 > calc_max_l1, "PFid"] = 0

	dcsec_l2 = determineSector(df_protonRec.PDc2Hitx, df_protonRec.PDc2Hity)
	x_rot_l2, y_rot_l2 = rotateDCHitPosition(df_protonRec.PDc2Hitx, df_protonRec.PDc2Hity, dcsec_l2)
	x_rot_l2, _ = rotateDCHitPosition_alongY(x_rot_l2, df_protonRec.PDc2Hitz)
	calc_min_l2 = -0.505 * (x_rot_l2 + 114 + adjustment_layer2)
	calc_max_l2 =  0.505 * (x_rot_l2 + 114 + adjustment_layer2)
	df_protonRec.loc[y_rot_l2 < calc_min_l2, "PFid"] = 0
	df_protonRec.loc[y_rot_l2 > calc_max_l2, "PFid"] = 0

	dcsec_l3 = determineSector(df_protonRec.PDc3Hitx, df_protonRec.PDc3Hity)
	x_rot_l3, y_rot_l3 = rotateDCHitPosition(df_protonRec.PDc3Hitx, df_protonRec.PDc3Hity, dcsec_l3)
	x_rot_l3, _ = rotateDCHitPosition_alongY(x_rot_l3, df_protonRec.PDc3Hitz)
	calc_min_l3 = -0.495 * (x_rot_l3 + 180 + adjustment_layer3)
	calc_max_l3 =  0.495 * (x_rot_l3 + 180 + adjustment_layer3)
	df_protonRec.loc[y_rot_l3 < calc_min_l3, "PFid"] = 0
	df_protonRec.loc[y_rot_l3 > calc_max_l3, "PFid"] = 0

	#proton CVT fiducial cut
	df_protonRec.loc[:, "PCvt12theta"] = -100000
	df_protonRec.loc[:, "PCvt12phi"] = -100000

	cut_CD = df_protonRec.Psector > 7

	df_protonRec.loc[cut_CD, "PCvt12theta"] = getTheta([df_protonRec.loc[cut_CD].PCvt12Hitx, df_protonRec.loc[cut_CD].PCvt12Hity, df_protonRec.loc[cut_CD].PCvt12Hitz])
	df_protonRec.loc[cut_CD, "PCvt12phi"] = getPhi([df_protonRec.loc[cut_CD].PCvt12Hitx, df_protonRec.loc[cut_CD].PCvt12Hity, df_protonRec.loc[cut_CD].PCvt12Hitz])

	df_protonRec.loc[cut_CD, "PFid"] = 0 #CD fid reset
	if fidlevel == 'mid':
		cut_right = cut_CD  & (df_protonRec.Ptheta      < 64.23)
		cut_bottom = cut_CD & (df_protonRec.PCvt12theta > 46.5)
		cut_sidel = cut_CD  & (df_protonRec.PCvt12theta < -2.942 + 1.274*df_protonRec.Ptheta)
		cut_sider = cut_CD  & (df_protonRec.PCvt12theta > -3.523 + 1.046*df_protonRec.Ptheta)

		cut_gaps1 = ~((df_protonRec.PCvt12phi>-95) & (df_protonRec.PCvt12phi<-80))
		cut_gaps2 = ~((df_protonRec.PCvt12phi>25) & (df_protonRec.PCvt12phi<40))
		cut_gaps3 = ~((df_protonRec.PCvt12phi>143) & (df_protonRec.PCvt12phi<158))

	elif fidlevel == 'loose':
		cut_right = cut_CD  & (df_protonRec.Ptheta      < 64.23 + 2.5)
		cut_bottom = cut_CD & (df_protonRec.PCvt12theta > 46.5  - 2.5)
		cut_sidel = cut_CD  & (df_protonRec.PCvt12theta < -2.942 + 1.274*df_protonRec.Ptheta + 2.5)
		cut_sider = cut_CD  & (df_protonRec.PCvt12theta > -3.523 + 1.046*df_protonRec.Ptheta - 2.5)

		cut_gaps1 = ~((df_protonRec.PCvt12phi>-95 +2.5) & (df_protonRec.PCvt12phi<-80-2.5))
		cut_gaps2 = ~((df_protonRec.PCvt12phi>25 +2.5) & (df_protonRec.PCvt12phi<40-2.5))
		cut_gaps3 = ~((df_protonRec.PCvt12phi>143 +2.5) & (df_protonRec.PCvt12phi<158-2.5))

	elif fidlevel == 'tight':
		cut_right = cut_CD  & (df_protonRec.Ptheta      < 64.23 - 2.5)
		cut_bottom = cut_CD & (df_protonRec.PCvt12theta > 46.5  + 2.5)
		cut_sidel = cut_CD  & (df_protonRec.PCvt12theta < -2.942 + 1.274*df_protonRec.Ptheta - 2.5)
		cut_sider = cut_CD  & (df_protonRec.PCvt12theta > -3.523 + 1.046*df_protonRec.Ptheta + 2.5)

		cut_gaps1 = ~((df_protonRec.PCvt12phi>-95-2.5) & (df_protonRec.PCvt12phi<-80+2.5))
		cut_gaps2 = ~((df_protonRec.PCvt12phi>25-2.5) & (df_protonRec.PCvt12phi<40+2.5))
		cut_gaps3 = ~((df_protonRec.PCvt12phi>143-2.5) & (df_protonRec.PCvt12phi<158+2.5))
	else:
		print("check fidlevel {}".format(fidlevel))


	cut_trapezoid = cut_CD & cut_right & cut_bottom & cut_sidel & cut_sider
	cut_gaps = cut_CD & cut_gaps1 & cut_gaps2 & cut_gaps3
	cut_total = cut_gaps & cut_trapezoid

	df_protonRec.loc[cut_total, "PFid"] = 1 #CD fid

	return df_protonRec.loc[df_protonRec.PFid==1, :]

def protonFiducialCVT(df_protonRec):
	df_protonRec.loc[:, "PCvt12theta"] = -100000
	df_protonRec.loc[:, "PCvt12phi"] = -100000

	pro = [df_protonRec['Ppx'], df_protonRec['Ppy'], df_protonRec['Ppz']]
	df_protonRec.loc[:, 'Pp'] = mag(pro)
	df_protonRec.loc[:, 'Pe'] = getEnergy(pro, M)
	df_protonRec.loc[:, 'Ptheta'] = getTheta(pro)
	df_protonRec.loc[:, 'Pphi'] = getPhi(pro)

	cut_CD = df_protonRec.Psector > 7

	df_protonRec.loc[cut_CD, "PCvt12theta"] = getTheta([df_protonRec.loc[cut_CD].PCvt12Hitx, df_protonRec.loc[cut_CD].PCvt12Hity, df_protonRec.loc[cut_CD].PCvt12Hitz])
	df_protonRec.loc[cut_CD, "PCvt12phi"] = getPhi([df_protonRec.loc[cut_CD].PCvt12Hitx, df_protonRec.loc[cut_CD].PCvt12Hity, df_protonRec.loc[cut_CD].PCvt12Hitz])

	df_protonRec.loc[cut_CD, "PFid"] = 0 #CD fid reset
	if fidlevel == 'mid':
	    cut_right = cut_CD & (df_protonRec.Ptheta<max_Ptheta)
	elif fidlevel == 'tight':
		cut_right = cut_CD & (df_protonRec.Ptheta<max_Ptheta-5)
	cut_bottom = cut_CD & (df_protonRec.PCvt12theta>44.5)
	cut_sidel = cut_CD & (df_protonRec.PCvt12theta<-2.942 + 1.274*df_protonRec.Ptheta)
	cut_sider = cut_CD & (df_protonRec.PCvt12theta>-3.523 + 1.046*df_protonRec.Ptheta)

	cut_trapezoid = cut_CD & cut_right & cut_bottom & cut_sidel & cut_sider

	cut_gaps1 = ~((df_protonRec.PCvt12phi>-95) & (df_protonRec.PCvt12phi<-80))
	cut_gaps2 = ~((df_protonRec.PCvt12phi>25) & (df_protonRec.PCvt12phi<40))
	cut_gaps3 = ~((df_protonRec.PCvt12phi>143) & (df_protonRec.PCvt12phi<158))
	cut_gaps = cut_CD & cut_gaps1 & cut_gaps2 & cut_gaps3
	cut_total = cut_gaps & cut_trapezoid

	df_protonRec.loc[cut_total, "PFid"] = 1 #CD fid
	return df_protonRec.loc[df_protonRec.PFid==1, :]

def protonFiducialChi2Cut(pol, df_protonRec):
	# proton fiducial cuts            
	if pol == "inbending":
		pchi2CD_lb,   pchi2CD_ub   = -5.000, 6.345
		pchi2FD_S1_lb, pchi2FD_S1_ub = -3.296, 3.508
		pchi2FD_S2_lb, pchi2FD_S2_ub = -3.552, 4.000
		pchi2FD_S3_lb, pchi2FD_S3_ub = -3.446, 3.937
		pchi2FD_S4_lb, pchi2FD_S4_ub = -2.747, 3.190
		pchi2FD_S5_lb, pchi2FD_S5_ub = -2.851, 3.418
		pchi2FD_S6_lb, pchi2FD_S6_ub = -3.174, 3.514
	elif pol == "outbending":
		pchi2CD_lb,   pchi2CD_ub   = -5.592,  6.785
		pchi2FD_S1_lb, pchi2FD_S1_ub = -3.905, 4.088
		pchi2FD_S2_lb, pchi2FD_S2_ub = -3.411, 3.939
		pchi2FD_S3_lb, pchi2FD_S3_ub = -4.042, 5.954
		pchi2FD_S4_lb, pchi2FD_S4_ub = -3.820, 5.065
		pchi2FD_S5_lb, pchi2FD_S5_ub = -3.384, 4.232
		pchi2FD_S6_lb, pchi2FD_S6_ub = -5.077, 5.100

	df_protonRec.loc[ (df_protonRec.Psector>4000) & ((df_protonRec.Pchi2pid<pchi2CD_lb)   | (df_protonRec.Pchi2pid>pchi2CD_ub)  ), "PFid"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==1)   & ((df_protonRec.Pchi2pid<pchi2FD_S1_lb) | (df_protonRec.Pchi2pid>pchi2FD_S1_ub)), "PFid"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==2)   & ((df_protonRec.Pchi2pid<pchi2FD_S2_lb) | (df_protonRec.Pchi2pid>pchi2FD_S2_ub)), "PFid"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==3)   & ((df_protonRec.Pchi2pid<pchi2FD_S3_lb) | (df_protonRec.Pchi2pid>pchi2FD_S3_ub)), "PFid"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==4)   & ((df_protonRec.Pchi2pid<pchi2FD_S4_lb) | (df_protonRec.Pchi2pid>pchi2FD_S4_ub)), "PFid"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==5)   & ((df_protonRec.Pchi2pid<pchi2FD_S5_lb) | (df_protonRec.Pchi2pid>pchi2FD_S5_ub)), "PFid"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==6)   & ((df_protonRec.Pchi2pid<pchi2FD_S6_lb) | (df_protonRec.Pchi2pid>pchi2FD_S6_ub)), "PFid"] = 0

	return df_protonRec.loc[df_protonRec.PFid==1, :]


def protonFiducialVzCut(pol, df_protonRec):
	# proton fiducial cuts            
	if pol == "inbending":
		vzdiffCD_lb,    vzdiffCD_ub    = -2.011, 2.314
		vzdiffFD_S1_lb, vzdiffFD_S1_ub = -3.209, 4.017
		vzdiffFD_S2_lb, vzdiffFD_S2_ub = -3.612, 4.139
		vzdiffFD_S3_lb, vzdiffFD_S3_ub = -3.328, 4.287
		vzdiffFD_S4_lb, vzdiffFD_S4_ub = -3.411, 4.108
		vzdiffFD_S5_lb, vzdiffFD_S5_ub = -3.607, 4.246
		vzdiffFD_S6_lb, vzdiffFD_S6_ub = -2.999, 3.927
	elif pol == "outbending":
		vzdiffCD_lb,    vzdiffCD_ub    = -2.737, 2.096
		vzdiffFD_S1_lb, vzdiffFD_S1_ub = -4.435, 3.429
		vzdiffFD_S2_lb, vzdiffFD_S2_ub = -4.646, 2.978
		vzdiffFD_S3_lb, vzdiffFD_S3_ub = -3.922, 3.040
		vzdiffFD_S4_lb, vzdiffFD_S4_ub = -4.646, 3.493
		vzdiffFD_S5_lb, vzdiffFD_S5_ub = -3.901, 3.750
		vzdiffFD_S6_lb, vzdiffFD_S6_ub = -3.846, 3.623

	df_protonRec.loc[ (df_protonRec.Psector>4000) & ((df_protonRec.vzdiff<vzdiffCD_lb)   | (df_protonRec.vzdiff>vzdiffCD_ub)  ), "PFid"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==1)   & ((df_protonRec.vzdiff<vzdiffFD_S1_lb) | (df_protonRec.vzdiff>vzdiffFD_S1_ub)), "PFid"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==2)   & ((df_protonRec.vzdiff<vzdiffFD_S2_lb) | (df_protonRec.vzdiff>vzdiffFD_S2_ub)), "PFid"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==3)   & ((df_protonRec.vzdiff<vzdiffFD_S3_lb) | (df_protonRec.vzdiff>vzdiffFD_S3_ub)), "PFid"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==4)   & ((df_protonRec.vzdiff<vzdiffFD_S4_lb) | (df_protonRec.vzdiff>vzdiffFD_S4_ub)), "PFid"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==5)   & ((df_protonRec.vzdiff<vzdiffFD_S5_lb) | (df_protonRec.vzdiff>vzdiffFD_S5_ub)), "PFid"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==6)   & ((df_protonRec.vzdiff<vzdiffFD_S6_lb) | (df_protonRec.vzdiff>vzdiffFD_S6_ub)), "PFid"] = 0

	return df_protonRec.loc[df_protonRec.PFid==1, :]


def electronFiducialCounting(df_electronRec, pol = "inbending", mc = False, fidlevel = 'mid'):
	df_electronRec.loc[:, "EFid_dw"] = 1
	df_electronRec.loc[:, "EFid_sf"] = 1
	df_electronRec.loc[:, "EFid_vz"] = 1
	df_electronRec.loc[:, "EFid_pcal1"] = 1
	df_electronRec.loc[:, "EFid_vz"] = 1
	df_electronRec.loc[:, "EFid_edep"] = 1
	df_electronRec.loc[:, "EFid_dc"] = 1
	df_electronRec.loc[:, "EFid_ap"] = 1

	#PCAL dead wires
	exclusion1_1 = (df_electronRec.EcalW1 > 74) & (df_electronRec.EcalW1 < 79.8)
	exclusion1_2 = (df_electronRec.EcalW1 > 83.6) & (df_electronRec.EcalW1 < 92.2)
	exclusion1_3 = (df_electronRec.EcalW1 > 212.5) & (df_electronRec.EcalW1 < 230)
	exclusion1 = exclusion1_1 | exclusion1_2 | exclusion1_3
	df_electronRec.loc[(df_electronRec.Esector == 1) & exclusion1, "EFid_dw"] = 0
	exclusion2_1 = (df_electronRec.EcalW1 < 14)
	exclusion2_2 = (df_electronRec.EcalU1 > 111.2) & (df_electronRec.EcalU1 < 119.3)
	exclusion2_3 = (df_electronRec.EcalV1 > 113) & (df_electronRec.EcalV1 < 118.7)
	exclusion2 = exclusion2_1 | exclusion2_2 | exclusion2_3
	df_electronRec.loc[(df_electronRec.Esector == 2) & exclusion2, "EFid_dw"] = 0
	exclusion3 = df_electronRec.EcalW1 < 14
	df_electronRec.loc[(df_electronRec.Esector == 3) & exclusion3, "EFid_dw"] = 0
	exclusion4_1 = (df_electronRec.EcalV1 < 14)
	exclusion4_2 = (df_electronRec.EcalV1 > 229.4) & (df_electronRec.EcalV1 < 240.7)
	exclusion4_3 = (df_electronRec.EcalW1 > 135) & (df_electronRec.EcalW1 < 150)
	exclusion4 = exclusion4_1 | exclusion4_2 | exclusion4_3
	df_electronRec.loc[(df_electronRec.Esector == 4) & exclusion4, "EFid_dw"] = 0
	exclusion6 = (df_electronRec.EcalW1 > 170) & (df_electronRec.EcalW1 < 192)
	df_electronRec.loc[(df_electronRec.Esector == 6) & exclusion6, "EFid_dw"] = 0

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
	if fidlevel == 'mid':
		df_electronRec.loc[df_electronRec.ESamplFrac < mean - e_sampl_sigma_range*sigma, "EFid_sf"]  = 0
		df_electronRec.loc[df_electronRec.ESamplFrac > mean + e_sampl_sigma_range*sigma, "EFid_sf"]  = 0
	elif fidlevel == 'tight':
		df_electronRec.loc[df_electronRec.ESamplFrac < mean - (e_sampl_sigma_range-0.5)*sigma, "EFid_sf"]  = 0
		df_electronRec.loc[df_electronRec.ESamplFrac > mean + (e_sampl_sigma_range-0.5)*sigma, "EFid_sf"]  = 0


	# #passElectronNpheCut
	# df_electronRec.loc[df_electronRec.Enphe <= min_nphe, "EFid"] = 0

	#passElectronVertexCut
	if pol == 'inbending':
		min_vz = vz_min_inb
		max_vz = vz_max_inb
	if pol == 'outbending':
		min_vz = vz_min_outb
		max_vz = vz_max_outb
	df_electronRec.loc[df_electronRec.Evz <= min_vz, "EFid_vz"] = 0
	df_electronRec.loc[df_electronRec.Evz >= max_vz, "EFid_vz"] = 0

	# passElectronPCALFiducialCut
	if fidlevel == 'mid':
		df_electronRec.loc[df_electronRec.EcalV1 <= min_v, "EFid_pcal1"] = 0
		df_electronRec.loc[df_electronRec.EcalW1 <= min_w, "EFid_pcal1"] = 0
	elif fidlevel == 'tight':
		df_electronRec.loc[df_electronRec.EcalV1 <= min_v+10, "EFid_pcal1"] = 0
		df_electronRec.loc[df_electronRec.EcalW1 <= min_w+10, "EFid_pcal1"] = 0

	#passElectronPCALEdepCut
	df_electronRec.loc[df_electronRec.Eedep1 <= min_pcal_dep, "EFid_edep"] = 0

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
	df_electronRec.loc[y_rot<=calc_min, "EFid_dc"] = 0
	df_electronRec.loc[y_rot>=calc_max, "EFid_dc"] = 0
	#passElectronDCR2
	dcsec = determineSector(df_electronRec.EDc2Hitx, df_electronRec.EDc2Hity)
	x_rot, y_rot = rotateDCHitPosition(df_electronRec.EDc2Hitx, df_electronRec.EDc2Hity, dcsec)
	calc_min, calc_max = e_DC_fiducial_cut_XY(x_rot, dcsec, 1, minparams, maxparams)
	df_electronRec.loc[y_rot<=calc_min, "EFid_dc"] = 0
	df_electronRec.loc[y_rot>=calc_max, "EFid_dc"] = 0

	#passElectronDCR3
	dcsec = determineSector(df_electronRec.EDc3Hitx, df_electronRec.EDc3Hity)
	x_rot, y_rot = rotateDCHitPosition(df_electronRec.EDc3Hitx, df_electronRec.EDc3Hity, dcsec)
	calc_min, calc_max = e_DC_fiducial_cut_XY(x_rot, dcsec, 2, minparams, maxparams)
	df_electronRec.loc[y_rot<=calc_min, "EFid_dc"] = 0
	df_electronRec.loc[y_rot>=calc_max, "EFid_dc"] = 0

	# #passElectronAntiPionCut
	df_electronRec.loc[(df_electronRec.Ep>4.5)&(-df_electronRec.Eedep1/df_electronRec.Ep + anti_pion_threshold > df_electronRec.Eedep2/df_electronRec.Ep), "EFid_ap"] = 0

	df_electronRec.loc[:, "EFid_pcal"] = df_electronRec.loc[:, "EFid_pcal1"] * df_electronRec.loc[:, "EFid_dw"]
	df_electronRec.loc[:, "EFid"] = df_electronRec.EFid_dw*df_electronRec.EFid_sf*df_electronRec.EFid_vz*df_electronRec.EFid_pcal1*df_electronRec.EFid_vz*df_electronRec.EFid_edep*df_electronRec.EFid_dc*df_electronRec.EFid_ap
	return df_electronRec

def gammaFiducialCounting(df_gammaRec):
	df_gammaRec.loc[:, "GFid_beta"] = 1
	df_gammaRec.loc[:, "GFid_Pcal1"] = 1
	df_gammaRec.loc[:, "GFid_Pcal2"] = 0
	df_gammaRec.loc[df_gammaRec.Gsector>7, "GFid_Pcal2"] = 1
	df_gammaRec.loc[:, "GFid_FT"] = 1
	#passGammaPCALFiducialCut
	df_gammaRec.loc[(df_gammaRec.GcalV1 <= g_min_v) & (df_gammaRec.Gsector<7), "GFid_Pcal1"] = 0
	df_gammaRec.loc[(df_gammaRec.GcalW1 <= g_min_w) & (df_gammaRec.Gsector<7), "GFid_Pcal1"] = 0
	#passGammaBetaCut
	df_gammaRec.loc[df_gammaRec.Gbeta <= min_Gbeta, "GFid_beta"] = 0
	df_gammaRec.loc[df_gammaRec.Gbeta >= max_Gbeta, "GFid_beta"] = 0

	#photon FD fiducial cuts by F.X. Girod
	#apply photon fiducial cuts
	sector_cond = [df_gammaRec.Gsector ==1, df_gammaRec.Gsector ==2, df_gammaRec.Gsector ==3, df_gammaRec.Gsector ==4, df_gammaRec.Gsector ==5, df_gammaRec.Gsector ==6]
	psplit = np.select(sector_cond, [87, 82, 85, 77, 78, 82])
	tleft = np.select(sector_cond, [58.7356, 62.8204, 62.2296, 53.7756, 58.2888, 54.5822])
	tright = np.select(sector_cond, [58.7477, 51.2589, 59.2357, 56.2415, 60.8219, 49.8914])
	sleft = np.select(sector_cond, [0.582053, 0.544976, 0.549788, 0.56899, 0.56414, 0.57343])
	sright = np.select(sector_cond, [-0.591876, -0.562926, -0.562246, -0.563726, -0.568902, -0.550729])
	rleft = np.select(sector_cond, [64.9348, 64.7541, 67.832, 55.9324, 55.9225, 60.0997])
	rright = np.select(sector_cond, [65.424, 54.6992, 63.6628, 57.8931, 56.5367, 56.4641])
	qleft = np.select(sector_cond, [0.745578, 0.606081, 0.729202, 0.627239, 0.503674, 0.717899])
	qright = np.select(sector_cond, [-0.775022, -0.633863, -0.678901, -0.612458, -0.455319, -0.692481])
	#first condition
	ang = np.radians((df_gammaRec.loc[df_gammaRec.Gsector<7, "Gsector"]-1) * 60)
	GcX_rot = df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] * np.sin(ang) + df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] * np.cos(ang)
	GcY_rot = df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] * np.cos(ang) - df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] * np.sin(ang)

	df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] = GcX_rot
	df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] = GcY_rot

	cond1_1 = df_gammaRec.GcX >= psplit
	cond1_2 = df_gammaRec.GcY < sleft * (df_gammaRec.GcX - tleft)
	cond1_3 = df_gammaRec.GcY > sright * (df_gammaRec.GcX - tright)
	cond1_4 = df_gammaRec.Gsector < 7
	cond1 = cond1_1 & cond1_2 & cond1_3 & cond1_4
	df_gammaRec.loc[cond1, "GFid_Pcal2"] = 1
	#second condition else if the first
	# cond2_0 = df_gammaRec.GFid == 0 # not necessary, because cond2_1 rules out the first (S. Lee)
	cond2_1 = df_gammaRec.GcX < psplit
	cond2_2 = df_gammaRec.GcY < qleft * (df_gammaRec.GcX - rleft)
	cond2_3 = df_gammaRec.GcY > qright * (df_gammaRec.GcX - rright)
	cond2_4 = df_gammaRec.Gsector < 7
	cond2 = cond2_1 & cond2_2 & cond2_3 & cond2_4
	df_gammaRec.loc[cond2, "GFid_Pcal2"] = 1

	#FT fiducial cuts
	circleCenterX1 = -8.419
	circleCenterY1 = 9.889
	circleRadius1 = 1.6

	circleCenterX2 = -9.89
	circleCenterY2 = -5.327
	circleRadius2 = 1.6

	circleCenterX3 = -6.15
	circleCenterY3 = -13
	circleRadius3 = 2.3

	circleCenterX4 = 3.7
	circleCenterY4 = -6.5
	circleRadius4 = 2

	circle1 = (df_gammaRec.GcX - circleCenterX1)**2 + (df_gammaRec.GcY - circleCenterY1)**2 < circleRadius1**2
	circle2 = (df_gammaRec.GcX - circleCenterX2)**2 + (df_gammaRec.GcY - circleCenterY2)**2 < circleRadius2**2
	circle3 = (df_gammaRec.GcX - circleCenterX3)**2 + (df_gammaRec.GcY - circleCenterY3)**2 < circleRadius3**2
	circle4 = (df_gammaRec.GcX - circleCenterX4)**2 + (df_gammaRec.GcY - circleCenterY4)**2 < circleRadius4**2

	df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle1, "GFid_FT"] = 0
	df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle2, "GFid_FT"] = 0
	df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle3, "GFid_FT"] = 0
	df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle4, "GFid_FT"] = 0

	exclusion1_1 = (df_gammaRec.GcalW1 > 74) & (df_gammaRec.GcalW1 < 79.8)
	exclusion1_2 = (df_gammaRec.GcalW1 > 83.6) & (df_gammaRec.GcalW1 < 92.2)
	exclusion1_3 = (df_gammaRec.GcalW1 > 212.5) & (df_gammaRec.GcalW1 < 230)
	exclusion1 = exclusion1_1 | exclusion1_2 | exclusion1_3
	df_gammaRec.loc[(df_gammaRec.Gsector == 1) & exclusion1, "GFid_FT"] = 0
	exclusion2_1 = (df_gammaRec.GcalW1 < 14)
	exclusion2_2 = (df_gammaRec.GcalU1 > 111.2) & (df_gammaRec.GcalU1 < 119.3)
	exclusion2_3 = (df_gammaRec.GcalV1 > 113) & (df_gammaRec.GcalV1 < 118.7)
	exclusion2 = exclusion2_1 | exclusion2_2 | exclusion2_3
	df_gammaRec.loc[(df_gammaRec.Gsector == 2) & exclusion2, "GFid_FT"] = 0
	exclusion3 = df_gammaRec.GcalW1 < 14
	df_gammaRec.loc[(df_gammaRec.Gsector == 3) & exclusion3, "GFid_FT"] = 0
	exclusion4_1 = (df_gammaRec.GcalV1 < 14)
	exclusion4_2 = (df_gammaRec.GcalV1 > 229.4) & (df_gammaRec.GcalV1 < 240.7)
	exclusion4_3 = (df_gammaRec.GcalW1 > 135) & (df_gammaRec.GcalW1 < 150)
	exclusion4 = exclusion4_1 | exclusion4_2 | exclusion4_3
	df_gammaRec.loc[(df_gammaRec.Gsector == 4) & exclusion4, "GFid_FT"] = 0
	exclusion6 = (df_gammaRec.GcalW1 > 170) & (df_gammaRec.GcalW1 < 192)
	df_gammaRec.loc[(df_gammaRec.Gsector == 6) & exclusion6, "GFid_FT"] = 0

	df_gammaRec.loc[:, "GFid_Pcal"] = 	df_gammaRec.loc[:, "GFid_Pcal1"] * 	df_gammaRec.loc[:, "GFid_Pcal2"]
	df_gammaRec.loc[:, "GFid"] = df_gammaRec.GFid_beta * df_gammaRec.GFid_Pcal1 * df_gammaRec.GFid_Pcal2 * df_gammaRec.GFid_FT 
	return df_gammaRec

def protonFiducialCounting(df_protonRec, pol = 'inbending', fidlevel = 'mid'):
	df_protonRec.loc[:, "PFid_dc"] = 1
	df_protonRec.loc[:, "PFid_cvt"] = 0
	df_protonRec.loc[df_protonRec.Psector<7, "PFid_cvt"] = 1 #FD fid done by previous pipeline
	df_protonRec.loc[:, "PFid_chi"] = 1
	df_protonRec.loc[:, "PFid_vz"] = 1

	dcsec = determineSector(df_protonRec.PDc1Hitx, df_protonRec.PDc1Hity)
	if pol == 'inbending':
		minparams = p_dc_minparams_in
		maxparams = p_dc_maxparams_in

		theta_DC, phi_DC = thetaphifromhit(df_protonRec.PDc1Hitx, df_protonRec.PDc1Hity, df_protonRec.PDc1Hitz)
		phi_DC_min, phi_DC_max = p_DC_fiducial_cut_thetaphi(theta_DC, dcsec, 0, minparams, maxparams)
		df_protonRec.loc[(phi_DC<=phi_DC_min) & (df_protonRec.Psector<7), "PFid_dc"] = 0
		df_protonRec.loc[(phi_DC>=phi_DC_max) & (df_protonRec.Psector<7), "PFid_dc"] = 0

		theta_DC, phi_DC = thetaphifromhit(df_protonRec.PDc2Hitx, df_protonRec.PDc2Hity, df_protonRec.PDc2Hitz)
		phi_DC_min, phi_DC_max = p_DC_fiducial_cut_thetaphi(theta_DC, dcsec, 1, minparams, maxparams)
		df_protonRec.loc[(phi_DC<=phi_DC_min) & (df_protonRec.Psector<7), "PFid_dc"] = 0
		df_protonRec.loc[(phi_DC>=phi_DC_max) & (df_protonRec.Psector<7), "PFid_dc"] = 0

		theta_DC, phi_DC = thetaphifromhit(df_protonRec.PDc3Hitx, df_protonRec.PDc3Hity, df_protonRec.PDc3Hitz)
		phi_DC_min, phi_DC_max = p_DC_fiducial_cut_thetaphi(theta_DC, dcsec, 2, minparams, maxparams)
		df_protonRec.loc[(phi_DC<=phi_DC_min) & (df_protonRec.Psector<7), "PFid_dc"] = 0
		df_protonRec.loc[(phi_DC>=phi_DC_max) & (df_protonRec.Psector<7), "PFid_dc"] = 0

		# pchi2CD_lb,   pchi2CD_ub   = -5.000, 6.345
		# pchi2FD_S1_lb, pchi2FD_S1_ub = -3.296, 3.508
		# pchi2FD_S2_lb, pchi2FD_S2_ub = -3.552, 4.000
		# pchi2FD_S3_lb, pchi2FD_S3_ub = -3.446, 3.937
		# pchi2FD_S4_lb, pchi2FD_S4_ub = -2.747, 3.190
		# pchi2FD_S5_lb, pchi2FD_S5_ub = -2.851, 3.418
		# pchi2FD_S6_lb, pchi2FD_S6_ub = -3.174, 3.514

		# vzdiffCD_lb,    vzdiffCD_ub    = -2.011, 2.314
		# vzdiffFD_S1_lb, vzdiffFD_S1_ub = -3.209, 4.017
		# vzdiffFD_S2_lb, vzdiffFD_S2_ub = -3.612, 4.139
		# vzdiffFD_S3_lb, vzdiffFD_S3_ub = -3.328, 4.287
		# vzdiffFD_S4_lb, vzdiffFD_S4_ub = -3.411, 4.108
		# vzdiffFD_S5_lb, vzdiffFD_S5_ub = -3.607, 4.246
		# vzdiffFD_S6_lb, vzdiffFD_S6_ub = -2.999, 3.927

	if pol == 'outbending':
		minparams = p_dc_minparams_out
		maxparams = p_dc_maxparams_out

		dcsec = determineSector(df_protonRec.PDc1Hitx, df_protonRec.PDc1Hity)
		x_rot, y_rot = rotateDCHitPosition(df_protonRec.PDc1Hitx, df_protonRec.PDc1Hity, dcsec)
		calc_min, calc_max = p_DC_fiducial_cut_XY(x_rot, dcsec, 0, minparams, maxparams)
		df_protonRec.loc[(y_rot<=calc_min) & (df_protonRec.Psector<7), "PFid_dc"] = 0
		df_protonRec.loc[(y_rot>=calc_max) & (df_protonRec.Psector<7), "PFid_dc"] = 0
		#passElectronDCR2
		dcsec = determineSector(df_protonRec.PDc2Hitx, df_protonRec.PDc2Hity)
		x_rot, y_rot = rotateDCHitPosition(df_protonRec.PDc2Hitx, df_protonRec.PDc2Hity, dcsec)
		calc_min, calc_max = p_DC_fiducial_cut_XY(x_rot, dcsec, 1, minparams, maxparams)
		df_protonRec.loc[(y_rot<=calc_min) & (df_protonRec.Psector<7), "PFid_dc"] = 0
		df_protonRec.loc[(y_rot>=calc_max) & (df_protonRec.Psector<7), "PFid_dc"] = 0

		#passElectronDCR3
		dcsec = determineSector(df_protonRec.PDc3Hitx, df_protonRec.PDc3Hity)
		x_rot, y_rot = rotateDCHitPosition(df_protonRec.PDc3Hitx, df_protonRec.PDc3Hity, dcsec)
		calc_min, calc_max = p_DC_fiducial_cut_XY(x_rot, dcsec, 2, minparams, maxparams)
		df_protonRec.loc[(y_rot<=calc_min) & (df_protonRec.Psector<7), "PFid_dc"] = 0
		df_protonRec.loc[(y_rot>=calc_max) & (df_protonRec.Psector<7), "PFid_dc"] = 0

		# pchi2CD_lb,   pchi2CD_ub   = -5.592,  6.785
		# pchi2FD_S1_lb, pchi2FD_S1_ub = -3.905, 4.088
		# pchi2FD_S2_lb, pchi2FD_S2_ub = -3.411, 3.939
		# pchi2FD_S3_lb, pchi2FD_S3_ub = -4.042, 5.954
		# pchi2FD_S4_lb, pchi2FD_S4_ub = -3.820, 5.065
		# pchi2FD_S5_lb, pchi2FD_S5_ub = -3.384, 4.232
		# pchi2FD_S6_lb, pchi2FD_S6_ub = -5.077, 5.100

		# vzdiffCD_lb,    vzdiffCD_ub    = -2.737, 2.096
		# vzdiffFD_S1_lb, vzdiffFD_S1_ub = -4.435, 3.429
		# vzdiffFD_S2_lb, vzdiffFD_S2_ub = -4.646, 2.978
		# vzdiffFD_S3_lb, vzdiffFD_S3_ub = -3.922, 3.040
		# vzdiffFD_S4_lb, vzdiffFD_S4_ub = -4.646, 3.493
		# vzdiffFD_S5_lb, vzdiffFD_S5_ub = -3.901, 3.750
		# vzdiffFD_S6_lb, vzdiffFD_S6_ub = -3.846, 3.623

	cut_CD = df_protonRec.Psector > 7
	if fidlevel == 'mid':
		cut_right = cut_CD & (df_protonRec.Ptheta<max_Ptheta)
	elif fidlevel == 'tight':
		cut_right = cut_CD & (df_protonRec.Ptheta<max_Ptheta-5)
	cut_bottom = cut_CD & (df_protonRec.PCvt12theta>44.5)
	cut_sidel = cut_CD & (df_protonRec.PCvt12theta<-2.942 + 1.274*df_protonRec.Ptheta)
	cut_sider = cut_CD & (df_protonRec.PCvt12theta>-3.523 + 1.046*df_protonRec.Ptheta)

	cut_trapezoid = cut_CD & cut_right & cut_bottom & cut_sidel & cut_sider

	cut_gaps1 = ~((df_protonRec.PCvt12phi>-95) & (df_protonRec.PCvt12phi<-80))
	cut_gaps2 = ~((df_protonRec.PCvt12phi>25) & (df_protonRec.PCvt12phi<40))
	cut_gaps3 = ~((df_protonRec.PCvt12phi>143) & (df_protonRec.PCvt12phi<158))
	cut_gaps = cut_CD & cut_gaps1 & cut_gaps2 & cut_gaps3
	cut_total = cut_gaps & cut_trapezoid

	df_protonRec.loc[cut_total, "PFid_cvt"] = 1 #CD fid

	# df_protonRec.loc[ (df_protonRec.Psector>4000) & ((df_protonRec.Pchi2pid<pchi2CD_lb)   | (df_protonRec.Pchi2pid>pchi2CD_ub)  ), "PFid_chi"] = 0
	# df_protonRec.loc[ (df_protonRec.Psector==1)   & ((df_protonRec.Pchi2pid<pchi2FD_S1_lb) | (df_protonRec.Pchi2pid>pchi2FD_S1_ub)), "PFid_chi"] = 0
	# df_protonRec.loc[ (df_protonRec.Psector==2)   & ((df_protonRec.Pchi2pid<pchi2FD_S2_lb) | (df_protonRec.Pchi2pid>pchi2FD_S2_ub)), "PFid_chi"] = 0
	# df_protonRec.loc[ (df_protonRec.Psector==3)   & ((df_protonRec.Pchi2pid<pchi2FD_S3_lb) | (df_protonRec.Pchi2pid>pchi2FD_S3_ub)), "PFid_chi"] = 0
	# df_protonRec.loc[ (df_protonRec.Psector==4)   & ((df_protonRec.Pchi2pid<pchi2FD_S4_lb) | (df_protonRec.Pchi2pid>pchi2FD_S4_ub)), "PFid_chi"] = 0
	# df_protonRec.loc[ (df_protonRec.Psector==5)   & ((df_protonRec.Pchi2pid<pchi2FD_S5_lb) | (df_protonRec.Pchi2pid>pchi2FD_S5_ub)), "PFid_chi"] = 0
	# df_protonRec.loc[ (df_protonRec.Psector==6)   & ((df_protonRec.Pchi2pid<pchi2FD_S6_lb) | (df_protonRec.Pchi2pid>pchi2FD_S6_ub)), "PFid_chi"] = 0


	# df_protonRec.loc[ (df_protonRec.Psector>4000) & ((df_protonRec.vzdiff<vzdiffCD_lb)   | (df_protonRec.vzdiff>vzdiffCD_ub)  ), "PFid_vz"] = 0
	# df_protonRec.loc[ (df_protonRec.Psector==1)   & ((df_protonRec.vzdiff<vzdiffFD_S1_lb) | (df_protonRec.vzdiff>vzdiffFD_S1_ub)), "PFid_vz"] = 0
	# df_protonRec.loc[ (df_protonRec.Psector==2)   & ((df_protonRec.vzdiff<vzdiffFD_S2_lb) | (df_protonRec.vzdiff>vzdiffFD_S2_ub)), "PFid_vz"] = 0
	# df_protonRec.loc[ (df_protonRec.Psector==3)   & ((df_protonRec.vzdiff<vzdiffFD_S3_lb) | (df_protonRec.vzdiff>vzdiffFD_S3_ub)), "PFid_vz"] = 0
	# df_protonRec.loc[ (df_protonRec.Psector==4)   & ((df_protonRec.vzdiff<vzdiffFD_S4_lb) | (df_protonRec.vzdiff>vzdiffFD_S4_ub)), "PFid_vz"] = 0
	# df_protonRec.loc[ (df_protonRec.Psector==5)   & ((df_protonRec.vzdiff<vzdiffFD_S5_lb) | (df_protonRec.vzdiff>vzdiffFD_S5_ub)), "PFid_vz"] = 0
	# df_protonRec.loc[ (df_protonRec.Psector==6)   & ((df_protonRec.vzdiff<vzdiffFD_S6_lb) | (df_protonRec.vzdiff>vzdiffFD_S6_ub)), "PFid_vz"] = 0

	df_protonRec.loc[:, "PFid"] = df_protonRec.PFid_dc * df_protonRec.PFid_cvt #* df_protonRec.PFid_chi * df_protonRec.PFid_vz 
	return df_protonRec
