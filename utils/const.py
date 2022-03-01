#!/bin/env python3
"""
Modules to define useful constants
"""

import numpy as np

M = 0.938272081 # target mass
me = 0.5109989461 * 0.001 # electron mass
ebeam = 10.604 # beam energy
pbeam = np.sqrt(ebeam * ebeam - me * me) # beam electron momentum
beam = [0, 0, pbeam] # beam vector
target = [0, 0, 0] # target vector


#reaction: dvcs, pi0
#polarity: inb, outb
#topology: CDFT, CD (CDFD), FD(FDFD)
#sigmas: 2, 3, 4
# total: 2*2*3*3=36

cuts_dvcs_CDFT_Inb_3sigma = {}
cuts_dvcs_CDFT_Inb_3sigma["MM2_ep_ub"] = 0.391
cuts_dvcs_CDFT_Inb_3sigma["MM2_ep_lb"] = -0.365
cuts_dvcs_CDFT_Inb_3sigma["MM2_eg_ub"] = 1.479
cuts_dvcs_CDFT_Inb_3sigma["MM2_eg_lb"] = 0.322
cuts_dvcs_CDFT_Inb_3sigma["ME_epg_ub"] = 0.360
cuts_dvcs_CDFT_Inb_3sigma["ME_epg_lb"] = -0.387
cuts_dvcs_CDFT_Inb_3sigma["MM2_epg_ub"] = 0.000629
cuts_dvcs_CDFT_Inb_3sigma["MM2_epg_lb"] = -0.00423
cuts_dvcs_CDFT_Inb_3sigma["coplanarity_ub"] = 5.993
cuts_dvcs_CDFT_Inb_3sigma["coneAngle_ub"] = [ 0.133, -0.27, 15.103]
cuts_dvcs_CDFT_Inb_3sigma["coneAngle_lb"] = [ 0.179, -0.410, 7.000]
cuts_dvcs_CDFT_Inb_3sigma["MPt_ub"] = 0.0742
cuts_dvcs_CDFT_Inb_3sigma["reconGam_ub"] = 0.663

cuts_dvcs_CDFT_Inb_2sigma = {}
cuts_dvcs_CDFT_Inb_2sigma["MM2_ep_ub"] = 0.265
cuts_dvcs_CDFT_Inb_2sigma["MM2_ep_lb"] = -0.239
cuts_dvcs_CDFT_Inb_2sigma["MM2_eg_ub"] = 1.272
cuts_dvcs_CDFT_Inb_2sigma["MM2_eg_lb"] = 0.527
cuts_dvcs_CDFT_Inb_2sigma["ME_epg_ub"] = 0.255
cuts_dvcs_CDFT_Inb_2sigma["ME_epg_lb"] = -0.2232
cuts_dvcs_CDFT_Inb_2sigma["MM2_epg_ub"] = 0.000405
cuts_dvcs_CDFT_Inb_2sigma["MM2_epg_lb"] = -0.00360
cuts_dvcs_CDFT_Inb_2sigma["coplanarity_ub"] = 4.075
cuts_dvcs_CDFT_Inb_2sigma["coneAngle_ub"] = [0.0943, 0.295, 11.847]
cuts_dvcs_CDFT_Inb_2sigma["coneAngle_lb"] = [0.182, -0.611, 9.293]
cuts_dvcs_CDFT_Inb_2sigma["MPt_ub"] = 0.0431
cuts_dvcs_CDFT_Inb_2sigma["reconGam_ub"] = 0.475

cuts_dvcs_CDFT_Inb_4sigma = {}
cuts_dvcs_CDFT_Inb_4sigma["MM2_ep_ub"] = 0.517
cuts_dvcs_CDFT_Inb_4sigma["MM2_ep_lb"] = -0.491
cuts_dvcs_CDFT_Inb_4sigma["MM2_eg_ub"] = 1.672
cuts_dvcs_CDFT_Inb_4sigma["MM2_eg_lb"] = 0.128
cuts_dvcs_CDFT_Inb_4sigma["ME_epg_ub"] = 0.551
cuts_dvcs_CDFT_Inb_4sigma["ME_epg_lb"] = -0.527
cuts_dvcs_CDFT_Inb_4sigma["MM2_epg_ub"] = 0.000941
cuts_dvcs_CDFT_Inb_4sigma["MM2_epg_lb"] = -0.00551
cuts_dvcs_CDFT_Inb_4sigma["coplanarity_ub"] = 7.801
cuts_dvcs_CDFT_Inb_4sigma["coneAngle_ub"] = [0.297, -2.514, 23.955]
cuts_dvcs_CDFT_Inb_4sigma["coneAngle_lb"] = [0.0246, 1.736, -1.543]
cuts_dvcs_CDFT_Inb_4sigma["MPt_ub"] = 0.0856
cuts_dvcs_CDFT_Inb_4sigma["reconGam_ub"] = 0.827

cuts_dvcs_CD_Inb_3sigma = {}
cuts_dvcs_CD_Inb_3sigma["MM2_ep_ub"] = 0.285
cuts_dvcs_CD_Inb_3sigma["MM2_ep_lb"] = -0.276
cuts_dvcs_CD_Inb_3sigma["MM2_eg_ub"] = 1.479
cuts_dvcs_CD_Inb_3sigma["MM2_eg_lb"] = 0.322
cuts_dvcs_CD_Inb_3sigma["ME_epg_ub"] = 0.360
cuts_dvcs_CD_Inb_3sigma["ME_epg_lb"] = -0.387
cuts_dvcs_CD_Inb_3sigma["MM2_epg_ub"] = 0.000629
cuts_dvcs_CD_Inb_3sigma["MM2_epg_lb"] = -0.00423
cuts_dvcs_CD_Inb_3sigma["coplanarity_ub"] = 5.993
cuts_dvcs_CD_Inb_3sigma["coneAngle_ub"] = [ 0.133, -0.27, 15.103]
cuts_dvcs_CD_Inb_3sigma["coneAngle_lb"] = [ 0.179, -0.410, 7.000]
cuts_dvcs_CD_Inb_3sigma["MPt_ub"] = 0.0742
cuts_dvcs_CD_Inb_3sigma["reconGam_ub"] = 0.663

cuts_dvcs_CD_Inb_2sigma = {}
cuts_dvcs_CD_Inb_2sigma["MM2_ep_ub"] = 0.265
cuts_dvcs_CD_Inb_2sigma["MM2_ep_lb"] = -0.239
cuts_dvcs_CD_Inb_2sigma["MM2_eg_ub"] = 1.272
cuts_dvcs_CD_Inb_2sigma["MM2_eg_lb"] = 0.527
cuts_dvcs_CD_Inb_2sigma["ME_epg_ub"] = 0.255
cuts_dvcs_CD_Inb_2sigma["ME_epg_lb"] = -0.2232
cuts_dvcs_CD_Inb_2sigma["MM2_epg_ub"] = 0.000405
cuts_dvcs_CD_Inb_2sigma["MM2_epg_lb"] = -0.00360
cuts_dvcs_CD_Inb_2sigma["coplanarity_ub"] = 4.075
cuts_dvcs_CD_Inb_2sigma["coneAngle_ub"] = [0.0943, 0.295, 11.847]
cuts_dvcs_CD_Inb_2sigma["coneAngle_lb"] = [0.182, -0.611, 9.293]
cuts_dvcs_CD_Inb_2sigma["MPt_ub"] = 0.0431
cuts_dvcs_CD_Inb_2sigma["reconGam_ub"] = 0.475

cuts_dvcs_CD_Inb_4sigma = {}
cuts_dvcs_CD_Inb_4sigma["MM2_ep_ub"] = 0.517
cuts_dvcs_CD_Inb_4sigma["MM2_ep_lb"] = -0.491
cuts_dvcs_CD_Inb_4sigma["MM2_eg_ub"] = 1.672
cuts_dvcs_CD_Inb_4sigma["MM2_eg_lb"] = 0.128
cuts_dvcs_CD_Inb_4sigma["ME_epg_ub"] = 0.551
cuts_dvcs_CD_Inb_4sigma["ME_epg_lb"] = -0.527
cuts_dvcs_CD_Inb_4sigma["MM2_epg_ub"] = 0.000941
cuts_dvcs_CD_Inb_4sigma["MM2_epg_lb"] = -0.00551
cuts_dvcs_CD_Inb_4sigma["coplanarity_ub"] = 7.801
cuts_dvcs_CD_Inb_4sigma["coneAngle_ub"] = [0.297, -2.514, 23.955]
cuts_dvcs_CD_Inb_4sigma["coneAngle_lb"] = [0.0246, 1.736, -1.543]
cuts_dvcs_CD_Inb_4sigma["MPt_ub"] = 0.0856
cuts_dvcs_CD_Inb_4sigma["reconGam_ub"] = 0.827


cuts_dvcs_FD_Inb_3sigma = {}
cuts_dvcs_FD_Inb_3sigma["MM2_ep_ub"] = 0.285
cuts_dvcs_FD_Inb_3sigma["MM2_ep_lb"] = -0.276
cuts_dvcs_FD_Inb_3sigma["MM2_eg_ub"] = 1.479
cuts_dvcs_FD_Inb_3sigma["MM2_eg_lb"] = 0.322
cuts_dvcs_FD_Inb_3sigma["ME_epg_ub"] = 0.360
cuts_dvcs_FD_Inb_3sigma["ME_epg_lb"] = -0.387
cuts_dvcs_FD_Inb_3sigma["MM2_epg_ub"] = 0.000629
cuts_dvcs_FD_Inb_3sigma["MM2_epg_lb"] = -0.00423
cuts_dvcs_FD_Inb_3sigma["coplanarity_ub"] = 5.993
cuts_dvcs_FD_Inb_3sigma["coneAngle_ub"] = [ 0.133, -0.27, 15.103]
cuts_dvcs_FD_Inb_3sigma["coneAngle_lb"] = [ 0.179, -0.410, 7.000]
cuts_dvcs_FD_Inb_3sigma["MPt_ub"] = 0.0742
cuts_dvcs_FD_Inb_3sigma["reconGam_ub"] = 0.663

cuts_dvcs_FD_Inb_2sigma = {}
cuts_dvcs_FD_Inb_2sigma["MM2_ep_ub"] = 0.265
cuts_dvcs_FD_Inb_2sigma["MM2_ep_lb"] = -0.239
cuts_dvcs_FD_Inb_2sigma["MM2_eg_ub"] = 1.272
cuts_dvcs_FD_Inb_2sigma["MM2_eg_lb"] = 0.527
cuts_dvcs_FD_Inb_2sigma["ME_epg_ub"] = 0.255
cuts_dvcs_FD_Inb_2sigma["ME_epg_lb"] = -0.2232
cuts_dvcs_FD_Inb_2sigma["MM2_epg_ub"] = 0.000405
cuts_dvcs_FD_Inb_2sigma["MM2_epg_lb"] = -0.00360
cuts_dvcs_FD_Inb_2sigma["coplanarity_ub"] = 4.075
cuts_dvcs_FD_Inb_2sigma["coneAngle_ub"] = [0.0943, 0.295, 11.847]
cuts_dvcs_FD_Inb_2sigma["coneAngle_lb"] = [0.182, -0.611, 9.293]
cuts_dvcs_FD_Inb_2sigma["MPt_ub"] = 0.0431
cuts_dvcs_FD_Inb_2sigma["reconGam_ub"] = 0.475

cuts_dvcs_FD_Inb_4sigma = {}
cuts_dvcs_FD_Inb_4sigma["MM2_ep_ub"] = 0.517
cuts_dvcs_FD_Inb_4sigma["MM2_ep_lb"] = -0.491
cuts_dvcs_FD_Inb_4sigma["MM2_eg_ub"] = 1.672
cuts_dvcs_FD_Inb_4sigma["MM2_eg_lb"] = 0.128
cuts_dvcs_FD_Inb_4sigma["ME_epg_ub"] = 0.551
cuts_dvcs_FD_Inb_4sigma["ME_epg_lb"] = -0.527
cuts_dvcs_FD_Inb_4sigma["MM2_epg_ub"] = 0.000941
cuts_dvcs_FD_Inb_4sigma["MM2_epg_lb"] = -0.00551
cuts_dvcs_FD_Inb_4sigma["coplanarity_ub"] = 7.801
cuts_dvcs_FD_Inb_4sigma["coneAngle_ub"] = [0.297, -2.514, 23.955]
cuts_dvcs_FD_Inb_4sigma["coneAngle_lb"] = [0.0246, 1.736, -1.543]
cuts_dvcs_FD_Inb_4sigma["MPt_ub"] = 0.0856
cuts_dvcs_FD_Inb_4sigma["reconGam_ub"] = 0.827


cuts_dvcs_CDFT_Outb_3sigma = {}
cuts_dvcs_CDFT_Outb_3sigma["MM2_ep_ub"] = 0.391
cuts_dvcs_CDFT_Outb_3sigma["MM2_ep_lb"] = -0.365
cuts_dvcs_CDFT_Outb_3sigma["MM2_eg_ub"] = 1.479
cuts_dvcs_CDFT_Outb_3sigma["MM2_eg_lb"] = 0.322
cuts_dvcs_CDFT_Outb_3sigma["ME_epg_ub"] = 0.360
cuts_dvcs_CDFT_Outb_3sigma["ME_epg_lb"] = -0.387
cuts_dvcs_CDFT_Outb_3sigma["MM2_epg_ub"] = 0.000629
cuts_dvcs_CDFT_Outb_3sigma["MM2_epg_lb"] = -0.00423
cuts_dvcs_CDFT_Outb_3sigma["coplanarity_ub"] = 5.993
cuts_dvcs_CDFT_Outb_3sigma["coneAngle_ub"] = [ 0.133, -0.27, 15.103]
cuts_dvcs_CDFT_Outb_3sigma["coneAngle_lb"] = [ 0.179, -0.410, 7.000]
cuts_dvcs_CDFT_Outb_3sigma["MPt_ub"] = 0.0742
cuts_dvcs_CDFT_Outb_3sigma["reconGam_ub"] = 0.663

cuts_dvcs_CDFT_Outb_2sigma = {}
cuts_dvcs_CDFT_Outb_2sigma["MM2_ep_ub"] = 0.265
cuts_dvcs_CDFT_Outb_2sigma["MM2_ep_lb"] = -0.239
cuts_dvcs_CDFT_Outb_2sigma["MM2_eg_ub"] = 1.272
cuts_dvcs_CDFT_Outb_2sigma["MM2_eg_lb"] = 0.527
cuts_dvcs_CDFT_Outb_2sigma["ME_epg_ub"] = 0.255
cuts_dvcs_CDFT_Outb_2sigma["ME_epg_lb"] = -0.2232
cuts_dvcs_CDFT_Outb_2sigma["MM2_epg_ub"] = 0.000405
cuts_dvcs_CDFT_Outb_2sigma["MM2_epg_lb"] = -0.00360
cuts_dvcs_CDFT_Outb_2sigma["coplanarity_ub"] = 4.075
cuts_dvcs_CDFT_Outb_2sigma["coneAngle_ub"] = [0.0943, 0.295, 11.847]
cuts_dvcs_CDFT_Outb_2sigma["coneAngle_lb"] = [0.182, -0.611, 9.293]
cuts_dvcs_CDFT_Outb_2sigma["MPt_ub"] = 0.0431
cuts_dvcs_CDFT_Outb_2sigma["reconGam_ub"] = 0.475

cuts_dvcs_CDFT_Outb_4sigma = {}
cuts_dvcs_CDFT_Outb_4sigma["MM2_ep_ub"] = 0.517
cuts_dvcs_CDFT_Outb_4sigma["MM2_ep_lb"] = -0.491
cuts_dvcs_CDFT_Outb_4sigma["MM2_eg_ub"] = 1.672
cuts_dvcs_CDFT_Outb_4sigma["MM2_eg_lb"] = 0.128
cuts_dvcs_CDFT_Outb_4sigma["ME_epg_ub"] = 0.551
cuts_dvcs_CDFT_Outb_4sigma["ME_epg_lb"] = -0.527
cuts_dvcs_CDFT_Outb_4sigma["MM2_epg_ub"] = 0.000941
cuts_dvcs_CDFT_Outb_4sigma["MM2_epg_lb"] = -0.00551
cuts_dvcs_CDFT_Outb_4sigma["coplanarity_ub"] = 7.801
cuts_dvcs_CDFT_Outb_4sigma["coneAngle_ub"] = [0.297, -2.514, 23.955]
cuts_dvcs_CDFT_Outb_4sigma["coneAngle_lb"] = [0.0246, 1.736, -1.543]
cuts_dvcs_CDFT_Outb_4sigma["MPt_ub"] = 0.0856
cuts_dvcs_CDFT_Outb_4sigma["reconGam_ub"] = 0.827

cuts_dvcs_CD_Outb_3sigma = {}
cuts_dvcs_CD_Outb_3sigma["MM2_ep_ub"] = 0.285
cuts_dvcs_CD_Outb_3sigma["MM2_ep_lb"] = -0.276
cuts_dvcs_CD_Outb_3sigma["MM2_eg_ub"] = 1.479
cuts_dvcs_CD_Outb_3sigma["MM2_eg_lb"] = 0.322
cuts_dvcs_CD_Outb_3sigma["ME_epg_ub"] = 0.360
cuts_dvcs_CD_Outb_3sigma["ME_epg_lb"] = -0.387
cuts_dvcs_CD_Outb_3sigma["MM2_epg_ub"] = 0.000629
cuts_dvcs_CD_Outb_3sigma["MM2_epg_lb"] = -0.00423
cuts_dvcs_CD_Outb_3sigma["coplanarity_ub"] = 5.993
cuts_dvcs_CD_Outb_3sigma["coneAngle_ub"] = [ 0.133, -0.27, 15.103]
cuts_dvcs_CD_Outb_3sigma["coneAngle_lb"] = [ 0.179, -0.410, 7.000]
cuts_dvcs_CD_Outb_3sigma["MPt_ub"] = 0.0742
cuts_dvcs_CD_Outb_3sigma["reconGam_ub"] = 0.663

cuts_dvcs_CD_Outb_2sigma = {}
cuts_dvcs_CD_Outb_2sigma["MM2_ep_ub"] = 0.265
cuts_dvcs_CD_Outb_2sigma["MM2_ep_lb"] = -0.239
cuts_dvcs_CD_Outb_2sigma["MM2_eg_ub"] = 1.272
cuts_dvcs_CD_Outb_2sigma["MM2_eg_lb"] = 0.527
cuts_dvcs_CD_Outb_2sigma["ME_epg_ub"] = 0.255
cuts_dvcs_CD_Outb_2sigma["ME_epg_lb"] = -0.2232
cuts_dvcs_CD_Outb_2sigma["MM2_epg_ub"] = 0.000405
cuts_dvcs_CD_Outb_2sigma["MM2_epg_lb"] = -0.00360
cuts_dvcs_CD_Outb_2sigma["coplanarity_ub"] = 4.075
cuts_dvcs_CD_Outb_2sigma["coneAngle_ub"] = [0.0943, 0.295, 11.847]
cuts_dvcs_CD_Outb_2sigma["coneAngle_lb"] = [0.182, -0.611, 9.293]
cuts_dvcs_CD_Outb_2sigma["MPt_ub"] = 0.0431
cuts_dvcs_CD_Outb_2sigma["reconGam_ub"] = 0.475

cuts_dvcs_CD_Outb_4sigma = {}
cuts_dvcs_CD_Outb_4sigma["MM2_ep_ub"] = 0.517
cuts_dvcs_CD_Outb_4sigma["MM2_ep_lb"] = -0.491
cuts_dvcs_CD_Outb_4sigma["MM2_eg_ub"] = 1.672
cuts_dvcs_CD_Outb_4sigma["MM2_eg_lb"] = 0.128
cuts_dvcs_CD_Outb_4sigma["ME_epg_ub"] = 0.551
cuts_dvcs_CD_Outb_4sigma["ME_epg_lb"] = -0.527
cuts_dvcs_CD_Outb_4sigma["MM2_epg_ub"] = 0.000941
cuts_dvcs_CD_Outb_4sigma["MM2_epg_lb"] = -0.00551
cuts_dvcs_CD_Outb_4sigma["coplanarity_ub"] = 7.801
cuts_dvcs_CD_Outb_4sigma["coneAngle_ub"] = [0.297, -2.514, 23.955]
cuts_dvcs_CD_Outb_4sigma["coneAngle_lb"] = [0.0246, 1.736, -1.543]
cuts_dvcs_CD_Outb_4sigma["MPt_ub"] = 0.0856
cuts_dvcs_CD_Outb_4sigma["reconGam_ub"] = 0.827


cuts_dvcs_FD_Outb_3sigma = {}
cuts_dvcs_FD_Outb_3sigma["MM2_ep_ub"] = 0.285
cuts_dvcs_FD_Outb_3sigma["MM2_ep_lb"] = -0.276
cuts_dvcs_FD_Outb_3sigma["MM2_eg_ub"] = 1.479
cuts_dvcs_FD_Outb_3sigma["MM2_eg_lb"] = 0.322
cuts_dvcs_FD_Outb_3sigma["ME_epg_ub"] = 0.360
cuts_dvcs_FD_Outb_3sigma["ME_epg_lb"] = -0.387
cuts_dvcs_FD_Outb_3sigma["MM2_epg_ub"] = 0.000629
cuts_dvcs_FD_Outb_3sigma["MM2_epg_lb"] = -0.00423
cuts_dvcs_FD_Outb_3sigma["coplanarity_ub"] = 5.993
cuts_dvcs_FD_Outb_3sigma["coneAngle_ub"] = [ 0.133, -0.27, 15.103]
cuts_dvcs_FD_Outb_3sigma["coneAngle_lb"] = [ 0.179, -0.410, 7.000]
cuts_dvcs_FD_Outb_3sigma["MPt_ub"] = 0.0742
cuts_dvcs_FD_Outb_3sigma["reconGam_ub"] = 0.663

cuts_dvcs_FD_Outb_2sigma = {}
cuts_dvcs_FD_Outb_2sigma["MM2_ep_ub"] = 0.265
cuts_dvcs_FD_Outb_2sigma["MM2_ep_lb"] = -0.239
cuts_dvcs_FD_Outb_2sigma["MM2_eg_ub"] = 1.272
cuts_dvcs_FD_Outb_2sigma["MM2_eg_lb"] = 0.527
cuts_dvcs_FD_Outb_2sigma["ME_epg_ub"] = 0.255
cuts_dvcs_FD_Outb_2sigma["ME_epg_lb"] = -0.2232
cuts_dvcs_FD_Outb_2sigma["MM2_epg_ub"] = 0.000405
cuts_dvcs_FD_Outb_2sigma["MM2_epg_lb"] = -0.00360
cuts_dvcs_FD_Outb_2sigma["coplanarity_ub"] = 4.075
cuts_dvcs_FD_Outb_2sigma["coneAngle_ub"] = [0.0943, 0.295, 11.847]
cuts_dvcs_FD_Outb_2sigma["coneAngle_lb"] = [0.182, -0.611, 9.293]
cuts_dvcs_FD_Outb_2sigma["MPt_ub"] = 0.0431
cuts_dvcs_FD_Outb_2sigma["reconGam_ub"] = 0.475

cuts_dvcs_FD_Outb_4sigma = {}
cuts_dvcs_FD_Outb_4sigma["MM2_ep_ub"] = 0.517
cuts_dvcs_FD_Outb_4sigma["MM2_ep_lb"] = -0.491
cuts_dvcs_FD_Outb_4sigma["MM2_eg_ub"] = 1.672
cuts_dvcs_FD_Outb_4sigma["MM2_eg_lb"] = 0.128
cuts_dvcs_FD_Outb_4sigma["ME_epg_ub"] = 0.551
cuts_dvcs_FD_Outb_4sigma["ME_epg_lb"] = -0.527
cuts_dvcs_FD_Outb_4sigma["MM2_epg_ub"] = 0.000941
cuts_dvcs_FD_Outb_4sigma["MM2_epg_lb"] = -0.00551
cuts_dvcs_FD_Outb_4sigma["coplanarity_ub"] = 7.801
cuts_dvcs_FD_Outb_4sigma["coneAngle_ub"] = [0.297, -2.514, 23.955]
cuts_dvcs_FD_Outb_4sigma["coneAngle_lb"] = [0.0246, 1.736, -1.543]
cuts_dvcs_FD_Outb_4sigma["MPt_ub"] = 0.0856
cuts_dvcs_FD_Outb_4sigma["reconGam_ub"] = 0.827
