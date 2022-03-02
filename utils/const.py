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
cuts_dvcs_CDFT_Inb_3sigma["coplanarity_ub"] = 5.633
cuts_dvcs_CDFT_Inb_3sigma["coneAngle_ub"] = [-0.106, 2.940, 5.527]
cuts_dvcs_CDFT_Inb_3sigma["coneAngle_lb"] = [0.434, -3.766, 16.994]
cuts_dvcs_CDFT_Inb_3sigma["MPt_ub"] = 0.0844
cuts_dvcs_CDFT_Inb_3sigma["reconGam_ub"] = 0.582
cuts_dvcs_CDFT_Inb_3sigma["MM2_epg_ub"] = 0.00800
cuts_dvcs_CDFT_Inb_3sigma["MM2_epg_lb"] = -0.0108

cuts_dvcs_CDFT_Inb_2sigma = {}
cuts_dvcs_CDFT_Inb_2sigma["MM2_ep_ub"] = 0.265
cuts_dvcs_CDFT_Inb_2sigma["MM2_ep_lb"] = -0.239
cuts_dvcs_CDFT_Inb_2sigma["MM2_eg_ub"] = 1.272
cuts_dvcs_CDFT_Inb_2sigma["MM2_eg_lb"] = 0.527
cuts_dvcs_CDFT_Inb_2sigma["ME_epg_ub"] = 0.255
cuts_dvcs_CDFT_Inb_2sigma["ME_epg_lb"] = -0.2232
cuts_dvcs_CDFT_Inb_2sigma["coplanarity_ub"] = 3.842
cuts_dvcs_CDFT_Inb_2sigma["coneAngle_ub"] = [0.0134, 1.399, 8.469]
cuts_dvcs_CDFT_Inb_2sigma["coneAngle_lb"] = [0.258, -1.632, 12.408]
cuts_dvcs_CDFT_Inb_2sigma["MPt_ub"] = 0.0458
cuts_dvcs_CDFT_Inb_2sigma["reconGam_ub"] = 0.574
cuts_dvcs_CDFT_Inb_2sigma["MM2_epg_ub"] = 0.00292
cuts_dvcs_CDFT_Inb_2sigma["MM2_epg_lb"] = -0.00439

cuts_dvcs_CDFT_Inb_4sigma = {}
cuts_dvcs_CDFT_Inb_4sigma["MM2_ep_ub"] = 0.517
cuts_dvcs_CDFT_Inb_4sigma["MM2_ep_lb"] = -0.491
cuts_dvcs_CDFT_Inb_4sigma["MM2_eg_ub"] = 1.672
cuts_dvcs_CDFT_Inb_4sigma["MM2_eg_lb"] = 0.128
cuts_dvcs_CDFT_Inb_4sigma["ME_epg_ub"] = 0.551
cuts_dvcs_CDFT_Inb_4sigma["ME_epg_lb"] = -0.527
cuts_dvcs_CDFT_Inb_4sigma["coplanarity_ub"] = 7.384
cuts_dvcs_CDFT_Inb_4sigma["coneAngle_ub"] = [-0.234, 4.460, 2.749]
cuts_dvcs_CDFT_Inb_4sigma["coneAngle_lb"] = [0.557, -5.221, 19.651]
cuts_dvcs_CDFT_Inb_4sigma["MPt_ub"] = 0.140
cuts_dvcs_CDFT_Inb_4sigma["reconGam_ub"] = 0.778
cuts_dvcs_CDFT_Inb_4sigma["MM2_epg_ub"] = 0.0157
cuts_dvcs_CDFT_Inb_4sigma["MM2_epg_lb"] = -0.0194

cuts_dvcs_CD_Inb_3sigma = {}
cuts_dvcs_CD_Inb_3sigma["MM2_ep_ub"] = 0.294
cuts_dvcs_CD_Inb_3sigma["MM2_ep_lb"] = -0.272
cuts_dvcs_CD_Inb_3sigma["MM2_eg_ub"] = 1.790
cuts_dvcs_CD_Inb_3sigma["MM2_eg_lb"] = 0.166
cuts_dvcs_CD_Inb_3sigma["ME_epg_ub"] = 0.672
cuts_dvcs_CD_Inb_3sigma["ME_epg_lb"] = -0.557
cuts_dvcs_CD_Inb_3sigma["coplanarity_ub"] = 5.358
cuts_dvcs_CD_Inb_3sigma["coneAngle_ub"] = [ 0.400, -5.193, 50.848]
cuts_dvcs_CD_Inb_3sigma["coneAngle_lb"] = [ 0.0674, -0.111, 9.876]
cuts_dvcs_CD_Inb_3sigma["MPt_ub"] = 0.201
cuts_dvcs_CD_Inb_3sigma["reconGam_ub"] = 0.669
cuts_dvcs_CD_Inb_3sigma["MM2_epg_ub"] = 0.00063
cuts_dvcs_CD_Inb_3sigma["MM2_epg_lb"] = -0.00431

cuts_dvcs_CD_Inb_2sigma = {}
cuts_dvcs_CD_Inb_2sigma["MM2_ep_ub"] = 0.200
cuts_dvcs_CD_Inb_2sigma["MM2_ep_lb"] = -0.178
cuts_dvcs_CD_Inb_2sigma["MM2_eg_ub"] = 1.542
cuts_dvcs_CD_Inb_2sigma["MM2_eg_lb"] = 0.402
cuts_dvcs_CD_Inb_2sigma["ME_epg_ub"] = 1.088
cuts_dvcs_CD_Inb_2sigma["ME_epg_lb"] = -0.933
cuts_dvcs_CD_Inb_2sigma["coplanarity_ub"] = 3.668
cuts_dvcs_CD_Inb_2sigma["coneAngle_ub"] = [0.178, -2.464, 37.165]
cuts_dvcs_CD_Inb_2sigma["coneAngle_lb"] = [-0.360, 4.472, 4.578]
cuts_dvcs_CD_Inb_2sigma["MPt_ub"] = 0.0749
cuts_dvcs_CD_Inb_2sigma["reconGam_ub"] = 0.477
cuts_dvcs_CD_Inb_2sigma["MM2_epg_ub"] = 0.000286
cuts_dvcs_CD_Inb_2sigma["MM2_epg_lb"] = -0.00308

cuts_dvcs_CD_Inb_4sigma = {}
cuts_dvcs_CD_Inb_4sigma["MM2_ep_ub"] = 0.388
cuts_dvcs_CD_Inb_4sigma["MM2_ep_lb"] = -0.366
cuts_dvcs_CD_Inb_4sigma["MM2_eg_ub"] = 2.074
cuts_dvcs_CD_Inb_4sigma["MM2_eg_lb"] = -0.120
cuts_dvcs_CD_Inb_4sigma["ME_epg_ub"] = 0.990
cuts_dvcs_CD_Inb_4sigma["ME_epg_lb"] = -0.845
cuts_dvcs_CD_Inb_4sigma["coplanarity_ub"] = 7.036
cuts_dvcs_CD_Inb_4sigma["coneAngle_ub"] = [0.779, -10.133, 66.560]
cuts_dvcs_CD_Inb_4sigma["coneAngle_lb"] = [-0.376, 5.057, -4.611]
cuts_dvcs_CD_Inb_4sigma["MPt_ub"] = 0.0952
cuts_dvcs_CD_Inb_4sigma["reconGam_ub"] = 0.899
cuts_dvcs_CD_Inb_4sigma["MM2_epg_ub"] = 0.000929
cuts_dvcs_CD_Inb_4sigma["MM2_epg_lb"] = -0.00560


cuts_dvcs_FD_Inb_3sigma = {}
cuts_dvcs_FD_Inb_3sigma["MM2_ep_ub"] = 0.285
cuts_dvcs_FD_Inb_3sigma["MM2_ep_lb"] = -0.276
cuts_dvcs_FD_Inb_3sigma["MM2_eg_ub"] = 1.479
cuts_dvcs_FD_Inb_3sigma["MM2_eg_lb"] = 0.322
cuts_dvcs_FD_Inb_3sigma["ME_epg_ub"] = 0.360
cuts_dvcs_FD_Inb_3sigma["ME_epg_lb"] = -0.387
cuts_dvcs_FD_Inb_3sigma["coplanarity_ub"] = 5.993
cuts_dvcs_FD_Inb_3sigma["coneAngle_ub"] = [ 0.133, -0.27, 15.103]
cuts_dvcs_FD_Inb_3sigma["coneAngle_lb"] = [ 0.179, -0.410, 7.000]
cuts_dvcs_FD_Inb_3sigma["MPt_ub"] = 0.0742
cuts_dvcs_FD_Inb_3sigma["reconGam_ub"] = 0.663
cuts_dvcs_FD_Inb_3sigma["MM2_epg_ub"] = 0.000629
cuts_dvcs_FD_Inb_3sigma["MM2_epg_lb"] = -0.00423

cuts_dvcs_FD_Inb_2sigma = {}
cuts_dvcs_FD_Inb_2sigma["MM2_ep_ub"] = 0.265
cuts_dvcs_FD_Inb_2sigma["MM2_ep_lb"] = -0.239
cuts_dvcs_FD_Inb_2sigma["MM2_eg_ub"] = 1.272
cuts_dvcs_FD_Inb_2sigma["MM2_eg_lb"] = 0.527
cuts_dvcs_FD_Inb_2sigma["ME_epg_ub"] = 0.255
cuts_dvcs_FD_Inb_2sigma["ME_epg_lb"] = -0.2232
cuts_dvcs_FD_Inb_2sigma["coplanarity_ub"] = 4.075
cuts_dvcs_FD_Inb_2sigma["coneAngle_ub"] = [0.0943, 0.295, 11.847]
cuts_dvcs_FD_Inb_2sigma["coneAngle_lb"] = [0.182, -0.611, 9.293]
cuts_dvcs_FD_Inb_2sigma["MPt_ub"] = 0.0431
cuts_dvcs_FD_Inb_2sigma["reconGam_ub"] = 0.475
cuts_dvcs_FD_Inb_2sigma["MM2_epg_ub"] = 0.000405
cuts_dvcs_FD_Inb_2sigma["MM2_epg_lb"] = -0.00360

cuts_dvcs_FD_Inb_4sigma = {}
cuts_dvcs_FD_Inb_4sigma["MM2_ep_ub"] = 0.517
cuts_dvcs_FD_Inb_4sigma["MM2_ep_lb"] = -0.491
cuts_dvcs_FD_Inb_4sigma["MM2_eg_ub"] = 1.672
cuts_dvcs_FD_Inb_4sigma["MM2_eg_lb"] = 0.128
cuts_dvcs_FD_Inb_4sigma["ME_epg_ub"] = 0.551
cuts_dvcs_FD_Inb_4sigma["ME_epg_lb"] = -0.527
cuts_dvcs_FD_Inb_4sigma["coplanarity_ub"] = 7.801
cuts_dvcs_FD_Inb_4sigma["coneAngle_ub"] = [0.297, -2.514, 23.955]
cuts_dvcs_FD_Inb_4sigma["coneAngle_lb"] = [0.0246, 1.736, -1.543]
cuts_dvcs_FD_Inb_4sigma["MPt_ub"] = 0.0856
cuts_dvcs_FD_Inb_4sigma["reconGam_ub"] = 0.827
cuts_dvcs_FD_Inb_4sigma["MM2_epg_ub"] = 0.000941
cuts_dvcs_FD_Inb_4sigma["MM2_epg_lb"] = -0.00551


cuts_dvcs_CDFT_Outb_3sigma = {}
cuts_dvcs_CDFT_Outb_3sigma["MM2_ep_ub"] = 0.391
cuts_dvcs_CDFT_Outb_3sigma["MM2_ep_lb"] = -0.365
cuts_dvcs_CDFT_Outb_3sigma["MM2_eg_ub"] = 1.479
cuts_dvcs_CDFT_Outb_3sigma["MM2_eg_lb"] = 0.322
cuts_dvcs_CDFT_Outb_3sigma["ME_epg_ub"] = 0.360
cuts_dvcs_CDFT_Outb_3sigma["ME_epg_lb"] = -0.387
cuts_dvcs_CDFT_Outb_3sigma["coplanarity_ub"] = 5.993
cuts_dvcs_CDFT_Outb_3sigma["coneAngle_ub"] = [ 0.133, -0.27, 15.103]
cuts_dvcs_CDFT_Outb_3sigma["coneAngle_lb"] = [ 0.179, -0.410, 7.000]
cuts_dvcs_CDFT_Outb_3sigma["MPt_ub"] = 0.0742
cuts_dvcs_CDFT_Outb_3sigma["reconGam_ub"] = 0.663
cuts_dvcs_CDFT_Outb_3sigma["MM2_epg_ub"] = 0.000629
cuts_dvcs_CDFT_Outb_3sigma["MM2_epg_lb"] = -0.00423

cuts_dvcs_CDFT_Outb_2sigma = {}
cuts_dvcs_CDFT_Outb_2sigma["MM2_ep_ub"] = 0.265
cuts_dvcs_CDFT_Outb_2sigma["MM2_ep_lb"] = -0.239
cuts_dvcs_CDFT_Outb_2sigma["MM2_eg_ub"] = 1.272
cuts_dvcs_CDFT_Outb_2sigma["MM2_eg_lb"] = 0.527
cuts_dvcs_CDFT_Outb_2sigma["ME_epg_ub"] = 0.255
cuts_dvcs_CDFT_Outb_2sigma["ME_epg_lb"] = -0.2232
cuts_dvcs_CDFT_Outb_2sigma["coplanarity_ub"] = 4.075
cuts_dvcs_CDFT_Outb_2sigma["coneAngle_ub"] = [0.0943, 0.295, 11.847]
cuts_dvcs_CDFT_Outb_2sigma["coneAngle_lb"] = [0.182, -0.611, 9.293]
cuts_dvcs_CDFT_Outb_2sigma["MPt_ub"] = 0.0431
cuts_dvcs_CDFT_Outb_2sigma["reconGam_ub"] = 0.475
cuts_dvcs_CDFT_Outb_2sigma["MM2_epg_ub"] = 0.000405
cuts_dvcs_CDFT_Outb_2sigma["MM2_epg_lb"] = -0.00360

cuts_dvcs_CDFT_Outb_4sigma = {}
cuts_dvcs_CDFT_Outb_4sigma["MM2_ep_ub"] = 0.517
cuts_dvcs_CDFT_Outb_4sigma["MM2_ep_lb"] = -0.491
cuts_dvcs_CDFT_Outb_4sigma["MM2_eg_ub"] = 1.672
cuts_dvcs_CDFT_Outb_4sigma["MM2_eg_lb"] = 0.128
cuts_dvcs_CDFT_Outb_4sigma["ME_epg_ub"] = 0.551
cuts_dvcs_CDFT_Outb_4sigma["ME_epg_lb"] = -0.527
cuts_dvcs_CDFT_Outb_4sigma["coplanarity_ub"] = 7.801
cuts_dvcs_CDFT_Outb_4sigma["coneAngle_ub"] = [0.297, -2.514, 23.955]
cuts_dvcs_CDFT_Outb_4sigma["coneAngle_lb"] = [0.0246, 1.736, -1.543]
cuts_dvcs_CDFT_Outb_4sigma["MPt_ub"] = 0.0856
cuts_dvcs_CDFT_Outb_4sigma["reconGam_ub"] = 0.827
cuts_dvcs_CDFT_Outb_4sigma["MM2_epg_ub"] = 0.000941
cuts_dvcs_CDFT_Outb_4sigma["MM2_epg_lb"] = -0.00551

cuts_dvcs_CD_Outb_3sigma = {}
cuts_dvcs_CD_Outb_3sigma["MM2_ep_ub"] = 0.285
cuts_dvcs_CD_Outb_3sigma["MM2_ep_lb"] = -0.276
cuts_dvcs_CD_Outb_3sigma["MM2_eg_ub"] = 1.479
cuts_dvcs_CD_Outb_3sigma["MM2_eg_lb"] = 0.322
cuts_dvcs_CD_Outb_3sigma["ME_epg_ub"] = 0.360
cuts_dvcs_CD_Outb_3sigma["ME_epg_lb"] = -0.387
cuts_dvcs_CD_Outb_3sigma["coplanarity_ub"] = 5.993
cuts_dvcs_CD_Outb_3sigma["coneAngle_ub"] = [ 0.133, -0.27, 15.103]
cuts_dvcs_CD_Outb_3sigma["coneAngle_lb"] = [ 0.179, -0.410, 7.000]
cuts_dvcs_CD_Outb_3sigma["MPt_ub"] = 0.0742
cuts_dvcs_CD_Outb_3sigma["reconGam_ub"] = 0.663
cuts_dvcs_CD_Outb_3sigma["MM2_epg_ub"] = 0.000629
cuts_dvcs_CD_Outb_3sigma["MM2_epg_lb"] = -0.00423

cuts_dvcs_CD_Outb_2sigma = {}
cuts_dvcs_CD_Outb_2sigma["MM2_ep_ub"] = 0.265
cuts_dvcs_CD_Outb_2sigma["MM2_ep_lb"] = -0.239
cuts_dvcs_CD_Outb_2sigma["MM2_eg_ub"] = 1.272
cuts_dvcs_CD_Outb_2sigma["MM2_eg_lb"] = 0.527
cuts_dvcs_CD_Outb_2sigma["ME_epg_ub"] = 0.255
cuts_dvcs_CD_Outb_2sigma["ME_epg_lb"] = -0.2232
cuts_dvcs_CD_Outb_2sigma["coplanarity_ub"] = 4.075
cuts_dvcs_CD_Outb_2sigma["coneAngle_ub"] = [0.0943, 0.295, 11.847]
cuts_dvcs_CD_Outb_2sigma["coneAngle_lb"] = [0.182, -0.611, 9.293]
cuts_dvcs_CD_Outb_2sigma["MPt_ub"] = 0.0431
cuts_dvcs_CD_Outb_2sigma["reconGam_ub"] = 0.475
cuts_dvcs_CD_Outb_2sigma["MM2_epg_ub"] = 0.000405
cuts_dvcs_CD_Outb_2sigma["MM2_epg_lb"] = -0.00360

cuts_dvcs_CD_Outb_4sigma = {}
cuts_dvcs_CD_Outb_4sigma["MM2_ep_ub"] = 0.517
cuts_dvcs_CD_Outb_4sigma["MM2_ep_lb"] = -0.491
cuts_dvcs_CD_Outb_4sigma["MM2_eg_ub"] = 1.672
cuts_dvcs_CD_Outb_4sigma["MM2_eg_lb"] = 0.128
cuts_dvcs_CD_Outb_4sigma["ME_epg_ub"] = 0.551
cuts_dvcs_CD_Outb_4sigma["ME_epg_lb"] = -0.527
cuts_dvcs_CD_Outb_4sigma["coplanarity_ub"] = 7.801
cuts_dvcs_CD_Outb_4sigma["coneAngle_ub"] = [0.297, -2.514, 23.955]
cuts_dvcs_CD_Outb_4sigma["coneAngle_lb"] = [0.0246, 1.736, -1.543]
cuts_dvcs_CD_Outb_4sigma["MPt_ub"] = 0.0856
cuts_dvcs_CD_Outb_4sigma["reconGam_ub"] = 0.827
cuts_dvcs_CD_Outb_4sigma["MM2_epg_ub"] = 0.000941
cuts_dvcs_CD_Outb_4sigma["MM2_epg_lb"] = -0.00551


cuts_dvcs_FD_Outb_3sigma = {}
cuts_dvcs_FD_Outb_3sigma["MM2_ep_ub"] = 0.285
cuts_dvcs_FD_Outb_3sigma["MM2_ep_lb"] = -0.276
cuts_dvcs_FD_Outb_3sigma["MM2_eg_ub"] = 1.479
cuts_dvcs_FD_Outb_3sigma["MM2_eg_lb"] = 0.322
cuts_dvcs_FD_Outb_3sigma["ME_epg_ub"] = 0.360
cuts_dvcs_FD_Outb_3sigma["ME_epg_lb"] = -0.387
cuts_dvcs_FD_Outb_3sigma["coplanarity_ub"] = 5.993
cuts_dvcs_FD_Outb_3sigma["coneAngle_ub"] = [ 0.133, -0.27, 15.103]
cuts_dvcs_FD_Outb_3sigma["coneAngle_lb"] = [ 0.179, -0.410, 7.000]
cuts_dvcs_FD_Outb_3sigma["MPt_ub"] = 0.0742
cuts_dvcs_FD_Outb_3sigma["reconGam_ub"] = 0.663
cuts_dvcs_FD_Outb_3sigma["MM2_epg_ub"] = 0.000629
cuts_dvcs_FD_Outb_3sigma["MM2_epg_lb"] = -0.00423

cuts_dvcs_FD_Outb_2sigma = {}
cuts_dvcs_FD_Outb_2sigma["MM2_ep_ub"] = 0.265
cuts_dvcs_FD_Outb_2sigma["MM2_ep_lb"] = -0.239
cuts_dvcs_FD_Outb_2sigma["MM2_eg_ub"] = 1.272
cuts_dvcs_FD_Outb_2sigma["MM2_eg_lb"] = 0.527
cuts_dvcs_FD_Outb_2sigma["ME_epg_ub"] = 0.255
cuts_dvcs_FD_Outb_2sigma["ME_epg_lb"] = -0.2232
cuts_dvcs_FD_Outb_2sigma["coplanarity_ub"] = 4.075
cuts_dvcs_FD_Outb_2sigma["coneAngle_ub"] = [0.0943, 0.295, 11.847]
cuts_dvcs_FD_Outb_2sigma["coneAngle_lb"] = [0.182, -0.611, 9.293]
cuts_dvcs_FD_Outb_2sigma["MPt_ub"] = 0.0431
cuts_dvcs_FD_Outb_2sigma["reconGam_ub"] = 0.475
cuts_dvcs_FD_Outb_2sigma["MM2_epg_ub"] = 0.000405
cuts_dvcs_FD_Outb_2sigma["MM2_epg_lb"] = -0.00360

cuts_dvcs_FD_Outb_4sigma = {}
cuts_dvcs_FD_Outb_4sigma["MM2_ep_ub"] = 0.517
cuts_dvcs_FD_Outb_4sigma["MM2_ep_lb"] = -0.491
cuts_dvcs_FD_Outb_4sigma["MM2_eg_ub"] = 1.672
cuts_dvcs_FD_Outb_4sigma["MM2_eg_lb"] = 0.128
cuts_dvcs_FD_Outb_4sigma["ME_epg_ub"] = 0.551
cuts_dvcs_FD_Outb_4sigma["ME_epg_lb"] = -0.527
cuts_dvcs_FD_Outb_4sigma["coplanarity_ub"] = 7.801
cuts_dvcs_FD_Outb_4sigma["coneAngle_ub"] = [0.297, -2.514, 23.955]
cuts_dvcs_FD_Outb_4sigma["coneAngle_lb"] = [0.0246, 1.736, -1.543]
cuts_dvcs_FD_Outb_4sigma["MPt_ub"] = 0.0856
cuts_dvcs_FD_Outb_4sigma["reconGam_ub"] = 0.827
cuts_dvcs_FD_Outb_4sigma["MM2_epg_ub"] = 0.000941
cuts_dvcs_FD_Outb_4sigma["MM2_epg_lb"] = -0.00551
