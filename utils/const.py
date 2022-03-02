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
cuts_dvcs_CD_Inb_3sigma["coplanarity_ub"] = 5.034
cuts_dvcs_CD_Inb_3sigma["coneAngle_ub"] = [ 0.541, -6.825, 49.430]
cuts_dvcs_CD_Inb_3sigma["coneAngle_lb"] = [ 0.161, -1.652, 20.548]
cuts_dvcs_CD_Inb_3sigma["MPt_ub"] = 0.0919
cuts_dvcs_CD_Inb_3sigma["reconGam_ub"] = 0.654
cuts_dvcs_CD_Inb_3sigma["MM2_epg_ub"] = 0.0139
cuts_dvcs_CD_Inb_3sigma["MM2_epg_lb"] = -0.0161

cuts_dvcs_CD_Inb_2sigma = {}
cuts_dvcs_CD_Inb_2sigma["MM2_ep_ub"] = 0.200
cuts_dvcs_CD_Inb_2sigma["MM2_ep_lb"] = -0.178
cuts_dvcs_CD_Inb_2sigma["MM2_eg_ub"] = 1.542
cuts_dvcs_CD_Inb_2sigma["MM2_eg_lb"] = 0.402
cuts_dvcs_CD_Inb_2sigma["ME_epg_ub"] = 1.088
cuts_dvcs_CD_Inb_2sigma["ME_epg_lb"] = -0.933
cuts_dvcs_CD_Inb_2sigma["coplanarity_ub"] = 3.452
cuts_dvcs_CD_Inb_2sigma["coneAngle_ub"] = [0.540, -6.842, 49.771]
cuts_dvcs_CD_Inb_2sigma["coneAngle_lb"] = [-0.0386, 0.622, 14.835]
cuts_dvcs_CD_Inb_2sigma["MPt_ub"] = 0.100
cuts_dvcs_CD_Inb_2sigma["reconGam_ub"] = 0.443
cuts_dvcs_CD_Inb_2sigma["MM2_epg_ub"] = 0.00589
cuts_dvcs_CD_Inb_2sigma["MM2_epg_lb"] = -0.00743

cuts_dvcs_CD_Inb_4sigma = {}
cuts_dvcs_CD_Inb_4sigma["MM2_ep_ub"] = 0.388
cuts_dvcs_CD_Inb_4sigma["MM2_ep_lb"] = -0.366
cuts_dvcs_CD_Inb_4sigma["MM2_eg_ub"] = 2.074
cuts_dvcs_CD_Inb_4sigma["MM2_eg_lb"] = -0.120
cuts_dvcs_CD_Inb_4sigma["ME_epg_ub"] = 0.990
cuts_dvcs_CD_Inb_4sigma["ME_epg_lb"] = -0.845
cuts_dvcs_CD_Inb_4sigma["coplanarity_ub"] = 6.536
cuts_dvcs_CD_Inb_4sigma["coneAngle_ub"] = [0.830, -10.426, 66.708]
cuts_dvcs_CD_Inb_4sigma["coneAngle_lb"] = [-0.000646, 0.352, 8.316]
cuts_dvcs_CD_Inb_4sigma["MPt_ub"] = 0.119
cuts_dvcs_CD_Inb_4sigma["reconGam_ub"] = 0.944
cuts_dvcs_CD_Inb_4sigma["MM2_epg_ub"] = 0.0270
cuts_dvcs_CD_Inb_4sigma["MM2_epg_lb"] = -0.0303

cuts_dvcs_FD_Inb_3sigma = {}
cuts_dvcs_FD_Inb_3sigma["MM2_ep_ub"] = 0.176
cuts_dvcs_FD_Inb_3sigma["MM2_ep_lb"] = -0.214
cuts_dvcs_FD_Inb_3sigma["MM2_eg_ub"] = 1.850
cuts_dvcs_FD_Inb_3sigma["MM2_eg_lb"] = 0.0892
cuts_dvcs_FD_Inb_3sigma["ME_epg_ub"] = 0.976
cuts_dvcs_FD_Inb_3sigma["ME_epg_lb"] = -0.809
cuts_dvcs_FD_Inb_3sigma["coplanarity_ub"] = 8.3
cuts_dvcs_FD_Inb_3sigma["coneAngle_ub"] = [ 0.710, -7.582, 54.386]
cuts_dvcs_FD_Inb_3sigma["coneAngle_lb"] = [ 0.589, -6.361, 42.476]
cuts_dvcs_FD_Inb_3sigma["MPt_ub"] = 0.219
cuts_dvcs_FD_Inb_3sigma["reconGam_ub"] = 1.359
cuts_dvcs_FD_Inb_3sigma["MM2_epg_ub"] = 0.0148
cuts_dvcs_FD_Inb_3sigma["MM2_epg_lb"] = -0.0179

cuts_dvcs_FD_Inb_2sigma = {}
cuts_dvcs_FD_Inb_2sigma["MM2_ep_ub"] = 0.149
cuts_dvcs_FD_Inb_2sigma["MM2_ep_lb"] = -0.111
cuts_dvcs_FD_Inb_2sigma["MM2_eg_ub"] = 1.481
cuts_dvcs_FD_Inb_2sigma["MM2_eg_lb"] = 0.446
cuts_dvcs_FD_Inb_2sigma["ME_epg_ub"] = 0.726
cuts_dvcs_FD_Inb_2sigma["ME_epg_lb"] = -0.511
cuts_dvcs_FD_Inb_2sigma["coplanarity_ub"] = 5.337
cuts_dvcs_FD_Inb_2sigma["coneAngle_ub"] = [0.585, -5.765, 49.165]
cuts_dvcs_FD_Inb_2sigma["coneAngle_lb"] = [0.513, -5.374, 39.871]
cuts_dvcs_FD_Inb_2sigma["MPt_ub"] = 0.129
cuts_dvcs_FD_Inb_2sigma["reconGam_ub"] = 0.856
cuts_dvcs_FD_Inb_2sigma["MM2_epg_ub"] = 0.00500
cuts_dvcs_FD_Inb_2sigma["MM2_epg_lb"] = -0.00693

cuts_dvcs_FD_Inb_4sigma = {}
cuts_dvcs_FD_Inb_4sigma["MM2_ep_ub"] = 0.279
cuts_dvcs_FD_Inb_4sigma["MM2_ep_lb"] = -0.241
cuts_dvcs_FD_Inb_4sigma["MM2_eg_ub"] = 2.127
cuts_dvcs_FD_Inb_4sigma["MM2_eg_lb"] = -0.177
cuts_dvcs_FD_Inb_4sigma["ME_epg_ub"] = 1.056
cuts_dvcs_FD_Inb_4sigma["ME_epg_lb"] = -0.881
cuts_dvcs_FD_Inb_4sigma["coplanarity_ub"] = 11.358
cuts_dvcs_FD_Inb_4sigma["coneAngle_ub"] = [0.608, -6.120, 54.752]
cuts_dvcs_FD_Inb_4sigma["coneAngle_lb"] = [0.559, -6.011, 36.742]
cuts_dvcs_FD_Inb_4sigma["MPt_ub"] = 0.251
cuts_dvcs_FD_Inb_4sigma["reconGam_ub"] = 1.736
cuts_dvcs_FD_Inb_4sigma["MM2_epg_ub"] = 0.0270
cuts_dvcs_FD_Inb_4sigma["MM2_epg_lb"] = -0.0312


cuts_dvcs_CDFT_Outb_3sigma = {}
cuts_dvcs_CDFT_Outb_3sigma["MM2_ep_ub"] = 0.321
cuts_dvcs_CDFT_Outb_3sigma["MM2_ep_lb"] = -0.244
cuts_dvcs_CDFT_Outb_3sigma["MM2_eg_ub"] = 1.352
cuts_dvcs_CDFT_Outb_3sigma["MM2_eg_lb"] = 0.418
cuts_dvcs_CDFT_Outb_3sigma["ME_epg_ub"] = 0.299
cuts_dvcs_CDFT_Outb_3sigma["ME_epg_lb"] = -0.275
cuts_dvcs_CDFT_Outb_3sigma["coplanarity_ub"] = 4.742
cuts_dvcs_CDFT_Outb_3sigma["coneAngle_ub"] = cuts_dvcs_CDFT_Inb_3sigma["coneAngle_ub"] 
cuts_dvcs_CDFT_Outb_3sigma["coneAngle_lb"] = cuts_dvcs_CDFT_Inb_3sigma["coneAngle_lb"] 
cuts_dvcs_CDFT_Outb_3sigma["MPt_ub"] = 0.0627
cuts_dvcs_CDFT_Outb_3sigma["reconGam_ub"] = 0.578
cuts_dvcs_CDFT_Outb_3sigma["MM2_epg_ub"] = 0.00653
cuts_dvcs_CDFT_Outb_3sigma["MM2_epg_lb"] = -0.00884

cuts_dvcs_CDFT_Outb_2sigma = {}
cuts_dvcs_CDFT_Outb_2sigma["MM2_ep_ub"] = 0.227
cuts_dvcs_CDFT_Outb_2sigma["MM2_ep_lb"] = -0.150
cuts_dvcs_CDFT_Outb_2sigma["MM2_eg_ub"] = 1.196
cuts_dvcs_CDFT_Outb_2sigma["MM2_eg_lb"] = 0.572
cuts_dvcs_CDFT_Outb_2sigma["ME_epg_ub"] = 0.272
cuts_dvcs_CDFT_Outb_2sigma["ME_epg_lb"] = -0.237
cuts_dvcs_CDFT_Outb_2sigma["coplanarity_ub"] = 3.000
cuts_dvcs_CDFT_Outb_2sigma["coneAngle_ub"] = cuts_dvcs_CDFT_Inb_2sigma["coneAngle_ub"] 
cuts_dvcs_CDFT_Outb_2sigma["coneAngle_lb"] = cuts_dvcs_CDFT_Inb_2sigma["coneAngle_lb"] 
cuts_dvcs_CDFT_Outb_2sigma["MPt_ub"] = 0.0487
cuts_dvcs_CDFT_Outb_2sigma["reconGam_ub"] = 0.264
cuts_dvcs_CDFT_Outb_2sigma["MM2_epg_ub"] = 0.00237
cuts_dvcs_CDFT_Outb_2sigma["MM2_epg_lb"] = -0.00315

cuts_dvcs_CDFT_Outb_4sigma = {}
cuts_dvcs_CDFT_Outb_4sigma["MM2_ep_ub"] = 0.415
cuts_dvcs_CDFT_Outb_4sigma["MM2_ep_lb"] = -0.338
cuts_dvcs_CDFT_Outb_4sigma["MM2_eg_ub"] = 1.512
cuts_dvcs_CDFT_Outb_4sigma["MM2_eg_lb"] = 0.257
cuts_dvcs_CDFT_Outb_4sigma["ME_epg_ub"] = 0.516
cuts_dvcs_CDFT_Outb_4sigma["ME_epg_lb"] = -0.494
cuts_dvcs_CDFT_Outb_4sigma["coplanarity_ub"] = 6.209
cuts_dvcs_CDFT_Outb_4sigma["coneAngle_ub"] = cuts_dvcs_CDFT_Inb_4sigma["coneAngle_ub"] 
cuts_dvcs_CDFT_Outb_4sigma["coneAngle_lb"] = cuts_dvcs_CDFT_Inb_4sigma["coneAngle_lb"] 
cuts_dvcs_CDFT_Outb_4sigma["MPt_ub"] = 0.0968
cuts_dvcs_CDFT_Outb_4sigma["reconGam_ub"] = 0.778
cuts_dvcs_CDFT_Outb_4sigma["MM2_epg_ub"] = 0.0132
cuts_dvcs_CDFT_Outb_4sigma["MM2_epg_lb"] = -0.0162

cuts_dvcs_CD_Outb_3sigma = {}
cuts_dvcs_CD_Outb_3sigma["MM2_ep_ub"] = 0.226
cuts_dvcs_CD_Outb_3sigma["MM2_ep_lb"] = -0.196
cuts_dvcs_CD_Outb_3sigma["MM2_eg_ub"] = 1.902
cuts_dvcs_CD_Outb_3sigma["MM2_eg_lb"] = 0.0356
cuts_dvcs_CD_Outb_3sigma["ME_epg_ub"] = 0.755
cuts_dvcs_CD_Outb_3sigma["ME_epg_lb"] = -0.631
cuts_dvcs_CD_Outb_3sigma["coplanarity_ub"] = 5.181
cuts_dvcs_CD_Outb_3sigma["coneAngle_ub"] = cuts_dvcs_CD_Inb_3sigma["coneAngle_ub"] 
cuts_dvcs_CD_Outb_3sigma["coneAngle_lb"] = cuts_dvcs_CD_Inb_3sigma["coneAngle_lb"] 
cuts_dvcs_CD_Outb_3sigma["MPt_ub"] = 0.147
cuts_dvcs_CD_Outb_3sigma["reconGam_ub"] = 0.839
cuts_dvcs_CD_Outb_3sigma["MM2_epg_ub"] = 0.0138
cuts_dvcs_CD_Outb_3sigma["MM2_epg_lb"] = -0.0163

cuts_dvcs_CD_Outb_2sigma = {}
cuts_dvcs_CD_Outb_2sigma["MM2_ep_ub"] = 0.156
cuts_dvcs_CD_Outb_2sigma["MM2_ep_lb"] = -0.126
cuts_dvcs_CD_Outb_2sigma["MM2_eg_ub"] = 1.524
cuts_dvcs_CD_Outb_2sigma["MM2_eg_lb"] = 0.407
cuts_dvcs_CD_Outb_2sigma["ME_epg_ub"] = 0.263
cuts_dvcs_CD_Outb_2sigma["ME_epg_lb"] = -0.145
cuts_dvcs_CD_Outb_2sigma["coplanarity_ub"] = 3.564
cuts_dvcs_CD_Outb_2sigma["coneAngle_ub"] = cuts_dvcs_CD_Inb_2sigma["coneAngle_ub"] 
cuts_dvcs_CD_Outb_2sigma["coneAngle_lb"] = cuts_dvcs_CD_Inb_2sigma["coneAngle_lb"] 
cuts_dvcs_CD_Outb_2sigma["MPt_ub"] = 0.0686
cuts_dvcs_CD_Outb_2sigma["reconGam_ub"] = 0.587
cuts_dvcs_CD_Outb_2sigma["MM2_epg_ub"] = 0.00324
cuts_dvcs_CD_Outb_2sigma["MM2_epg_lb"] = -0.00456

cuts_dvcs_CD_Outb_4sigma = {}
cuts_dvcs_CD_Outb_4sigma["MM2_ep_ub"] = 0.297
cuts_dvcs_CD_Outb_4sigma["MM2_ep_lb"] = -0.267
cuts_dvcs_CD_Outb_4sigma["MM2_eg_ub"] = 2.220
cuts_dvcs_CD_Outb_4sigma["MM2_eg_lb"] = -0.277
cuts_dvcs_CD_Outb_4sigma["ME_epg_ub"] = 1.023
cuts_dvcs_CD_Outb_4sigma["ME_epg_lb"] = -0.898
cuts_dvcs_CD_Outb_4sigma["coplanarity_ub"] = 6.695
cuts_dvcs_CD_Outb_4sigma["coneAngle_ub"] = cuts_dvcs_CD_Inb_4sigma["coneAngle_ub"] 
cuts_dvcs_CD_Outb_4sigma["coneAngle_lb"] = cuts_dvcs_CD_Inb_4sigma["coneAngle_lb"] 
cuts_dvcs_CD_Outb_4sigma["MPt_ub"] = 0.131
cuts_dvcs_CD_Outb_4sigma["reconGam_ub"] = 1.148
cuts_dvcs_CD_Outb_4sigma["MM2_epg_ub"] = 0.0243
cuts_dvcs_CD_Outb_4sigma["MM2_epg_lb"] = -0.0273


cuts_dvcs_FD_Outb_3sigma = {}
cuts_dvcs_FD_Outb_3sigma["MM2_ep_ub"] = 0.285
cuts_dvcs_FD_Outb_3sigma["MM2_ep_lb"] = -0.276
cuts_dvcs_FD_Outb_3sigma["MM2_eg_ub"] = 1.479
cuts_dvcs_FD_Outb_3sigma["MM2_eg_lb"] = 0.322
cuts_dvcs_FD_Outb_3sigma["ME_epg_ub"] = 0.360
cuts_dvcs_FD_Outb_3sigma["ME_epg_lb"] = -0.387
cuts_dvcs_FD_Outb_3sigma["coplanarity_ub"] = 5.993
cuts_dvcs_FD_Outb_3sigma["coneAngle_ub"] = cuts_dvcs_FD_Inb_3sigma["coneAngle_ub"] 
cuts_dvcs_FD_Outb_3sigma["coneAngle_lb"] = cuts_dvcs_FD_Inb_3sigma["coneAngle_lb"] 
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
cuts_dvcs_FD_Outb_2sigma["coneAngle_ub"] = cuts_dvcs_FD_Inb_2sigma["coneAngle_ub"] 
cuts_dvcs_FD_Outb_2sigma["coneAngle_lb"] = cuts_dvcs_FD_Inb_2sigma["coneAngle_lb"] 
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
cuts_dvcs_FD_Outb_4sigma["coneAngle_ub"] = cuts_dvcs_FD_Inb_4sigma["coneAngle_ub"] 
cuts_dvcs_FD_Outb_4sigma["coneAngle_lb"] = cuts_dvcs_FD_Inb_4sigma["coneAngle_lb"] 
cuts_dvcs_FD_Outb_4sigma["MPt_ub"] = 0.0856
cuts_dvcs_FD_Outb_4sigma["reconGam_ub"] = 0.827
cuts_dvcs_FD_Outb_4sigma["MM2_epg_ub"] = 0.000941
cuts_dvcs_FD_Outb_4sigma["MM2_epg_lb"] = -0.00551
