#!/bin/env python3
"""
Modules to define useful constants
"""

import numpy as np
from copy import copy

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

#binning variables
x1 = 1/2/M/8.604
x2 = 1/(5-M**2)
x3 = (10.604/8.604-1)/M*10.604* (1-np.cos(np.radians(35)))
x4 = (1-(4-M**2)/2/10.604/M)/(1+(4-M**2)/2/10.604**2/(1-np.cos(np.radians(35))))

y1 = 1
y2 = 1.456
y3 = 2.510
y4 = 4.326
y5 = 7.671

c0 = y2/2/M/8.604
d0 = 1/(1+(4-M*M)/y2)
c1 = np.sqrt(y2*y3)/2/M/8.604
d1 = 1/(1+(4-M*M)/np.sqrt(y2*y3))
c2 = y3/2/M/8.604
d2 =  1/(1+(4-M*M)/y3)
c3 = np.sqrt(y3*y4)/2/M/8.604
d3 = 1/(1+(4-M*M)/np.sqrt(y3*y4))
c4 = y4/2/M/8.604
d4 = 1/(1+(4-M*M)/y4)
c5 = np.sqrt(y4*y5)/2/M/8.604
d5 = 1/(1+(4-M*M)/np.sqrt(y4*y5))
c6 = y5/2/M/8.604
d6 = 1/(1+(4-M*M)/y5)

Q2bins = [y1, y2, np.sqrt(y2*y3), y3, np.sqrt(y3*y4), y4, np.sqrt(y4*y5), y5]
Q2bin_i = Q2bins[:-1]
Q2bin_f = Q2bins[1:]
xBbin_i = {0:[[x1, c0], c0, c1, c2, c3], 1: [[c0, c1], c1, c2, c3, c4], 2: [[c1, c2], c2, c3, c4, d1], 3: [[c2, c3], c3, c4, d1], 4: [[c3, c4], c4, d1], 5: [[c4, c5], d1], 6: [[c5, c6]]}
xBbin_f = {0: [c0, c1, c2, c3, [x2, d0]], 1: [c1, c2, c3, c4, [d0, d1]], 2: [c2, c3, c4, d1, [d1, d2]], 3: [c3, c4, d1, [d2, d3]], 4: [c4, d1, [d3, d4]], 5: [d1, [d4, d5]], 6: [[d5, d6]]}
tbins = [0.088, 0.177, 0.321, 0.523, 0.813, 1.187, 1.46, 1.72]
tbin_i = tbins[:-1]
tbin_f = tbins[1:]
phibins = np.linspace(0, 360, 25)
phibin_i = phibins[:-1]
phibin_f = phibins[1:]

badBins = ['004', '005', '006', '015', '016', '026', '045', '046', '104', '105', '106', '115', '116', '126', '146', '205', '206', '216', '240', '245', '246', '306', '320', '330', '406', '410', '420', '500', '505', '506', '510', '600', '601']
goodBins = ['000', '001', '002', '003', '010', '011', '012', '013', '014', '020', '021', '022', '023', '024', '025', '030', '031', '032', '033', '034', '035', '036', '040', '041', '042', '043', '044', '100', '101', '102', '103', '110', '111', '112', '113', '114', '120', '121', '122', '123', '124', '125', '130', '131', '132', '133', '134', '135', '136', '140', '141', '142', '143', '144', '145', '200', '201', '202', '203', '204', '210', '211', '212', '213', '214', '215', '220', '221', '222', '223', '224', '225', '226', '230', '231', '232', '233', '234', '235', '236', '241', '242', '243', '244', '300', '301', '302', '303', '304', '305', '310', '311', '312', '313', '314', '315', '316', '321', '322', '323', '324', '325', '326', '331', '332', '333', '334', '335', '336', '400', '401', '402', '403', '404', '405', '411', '412', '413', '414', '415', '416', '421', '422', '423', '424', '425', '426', '501', '502', '503', '504', '511', '512', '513', '514', '515', '516', '602', '603', '604', '605', '606']


inactivebins = np.array([   0.,   22.,   23.,   24.,   25.,   46.,   47.,   48.,   71.,
        122.,  141.,  145.,  146.,  147.,  152.,  156.,  157.,  163.,
        165.,  166.,  167.,  313.,  334.,  481.,  680.,  681.,  682.,
        683.,  684.,  685.,  686.,  687.,  707.,  819.,  827.,  828.,
        888.,  911.,  912.,  974.,  986.,  995., 1000., 1004., 1005.,
       1352., 1353., 1354., 1355., 1356., 1357., 1358., 1359., 1518.,
       1519., 1520., 1521., 1522., 1523., 1524., 1525., 1526., 1527.,
       1528., 1529., 1642., 1645., 1666., 1667., 1668., 1669., 1675.,
       1752., 1806., 1826., 1833., 1835., 1839., 1999., 2025., 2026.,
       2027., 2028., 2029., 2030., 2190., 2191., 2193., 2194., 2195.,
       2196., 2197., 2198., 2199., 2200., 2201., 2202., 2352., 2353.,
       2354., 2355., 2356., 2357., 2358., 2359., 2360., 2361., 2362.,
       2363., 2364., 2365., 2366., 2367., 2368., 2369., 2370., 2371.,
       2372., 2373., 2374., 2375., 2387., 2388., 2459., 2476., 2482.,
       2483., 2484., 2485., 2505., 2506., 2507., 2508., 2509., 2510.,
       2511., 2616., 2648., 2654., 2671., 2672., 2677., 2678., 2679.,
       2700., 2861., 2862., 2863., 2864., 2865., 2866., 2867., 2868.,
       2869., 2870., 2871., 2872., 2873., 2874., 2875., 3024., 3025.,
       3026., 3027., 3028., 3029., 3030., 3031., 3032., 3033., 3034.,
       3035., 3036., 3037., 3038., 3039., 3040., 3041., 3042., 3043.,
       3044., 3045., 3046., 3047., 3179., 3180., 3202., 3203., 3204.,
       3205., 3207., 3208., 3343., 3363., 3364., 3365., 3366., 3367.,
       3368., 3369., 3370., 3371., 3372., 3373., 3374., 3375., 3376.,
       3377., 3378., 3379., 3380., 3528., 3529., 3530., 3531., 3532.,
       3533., 3534., 3535., 3536., 3537., 3538., 3539., 3540., 3541.,
       3542., 3543., 3544., 3545., 3546., 3547., 3548., 3549., 3550.,
       3551., 3696., 3697., 3698., 3699., 3700., 3701., 3702., 3703.,
       3704., 3705., 3706., 3707., 3708., 3709., 3710., 3711., 3712.,
       3713., 3714., 3715., 3716., 3717., 3718., 3719., 3864., 3865.,
       3866., 3867., 3868., 3869., 3870., 3871., 3872., 3873., 3874.,
       3875., 3876., 3877., 3878., 3879., 3880., 3881., 3882., 3883.,
       3884., 3885., 3886., 3887., 3899., 4032., 4033., 4034., 4035.,
       4036., 4037., 4038., 4039., 4040., 4041., 4042., 4043., 4044.,
       4045., 4046., 4047., 4048., 4049., 4050., 4051., 4052., 4053.,
       4054., 4055.])
activebins = np.linspace(0, 4199, 4200, dtype=int)[~np.isin(np.linspace(0, 4199, 4200), inactivebins)]

inactiveQ2xBtbins = np.array([98, 126, 147, 154, 161, 168])
activeQ2xBtbins = np.linspace(0, 174, 175, dtype=int)[~np.isin(np.linspace(0, 174, 175, dtype=int), inactiveQ2xBtbins)]

newxBbins = [x1, c0, c1, c2, c3, c4, d2]
newQ2bins = [y1, y2, np.sqrt(y2*y3), y3, np.sqrt(y3*y4), y4, np.sqrt(y4*y5)]
newtbins = [0.11, 0.15, 0.25, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 1.79]

collection_cont_xBbins  = [np.linspace(0.05, 0.85, 6), newxBbins, np.linspace(0.05, 0.85, 2)]
collection_cont_Q2bins  = [np.array([1, 1.5, 2, 2.5, 3.5, 4.5, 6, 7.5, 10.3]), newQ2bins, np.linspace(1, 10.3, 2)]
collection_cont_tbins   = [np.array([0.09, 0.2, 0.4, 0.8, 1.8]), newtbins, np.linspace(0.09, 1.8, 2)]
collection_cont_phibins = [np.linspace(0, 360, 6), phibins, np.linspace(0, 360, 2)]

collection_xBbins = [np.linspace(0.05, 0.85, 6), newxBbins, newxBbins]
collection_Q2bins = [np.array([1, 1.5, 2, 2.5, 3.5, 4.5, 6, 7.5, 10.3]), newQ2bins, newQ2bins]
collection_tbins = [np.array([0.09, 0.2, 0.4, 0.8, 1.8]), tbins, newtbins]
collection_phibins = [phibins, phibins, phibins]

# simulation run numbers
runs_inb_vgg50nA = [3987, 4124, 4139, 4181, 4182, 4397, 4528, 4529, 4535, 4539]
runs_inb_vgg55nA = [4186, 4545]
runs_inb_vgg45nA = [4188, 4547, 4893]
runs_inb_vgg0nA = [4192, 4561]
runs_inb_bh50nA = [4238, 4542]
# runs_inb_bh45nA = [4740, 4742, 4745, 4751, 4760, 4892]
runs_inb_bh45nA = [4892]
runs_inb_bkg50nA = [4076, 4202, 4209]
runs_inb_bkg55nA = [4212]
runs_inb_bkg45nA = [4217]
runs_inb_bkg0nA = [4231]

runs_outb_vgg50nA = [4240, 4250, 4251, 4252, 4255, 4398, 4532, 4534, 4540, 4541, 4717]
runs_outb_vgg40nA = [4263, 4546]
runs_outb_vgg0nA = [4262, 4554]
runs_outb_vgg40nAT = [4266, 4562]
runs_outb_bh50nA = [4249, 4544, 4780, 4808, 4812, 4818, 4902, 4903, 4905]
runs_outb_bkg50nA = [4243, 4271, 4290]
runs_outb_bkg40nA = [4293]
runs_outb_bkg0nA = [4304]
runs_outb_bkg40nAT = [4306]

def numberingDF(total, Q2bin_i=Q2bin_i, Q2bin_f=Q2bin_f, xBbin_i=xBbin_i, xBbin_f=xBbin_f, tbin_i=tbin_i, tbin_f=tbin_f, goodBins=goodBins, badBins=badBins):
    df_allBins = {}
    Q2xBtphi = 0

    for Q2bin in range(0, len(Q2bin_i)):#Q2 bin
        for xBbin in range(0, len(xBbin_i[Q2bin])):
            for tbin in range(0, len(tbin_i)):
                local = total
                Q2_i = Q2bin_i[Q2bin]
                Q2_f = Q2bin_f[Q2bin]
                xB_i = xBbin_i[Q2bin][xBbin]
                xB_f = xBbin_f[Q2bin][xBbin]
                t_i = tbin_i[tbin]
                t_f = tbin_f[tbin]
                #cut by Q2
                local = local.loc[(local.Q2>=Q2_i) & (local.Q2<Q2_f)]
                #cut by xB
                #xB lower bound
                if xBbin == 0:
                    local = local.loc[local.Q2<=2*M*(10.604-2)*local.xB, :]
                else:
                    local = local.loc[local.xB>=xB_i] 
                #xB upper bound
                if (xBbin == len(xBbin_i[Q2bin])-1):
                    local = local.loc[local.Q2>=(4 - M*M)*local.xB/(1 - local.xB)]
                else:
                    local = local.loc[local.xB<xB_f]
                #cut by t
                local = local.loc[(local.t1>=t_i) & (local.t1<t_f)]
                Q2xBtbin = "{}{}{}".format(Q2bin,xBbin,tbin)

                # if Q2xBtbin in badBins:
                #     continue
                if (len(local)==0):
                    print(Q2xBtbin, Q2xBtphi, xB_i)#continue
                    continue
                    
                for phi_ind in range(0, len(phibin_i)):
                    local.loc[:, "xBbin"] = xBbin
                    local.loc[:, "Q2bin"] = Q2bin
                    local.loc[:, "tbin"] = tbin
                    local.loc[:, "phibin"] = phi_ind
                    local.loc[:, "Q2xBtbin"] = Q2xBtbin
                    local.loc[:, "Q2xBtphi"] = Q2xBtphi
                    df_allBins[Q2xBtphi] = local.loc[(local.phi1>=phibin_i[phi_ind])&(local.phi1<phibin_f[phi_ind])]
                    Q2xBtphi += 1

    total = pd.concat(df_allBins.values()).sort_values( by = 'event')
    return total

CD_Ptheta_ub = 80
CD_Ptheta_lb = 40
FD_Ptheta_inb_ub = 40
FD_Ptheta_outb_ub = 40
FD_Ptheta_lb = 5

Ge2Threshold_default = 0
Ge2Threshold_tight = 0.6
Ge2Threshold_mid = 0.4
Ge2Threshold_loose = 0.3

#the default compatible with F. X. Girod's dvcs wagon
cuts_dvcs_default_CDFT = {}
cuts_dvcs_default_CDFT["MM2_ep_ub"] = 0.6
cuts_dvcs_default_CDFT["MM2_ep_lb"] = -0.6
cuts_dvcs_default_CDFT["MM2_eg_ub"] = 2
cuts_dvcs_default_CDFT["MM2_eg_lb"] = 0.25
cuts_dvcs_default_CDFT["ME_epg_ub"] = 2
cuts_dvcs_default_CDFT["ME_epg_lb"] = -1
cuts_dvcs_default_CDFT["coplanarity_ub"] = 360
cuts_dvcs_default_CDFT["coneAngle_ub"] = [-0.00221, 0.863, 10.287]
cuts_dvcs_default_CDFT["coneAngle_lb"] = [0.0267, -0.0625, 7.730]
cuts_dvcs_default_CDFT["coneAngleCR_ub"] = [-0.000382, 0.777, 0.867]
cuts_dvcs_default_CDFT["coneAngleCR_lb"] = [0.0510, -0.0470, -0.492]
cuts_dvcs_default_CDFT["MPt_ub"] = 0.75
cuts_dvcs_default_CDFT["reconGam_ub"] = 7.5
cuts_dvcs_default_CDFT["MM2_epg_ub"] = 0.1
cuts_dvcs_default_CDFT["MM2_epg_lb"] = -0.1

#the default compatible with F. X. Girod's dvcs wagon
cuts_dvcs_default_CD = {}
cuts_dvcs_default_CD["MM2_ep_ub"] = 0.6
cuts_dvcs_default_CD["MM2_ep_lb"] = -0.6
cuts_dvcs_default_CD["MM2_eg_ub"] = 2
cuts_dvcs_default_CD["MM2_eg_lb"] = 0.25
cuts_dvcs_default_CD["ME_epg_ub"] = 2
cuts_dvcs_default_CD["ME_epg_lb"] = -1
cuts_dvcs_default_CD["coplanarity_ub"] = 360
cuts_dvcs_default_CD["coneAngle_ub"] = [0.0470, -1.677, 46.014]
cuts_dvcs_default_CD["coneAngle_lb"] = [0.0164, 0.408, 4.901]
cuts_dvcs_default_CD["MPt_ub"] = 0.75
cuts_dvcs_default_CD["reconGam_ub"] = 7.5
cuts_dvcs_default_CD["MM2_epg_ub"] = 0.1
cuts_dvcs_default_CD["MM2_epg_lb"] = -0.1

#the default compatible with F. X. Girod's dvcs wagon
cuts_dvcs_default_FD = {}
cuts_dvcs_default_FD["MM2_ep_ub"] = 0.6
cuts_dvcs_default_FD["MM2_ep_lb"] = -0.6
cuts_dvcs_default_FD["MM2_eg_ub"] = 2
cuts_dvcs_default_FD["MM2_eg_lb"] = 0.25
cuts_dvcs_default_FD["ME_epg_ub"] = 2
cuts_dvcs_default_FD["ME_epg_lb"] = -1
cuts_dvcs_default_FD["coplanarity_ub"] = 360
cuts_dvcs_default_FD["coneAngle_ub"] = [0.0280, -1.001, 49.895]
cuts_dvcs_default_FD["coneAngle_lb"] = [0.0214, -0.379, 21.998]
cuts_dvcs_default_FD["MPt_ub"] = 0.75
cuts_dvcs_default_FD["reconGam_ub"] = 7.5
cuts_dvcs_default_FD["MM2_epg_ub"] = 0.1
cuts_dvcs_default_FD["MM2_epg_lb"] = -0.1

cuts_dvcs_CDFT_Inb_3sigma = {}
cuts_dvcs_CDFT_Inb_3sigma["MM2_ep_ub"] = 0.391
cuts_dvcs_CDFT_Inb_3sigma["MM2_ep_lb"] = -0.365
cuts_dvcs_CDFT_Inb_3sigma["MM2_eg_ub"] = 1.479
cuts_dvcs_CDFT_Inb_3sigma["MM2_eg_lb"] = 0.322
cuts_dvcs_CDFT_Inb_3sigma["ME_epg_ub"] = 0.360
cuts_dvcs_CDFT_Inb_3sigma["ME_epg_lb"] = -0.387
cuts_dvcs_CDFT_Inb_3sigma["coplanarity_ub"] = 5.633
cuts_dvcs_CDFT_Inb_3sigma["coneAngle_ub"] = [-0.00221, 0.863, 10.287]
cuts_dvcs_CDFT_Inb_3sigma["coneAngle_lb"] = [0.0267, -0.0625, 7.730]
cuts_dvcs_CDFT_Inb_3sigma["coneAngleCR_ub"] = [-0.000382, 0.777, 0.867]
cuts_dvcs_CDFT_Inb_3sigma["coneAngleCR_lb"] = [0.0510, -0.0470, -0.492]
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
cuts_dvcs_CDFT_Inb_2sigma["coneAngle_ub"] = cuts_dvcs_CDFT_Inb_3sigma["coneAngle_ub"]
cuts_dvcs_CDFT_Inb_2sigma["coneAngle_lb"] = cuts_dvcs_CDFT_Inb_3sigma["coneAngle_lb"]
cuts_dvcs_CDFT_Inb_2sigma["coneAngleCR_ub"] = cuts_dvcs_CDFT_Inb_3sigma["coneAngleCR_ub"]
cuts_dvcs_CDFT_Inb_2sigma["coneAngleCR_lb"] = cuts_dvcs_CDFT_Inb_3sigma["coneAngleCR_lb"]
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
cuts_dvcs_CDFT_Inb_4sigma["coneAngle_ub"] = cuts_dvcs_CDFT_Inb_3sigma["coneAngle_ub"]
cuts_dvcs_CDFT_Inb_4sigma["coneAngle_lb"] = cuts_dvcs_CDFT_Inb_3sigma["coneAngle_lb"]
cuts_dvcs_CDFT_Inb_4sigma["coneAngleCR_ub"] = cuts_dvcs_CDFT_Inb_3sigma["coneAngleCR_ub"]
cuts_dvcs_CDFT_Inb_4sigma["coneAngleCR_lb"] = cuts_dvcs_CDFT_Inb_3sigma["coneAngleCR_lb"]
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
cuts_dvcs_CD_Inb_3sigma["coneAngle_ub"] = [0.0470, -1.677, 46.014]
cuts_dvcs_CD_Inb_3sigma["coneAngle_lb"] = [0.0164, 0.408, 4.901]
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
cuts_dvcs_CD_Inb_2sigma["coneAngle_ub"] = cuts_dvcs_CD_Inb_3sigma["coneAngle_ub"]
cuts_dvcs_CD_Inb_2sigma["coneAngle_lb"] = cuts_dvcs_CD_Inb_3sigma["coneAngle_lb"]
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
cuts_dvcs_CD_Inb_4sigma["coneAngle_ub"] = cuts_dvcs_CD_Inb_3sigma["coneAngle_ub"]
cuts_dvcs_CD_Inb_4sigma["coneAngle_lb"] = cuts_dvcs_CD_Inb_3sigma["coneAngle_lb"]
cuts_dvcs_CD_Inb_4sigma["MPt_ub"] = 0.119
cuts_dvcs_CD_Inb_4sigma["reconGam_ub"] = 0.944
cuts_dvcs_CD_Inb_4sigma["MM2_epg_ub"] = 0.0270
cuts_dvcs_CD_Inb_4sigma["MM2_epg_lb"] = -0.0303

cuts_dvcs_FD_Inb_3sigma = {}
cuts_dvcs_FD_Inb_3sigma["MM2_ep_ub"] = 0.190
cuts_dvcs_FD_Inb_3sigma["MM2_ep_lb"] = -0.144
cuts_dvcs_FD_Inb_3sigma["MM2_eg_ub"] = 1.989
cuts_dvcs_FD_Inb_3sigma["MM2_eg_lb"] = 0.0505
cuts_dvcs_FD_Inb_3sigma["ME_epg_ub"] = 0.976
cuts_dvcs_FD_Inb_3sigma["ME_epg_lb"] = -0.740
cuts_dvcs_FD_Inb_3sigma["coplanarity_ub"] = 8.36
cuts_dvcs_FD_Inb_3sigma["coneAngle_ub"] = [0.0280, -1.001, 49.895]
cuts_dvcs_FD_Inb_3sigma["coneAngle_lb"] = [0.0214, -0.379, 21.998]
cuts_dvcs_FD_Inb_3sigma["MPt_ub"] = 0.291
cuts_dvcs_FD_Inb_3sigma["reconGam_ub"] = 1.487
cuts_dvcs_FD_Inb_3sigma["MM2_epg_ub"] = 0.0142
cuts_dvcs_FD_Inb_3sigma["MM2_epg_lb"] = -0.0177

cuts_dvcs_FD_Inb_2sigma = {}
cuts_dvcs_FD_Inb_2sigma["MM2_ep_ub"] = 0.149
cuts_dvcs_FD_Inb_2sigma["MM2_ep_lb"] = -0.111
cuts_dvcs_FD_Inb_2sigma["MM2_eg_ub"] = 1.481
cuts_dvcs_FD_Inb_2sigma["MM2_eg_lb"] = 0.446
cuts_dvcs_FD_Inb_2sigma["ME_epg_ub"] = 0.726
cuts_dvcs_FD_Inb_2sigma["ME_epg_lb"] = -0.511
cuts_dvcs_FD_Inb_2sigma["coplanarity_ub"] = 5.337
cuts_dvcs_FD_Inb_2sigma["coneAngle_ub"] = cuts_dvcs_FD_Inb_3sigma["coneAngle_ub"]
cuts_dvcs_FD_Inb_2sigma["coneAngle_lb"] = cuts_dvcs_FD_Inb_3sigma["coneAngle_lb"]
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
cuts_dvcs_FD_Inb_4sigma["coneAngle_ub"] = cuts_dvcs_FD_Inb_3sigma["coneAngle_ub"]
cuts_dvcs_FD_Inb_4sigma["coneAngle_lb"] = cuts_dvcs_FD_Inb_3sigma["coneAngle_lb"]
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
cuts_dvcs_CDFT_Outb_3sigma["coneAngleCR_ub"] = cuts_dvcs_CDFT_Inb_3sigma["coneAngleCR_ub"] 
cuts_dvcs_CDFT_Outb_3sigma["coneAngleCR_lb"] = cuts_dvcs_CDFT_Inb_3sigma["coneAngleCR_lb"] 
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
cuts_dvcs_CDFT_Outb_2sigma["coneAngleCR_ub"] = cuts_dvcs_CDFT_Inb_2sigma["coneAngleCR_ub"] 
cuts_dvcs_CDFT_Outb_2sigma["coneAngleCR_lb"] = cuts_dvcs_CDFT_Inb_2sigma["coneAngleCR_lb"] 
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
cuts_dvcs_CDFT_Outb_4sigma["coneAngleCR_ub"] = cuts_dvcs_CDFT_Inb_4sigma["coneAngleCR_ub"] 
cuts_dvcs_CDFT_Outb_4sigma["coneAngleCR_lb"] = cuts_dvcs_CDFT_Inb_4sigma["coneAngleCR_lb"] 
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
cuts_dvcs_FD_Outb_3sigma["MM2_ep_ub"] = 0.225
cuts_dvcs_FD_Outb_3sigma["MM2_ep_lb"] = -0.174
cuts_dvcs_FD_Outb_3sigma["MM2_eg_ub"] = 1.940
cuts_dvcs_FD_Outb_3sigma["MM2_eg_lb"] = 0.160
cuts_dvcs_FD_Outb_3sigma["ME_epg_ub"] = 0.772
cuts_dvcs_FD_Outb_3sigma["ME_epg_lb"] = -0.519
cuts_dvcs_FD_Outb_3sigma["coplanarity_ub"] = 6.342
cuts_dvcs_FD_Outb_3sigma["coneAngle_ub"] = cuts_dvcs_FD_Inb_3sigma["coneAngle_ub"] 
cuts_dvcs_FD_Outb_3sigma["coneAngle_lb"] = cuts_dvcs_FD_Inb_3sigma["coneAngle_lb"] 
cuts_dvcs_FD_Outb_3sigma["MPt_ub"] = 0.340
cuts_dvcs_FD_Outb_3sigma["reconGam_ub"] = 1.645
cuts_dvcs_FD_Outb_3sigma["MM2_epg_ub"] = 0.0162
cuts_dvcs_FD_Outb_3sigma["MM2_epg_lb"] = -0.0202

cuts_dvcs_FD_Outb_2sigma = {}
cuts_dvcs_FD_Outb_2sigma["MM2_ep_ub"] = 0.169
cuts_dvcs_FD_Outb_2sigma["MM2_ep_lb"] = -0.116
cuts_dvcs_FD_Outb_2sigma["MM2_eg_ub"] = 1.678
cuts_dvcs_FD_Outb_2sigma["MM2_eg_lb"] = 0.349
cuts_dvcs_FD_Outb_2sigma["ME_epg_ub"] = 1.088
cuts_dvcs_FD_Outb_2sigma["ME_epg_lb"] = -0.828
cuts_dvcs_FD_Outb_2sigma["coplanarity_ub"] = 3.695
cuts_dvcs_FD_Outb_2sigma["coneAngle_ub"] = cuts_dvcs_FD_Inb_2sigma["coneAngle_ub"] 
cuts_dvcs_FD_Outb_2sigma["coneAngle_lb"] = cuts_dvcs_FD_Inb_2sigma["coneAngle_lb"] 
cuts_dvcs_FD_Outb_2sigma["MPt_ub"] = 0.110
cuts_dvcs_FD_Outb_2sigma["reconGam_ub"] = 1.125
cuts_dvcs_FD_Outb_2sigma["MM2_epg_ub"] = 0.00553
cuts_dvcs_FD_Outb_2sigma["MM2_epg_lb"] = -0.00745

cuts_dvcs_FD_Outb_4sigma = {}
cuts_dvcs_FD_Outb_4sigma["MM2_ep_ub"] = 0.312
cuts_dvcs_FD_Outb_4sigma["MM2_ep_lb"] = -0.259
cuts_dvcs_FD_Outb_4sigma["MM2_eg_ub"] = 2.450
cuts_dvcs_FD_Outb_4sigma["MM2_eg_lb"] = -0.397
cuts_dvcs_FD_Outb_4sigma["ME_epg_ub"] = 1.192
cuts_dvcs_FD_Outb_4sigma["ME_epg_lb"] = -0.960
cuts_dvcs_FD_Outb_4sigma["coplanarity_ub"] = 7.762
cuts_dvcs_FD_Outb_4sigma["coneAngle_ub"] = cuts_dvcs_FD_Inb_4sigma["coneAngle_ub"] 
cuts_dvcs_FD_Outb_4sigma["coneAngle_lb"] = cuts_dvcs_FD_Inb_4sigma["coneAngle_lb"] 
cuts_dvcs_FD_Outb_4sigma["MPt_ub"] = 0.341
cuts_dvcs_FD_Outb_4sigma["reconGam_ub"] = 2.725
cuts_dvcs_FD_Outb_4sigma["MM2_epg_ub"] = 0.0374
cuts_dvcs_FD_Outb_4sigma["MM2_epg_lb"] = -0.0444
