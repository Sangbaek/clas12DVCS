import pandas as pd
import numpy as np
import argparse
M = 0.938272081 # target mass

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
phibin = [0, 12, 24, 36, 48, 60, 72, 96, 120, 144, 168, 192, 216, 240, 264, 288, 300, 312, 324, 336, 348, 360]
phibin_i = phibin[:-1]
phibin_f = phibin[1:]

badBins = ['004', '005', '006', '015', '016', '026', '045', '046', '104', '105', '106', '115', '116', '126', '146', '205', '206', '216', '240', '245', '246', '306', '320', '330', '406', '410', '420', '500', '505', '506', '510', '600', '601']
goodBins = ['000', '001', '002', '003', '010', '011', '012', '013', '014', '020', '021', '022', '023', '024', '025', '030', '031', '032', '033', '034', '035', '036', '040', '041', '042', '043', '044', '100', '101', '102', '103', '110', '111', '112', '113', '114', '120', '121', '122', '123', '124', '125', '130', '131', '132', '133', '134', '135', '136', '140', '141', '142', '143', '144', '145', '200', '201', '202', '203', '204', '210', '211', '212', '213', '214', '215', '220', '221', '222', '223', '224', '225', '226', '230', '231', '232', '233', '234', '235', '236', '241', '242', '243', '244', '300', '301', '302', '303', '304', '305', '310', '311', '312', '313', '314', '315', '316', '321', '322', '323', '324', '325', '326', '331', '332', '333', '334', '335', '336', '400', '401', '402', '403', '404', '405', '411', '412', '413', '414', '415', '416', '421', '422', '423', '424', '425', '426', '501', '502', '503', '504', '511', '512', '513', '514', '515', '516', '602', '603', '604', '605', '606']

df_global = pd.read_pickle("/volatile/clas12/sangbaek/clas12DVCS/df_global_Feb.pkl")

def TruebinVol(Q2bin, xBbin, tbin, phibin, Q2xBtphi, df, N1=10, N2=10, N3=10, N4=10):
    
    count = 0 
    
    Q2_i = Q2bin_i[Q2bin]
    Q2_f = Q2bin_f[Q2bin]
    xB_i = xBbin_i[Q2bin][xBbin]
    xB_f = xBbin_f[Q2bin][xBbin]
    t_i = tbin_i[tbin]
    t_f = tbin_f[tbin]
    phi_i = phibin_i[phibin]
    phi_f = phibin_f[phibin]
    
    local = df.loc[df.Q2xBtphi == Q2xBtphi]
    print(Q2xBtphi, len(local))
                
    if isinstance(xB_i, list):
        xB_i = min(xB_i)
    if isinstance(xB_f, list):
        xB_f = max(xB_f)
    for q2ind in range(N1):
        for xind in range(N2):
            for tind in range(N3):
                for phiind in range(N4):
                    qmin = Q2_i + (Q2_f - Q2_i)*q2ind/N1
                    qmax = Q2_i + (Q2_f - Q2_i)*(q2ind+1)/N1
                    xmin = xB_i + (xB_f - xB_i)*xind/N2
                    xmax = xB_i + (xB_f - xB_i)*(xind+1)/N2
                    tmin = t_i + (t_f - t_i)*tind/N3
                    tmax = t_i + (t_f - t_i)*(tind+1)/N3
                    pmin = phi_i + (phi_f - phi_i)*phiind/N4
                    pmax = phi_i + (phi_f - phi_i)*(phiind+1)/N4
                    local2 = local.loc[(local.Q2>=qmin) & (local.Q2<qmax)]
                    local2 = local2.loc[(local2.xB>=xmin) & (local2.xB<xmax)]
                    local2 = local2.loc[(local2.t1>=tmin) & (local2.t1<tmax)]
                    local2 = local2.loc[(local2.phi1>=pmin) & (local2.phi1<pmax)]
                    if(len(local2)):
                        count += 1
    return count/N1/N2/N3/N4*(Q2_f - Q2_i)*(xB_f - xB_i)*(t_f - t_i)*np.radians(phi_f - phi_i)


import argparse
parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-c","--config", help="config", default="1")
parser.add_argument("-p","--polarity", help="polarity", default="inbending")

if args.polarity == "inbending":
    print("inbending")
    parent_MC = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/dvcs/"
    parent_bhMC = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/bh/"


    #dvcs inb 50 nA
    print("reading dvcs inb 50 nA")
    df_3987_corr = pd.read_pickle(parent_MC + "3987.pkl")
    df_4124_corr = pd.read_pickle(parent_MC + "4124.pkl")
    df_4139_corr = pd.read_pickle(parent_MC + "4139.pkl")
    df_4181_corr = pd.read_pickle(parent_MC + "4181.pkl")
    df_4182_corr = pd.read_pickle(parent_MC + "4182.pkl")
    df_4397_corr = pd.read_pickle(parent_MC + "4397.pkl")

    #dvcs inb 55 nA
    print("reading dvcs inb 55 nA")
    df_4186_corr = pd.read_pickle(parent_MC + "4186.pkl")

    #dvcs inb 45 nA
    print("reading dvcs inb 45 nA")
    df_4188_corr = pd.read_pickle(parent_MC + "4188.pkl")

    #dvcs inb 0 nA
    print("reading dvcs inb 0 nA")
    df_4192_corr = pd.read_pickle(parent_MC + "4192.pkl")

    #bh inb 50 nA
    print("reading bh inb 50 nA")
    df_4238_corr = pd.read_pickle(parent_bhMC + "4238.pkl")

    dvcsSimInb50nA = pd.concat([df_3987_corr, df_4124_corr, df_4139_corr, df_4181_corr, df_4182_corr, df_4397_corr])
    dvcsSimInb55nA = df_4186_corr
    dvcsSimInb45nA = df_4188_corr
    dvcsSimInb0nA = df_4192_corr
    bhSimInb50nA = df_4238_corr

    dvcsBHSim = pd.concat([dvcsSimInb50nA, dvcsSimInb55nA, dvcsSimInb45nA, dvcsSimInb0nA, bhSimInb50nA])

else:
    print("outbending")
    parent_MC = "/volatile/clas12/sangbaek/nov2021/convPkl_full/outb/dvcs/"
    parent_bhMC = "/volatile/clas12/sangbaek/nov2021/convPkl_full/outb/bh/"

    df_4240_corr = pd.read_pickle(parent_MC + "4240.pkl")
    df_4250_corr = pd.read_pickle(parent_MC + "4250.pkl")
    df_4251_corr = pd.read_pickle(parent_MC + "4251.pkl")
    df_4252_corr = pd.read_pickle(parent_MC + "4252.pkl")
    df_4255_corr = pd.read_pickle(parent_MC + "4255.pkl")
    df_4263_corr = pd.read_pickle(parent_MC + "4263.pkl")
    df_4262_corr = pd.read_pickle(parent_MC + "4262.pkl")
    df_4266_corr = pd.read_pickle(parent_MC + "4266.pkl")
    df_4249_corr = pd.read_pickle(parent_bhMC + "4249.pkl")

    dvcsSimOutb50nA = pd.concat([df_4240_corr, df_4250_corr, df_4251_corr, df_4252_corr, df_4255_corr])
    dvcsSimOutb40nA = df_4263_corr
    dvcsSimOutb0nA = df_4262_corr
    dvcsSimOutb40nAT = df_4266_corr
    bhSimOutb50nA = df_4249_corr

    dvcsBHSim = pd.concat([dvcsSimOutb50nA, dvcsSimOutb40nA, dvcsSimOutb0nA, dvcsSimOutb40nAT, bhSimOutb50nA])

if args.config == "0":
    print("counting inb bin volumes for All configs")
if args.config == "1":
    print("counting inb bin volumes for FD")
    dvcsBHSim = dvcsBHSim.loc[dvcsBHSim.config == 1]
if args.config == "2":
    print("counting inb bin volumes for CD")
    dvcsBHSim = dvcsBHSim.loc[dvcsBHSim.config == 2]
if args.config == "3":
    print("counting inb bin volumes for CDFT")
    dvcsBHSim = dvcsBHSim.loc[dvcsBHSim.config == 3]

TrueVol = []
for i in range(len(df_global)):
    TrueVol.append(TruebinVol(df_global.Q2[i], df_global.xB[i], df_global.t[i], df_global.phi[i], df_global.Q2xBtphi[i], dvcsBHSim, 6, 6, 6, 6))
df_global.loc[:, "binVolInb"] = TrueVolInb
df_global.to_pickle("/volatile/clas12/sangbaek/clas12DVCS/results/truebinVol_{}{}.pkl".format(args.polarity, args.config))