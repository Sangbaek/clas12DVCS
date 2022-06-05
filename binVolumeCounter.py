import pandas as pd
import numpy as np
import argparse
from utils.const import *

df_global = pd.read_pickle("/volatile/clas12/sangbaek/clas12DVCS/df_global_May.pkl")

def TruebinVol(df_global, df1, df2, df3, df4, N1=10, N2=10, N3=10, N4=10):
    

    df = pd.concat([df1, df2, df3, df4])

    count_by_Q2xBt, _  = np.histogram(df.Q2xBtbin, bins = np.linspace(-0.5, 174.5, 176))
    Q2s_by_Q2xBt, _  = np.histogram(df.Q2xBtbin, bins = np.linspace(-0.5, 174.5, 176), weights = df.Q2)
    xBs_by_Q2xBt, _  = np.histogram(df.Q2xBtbin, bins = np.linspace(-0.5, 174.5, 176), weights = df.xB)
    ts_by_Q2xBt, _  = np.histogram(df.Q2xBtbin, bins = np.linspace(-0.5, 174.5, 176), weights = df.t)
    
    xBavg1 = np.divide(xBs_by_Q2xBt, count_by_Q2xBt, where  = count_by_Q2xBt>0, out = np.zeros(count_by_Q2xBt.shape))
    Q2avg1 = np.divide(Q2s_by_Q2xBt, count_by_Q2xBt, where  = count_by_Q2xBt>0, out = np.zeros(count_by_Q2xBt.shape))
    tavg1  = np.divide(ts_by_Q2xBt, count_by_Q2xBt, where  = count_by_Q2xBt>0, out = np.zeros(count_by_Q2xBt.shape))

    xBavg1 = xBavg1[activeQ2xBtbins]
    Q2avg1 = Q2avg1[activeQ2xBtbins]
    tavg1 = tavg1[activeQ2xBtbins]

    dict_xB1 = {activeQ2xBtbins[i]:xBavg1[i] for i in range(len(xBavg1))}
    dict_Q21 = {activeQ2xBtbins[i]:Q2avg1[i] for i in range(len(Q2avg1))}
    dict_t1 = {activeQ2xBtbins[i]:tavg1[i] for i in range(len(tavg1))}

    df.loc[:, "xBavg1"] = df_global.Q2xBtbin.map(dict_xB1)
    df.loc[:, "Q2avg1"] = df_global.Q2xBtbin.map(dict_Q21)
    df.loc[:, "tavg1"] = df_global.Q2xBtbin.map(dict_t1)


    count_by_Q2xBtphi, _  = np.histogram(df.Q2xBtphibin, bins = np.linspace(-0.5, 4199.5, 4201))
    Q2s_by_Q2xBtphi, _  = np.histogram(df.Q2xBtphibin, bins = np.linspace(-0.5, 4199.5, 4201), weights = df.Q2)
    xBs_by_Q2xBtphi, _  = np.histogram(df.Q2xBtphibin, bins = np.linspace(-0.5, 4199.5, 4201), weights = df.xB)
    ts_by_Q2xBtphi, _  = np.histogram(df.Q2xBtphibin, bins = np.linspace(-0.5, 4199.5, 4201), weights = df.t)
    phis_by_Q2xBtphi, _  = np.histogram(df.Q2xBtphibin, bins = np.linspace(-0.5, 4199.5, 4201), weights = df.phi)
    
    xBavg2 = np.divide(xBs_by_Q2xBtphi, count_by_Q2xBtphi, where  = count_by_Q2xBtphi>0, out = np.zeros(count_by_Q2xBtphi.shape))
    Q2avg2 = np.divide(Q2s_by_Q2xBtphi, count_by_Q2xBtphi, where  = count_by_Q2xBtphi>0, out = np.zeros(count_by_Q2xBtphi.shape))
    tavg2  = np.divide(ts_by_Q2xBtphi, count_by_Q2xBtphi, where  = count_by_Q2xBtphi>0, out = np.zeros(count_by_Q2xBtphi.shape))
    phiavg2  = np.divide(phis_by_Q2xBtphi, count_by_Q2xBtphi, where  = count_by_Q2xBtphi>0, out = np.zeros(count_by_Q2xBtphi.shape))

    xBavg2 = xBavg2[activebins]
    Q2avg2 = Q2avg2[activebins]
    tavg2 = tavg2[activebins]
    phiavg2 = phiavg2[activebins]

    dict_xB2 = {activebins[i]:xBavg2[i] for i in range(len(xBavg2))}
    dict_Q22 = {activebins[i]:Q2avg2[i] for i in range(len(Q2avg2))}
    dict_t2 = {activebins[i]:tavg2[i] for i in range(len(tavg2))}
    dict_phi2 = {activebins[i]:phiavg2[i] for i in range(len(phiavg2))}

    df.loc[:, "xBavg2"] = df_global.Q2xBtphibin.map(dict_xB2)
    df.loc[:, "Q2avg2"] = df_global.Q2xBtphibin.map(dict_Q22)
    df.loc[:, "tavg2"] = df_global.Q2xBtphibin.map(dict_t2)
    df.loc[:, "phiavg2"] = df_global.Q2xBtphibin.map(dict_phi2)

    local1 = df1.loc[df1.Q2xBtphibin == Q2xBtphibin]
    local2 = df2.loc[df2.Q2xBtphibin == Q2xBtphibin]
    local3 = df3.loc[df3.Q2xBtphibin == Q2xBtphibin]
    local4 = df4.loc[df4.Q2xBtphibin == Q2xBtphibin]
    local = pd.concat([local1, local2, local3, local4])
    if len(local)==0:
        return 0, 0, 0, 0, 0, 0, 0, 0, 0

    print(Q2xBtphibin, len(local))
    xBavg2 = local.xB.mean()
    Q2avg2 = local.Q2.mean()
    tavg2 =  local.t1.mean()
    phiavg2 = local.phi1.mean()
    local2 = local.loc[(local.GenWeight>0)&(local.BornWeight>0), :]
    xsecavg = 1/((1/local2.GenWeight).mean())
    xsecavg_born = 1/((1/local2.BornWeight).mean())

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
                    local2 = local.loc[(local.Q2>=qmin) & (local.Q2<qmax) & (local.xB>=xmin) & (local.xB<xmax) & (local.t1>=tmin) & (local.t1<tmax) & (local.phi1>=pmin) & (local.phi1<pmax)]
                    if(len(local2)):
                        count += 1
    local2 = 0                                                                                     
    return count/N1/N2/N3/N4*(Q2_f - Q2_i)*(xB_f - xB_i)*(t_f - t_i)*np.radians(phi_f - phi_i), xBavg1, Q2avg1, tavg1, xBavg2, Q2avg2, tavg2, phiavg2, xsecavg, xsecavg_born


# import argparse
# parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# parser.add_argument("-p","--polarity", help="polarity", default="inbending")
# args = parser.parse_args()

print("inbending")
parent_MC = "/volatile/clas12/sangbaek/nov2021/convPkl_Gen/inb/dvcs/"
parent_bhMC = "/volatile/clas12/sangbaek/nov2021/convPkl_Gen/inb/bh/"

# #dvcs inb 50 nA
# print("reading dvcs inb 50 nA")
# df_3987_corr = pd.read_pickle(parent_MC + "3987.pkl")
# df_3987_corr = df_3987_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4124_corr = pd.read_pickle(parent_MC + "4124.pkl")
# df_4124_corr = df_4124_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4139_corr = pd.read_pickle(parent_MC + "4139.pkl")
# df_4139_corr = df_4139_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4181_corr = pd.read_pickle(parent_MC + "4181.pkl")       
# df_4181_corr = df_4181_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4182_corr = pd.read_pickle(parent_MC + "4182.pkl")
# df_4182_corr = df_4182_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4397_corr = pd.read_pickle(parent_MC + "4397.pkl")
# df_4397_corr = df_4397_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4528_corr = pd.read_pickle(parent_MC + "4528.pkl")
# df_4528_corr = df_4528_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4529_corr = pd.read_pickle(parent_MC + "4529.pkl")
# df_4529_corr = df_4529_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4535_corr = pd.read_pickle(parent_MC + "4535.pkl")
# df_4535_corr = df_4535_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4539_corr = pd.read_pickle(parent_MC + "4539.pkl")
# df_4539_corr = df_4539_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]

# #dvcs inb 55 nA
# print("reading dvcs inb 55 nA")
# df_4186_corr = pd.read_pickle(parent_MC + "4186.pkl")
# df_4186_corr = df_4186_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4545_corr = pd.read_pickle(parent_MC + "4545.pkl")
# df_4545_corr = df_4545_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]

# #dvcs inb 45 nA
# print("reading dvcs inb 45 nA")
# df_4188_corr = pd.read_pickle(parent_MC + "4188.pkl")
# df_4188_corr = df_4188_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4547_corr = pd.read_pickle(parent_MC + "4547.pkl")
# df_4547_corr = df_4547_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]

# #dvcs inb 0 nA
# print("reading dvcs inb 0 nA")
# df_4192_corr = pd.read_pickle(parent_MC + "4192.pkl")
# df_4192_corr = df_4192_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4561_corr = pd.read_pickle(parent_MC + "4561.pkl")
# df_4561_corr = df_4561_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]

#bh inb 50 nA
print("reading bh inb 50 nA")
df_4238_corr = pd.read_pickle(parent_bhMC + "4238.pkl")
df_4238_corr = df_4238_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtbin", "Q2xBtphibin", "GenWeight", "BornWeight"]]
df_4542_corr = pd.read_pickle(parent_bhMC + "4542.pkl")
df_4542_corr = df_4542_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtbin", "Q2xBtphibin", "GenWeight", "BornWeight"]]

# dvcsSimInb50nA = pd.concat([df_3987_corr, df_4124_corr, df_4139_corr, df_4181_corr, df_4182_corr, df_4397_corr, df_4528_corr, df_4529_corr, df_4535_corr, df_4539_corr])
# dvcsSimInb55nA = pd.concat([df_4186_corr, df_4545_corr])
# dvcsSimInb45nA = pd.concat([df_4188_corr, df_4547_corr])
# dvcsSimInb0nA = pd.concat([df_4192_corr, df_4561_corr])
# bhSimInb50nA = pd.concat([df_4238_corr, df_4542_corr])

print("outbending")
parent_MC = "/volatile/clas12/sangbaek/nov2021/convPkl_Gen/outb/dvcs/"
parent_bhMC = "/volatile/clas12/sangbaek/nov2021/convPkl_Gen/outb/bh/"

# #dvcs outb 50 nA
# print("reading dvcs outb 50 nA")
# df_4240_corr = pd.read_pickle(parent_MC + "4240.pkl")
# df_4240_corr = df_4240_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4250_corr = pd.read_pickle(parent_MC + "4250.pkl")
# df_4250_corr = df_4250_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4251_corr = pd.read_pickle(parent_MC + "4251.pkl")
# df_4251_corr = df_4251_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4252_corr = pd.read_pickle(parent_MC + "4252.pkl")
# df_4252_corr = df_4252_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4255_corr = pd.read_pickle(parent_MC + "4255.pkl")
# df_4255_corr = df_4255_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4398_corr = pd.read_pickle(parent_MC + "4398.pkl")
# df_4398_corr = df_4398_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4532_corr = pd.read_pickle(parent_MC + "4532.pkl")
# df_4532_corr = df_4532_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4534_corr = pd.read_pickle(parent_MC + "4534.pkl")
# df_4534_corr = df_4534_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4540_corr = pd.read_pickle(parent_MC + "4540.pkl")
# df_4540_corr = df_4540_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4541_corr = pd.read_pickle(parent_MC + "4541.pkl")
# df_4541_corr = df_4541_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]

# #dvcs outb 40 nA
# print("reading dvcs outb 40 nA")
# df_4263_corr = pd.read_pickle(parent_MC + "4263.pkl")
# df_4263_corr = df_4263_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4546_corr = pd.read_pickle(parent_MC + "4546.pkl")
# df_4546_corr = df_4544_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]

# #dvcs outb 0 nA
# print("reading dvcs outb 40 nA")
# df_4262_corr = pd.read_pickle(parent_MC + "4262.pkl")
# df_4262_corr = df_4262_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4554_corr = pd.read_pickle(parent_MC + "4554.pkl")
# df_4554_corr = df_4546_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]

# #dvcs outb +1.01, 40 nA
# print("reading dvcs outb +1.01 40 nA")
# df_4266_corr = pd.read_pickle(parent_MC + "4266.pkl")
# df_4266_corr = df_4266_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
# df_4562_corr = pd.read_pickle(parent_MC + "4562.pkl")
# df_4562_corr = df_4562_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]

#bh outb 50 nA
print("reading bh outb 50 nA")
df_4249_corr = pd.read_pickle(parent_bhMC + "4249.pkl")
df_4249_corr = df_4249_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtbin", "Q2xBtphibin", "GenWeight", "BornWeight"]]
df_4544_corr = pd.read_pickle(parent_bhMC + "4544.pkl")
df_4544_corr = df_4544_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtbin", "Q2xBtphibin", "GenWeight", "BornWeight"]]

# dvcsSimOutb50nA = pd.concat([df_4240_corr, df_4250_corr, df_4251_corr, df_4252_corr, df_4255_corr, df_4398_corr, df_4532_corr, df_4534_corr, df_4540_corr, df_4541_corr])
# dvcsSimOutb40nA = pd.concat([df_4263_corr, df_4546_corr])
# dvcsSimOutb0nA = pd.concat([df_4262_corr, df_4554_corr])
# dvcsSimOutb40nAT = pd.concat([df_4266_corr, df_4562_corr])
# bhSimOutb50nA = pd.concat([df_4249_corr, df_4544_corr])

# dvcsBHSim = pd.concat([dvcsSimOutb50nA, dvcsSimOutb40nA, dvcsSimOutb0nA, dvcsSimOutb40nAT, bhSimOutb50nA, dvcsSimInb50nA, dvcsSimInb55nA, dvcsSimInb45nA, dvcsSimInb0nA, bhSimInb50nA])
# dvcsBHSim = pd.concat([bhSimOutb50nA, bhSimInb50nA])


# dvcsBHSim = dvcsBHSim.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "config"]]

TrueVols = []
xBavgs1 = []
Q2avgs1 = []
tavgs1 = []
xBavgs2 = []
Q2avgs2 = []
tavgs2 = []
phiavgs2 = []
xsecavgs = []
xsecavgs_born = []

df1 = df_4238_corr
df2 = df_4542_corr
df3 = df_4249_corr
df4 = df_4544_corr
# for i in range(len(df_global)):
#     if i%50==0:
#         print("{}th event".format(i))
#     TrueVol, xBavg1, Q2avg1, tavg1, xBavg2, Q2avg2, tavg2, phiavg2, xsecavg, xsecavg_born = TruebinVol(df_global.Q2bin[i], df_global.xBbin[i], df_global.tbin[i], df_global.phibin[i], df_global.Q2xBtbin[i], df_global.Q2xBtphibin[i], df1, df2, df3, df4, 6, 6, 6, 6)
#     # TrueVol = TruebinVol(df_global.Q2bin[i], df_global.xBbin[i], df_global.tbin[i], df_global.phibin[i], df_global.Q2xBtphibin[i], dvcsBHSim, 6, 6, 6, 6)
#     TrueVols.append(TrueVol)
#     xBavgs1.append(xBavg1)
#     Q2avgs1.append(Q2avg1)
#     tavgs1.append(tavg1)
#     xBavgs2.append(xBavg2)
#     Q2avgs2.append(Q2avg2)
#     tavgs2.append(tavg2)
#     phiavgs2.append(phiavg2)
#     xsecavgs.append(xsecavg)
#     xsecavgs_born.append(xsecavg_born)
# df_global.loc[:, "TruebinVol"] = TrueVols
# df_global.loc[:, "xBavg1"] = xBavgs1
# df_global.loc[:, "Q2avg1"] = Q2avgs1
# df_global.loc[:, "tavg1"] = tavgs1
# df_global.loc[:, "xBavg2"] = xBavgs2
# df_global.loc[:, "Q2avg2"] = Q2avgs2
# df_global.loc[:, "tavg2"] = tavgs2
# df_global.loc[:, "phiavg2"] = phiavgs2
# df_global.loc[:, "xsecavg"] = xsecavgs
# df_global.loc[:, "xsecavg_born"] = xsecavgs_born

df_global = TruebinVol(df_global, df1, df2, df3, df4, 6, 6, 6, 6)
df_global.to_pickle("/volatile/clas12/sangbaek/clas12DVCS/results/truebinVolMay_Gen3.pkl")