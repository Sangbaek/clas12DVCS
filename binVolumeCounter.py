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

df_global = pd.read_pickle("/volatile/clas12/sangbaek/clas12DVCS/df_global_Mar.pkl")

def TruebinVol(Q2bin, xBbin, tbin, phibin, Q2xBtbin, Q2xBtphibin, df1, df2, df3, df4, N1=10, N2=10, N3=10, N4=10):
    
    count = 0 
    
    Q2_i = Q2bin_i[Q2bin]
    Q2_f = Q2bin_f[Q2bin]
    xB_i = xBbin_i[Q2bin][xBbin]
    xB_f = xBbin_f[Q2bin][xBbin]
    t_i = tbin_i[tbin]
    t_f = tbin_f[tbin]
    phi_i = phibin_i[phibin]
    phi_f = phibin_f[phibin]
    
    local1 = df1.loc[df1.Q2xBtbin == Q2xBtbin]
    local2 = df2.loc[df2.Q2xBtbin == Q2xBtbin]
    local3 = df3.loc[df3.Q2xBtbin == Q2xBtbin]
    local4 = df4.loc[df4.Q2xBtbin == Q2xBtbin]
    local = pd.concat([local1, local2, local3, local4])

    xBavg = sum(local.xB * local.GenWeight)/sum(local.GenWeight)
    Q2avg = sum(local.Q2 * local.GenWeight)/sum(local.GenWeight)
    tavg = sum(local.t * local.GenWeight)/sum(local.GenWeight)
    if len(local)==0:
        return 0, 0, 0, 0, 0, 0

    local1 = df1.loc[df1.Q2xBtphibin == Q2xBtphibin]
    local2 = df2.loc[df2.Q2xBtphibin == Q2xBtphibin]
    local3 = df3.loc[df3.Q2xBtphibin == Q2xBtphibin]
    local4 = df4.loc[df4.Q2xBtphibin == Q2xBtphibin]
    local = pd.concat([local1, local2, local3, local4])
    if len(local)==0:
        return 0, 0, 0, 0, 0, 0

    print(Q2xBtphibin, len(local))
    phiavg = sum(local.phi * local.GenWeight)/sum(local.GenWeight)
    xsecavg = local.GenWeight.mean()

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
                    if(sum((local.Q2>=qmin) & (local.Q2<qmax) & (local.xB>=xmin) & (local.xB<xmax) & (local.t1>=tmin) & (local.t1<tmax) & (local.phi1>=pmin) & (local.phi1<pmax)) > 0):
                    # if(len(local2)):
                        count += 1
    local2 = 0                                                                                     
    return count/N1/N2/N3/N4*(Q2_f - Q2_i)*(xB_f - xB_i)*(t_f - t_i)*np.radians(phi_f - phi_i), xBavg, Q2avg, tavg, phiavg, xsecavg


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
df_4238_corr = df_4238_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
df_4542_corr = pd.read_pickle(parent_bhMC + "4542.pkl")
df_4542_corr = df_4542_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]

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
df_4249_corr = df_4249_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]
df_4544_corr = pd.read_pickle(parent_bhMC + "4544.pkl")
df_4544_corr = df_4544_corr.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "GenWeight"]]

# dvcsSimOutb50nA = pd.concat([df_4240_corr, df_4250_corr, df_4251_corr, df_4252_corr, df_4255_corr, df_4398_corr, df_4532_corr, df_4534_corr, df_4540_corr, df_4541_corr])
# dvcsSimOutb40nA = pd.concat([df_4263_corr, df_4546_corr])
# dvcsSimOutb0nA = pd.concat([df_4262_corr, df_4554_corr])
# dvcsSimOutb40nAT = pd.concat([df_4266_corr, df_4562_corr])
# bhSimOutb50nA = pd.concat([df_4249_corr, df_4544_corr])

# dvcsBHSim = pd.concat([dvcsSimOutb50nA, dvcsSimOutb40nA, dvcsSimOutb0nA, dvcsSimOutb40nAT, bhSimOutb50nA, dvcsSimInb50nA, dvcsSimInb55nA, dvcsSimInb45nA, dvcsSimInb0nA, bhSimInb50nA])
# dvcsBHSim = pd.concat([bhSimOutb50nA, bhSimInb50nA])


# dvcsBHSim = dvcsBHSim.loc[:, ["xB", "Q2", "t1", "phi1", "Q2xBtphibin", "config"]]

TrueVols = []
xBavgs = []
Q2avgs = []
tavgs = []
phiavgs = []
xsecavgs = []

df1 = df_4238_corr
df2 = df_4542_corr
df3 = df_4249_corr
df4 = df_4544_corr
for i in range(len(df_global)):
    if i%50==0:
        print("{}th event".format(i))
    TrueVol, xBavg, Q2avg, tavg, phiavg, xsecavg = TruebinVol(df_global.Q2bin[i], df_global.xBbin[i], df_global.tbin[i], df_global.phibin[i], df_global.Q2xBtbin[i], df_global.Q2xBtphibin[i], df1, df2, df3, df4, 6, 6, 6, 6)
    # TrueVol = TruebinVol(df_global.Q2bin[i], df_global.xBbin[i], df_global.tbin[i], df_global.phibin[i], df_global.Q2xBtphibin[i], dvcsBHSim, 6, 6, 6, 6)
    TrueVols.append(TrueVol)
    xBavgs.append(xBavg)
    Q2avgs.append(Q2avg)
    tavgs.append(tavg)
    phiavgs.append(phiavg)
    xsecavgs.append(xsecavg)
df_global.loc[:, "TruebinVol"] = TrueVol
df_global.loc[:, "xBavg"] = xBavgs
df_global.loc[:, "Q2avg"] = Q2avgs
df_global.loc[:, "tavg"] = tavgs
df_global.loc[:, "phiavg"] = phiavgs
df_global.loc[:, "xsecavg"] = xsecavgs
df_global.to_pickle("/volatile/clas12/sangbaek/clas12DVCS/results/truebinVol_Gen.pkl")