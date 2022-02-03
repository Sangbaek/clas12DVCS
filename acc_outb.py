import pandas as pd
import numpy as np

M = 0.938272081 # target mass

k= 2*M*(np.sqrt(0.8**2+M**2)-M)
x1 = 1/2/M/8.604
x2 = 1/2/M/3
x3 = 1-(4-M**2)/2/M/3
x4 = 2*(1-np.sqrt((k-M**4+4*M**2)/k))/(M**4/k - 4*M**2/k)
x5 = (-1+np.sqrt(1+4*M**2/k*(1+M/8.604/2)))/2/M**2*k
y1 = 1
y2 = 1.456
y3 = 6*M*x3
y4 = (4-M**2)*x4/(1-x4)
y5 = 2*M*8.604*x5

c0 = y2/2/M/8.604
d0 = y2/2/M/3
c1 = np.sqrt(y2*y3)/2/M/8.604
d1 = np.sqrt(y2*y3)/2/M/3
c2 = y3/2/M/8.604
d2 = x3
c3 = np.sqrt(y3*y4)/2/M/8.604
d3 = 1/(1+(4-M*M)/np.sqrt(y3*y4))
c4 = y4/2/M/8.604
d4 = x5
c5 = np.sqrt(y4*y5)/2/M/8.604
d5 = x5

tbin0 = 0.088#2*M*(np.sqrt(0.3**2+M**2)-M)
tbin1 = 0.168#2*M*(np.sqrt(0.42**2+M**2)-M)
tbin2 = 0.234#2*M*(np.sqrt(0.5**2+M**2)-M)
tbin3 = 0.414#2*M*(np.sqrt(0.68**2+M**2)-M)
tbin4 = 0.553#2*M*(np.sqrt(0.8**2+M**2)-M)

Q2bin_i = [y1, y2, np.sqrt(y2*y3), y3, np.sqrt(y3*y4), y4, np.sqrt(y4*y5)]
Q2bin_f = [y2, np.sqrt(y2*y3), y3, np.sqrt(y3*y4), y4, np.sqrt(y4*y5), y5]
xBbin_i = {0:[[x1, c0], c1, c2], 1: [[c0, c1], c2, c3], 2: [[c1, c2], c2, c3, c4], 3: [[c2, c3], c3, c4, c5], 4: [[c3, c4], c4, c5], 5: [[c4, c5], c5], 6: [[c5, x5]]}
xBbin_f = {0: [c1, c2, [x2, d0]], 1: [c2, c3, [x2, d1]], 2: [c2, c3, c4, [d1, d2]], 3: [c3, c4, c5, [d2, d3]], 4: [c4, c5, [d3, d4]], 5: [c5, x5], 6: [x5]}
tbin_i = [tbin0, tbin1, tbin2, tbin3]
tbin_f = [tbin1, tbin2, tbin3, tbin4]
phibin = [0, 12, 24, 36, 48, 60, 72, 96, 120, 144, 168, 192, 216, 240, 264, 288, 300, 312, 324, 336, 348, 360]
phibin_i = phibin[:-1]
phibin_f = phibin[1:]

goodBins = ['000', '001', '002', '003', '010', '011', '012', '013', '020', '021', '022', '023', '100', '101', '102', '103', '110', '111', '112', '113', '120', '121', '122', '123', '200', '201', '202', '203', '210', '211', '212', '213', '220', '221', '222', '223', '231', '232', '233', '300', '301', '302', '303', '310', '311', '312', '313', '321', '322', '323', '332', '333', '400', '401', '402', '403', '411', '412', '413', '422', '423', '502', '503', '512', '513']
badBins = ['230', '320', '330', '331', '410', '420', '421', '500', '501', '510', '511', '600', '601', '602', '603']

df_global = pd.read_pickle("df_globalFeb.pkl")

def numberingDF(total, Q2bin_i=Q2bin_i, Q2bin_f=Q2bin_f, xBbin_i=xBbin_i, xBbin_f=xBbin_f, tbin_i=tbin_i, tbin_f=tbin_f, goodBins=goodBins, badBins=badBins, df_global=df_global):
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
                if Q2bin == len(Q2bin_i)-1:
                    local = local.loc[(local.Q2>=Q2_i) & (local.Q2<=Q2_f)]
                else:
                    local = local.loc[(local.Q2>=Q2_i) & (local.Q2<Q2_f)]
                #cut by xB
                #xB lower bound
                if xBbin == 0:
                    local = local.loc[local.Q2<=2*M*(10.604-2)*local.xB, :]
                else:
                    local = local.loc[local.xB>=xB_i] 
                #xB upper bound
                if (xBbin == len(xBbin_i[Q2bin])-1) & (Q2bin < 3):
                    local = local.loc[local.Q2>=2*M*3*local.xB]
                elif (xBbin == len(xBbin_i[Q2bin])-1) & (Q2bin < 5):
                    local = local.loc[local.Q2>=(4 - M*M)*local.xB/(1 - local.xB)]
                else:
                    local = local.loc[local.xB<xB_f]
                #cut by t
                if tbin == len(tbin_i)-1:
                    local = local.loc[(local.t1>=t_i) & (local.t1<=t_f)]
                else:
                    local = local.loc[(local.t1>=t_i) & (local.t1<t_f)]
                Q2xBtbin = "{}{}{}".format(Q2bin,xBbin,tbin)

                if Q2xBtbin in badBins:
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

parent_exp = "/volatile/clas12/sangbaek/nov2021/convPkl_full/outb/exp/"
parent_MC = "/volatile/clas12/sangbaek/nov2021/convPkl/dvcs/"
parent_bhMC = "/volatile/clas12/sangbaek/nov2021/convPkl/bh/"
parent_MC_bkg1g = "/volatile/clas12/sangbaek/nov2021/convPkl/bkg_1g/"
parent_MC_bkg2g = "/volatile/clas12/sangbaek/nov2021/convPkl/bkg_2g/"
parent_Gen = "/volatile/clas12/sangbaek/nov2021/convPkl_Gen/dvcs/outb/"
parent_bhGen = "/volatile/clas12/sangbaek/nov2021/convPkl_Gen/bh/outb/"



#dvcs outb 50 nA
print("reading dvcs outb 50nA")
df_4240_corr = pd.read_pickle(parent_MC + "4240.pkl")
df_4250_corr = pd.read_pickle(parent_MC + "4250.pkl")
df_4251_corr = pd.read_pickle(parent_MC + "4251.pkl")
df_4252_corr = pd.read_pickle(parent_MC + "4252.pkl")
df_4255_corr = pd.read_pickle(parent_MC + "4255.pkl")
df_4240_Gen = pd.read_pickle(parent_Gen + "4240.pkl")
df_4250_Gen = pd.read_pickle(parent_Gen + "4250.pkl")
df_4251_Gen = pd.read_pickle(parent_Gen + "4251.pkl")
df_4252_Gen = pd.read_pickle(parent_Gen + "4252.pkl")
df_4255_Gen = pd.read_pickle(parent_Gen + "4255.pkl")

df_4240_Gen.loc[:, "event"] = df_4240_Gen.index
df_4250_Gen.loc[:, "event"] = df_4250_Gen.index
df_4251_Gen.loc[:, "event"] = df_4251_Gen.index
df_4252_Gen.loc[:, "event"] = df_4252_Gen.index
df_4255_Gen.loc[:, "event"] = df_4255_Gen.index

#dvcs outb 40 nA
print("reading dvcs outb 40nA")
df_4263_corr = pd.read_pickle(parent_MC + "4263.pkl")
df_4263_Gen = pd.read_pickle(parent_Gen + "4263.pkl")

df_4263_Gen.loc[:, "event"] = df_4263_Gen.index

#dvcs outb 0 nA
print("reading dvcs outb 0 nA")
df_4262_corr = pd.read_pickle(parent_MC + "4262.pkl")
df_4262_Gen = pd.read_pickle(parent_Gen + "4262.pkl")

df_4262_Gen.loc[:, "event"] = df_4262_Gen.index

#dvcs outb 40 nA torus+1.01
print("reading dvcs outb 40nA torus+1.01")
df_4266_corr = pd.read_pickle(parent_MC + "4266.pkl")
df_4266_Gen = pd.read_pickle(parent_Gen + "4266.pkl")

df_4266_Gen.loc[:, "event"] = df_4266_Gen.index

#bh outb 50 nA
print("reading bh outb 50 nA")
df_4249_corr = pd.read_pickle(parent_bhMC + "4249.pkl")
df_4249_Gen = pd.read_pickle(parent_bhGen + "4249.pkl")

df_4249_Gen.loc[:, "event"] = df_4249_Gen.index

#bkg1g outb 50 nA
df_4076_1g_corr = pd.read_pickle(parent_MC_bkg1g + "4076.pkl")
df_4202_1g_corr = pd.read_pickle(parent_MC_bkg1g + "4202.pkl")
df_4209_1g_corr = pd.read_pickle(parent_MC_bkg1g + "4209.pkl")
#bkg2g outb 50 nA
df_4076_2g_corr = pd.read_pickle(parent_MC_bkg2g + "4076.pkl")
df_4202_2g_corr = pd.read_pickle(parent_MC_bkg2g + "4202.pkl")
df_4209_2g_corr = pd.read_pickle(parent_MC_bkg2g + "4209.pkl")

#Exp
epgExpOutb = pd.read_pickle(parent_exp + "dvcs.pkl")
pi0ExpOutb = pd.read_pickle(parent_exp + "pi0.pkl")

print("done with reading outbending")

print("4240")
df_4240_corr = numberingDF(df_4240_corr)
df_4240_Gen = df_4240_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
df_4240_Gen = numberingDF(df_4240_Gen)
print("4250")
df_4250_corr = numberingDF(df_4250_corr)
df_4250_Gen = df_4250_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
df_4250_Gen = numberingDF(df_4250_Gen)
print("4251")
df_4251_corr = numberingDF(df_4251_corr)
df_4251_Gen = df_4251_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
df_4251_Gen = numberingDF(df_4251_Gen)
print("4252")
df_4252_corr = numberingDF(df_4252_corr)
df_4252_Gen = df_4252_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
df_4252_Gen = numberingDF(df_4252_Gen)
print("4255")
df_4255_corr = numberingDF(df_4255_corr)
df_4255_Gen = df_4255_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
df_4255_Gen = numberingDF(df_4255_Gen)
print("4263")
df_4263_corr = numberingDF(df_4263_corr)
df_4263_Gen = df_4263_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
df_4263_Gen = numberingDF(df_4263_Gen)
print("4262")
df_4262_corr = numberingDF(df_4262_corr)
df_4262_Gen = df_4262_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
df_4262_Gen = numberingDF(df_4262_Gen)
print("4266")
df_4266_corr = numberingDF(df_4266_corr)
df_4266_Gen = df_4266_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
df_4266_Gen = numberingDF(df_4266_Gen)

print("4249")
df_4249_corr = numberingDF(df_4249_corr)
df_4249_Gen = df_4249_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
df_4249_Gen = numberingDF(df_4249_Gen)

print("4076")
df_4076_1g_corr = numberingDF(df_4076_1g_corr)
df_4076_2g_corr = numberingDF(df_4076_2g_corr)
print("4202")
df_4202_1g_corr = numberingDF(df_4202_1g_corr)
df_4202_2g_corr = numberingDF(df_4202_2g_corr)
print("4209")
df_4209_1g_corr = numberingDF(df_4209_1g_corr)
df_4209_2g_corr = numberingDF(df_4209_2g_corr)

print("exp dvcs")
epgExpOutb = numberingDF(epgExpOutb)
print("exp pi0")
pi0ExpOutb = numberingDF(pi0ExpOutb)

print("done with numbering")


dvcsGenOutb50nA = pd.concat([df_4240_Gen, df_4250_Gen, df_4251_Gen, df_4252_Gen, df_4255_Gen])
dvcsGenOutb40nA = df_4263_Gen
dvcsGenOutb0nA = df_4262_Gen
dvcsGenOutb40nAT = df_4266_Gen
bhGenOutb50nA = df_4249_Gen

dvcsSimOutb50nA = pd.concat([df_4240_corr, df_4250_corr, df_4251_corr, df_4252_corr, df_4255_corr])
dvcsSimOutb40nA = df_4263_corr
dvcsSimOutb0nA = df_4262_corr
dvcsSimOutb40nAT = df_4266_corr
bhSimOutb50nA = df_4249_corr

bkgSimOutb50nA = pd.concat([df_4076_1g_corr, df_4202_1g_corr, df_4209_1g_corr])
pi0SimOutb50nA = pd.concat([df_4076_2g_corr, df_4202_2g_corr, df_4209_2g_corr])

def countDF(total, df_global, colName = "new"):
    numbers1 = []
    numbers2 = []
    numbers3 = []
    if 'Q2xBtphi' not in total:
        total = numberingDF(total, Q2bin_i, Q2bin_f, xBbin_i, xBbin_f, tbin_i, tbin_f, goodBins, badBins, df_global)
    for i in range(len(df_global)):
        if i%50 == 0:
            print(i)
        number1 = sum((total.Q2xBtphi == i) & (total.config == 1))
        number2 = sum((total.Q2xBtphi == i) & (total.config == 2))
        number3 = sum((total.Q2xBtphi == i) & (total.config == 3))
        numbers1.append(number1)
        numbers2.append(number2)
        numbers3.append(number3)
    df_global.loc[:, colName+"1"] = numbers1
    df_global.loc[:, colName+"2"] = numbers2
    df_global.loc[:, colName+"3"] = numbers3
    return df_global

def countGenDF(total, df_global, colName = "new"):
    numbers = []
    if 'Q2xBtphi' not in total:
        print("invalid df input.")
        return df_global
    for i in range(len(df_global)):
        if i%10==0:
            print(i)
        number = sum((total.Q2xBtphi == i))
        numbers.append(number)
    df_global.loc[:, colName] = numbers
    return df_global

#total acc.
df_global = countGenDF(dvcsGenOutb50nA, df_global, "dvcsGenOutb50nA")
df_global = countDF(dvcsSimOutb50nA, df_global, "dvcsSimOutb50nA")

df_global = countGenDF(dvcsGenOutb40nA, df_global, "dvcsGenOutb40nA")
df_global = countDF(dvcsSimOutb40nA, df_global, "dvcsSimOutb40nA")

df_global = countGenDF(dvcsGenOutb0nA, df_global, "dvcsGenOutb0nA")
df_global = countDF(dvcsSimOutb0nA, df_global, "dvcsSimOutb0nA")

df_global = countGenDF(dvcsGenOutb40nAT, df_global, "dvcsGenOutb40nAT")
df_global = countDF(dvcsSimOutb40nAT, df_global, "dvcsSimOutb40nAT")

df_global = countGenDF(bhGenOutb50nA, df_global, "bhGenOutb50nA")
df_global = countDF(bhSimOutb50nA, df_global, "bhSimOutb50nA")

#fractional acc.
df_global = countGenDF(dvcsGenOutb50nA.loc[ radMode == 1], df_global, "dvcsGenOutb50nA_non")
df_global = countDF(dvcsSimOutb50nA.loc[ radMode == 1], df_global, "dvcsSimOutb50nA_non")

df_global = countGenDF(dvcsGenOutb40nA.loc[ radMode == 1], df_global, "dvcsGenOutb40nA_non")
df_global = countDF(dvcsSimOutb40nA.loc[ radMode == 1], df_global, "dvcsSimOutb40nA_non")

df_global = countGenDF(dvcsGenOutb0nA.loc[ radMode == 1], df_global, "dvcsGenOutb0nA_non")
df_global = countDF(dvcsSimOutb0nA.loc[ radMode == 1], df_global, "dvcsSimOutb0nA_non")

df_global = countGenDF(dvcsGenOutb40nAT.loc[ radMode == 1], df_global, "dvcsGenOutb40nAT_non")
df_global = countDF(dvcsSimOutb40nAT.loc[ radMode == 1], df_global, "dvcsSimOutb40nAT_non")

df_global = countGenDF(bhGenOutb50nA.loc[ radMode == 1], df_global, "bhGenOutb50nA_non")
df_global = countDF(bhSimOutb50nA.loc[ radMode == 1], df_global, "bhSimOutb50nA_non")

#fractional acc. s peak
df_global = countGenDF(dvcsGenOutb50nA.loc[ radMode == 2], df_global, "dvcsGenOutb50nA_sPeak")
df_global = countDF(dvcsSimOutb50nA.loc[ radMode == 2], df_global, "dvcsSimOutb50nA_sPeak")

df_global = countGenDF(dvcsGenOutb40nA.loc[ radMode == 2], df_global, "dvcsGenOutb40nA_sPeak")
df_global = countDF(dvcsSimOutb40nA.loc[ radMode == 2], df_global, "dvcsSimOutb40nA_sPeak")

df_global = countGenDF(dvcsGenOutb0nA.loc[ radMode == 2], df_global, "dvcsGenOutb0nA_sPeak")
df_global = countDF(dvcsSimOutb0nA.loc[ radMode == 2], df_global, "dvcsSimOutb0nA_sPeak")

df_global = countGenDF(dvcsGenOutb40nAT.loc[ radMode == 2], df_global, "dvcsGenOutb40nAT_sPeak")
df_global = countDF(dvcsSimOutb40nAT.loc[ radMode == 2], df_global, "dvcsSimOutb40nAT_sPeak")

df_global = countGenDF(bhGenOutb50nA.loc[ radMode == 2], df_global, "bhGenOutb50nA_sPeak")
df_global = countDF(bhSimOutb50nA.loc[ radMode == 2], df_global, "bhSimOutb50nA_sPeak")

#fractional acc. p peak
df_global = countGenDF(dvcsGenOutb50nA.loc[ radMode == 3], df_global, "dvcsGenOutb50nA_pPeak")
df_global = countDF(dvcsSimOutb50nA.loc[ radMode == 3], df_global, "dvcsSimOutb50nA_pPeak")

df_global = countGenDF(dvcsGenOutb40nA.loc[ radMode == 3], df_global, "dvcsGenOutb40nA_pPeak")
df_global = countDF(dvcsSimOutb40nA.loc[ radMode == 3], df_global, "dvcsSimOutb40nA_pPeak")

df_global = countGenDF(dvcsGenOutb0nA.loc[ radMode == 3], df_global, "dvcsGenOutb0nA_pPeak")
df_global = countDF(dvcsSimOutb0nA.loc[ radMode == 3], df_global, "dvcsSimOutb0nA_pPeak")

df_global = countGenDF(dvcsGenOutb40nT.loc[ radMode == 3], df_global, "dvcsGenOutb40nT_pPeak")
df_global = countDF(dvcsSimOutb40nT.loc[ radMode == 3], df_global, "dvcsSimOutb40nT_pPeak")

df_global = countGenDF(bhGenOutb50nA.loc[ radMode == 3], df_global, "bhGenOutb50nA_pPeak")
df_global = countDF(bhSimOutb50nA.loc[ radMode == 3], df_global, "bhSimOutb50nA_pPeak")

df_global = countDF(epgExpOutb, df_global, "epgExpOutb")
df_global = countDF(pi0ExpOutb, df_global, "pi0ExpOutb")
df_global = countDF(bkgSimOutb50nA, df_global, "bkgSimOutb50nA")
df_global = countDF(pi0SimOutb50nA, df_global, "pi0SimOutb50nA")


df_global.to_pickle("df_globalFeb_outb.pkl")