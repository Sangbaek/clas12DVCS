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

local502 = pd.read_pickle("df_globalFeb.pkl")
local502 = local502.loc[local502.Q2xBt == "502"]

def numberingDF(total, Q2bin_i=Q2bin_i, Q2bin_f=Q2bin_f, xBbin_i=xBbin_i, xBbin_f=xBbin_f, tbin_i=tbin_i, tbin_f=tbin_f, goodBins=goodBins, badBins=badBins):
    df_allBins = {}
    Q2xBtphi = 1281#0

    for Q2bin in range(5, 6):#Q2 bin
        for xBbin in range(0, 1):
            for tbin in range(2, 3):
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

parent_exp = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/exp/"
parent_MC = "/volatile/clas12/sangbaek/nov2021/convPkl/dvcs/"
parent_bhMC = "/volatile/clas12/sangbaek/nov2021/convPkl/bh/"
parent_MC_bkg1g = "/volatile/clas12/sangbaek/nov2021/convPkl/bkg_1g/"
parent_MC_bkg2g = "/volatile/clas12/sangbaek/nov2021/convPkl/bkg_2g/"
parent_Gen = "/volatile/clas12/sangbaek/nov2021/convPkl_Gen/dvcs/inb/"
parent_bhGen = "/volatile/clas12/sangbaek/nov2021/convPkl_Gen/bh/inb/"



#dvcs inb 50 nA
print("reading dvcs inb 50 nA")
df_3987_corr = pd.read_pickle(parent_MC + "3987.pkl")
df_4124_corr = pd.read_pickle(parent_MC + "4124.pkl")
df_4139_corr = pd.read_pickle(parent_MC + "4139.pkl")
df_4181_corr = pd.read_pickle(parent_MC + "4181.pkl")
df_4182_corr = pd.read_pickle(parent_MC + "4182.pkl")
df_3987_Gen = pd.read_pickle(parent_Gen + "3987.pkl")
df_4124_Gen = pd.read_pickle(parent_Gen + "4124.pkl")
df_4139_Gen = pd.read_pickle(parent_Gen + "4139.pkl")
df_4181_Gen = pd.read_pickle(parent_Gen + "4181.pkl")
df_4182_Gen = pd.read_pickle(parent_Gen + "4182.pkl")

df_3987_Gen.loc[:, "event"] = df_3987_Gen.index
df_4124_Gen.loc[:, "event"] = df_4124_Gen.index
df_4139_Gen.loc[:, "event"] = df_4139_Gen.index
df_4181_Gen.loc[:, "event"] = df_4181_Gen.index
df_4182_Gen.loc[:, "event"] = df_4182_Gen.index

#dvcs inb 55 nA
print("reading dvcs inb 55 nA")
df_4186_corr = pd.read_pickle(parent_MC + "4186.pkl")
df_4186_Gen = pd.read_pickle(parent_Gen + "4186.pkl")

df_4186_Gen.loc[:, "event"] = df_4186_Gen.index

#dvcs inb 45 nA
print("reading dvcs inb 45 nA")
df_4188_corr = pd.read_pickle(parent_MC + "4188.pkl")
df_4188_Gen = pd.read_pickle(parent_Gen + "4188.pkl")

df_4188_Gen.loc[:, "event"] = df_4188_Gen.index

#dvcs inb 0 nA
print("reading dvcs inb 0 nA")
df_4192_corr = pd.read_pickle(parent_MC + "4192.pkl")
df_4192_Gen = pd.read_pickle(parent_Gen + "4192.pkl")

df_4192_Gen.loc[:, "event"] = df_4192_Gen.index

#bh inb 50 nA
print("reading bh inb 50 nA")
df_4238_corr = pd.read_pickle(parent_bhMC + "4238.pkl")
df_4238_Gen = pd.read_pickle(parent_bhGen + "4238.pkl")

df_4238_Gen.loc[:, "event"] = df_4238_Gen.index

#bkg1g inb 50 nA
df_4076_1g_corr = pd.read_pickle(parent_MC_bkg1g + "4076.pkl")
df_4202_1g_corr = pd.read_pickle(parent_MC_bkg1g + "4202.pkl")
df_4209_1g_corr = pd.read_pickle(parent_MC_bkg1g + "4209.pkl")
#bkg2g inb 50 nA
df_4076_2g_corr = pd.read_pickle(parent_MC_bkg2g + "4076.pkl")
df_4202_2g_corr = pd.read_pickle(parent_MC_bkg2g + "4202.pkl")
df_4209_2g_corr = pd.read_pickle(parent_MC_bkg2g + "4209.pkl")

#Exp
epgExpInb = pd.read_pickle(parent_exp + "dvcs.pkl")
pi0ExpInb = pd.read_pickle(parent_exp + "pi0.pkl")

print("done with reading inbending")

print("3987")
df_3987_corr = numberingDF(df_3987_corr)
df_3987_Gen = df_3987_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
df_3987_Gen = numberingDF(df_3987_Gen)
print("4124")
df_4124_corr = numberingDF(df_4124_corr)
df_4124_Gen = df_4124_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
df_4124_Gen = numberingDF(df_4124_Gen)
print("4139")
df_4139_corr = numberingDF(df_4139_corr)
df_4139_Gen = df_4139_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
df_4139_Gen = numberingDF(df_4139_Gen)
print("4181")
df_4181_corr = numberingDF(df_4181_corr)
df_4181_Gen = df_4181_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
df_4181_Gen = numberingDF(df_4181_Gen)
print("4182")
df_4182_corr = numberingDF(df_4182_corr)
df_4182_Gen = df_4182_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
df_4182_Gen = numberingDF(df_4182_Gen)
print("4186")
df_4186_corr = numberingDF(df_4186_corr)
df_4186_Gen = df_4186_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
df_4186_Gen = numberingDF(df_4186_Gen)
print("4188")
df_4188_corr = numberingDF(df_4188_corr)
df_4188_Gen = df_4188_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
df_4188_Gen = numberingDF(df_4188_Gen)
print("4192")
df_4192_corr = numberingDF(df_4192_corr)
df_4192_Gen = df_4192_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
df_4192_Gen = numberingDF(df_4192_Gen)
print("4238")
df_4238_corr = numberingDF(df_4238_corr)
df_4238_Gen = df_4238_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
df_4238_Gen = numberingDF(df_4238_Gen)

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
epgExpInb = numberingDF(epgExpInb)
print("exp pi0")
pi0ExpInb = numberingDF(pi0ExpInb)

print("done with numbering")


dvcsGenInb50nA = pd.concat([df_3987_Gen, df_4124_Gen, df_4139_Gen, df_4181_Gen, df_4182_Gen])
dvcsGenInb55nA = df_4186_Gen
dvcsGenInb45nA = df_4188_Gen
dvcsGenInb0nA = df_4192_Gen
bhGenInb50nA = df_4238_Gen

dvcsSimInb50nA = pd.concat([df_3987_corr, df_4124_corr, df_4139_corr, df_4181_corr, df_4182_corr])
dvcsSimInb55nA = df_4186_corr
dvcsSimInb45nA = df_4188_corr
dvcsSimInb0nA = df_4192_corr
bhSimInb50nA = df_4238_corr

bkgSimInb50nA = pd.concat([df_4076_1g_corr, df_4202_1g_corr, df_4209_1g_corr])
pi0SimInb50nA = pd.concat([df_4076_2g_corr, df_4202_2g_corr, df_4209_2g_corr])

def countDF(total, df_global, colName = "new"):
    numbers1 = []
    numbers2 = []
    numbers3 = []
    if 'Q2xBtphi' not in total:
        total = numberingDF(total, Q2bin_i, Q2bin_f, xBbin_i, xBbin_f, tbin_i, tbin_f, goodBins, badBins, df_global)
    for i in range(1281, 1302):
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
    for i in range(1281, 1302):
        if i%10==0:
            print(i)
        number = sum((total.Q2xBtphi == i))
        numbers.append(number)
    df_global.loc[:, colName] = numbers
    return df_global

#total acc.
local502 = countGenDF(dvcsGenInb50nA, local502, "dvcsGenInb50nA")
local502 = countDF(dvcsSimInb50nA, local502, "dvcsSimInb50nA")

local502 = countGenDF(dvcsGenInb55nA, local502, "dvcsGenInb55nA")
local502 = countDF(dvcsSimInb55nA, local502, "dvcsSimInb55nA")

local502 = countGenDF(dvcsGenInb45nA, local502, "dvcsGenInb45nA")
local502 = countDF(dvcsSimInb45nA, local502, "dvcsSimInb45nA")

local502 = countGenDF(dvcsGenInb0nA, local502, "dvcsGenInb0nA")
local502 = countDF(dvcsSimInb0nA, local502, "dvcsSimInb0nA")

local502 = countGenDF(bhGenInb50nA, local502, "bhGenInb50nA")
local502 = countDF(bhSimInb50nA, local502, "bhSimInb50nA")

#fractional acc.
local502 = countGenDF(dvcsGenInb50nA.loc[ dvcsGenInb50nA.radMode == 1], local502, "dvcsGenInb50nA_non")
local502 = countDF(dvcsSimInb50nA.loc[ dvcsSimInb50nA.radMode == 1], local502, "dvcsSimInb50nA_non")

local502 = countGenDF(dvcsGenInb55nA.loc[ dvcsGenInb55nA.radMode == 1], local502, "dvcsGenInb55nA_non")
local502 = countDF(dvcsSimInb55nA.loc[ dvcsSimInb55nA.radMode == 1], local502, "dvcsSimInb55nA_non")

local502 = countGenDF(dvcsGenInb45nA.loc[ dvcsGenInb45nA.radMode == 1], local502, "dvcsGenInb45nA_non")
local502 = countDF(dvcsSimInb45nA.loc[ dvcsSimInb45nA.radMode == 1], local502, "dvcsSimInb45nA_non")

local502 = countGenDF(dvcsGenInb0nA.loc[ dvcsGenInb0nA.radMode == 1], local502, "dvcsGenInb0nA_non")
local502 = countDF(dvcsSimInb0nA.loc[ dvcsSimInb0nA.radMode == 1], local502, "dvcsSimInb0nA_non")

local502 = countGenDF(bhGenInb50nA.loc[ bhGenInb50nA.radMode == 1], local502, "bhGenInb50nA_non")
local502 = countDF(bhSimInb50nA.loc[ bhSimInb50nA.radMode == 1], local502, "bhSimInb50nA_non")

#fractional acc. s peak
local502 = countGenDF(dvcsGenInb50nA.loc[ dvcsGenInb50nA.radMode == 2], local502, "dvcsGenInb50nA_sPeak")
local502 = countDF(dvcsSimInb50nA.loc[ dvcsSimInb50nA.radMode == 2], local502, "dvcsSimInb50nA_sPeak")

local502 = countGenDF(dvcsGenInb55nA.loc[ dvcsGenInb55nA.radMode == 2], local502, "dvcsGenInb55nA_sPeak")
local502 = countDF(dvcsSimInb55nA.loc[ dvcsSimInb55nA.radMode == 2], local502, "dvcsSimInb55nA_sPeak")

local502 = countGenDF(dvcsGenInb45nA.loc[ dvcsGenInb45nA.radMode == 2], local502, "dvcsGenInb45nA_sPeak")
local502 = countDF(dvcsSimInb45nA.loc[ dvcsSimInb45nA.radMode == 2], local502, "dvcsSimInb45nA_sPeak")

local502 = countGenDF(dvcsGenInb0nA.loc[ dvcsGenInb0nA.radMode == 2], local502, "dvcsGenInb0nA_sPeak")
local502 = countDF(dvcsSimInb0nA.loc[ dvcsSimInb0nA.radMode == 2], local502, "dvcsSimInb0nA_sPeak")

local502 = countGenDF(bhGenInb50nA.loc[ bhGenInb50nA.radMode == 2], local502, "bhGenInb50nA_sPeak")
local502 = countDF(bhSimInb50nA.loc[ bhSimInb50nA.radMode == 2], local502, "bhSimInb50nA_sPeak")

#fractional acc. p peak
local502 = countGenDF(dvcsGenInb50nA.loc[ dvcsGenInb50nA.radMode == 3], local502, "dvcsGenInb50nA_pPeak")
local502 = countDF(dvcsSimInb50nA.loc[ dvcsSimInb50nA.radMode == 3], local502, "dvcsSimInb50nA_pPeak")

local502 = countGenDF(dvcsGenInb55nA.loc[ dvcsGenInb55nA.radMode == 3], local502, "dvcsGenInb55nA_pPeak")
local502 = countDF(dvcsSimInb55nA.loc[ dvcsSimInb55nA.radMode == 3], local502, "dvcsSimInb55nA_pPeak")

local502 = countGenDF(dvcsGenInb45nA.loc[ dvcsGenInb45nA.radMode == 3], local502, "dvcsGenInb45nA_pPeak")
local502 = countDF(dvcsSimInb45nA.loc[ dvcsSimInb45nA.radMode == 3], local502, "dvcsSimInb45nA_pPeak")

local502 = countGenDF(dvcsGenInb0nA.loc[ dvcsGenInb0nA.radMode == 3], local502, "dvcsGenInb0nA_pPeak")
local502 = countDF(dvcsSimInb0nA.loc[ dvcsSimInb0nA.radMode == 3], local502, "dvcsSimInb0nA_pPeak")

local502 = countGenDF(bhGenInb50nA.loc[ bhGenInb50nA.radMode == 3], local502, "bhGenInb50nA_pPeak")
local502 = countDF(bhSimInb50nA.loc[ bhSimInb50nA.radMode == 3], local502, "bhSimInb50nA_pPeak")

local502 = countDF(epgExpInb, local502, "epgExpInb")
local502 = countDF(pi0ExpInb, local502, "pi0ExpInb")
local502 = countDF(bkgSimInb50nA, local502, "bkgSimInb50nA")
local502 = countDF(pi0SimInb50nA, local502, "pi0SimInb50nA")


local502.to_pickle("local502_inb.pkl")