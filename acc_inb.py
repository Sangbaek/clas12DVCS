import pandas as pd
import numpy as np

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

df_global = pd.read_pickle("df_global_Feb.pkl")

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
parent_MC = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/dvcs/"
parent_bhMC = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/bh/"
parent_MC_bkg1g = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/bkg_1g/"
parent_MC_bkg2g = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/bkg_2g/"
parent_Gen = "/volatile/clas12/sangbaek/nov2021/convPkl_Gen/dvcs/inb/"
parent_bhGen = "/volatile/clas12/sangbaek/nov2021/convPkl_Gen/bh/inb/"



#dvcs inb 50 nA
print("reading dvcs inb 50 nA")
df_3987_corr = pd.read_pickle(parent_MC + "3987.pkl")
df_4124_corr = pd.read_pickle(parent_MC + "4124.pkl")
df_4139_corr = pd.read_pickle(parent_MC + "4139.pkl")
df_4181_corr = pd.read_pickle(parent_MC + "4181.pkl")
df_4182_corr = pd.read_pickle(parent_MC + "4182.pkl")
df_4397_corr = pd.read_pickle(parent_MC + "4397.pkl")
df_3987_Gen = pd.read_pickle(parent_Gen + "3987.pkl")
df_4124_Gen = pd.read_pickle(parent_Gen + "4124.pkl")
df_4139_Gen = pd.read_pickle(parent_Gen + "4139.pkl")
df_4181_Gen = pd.read_pickle(parent_Gen + "4181.pkl")
df_4182_Gen = pd.read_pickle(parent_Gen + "4182.pkl")
df_4397_Gen = pd.read_pickle(parent_Gen + "4397.pkl")

df_3987_Gen.loc[:, "event"] = df_3987_Gen.index
df_4124_Gen.loc[:, "event"] = df_4124_Gen.index
df_4139_Gen.loc[:, "event"] = df_4139_Gen.index
df_4181_Gen.loc[:, "event"] = df_4181_Gen.index
df_4182_Gen.loc[:, "event"] = df_4182_Gen.index
df_4397_Gen.loc[:, "event"] = df_4397_Gen.index

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
print("reading pi0->dvcs inb 50 nA")
df_4076_1g_corr = pd.read_pickle(parent_MC_bkg1g + "4076.pkl")
df_4202_1g_corr = pd.read_pickle(parent_MC_bkg1g + "4202.pkl")
df_4209_1g_corr = pd.read_pickle(parent_MC_bkg1g + "4209.pkl")
#bkg2g inb 50 nA
print("reading pi0->pi0 inb 50 nA")
df_4076_2g_corr = pd.read_pickle(parent_MC_bkg2g + "4076.pkl")
df_4202_2g_corr = pd.read_pickle(parent_MC_bkg2g + "4202.pkl")
df_4209_2g_corr = pd.read_pickle(parent_MC_bkg2g + "4209.pkl")

#bkg1g inb 55 nA
print("reading pi0->dvcs inb 50 nA")
df_4212_1g_corr = pd.read_pickle(parent_MC_bkg1g + "4212.pkl")
#bkg2g inb 55 nA
print("reading pi0->pi0 inb 50 nA")
df_4212_2g_corr = pd.read_pickle(parent_MC_bkg2g + "4212.pkl")

#bkg1g inb 45 nA
print("reading pi0->dvcs inb 50 nA")
df_4217_1g_corr = pd.read_pickle(parent_MC_bkg1g + "4217.pkl")
#bkg2g inb 45 nA
print("reading pi0->pi0 inb 50 nA")
df_4217_2g_corr = pd.read_pickle(parent_MC_bkg2g + "4217.pkl")

#bkg1g inb 0 nA
print("reading pi0->dvcs inb 50 nA")
df_4231_1g_corr = pd.read_pickle(parent_MC_bkg1g + "4231.pkl")
#bkg2g inb 0 nA
print("reading pi0->pi0 inb 50 nA")
df_4231_2g_corr = pd.read_pickle(parent_MC_bkg2g + "4231.pkl")

#Exp
print("reading pi0 exp")
# epgExpInb = pd.read_pickle(parent_exp + "dvcs.pkl")
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
print("4397")
df_4397_corr = numberingDF(df_4397_corr)
df_4397_Gen = df_4397_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
df_4397_Gen = numberingDF(df_4397_Gen)
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

print("4212")
df_4212_1g_corr = numberingDF(df_4212_1g_corr)
df_4212_2g_corr = numberingDF(df_4212_2g_corr)
print("4217")
df_4217_1g_corr = numberingDF(df_4217_1g_corr)
df_4217_2g_corr = numberingDF(df_4217_2g_corr)
print("4231")
df_4231_1g_corr = numberingDF(df_4231_1g_corr)
df_4231_2g_corr = numberingDF(df_4231_2g_corr)

# print("exp dvcs")
# epgExpInb = numberingDF(epgExpInb)
print("exp pi0")
pi0ExpInb = numberingDF(pi0ExpInb)

print("done with numbering")


dvcsGenInb50nA = pd.concat([df_3987_Gen, df_4124_Gen, df_4139_Gen, df_4181_Gen, df_4182_Gen, df_4397_Gen])
dvcsGenInb55nA = df_4186_Gen
dvcsGenInb45nA = df_4188_Gen
dvcsGenInb0nA = df_4192_Gen
bhGenInb50nA = df_4238_Gen

dvcsSimInb50nA = pd.concat([df_3987_corr, df_4124_corr, df_4139_corr, df_4181_corr, df_4182_corr, df_4397_corr])
dvcsSimInb55nA = df_4186_corr
dvcsSimInb45nA = df_4188_corr
dvcsSimInb0nA = df_4192_corr
bhSimInb50nA = df_4238_corr

bkgSimInb50nA = pd.concat([df_4076_1g_corr, df_4202_1g_corr, df_4209_1g_corr])
bkgSimInb55nA = df_4212_1g_corr
bkgSimInb45nA = df_4217_1g_corr
bkgSimInb0nA = df_4231_1g_corr
pi0SimInb50nA = pd.concat([df_4076_2g_corr, df_4202_2g_corr, df_4209_2g_corr])
pi0SimInb55nA = df_4212_2g_corr
pi0SimInb45nA = df_4217_2g_corr
pi0SimInb0nA = df_4231_2g_corr

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
df_global = countGenDF(dvcsGenInb50nA, df_global, "dvcsGenInb50nA")
df_global = countDF(dvcsSimInb50nA, df_global, "dvcsSimInb50nA")

df_global = countGenDF(dvcsGenInb55nA, df_global, "dvcsGenInb55nA")
df_global = countDF(dvcsSimInb55nA, df_global, "dvcsSimInb55nA")

df_global = countGenDF(dvcsGenInb45nA, df_global, "dvcsGenInb45nA")
df_global = countDF(dvcsSimInb45nA, df_global, "dvcsSimInb45nA")

df_global = countGenDF(dvcsGenInb0nA, df_global, "dvcsGenInb0nA")
df_global = countDF(dvcsSimInb0nA, df_global, "dvcsSimInb0nA")

df_global = countGenDF(bhGenInb50nA, df_global, "bhGenInb50nA")
df_global = countDF(bhSimInb50nA, df_global, "bhSimInb50nA")

#fractional acc.
df_global = countGenDF(dvcsGenInb50nA.loc[ dvcsGenInb50nA.radMode == 1], df_global, "dvcsGenInb50nA_non")
df_global = countDF(dvcsSimInb50nA.loc[ dvcsSimInb50nA.radMode == 1], df_global, "dvcsSimInb50nA_non")

df_global = countGenDF(dvcsGenInb55nA.loc[ dvcsGenInb55nA.radMode == 1], df_global, "dvcsGenInb55nA_non")
df_global = countDF(dvcsSimInb55nA.loc[ dvcsSimInb55nA.radMode == 1], df_global, "dvcsSimInb55nA_non")

df_global = countGenDF(dvcsGenInb45nA.loc[ dvcsGenInb45nA.radMode == 1], df_global, "dvcsGenInb45nA_non")
df_global = countDF(dvcsSimInb45nA.loc[ dvcsSimInb45nA.radMode == 1], df_global, "dvcsSimInb45nA_non")

df_global = countGenDF(dvcsGenInb0nA.loc[ dvcsGenInb0nA.radMode == 1], df_global, "dvcsGenInb0nA_non")
df_global = countDF(dvcsSimInb0nA.loc[ dvcsSimInb0nA.radMode == 1], df_global, "dvcsSimInb0nA_non")

df_global = countGenDF(bhGenInb50nA.loc[ bhGenInb50nA.radMode == 1], df_global, "bhGenInb50nA_non")
df_global = countDF(bhSimInb50nA.loc[ bhSimInb50nA.radMode == 1], df_global, "bhSimInb50nA_non")

#fractional acc. s peak
df_global = countGenDF(dvcsGenInb50nA.loc[ dvcsGenInb50nA.radMode == 2], df_global, "dvcsGenInb50nA_sPeak")
df_global = countDF(dvcsSimInb50nA.loc[ dvcsSimInb50nA.radMode == 2], df_global, "dvcsSimInb50nA_sPeak")

df_global = countGenDF(dvcsGenInb55nA.loc[ dvcsGenInb55nA.radMode == 2], df_global, "dvcsGenInb55nA_sPeak")
df_global = countDF(dvcsSimInb55nA.loc[ dvcsSimInb55nA.radMode == 2], df_global, "dvcsSimInb55nA_sPeak")

df_global = countGenDF(dvcsGenInb45nA.loc[ dvcsGenInb45nA.radMode == 2], df_global, "dvcsGenInb45nA_sPeak")
df_global = countDF(dvcsSimInb45nA.loc[ dvcsSimInb45nA.radMode == 2], df_global, "dvcsSimInb45nA_sPeak")

df_global = countGenDF(dvcsGenInb0nA.loc[ dvcsGenInb0nA.radMode == 2], df_global, "dvcsGenInb0nA_sPeak")
df_global = countDF(dvcsSimInb0nA.loc[ dvcsSimInb0nA.radMode == 2], df_global, "dvcsSimInb0nA_sPeak")

df_global = countGenDF(bhGenInb50nA.loc[ bhGenInb50nA.radMode == 2], df_global, "bhGenInb50nA_sPeak")
df_global = countDF(bhSimInb50nA.loc[ bhSimInb50nA.radMode == 2], df_global, "bhSimInb50nA_sPeak")

#fractional acc. p peak
df_global = countGenDF(dvcsGenInb50nA.loc[ dvcsGenInb50nA.radMode == 3], df_global, "dvcsGenInb50nA_pPeak")
df_global = countDF(dvcsSimInb50nA.loc[ dvcsSimInb50nA.radMode == 3], df_global, "dvcsSimInb50nA_pPeak")

df_global = countGenDF(dvcsGenInb55nA.loc[ dvcsGenInb55nA.radMode == 3], df_global, "dvcsGenInb55nA_pPeak")
df_global = countDF(dvcsSimInb55nA.loc[ dvcsSimInb55nA.radMode == 3], df_global, "dvcsSimInb55nA_pPeak")

df_global = countGenDF(dvcsGenInb45nA.loc[ dvcsGenInb45nA.radMode == 3], df_global, "dvcsGenInb45nA_pPeak")
df_global = countDF(dvcsSimInb45nA.loc[ dvcsSimInb45nA.radMode == 3], df_global, "dvcsSimInb45nA_pPeak")

df_global = countGenDF(dvcsGenInb0nA.loc[ dvcsGenInb0nA.radMode == 3], df_global, "dvcsGenInb0nA_pPeak")
df_global = countDF(dvcsSimInb0nA.loc[ dvcsSimInb0nA.radMode == 3], df_global, "dvcsSimInb0nA_pPeak")

df_global = countGenDF(bhGenInb50nA.loc[ bhGenInb50nA.radMode == 3], df_global, "bhGenInb50nA_pPeak")
df_global = countDF(bhSimInb50nA.loc[ bhSimInb50nA.radMode == 3], df_global, "bhSimInb50nA_pPeak")

# df_global = countDF(epgExpInb, df_global, "epgExpInb")
df_global = countDF(pi0ExpInb, df_global, "pi0ExpInb")

df_global = countDF(bkgSimInb50nA, df_global, "bkgSimInb50nA")
df_global = countDF(pi0SimInb50nA, df_global, "pi0SimInb50nA")

df_global = countDF(bkgSimInb45nA, df_global, "bkgSimInb45nA")
df_global = countDF(pi0SimInb45nA, df_global, "pi0SimInb45nA")

df_global = countDF(bkgSimInb55nA, df_global, "bkgSimInb55nA")
df_global = countDF(pi0SimInb55nA, df_global, "pi0SimInb55nA")

df_global = countDF(bkgSimInb0nA, df_global, "bkgSimInb0nA")
df_global = countDF(pi0SimInb0nA, df_global, "pi0SimInb0nA")

df_global.to_pickle("df_global_inb.pkl")