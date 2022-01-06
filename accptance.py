import pandas as pd
import numpy as np

k= 2*M*(np.sqrt(0.8**2+M**2)-M)
x1 = 1.456/2/M/8.604
x2 = 1.456/2/M/3
x3 = 1-(4-M**2)/2/M/3
x4 = 2*(1-np.sqrt((k-M**4+4*M**2)/k))/(M**4/k - 4*M**2/k)
x5 = (-1+np.sqrt(1+4*M**2/k*(1+M/8.604/2)))/2/M**2*k
y1 = 1.456
y2 = 1.456
y3 = 6*M*x3
y4 = (4-M**2)*x4/(1-x4)
y5 = 2*M*8.604*x5

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

Q2bin_i = [y1, np.sqrt(y1*y3), y3, np.sqrt(y3*y4), y4, np.sqrt(y4*y5)]
Q2bin_f = [np.sqrt(y1*y3), y3, np.sqrt(y3*y4), y4, np.sqrt(y4*y5), y5]
xBbin_i = {1: [[x1, c1], c2, c3], 2: [[c1, c2], c2, c3, c4], 3: [[c2, c3], c3, c4, c5], 4: [[c3, c4], c4, c5], 5: [[c4, c5], c5], 6: [[c5, x5]]}
xBbin_f = {1: [c2, c3, [x2, d1]], 2: [c2, c3, c4, [d1, d2]], 3: [c3, c4, c5, [d2, d3]], 4: [c4, c5, [d3, d4]], 5: [c5, x5], 6: [x5]}
tbin_i = [10**(-1.1), 10**(-0.8), 10**(-0.6), 10**(-0.4)]
tbin_f = [10**(-0.8), 10**(-0.6), 10**(-0.4), 10**(-0.2)]
goodBins = [111, 112, 113, 114, 121, 122, 123, 124, 131, 132, 133, 134, 211,
       212, 213, 214, 221, 222, 223, 224, 231, 232, 233, 234, 242, 243,
       244, 311, 312, 313, 314, 321, 322, 323, 324, 332, 333, 334, 343,
       344, 411, 412, 413, 414, 422, 423, 424, 433, 434, 512, 513, 514,
       523, 524, 613, 614]
badBins = [241, 331, 341, 342, 421, 431, 432, 511, 521, 522, 611, 612]

df_global = pd.read_pickle("df_global2.pkl")

def numberingDF(total, Q2bin_i=Q2bin_i, Q2bin_f=Q2bin_f, xBbin_i=xBbin_i, xBbin_f=xBbin_f, tbin_i=tbin_i, tbin_f=tbin_f, goodBins=goodBins, badBins=badBins, df_global=df_global):
    df_allBins = {}
    Q2xBtphi = 0

    for Q2bin in range(0, len(Q2bin_i)):#Q2 bin
        for xBbin in range(0, len(xBbin_i[Q2bin+1])):
            for tbin in range(0, len(tbin_i)):
                local = total
                Q2_i = Q2bin_i[Q2bin]
                Q2_f = Q2bin_f[Q2bin]
                xB_i = xBbin_i[Q2bin+1][xBbin]
                xB_f = xBbin_f[Q2bin+1][xBbin]
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
                if (xBbin == len(xBbin_i[Q2bin+1])-1) & (Q2bin < 2):
                    local = local.loc[local.Q2>=2*M*3*local.xB]
                elif (xBbin == len(xBbin_i[Q2bin+1])-1) & (Q2bin < 4):
                    local = local.loc[local.Q2>=(4 - M*M)*local.xB/(1 - local.xB)]
                else:
                    local = local.loc[local.xB<xB_f]
                #cut by t
                if tbin == len(tbin_i)-1:
                    local = local.loc[(local.t2>=t_i) & (local.t2<=t_f)]
                else:
                    local = local.loc[(local.t2>=t_i) & (local.t2<t_f)]
                Q2xBtbin = 100*(Q2bin+1) + 10*(xBbin+1) + (tbin+1)

                if Q2xBtbin in badBins:
                    continue

                phibin_i = df_global.loc[(df_global.Q2xBt == Q2xBtbin), "phi_i"].to_numpy()
                phibin_f = df_global.loc[(df_global.Q2xBt == Q2xBtbin), "phi_f"].to_numpy()
                phibin_f[-1] = 360
                for phi_ind in range(0, len(phibin_i)):
                    local.loc[:, "xBbin"] = xBbin
                    local.loc[:, "Q2bin"] = Q2bin
                    local.loc[:, "tbin"] = tbin
                    local.loc[:, "phibin"] = phi_ind
                    local.loc[:, "Q2xBtbin"] = Q2xBtbin
                    local.loc[:, "Q2xBtphi"] = Q2xBtphi
                    df_allBins[Q2xBtphi] = local.loc[(local.phi2>=phibin_i[phi_ind])&(local.phi2<phibin_f[phi_ind])]
                    Q2xBtphi += 1

    total = pd.concat(df_allBins.values()).sort_values( by = 'event')
    return total

parent_Gen = "/volatile/clas12/sangbaek/nov2021/convPkl_Gen/dvcs/inb"
parent_bhGen = "/volatile/clas12/sangbaek/nov2021/convPkl_Gen/bh/inb/"

df_3987_Gen = pd.read_pickle(parent_Gen + "3987.pkl")
df_4124_Gen = pd.read_pickle(parent_Gen + "4124.pkl")
df_4139_Gen = pd.read_pickle(parent_Gen + "4139.pkl")
df_4181_Gen = pd.read_pickle(parent_Gen + "4181.pkl")
df_4182_Gen = pd.read_pickle(parent_Gen + "4182.pkl")

df_4238_Gen = pd.read_pickle(parent_bhGen + "4238.pkl")

df_3987_Gen.loc[:, "event"] = df_3987_Gen.index
df_4124_Gen.loc[:, "event"] = df_4124_Gen.index
df_4139_Gen.loc[:, "event"] = df_4139_Gen.index
df_4181_Gen.loc[:, "event"] = df_4181_Gen.index
df_4182_Gen.loc[:, "event"] = df_4182_Gen.index
df_4238_Gen.loc[:, "event"] = df_4238_Gen.index

print("3987")
df_3987_Gen = numberingDF(df_3987_Gen)
print("4124")
df_4124_Gen = numberingDF(df_4124_Gen)
print("4139")
df_4139_Gen = numberingDF(df_4139_Gen)
print("4181")
df_4181_Gen = numberingDF(df_4181_Gen)
print("4182")
df_4182_Gen = numberingDF(df_4182_Gen)
print("4238")
df_4238_Gen = numberingDF(df_4238_Gen)

dvcsGenInb = pd.concat([df_3987_Gen, df_4124_Gen, df_4139_Gen, df_4181_Gen, df_4182_Gen])
bhGenInb = df_4238_Gen

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

df_global = countGenDF(dvcsGenInb, df_global, "dvcsGenInb")
df_global = countGenDF(bhGenInb, df_global, "bhGenInb")
df_global.to_pickle("df_global3.pkl")