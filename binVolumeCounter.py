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

local502 = pd.read_pickle("local502_outb.pkl")

def TruebinVol(Q2bin, xBbin, tbin, phibin, df, N1=10, N2=10, N3=10, N4=10):
    
    count = 0 
    
    Q2_i = Q2bin_i[Q2bin]
    Q2_f = Q2bin_f[Q2bin]
    xB_i = xBbin_i[Q2bin][xBbin]
    xB_f = xBbin_f[Q2bin][xBbin]
    t_i = tbin_i[tbin]
    t_f = tbin_f[tbin]
    phi_i = phibin_i[phibin]
    phi_f = phibin_f[phibin]
    
    local = df
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

parent_MC = "/volatile/clas12/sangbaek/nov2021/convPkl/dvcs/"
parent_bhMC = "/volatile/clas12/sangbaek/nov2021/convPkl/bh/"


#dvcs inb 50 nA
print("reading dvcs inb 50 nA")
df_3987_corr = pd.read_pickle(parent_MC + "3987.pkl")
df_4124_corr = pd.read_pickle(parent_MC + "4124.pkl")
df_4139_corr = pd.read_pickle(parent_MC + "4139.pkl")
df_4181_corr = pd.read_pickle(parent_MC + "4181.pkl")
df_4182_corr = pd.read_pickle(parent_MC + "4182.pkl")

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

dvcsSimInb50nA = pd.concat([df_3987_corr, df_4124_corr, df_4139_corr, df_4181_corr, df_4182_corr])
dvcsSimInb55nA = df_4186_corr
dvcsSimInb45nA = df_4188_corr
dvcsSimInb0nA = df_4192_corr
bhSimInb50nA = df_4238_corr

dvcsBHSimInb = pd.concat([dvcsSimInb50nA, dvcsSimInb55nA, dvcsSimInb45nA, dvcsSimInb0nA, bhSimInb50nA])

print("counting inb bin volumes")
TrueVolInb = []
for i in range(len(phibin_i)):
    print(i)
    TrueVolInb.append(TruebinVol(5, 0, 2, i, dvcsBHSimInb, 10, 10, 10, 10))
local502.loc[:, "binVolInb"] = TrueVolInb
local502.to_pickle("truebinVol_inb.pkl")

parent_MC = "/volatile/clas12/sangbaek/nov2021/convPkl_outb/dvcs/"
parent_bhMC = "/volatile/clas12/sangbaek/nov2021/convPkl_outb/bh/"

df_4240_corr = pd.read_pickle(parent_MC + "4240.pkl")
df_4250_corr = pd.read_pickle(parent_MC + "4250.pkl")
df_4251_corr = pd.read_pickle(parent_MC + "4251.pkl")
df_4252_corr = pd.read_pickle(parent_MC + "4252.pkl")
df_4255_corr = pd.read_pickle(parent_MC + "4255.pkl")
df_4263_corr = pd.read_pickle(parent_MC + "4263.pkl")
df_4262_corr = pd.read_pickle(parent_MC + "4262.pkl")
df_4266_corr = pd.read_pickle(parent_MC + "4266.pkl")
df_4249_Gen = pd.read_pickle(parent_bhGen + "4249.pkl")

dvcsSimOutb50nA = pd.concat([df_4240_corr, df_4250_corr, df_4251_corr, df_4252_corr, df_4255_corr])
dvcsSimOutb40nA = df_4263_corr
dvcsSimOutb0nA = df_4262_corr
dvcsSimOutb40nAT = df_4266_corr
bhSimOutb50nA = df_4249_corr

dvcsBHSimOutb = pd.concat([dvcsSimOutb50nA, dvcsSimOutb40nA, dvcsSimOutb0nA, dvcsSimOutb40nAT, bhSimOutb50nA])

print("counting outb bin volumes")
TrueVolOutb = []
for i in range(len(phibin_i)):
    print(i)
    TrueVolOutb.append(TruebinVol(5, 0, 2, i, dvcsBHSimOutb, 10, 10, 10, 10))
local502.loc[:, "binVolInb"] = TrueVolOutb
local502.to_pickle("truebinVol_outb.pkl")
