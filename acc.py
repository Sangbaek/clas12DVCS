import pandas as pd
import numpy as np
import argparse
from utils.const import *
from utils.physics import *

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

def countDF(total, df_global, colName = "new"):
    numbers1 = []
    numbers2 = []
    numbers3 = []
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
    for i in range(len(df_global)):
        if i%10==0:
            print(i)
        number = sum((total.Q2xBtphi == i))
        numbers.append(number)
    df_global.loc[:, colName] = numbers
    return df_global

# parent_exp = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/exp/"
# parent_MC = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/dvcs/"
# parent_bhMC = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/bh/"
# parent_MC_bkg1g = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/bkg_1g/"
# parent_MC_bkg2g = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/bkg_2g/"
# parent_Gen = "/volatile/clas12/sangbaek/nov2021/convPkl_Gen/dvcs/inb/"
# parent_bhGen = "/volatile/clas12/sangbaek/nov2021/convPkl_Gen/bh/inb/"



# #dvcs inb 50 nA
# print("reading dvcs inb 50 nA")
# df_3987_corr = pd.read_pickle(parent_MC + "3987.pkl")
# df_4124_corr = pd.read_pickle(parent_MC + "4124.pkl")
# df_4139_corr = pd.read_pickle(parent_MC + "4139.pkl")
# df_4181_corr = pd.read_pickle(parent_MC + "4181.pkl")
# df_4182_corr = pd.read_pickle(parent_MC + "4182.pkl")
# df_4397_corr = pd.read_pickle(parent_MC + "4397.pkl")
# df_3987_Gen = pd.read_pickle(parent_Gen + "3987.pkl")
# df_4124_Gen = pd.read_pickle(parent_Gen + "4124.pkl")
# df_4139_Gen = pd.read_pickle(parent_Gen + "4139.pkl")
# df_4181_Gen = pd.read_pickle(parent_Gen + "4181.pkl")
# df_4182_Gen = pd.read_pickle(parent_Gen + "4182.pkl")
# df_4397_Gen = pd.read_pickle(parent_Gen + "4397.pkl")

# df_3987_Gen.loc[:, "event"] = df_3987_Gen.index
# df_4124_Gen.loc[:, "event"] = df_4124_Gen.index
# df_4139_Gen.loc[:, "event"] = df_4139_Gen.index
# df_4181_Gen.loc[:, "event"] = df_4181_Gen.index
# df_4182_Gen.loc[:, "event"] = df_4182_Gen.index
# df_4397_Gen.loc[:, "event"] = df_4397_Gen.index

# #dvcs inb 55 nA
# print("reading dvcs inb 55 nA")
# df_4186_corr = pd.read_pickle(parent_MC + "4186.pkl")
# df_4186_Gen = pd.read_pickle(parent_Gen + "4186.pkl")

# df_4186_Gen.loc[:, "event"] = df_4186_Gen.index

# #dvcs inb 45 nA
# print("reading dvcs inb 45 nA")
# df_4188_corr = pd.read_pickle(parent_MC + "4188.pkl")
# df_4188_Gen = pd.read_pickle(parent_Gen + "4188.pkl")

# df_4188_Gen.loc[:, "event"] = df_4188_Gen.index

# #dvcs inb 0 nA
# print("reading dvcs inb 0 nA")
# df_4192_corr = pd.read_pickle(parent_MC + "4192.pkl")
# df_4192_Gen = pd.read_pickle(parent_Gen + "4192.pkl")

# df_4192_Gen.loc[:, "event"] = df_4192_Gen.index

# #bh inb 50 nA
# print("reading bh inb 50 nA")
# df_4238_corr = pd.read_pickle(parent_bhMC + "4238.pkl")
# df_4238_Gen = pd.read_pickle(parent_bhGen + "4238.pkl")

# df_4238_Gen.loc[:, "event"] = df_4238_Gen.index

# #bkg1g inb 50 nA
# print("reading pi0->dvcs inb 50 nA")
# df_4076_1g_corr = pd.read_pickle(parent_MC_bkg1g + "4076.pkl")
# df_4202_1g_corr = pd.read_pickle(parent_MC_bkg1g + "4202.pkl")
# df_4209_1g_corr = pd.read_pickle(parent_MC_bkg1g + "4209.pkl")
# #bkg2g inb 50 nA
# print("reading pi0->pi0 inb 50 nA")
# df_4076_2g_corr = pd.read_pickle(parent_MC_bkg2g + "4076.pkl")
# df_4202_2g_corr = pd.read_pickle(parent_MC_bkg2g + "4202.pkl")
# df_4209_2g_corr = pd.read_pickle(parent_MC_bkg2g + "4209.pkl")

# #bkg1g inb 55 nA
# print("reading pi0->dvcs inb 55 nA")
# df_4212_1g_corr = pd.read_pickle(parent_MC_bkg1g + "4212.pkl")
# #bkg2g inb 55 nA
# print("reading pi0->pi0 inb 55 nA")
# df_4212_2g_corr = pd.read_pickle(parent_MC_bkg2g + "4212.pkl")

# #bkg1g inb 45 nA
# print("reading pi0->dvcs inb 45 nA")
# df_4217_1g_corr = pd.read_pickle(parent_MC_bkg1g + "4217.pkl")
# #bkg2g inb 45 nA
# print("reading pi0->pi0 inb 45 nA")
# df_4217_2g_corr = pd.read_pickle(parent_MC_bkg2g + "4217.pkl")

# #bkg1g inb 0 nA
# print("reading pi0->dvcs inb 0 nA")
# df_4231_1g_corr = pd.read_pickle(parent_MC_bkg1g + "4231.pkl")
# #bkg2g inb 0 nA
# print("reading pi0->pi0 inb 0 nA")
# df_4231_2g_corr = pd.read_pickle(parent_MC_bkg2g + "4231.pkl")

# #Exp
# print("reading pi0 exp")
# # epgExpInb = pd.read_pickle(parent_exp + "dvcs.pkl")
# pi0ExpInb = pd.read_pickle(parent_exp + "pi0.pkl")

# print("done with reading inbending")

# print("3987")
# df_3987_corr = numberingDF(df_3987_corr)
# df_3987_Gen = df_3987_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
# df_3987_Gen = numberingDF(df_3987_Gen)
# print("4124")
# df_4124_corr = numberingDF(df_4124_corr)
# df_4124_Gen = df_4124_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
# df_4124_Gen = numberingDF(df_4124_Gen)
# print("4139")
# df_4139_corr = numberingDF(df_4139_corr)
# df_4139_Gen = df_4139_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
# df_4139_Gen = numberingDF(df_4139_Gen)
# print("4181")
# df_4181_corr = numberingDF(df_4181_corr)
# df_4181_Gen = df_4181_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
# df_4181_Gen = numberingDF(df_4181_Gen)
# print("4182")
# df_4182_corr = numberingDF(df_4182_corr)
# df_4182_Gen = df_4182_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
# df_4182_Gen = numberingDF(df_4182_Gen)
# print("4397")
# df_4397_corr = numberingDF(df_4397_corr)
# df_4397_Gen = df_4397_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
# df_4397_Gen = numberingDF(df_4397_Gen)
# print("4186")
# df_4186_corr = numberingDF(df_4186_corr)
# df_4186_Gen = df_4186_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
# df_4186_Gen = numberingDF(df_4186_Gen)
# print("4188")
# df_4188_corr = numberingDF(df_4188_corr)
# df_4188_Gen = df_4188_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
# df_4188_Gen = numberingDF(df_4188_Gen)
# print("4192")
# df_4192_corr = numberingDF(df_4192_corr)
# df_4192_Gen = df_4192_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
# df_4192_Gen = numberingDF(df_4192_Gen)
# print("4238")
# df_4238_corr = numberingDF(df_4238_corr)
# df_4238_Gen = df_4238_Gen.rename(columns ={"phi2": "phi1", "t2": "t1"})
# df_4238_Gen = numberingDF(df_4238_Gen)

# print("4076")
# df_4076_1g_corr = numberingDF(df_4076_1g_corr)
# df_4076_2g_corr = numberingDF(df_4076_2g_corr)
# print("4202")
# df_4202_1g_corr = numberingDF(df_4202_1g_corr)
# df_4202_2g_corr = numberingDF(df_4202_2g_corr)
# print("4209")
# df_4209_1g_corr = numberingDF(df_4209_1g_corr)
# df_4209_2g_corr = numberingDF(df_4209_2g_corr)

# print("4212")
# df_4212_1g_corr = numberingDF(df_4212_1g_corr)
# df_4212_2g_corr = numberingDF(df_4212_2g_corr)
# print("4217")
# df_4217_1g_corr = numberingDF(df_4217_1g_corr)
# df_4217_2g_corr = numberingDF(df_4217_2g_corr)
# print("4231")
# df_4231_1g_corr = numberingDF(df_4231_1g_corr)
# df_4231_2g_corr = numberingDF(df_4231_2g_corr)

# # print("exp dvcs")
# # epgExpInb = numberingDF(epgExpInb)
# print("exp pi0")
# pi0ExpInb = numberingDF(pi0ExpInb)

# print("done with numbering")


# dvcsGenInb50nA = pd.concat([df_3987_Gen, df_4124_Gen, df_4139_Gen, df_4181_Gen, df_4182_Gen, df_4397_Gen])
# dvcsGenInb55nA = df_4186_Gen
# dvcsGenInb45nA = df_4188_Gen
# dvcsGenInb0nA = df_4192_Gen
# bhGenInb50nA = df_4238_Gen

# dvcsSimInb50nA = pd.concat([df_3987_corr, df_4124_corr, df_4139_corr, df_4181_corr, df_4182_corr, df_4397_corr])
# dvcsSimInb55nA = df_4186_corr
# dvcsSimInb45nA = df_4188_corr
# dvcsSimInb0nA = df_4192_corr
# bhSimInb50nA = df_4238_corr

# bkgSimInb50nA = pd.concat([df_4076_1g_corr, df_4202_1g_corr, df_4209_1g_corr])
# bkgSimInb55nA = df_4212_1g_corr
# bkgSimInb45nA = df_4217_1g_corr
# bkgSimInb0nA = df_4231_1g_corr
# pi0SimInb50nA = pd.concat([df_4076_2g_corr, df_4202_2g_corr, df_4209_2g_corr])
# pi0SimInb55nA = df_4212_2g_corr
# pi0SimInb45nA = df_4217_2g_corr
# pi0SimInb0nA = df_4231_2g_corr



# #total acc.
# df_global = countGenDF(dvcsGenInb50nA, df_global, "dvcsGenInb50nA")
# df_global = countDF(dvcsSimInb50nA, df_global, "dvcsSimInb50nA")

# df_global = countGenDF(dvcsGenInb55nA, df_global, "dvcsGenInb55nA")
# df_global = countDF(dvcsSimInb55nA, df_global, "dvcsSimInb55nA")

# df_global = countGenDF(dvcsGenInb45nA, df_global, "dvcsGenInb45nA")
# df_global = countDF(dvcsSimInb45nA, df_global, "dvcsSimInb45nA")

# df_global = countGenDF(dvcsGenInb0nA, df_global, "dvcsGenInb0nA")
# df_global = countDF(dvcsSimInb0nA, df_global, "dvcsSimInb0nA")

# df_global = countGenDF(bhGenInb50nA, df_global, "bhGenInb50nA")
# df_global = countDF(bhSimInb50nA, df_global, "bhSimInb50nA")

# #fractional acc.
# df_global = countGenDF(dvcsGenInb50nA.loc[ dvcsGenInb50nA.radMode == 1], df_global, "dvcsGenInb50nA_non")
# df_global = countDF(dvcsSimInb50nA.loc[ dvcsSimInb50nA.radMode == 1], df_global, "dvcsSimInb50nA_non")

# df_global = countGenDF(dvcsGenInb55nA.loc[ dvcsGenInb55nA.radMode == 1], df_global, "dvcsGenInb55nA_non")
# df_global = countDF(dvcsSimInb55nA.loc[ dvcsSimInb55nA.radMode == 1], df_global, "dvcsSimInb55nA_non")

# df_global = countGenDF(dvcsGenInb45nA.loc[ dvcsGenInb45nA.radMode == 1], df_global, "dvcsGenInb45nA_non")
# df_global = countDF(dvcsSimInb45nA.loc[ dvcsSimInb45nA.radMode == 1], df_global, "dvcsSimInb45nA_non")

# df_global = countGenDF(dvcsGenInb0nA.loc[ dvcsGenInb0nA.radMode == 1], df_global, "dvcsGenInb0nA_non")
# df_global = countDF(dvcsSimInb0nA.loc[ dvcsSimInb0nA.radMode == 1], df_global, "dvcsSimInb0nA_non")

# df_global = countGenDF(bhGenInb50nA.loc[ bhGenInb50nA.radMode == 1], df_global, "bhGenInb50nA_non")
# df_global = countDF(bhSimInb50nA.loc[ bhSimInb50nA.radMode == 1], df_global, "bhSimInb50nA_non")

# #fractional acc. s peak
# df_global = countGenDF(dvcsGenInb50nA.loc[ dvcsGenInb50nA.radMode == 2], df_global, "dvcsGenInb50nA_sPeak")
# df_global = countDF(dvcsSimInb50nA.loc[ dvcsSimInb50nA.radMode == 2], df_global, "dvcsSimInb50nA_sPeak")

# df_global = countGenDF(dvcsGenInb55nA.loc[ dvcsGenInb55nA.radMode == 2], df_global, "dvcsGenInb55nA_sPeak")
# df_global = countDF(dvcsSimInb55nA.loc[ dvcsSimInb55nA.radMode == 2], df_global, "dvcsSimInb55nA_sPeak")

# df_global = countGenDF(dvcsGenInb45nA.loc[ dvcsGenInb45nA.radMode == 2], df_global, "dvcsGenInb45nA_sPeak")
# df_global = countDF(dvcsSimInb45nA.loc[ dvcsSimInb45nA.radMode == 2], df_global, "dvcsSimInb45nA_sPeak")

# df_global = countGenDF(dvcsGenInb0nA.loc[ dvcsGenInb0nA.radMode == 2], df_global, "dvcsGenInb0nA_sPeak")
# df_global = countDF(dvcsSimInb0nA.loc[ dvcsSimInb0nA.radMode == 2], df_global, "dvcsSimInb0nA_sPeak")

# df_global = countGenDF(bhGenInb50nA.loc[ bhGenInb50nA.radMode == 2], df_global, "bhGenInb50nA_sPeak")
# df_global = countDF(bhSimInb50nA.loc[ bhSimInb50nA.radMode == 2], df_global, "bhSimInb50nA_sPeak")

# #fractional acc. p peak
# df_global = countGenDF(dvcsGenInb50nA.loc[ dvcsGenInb50nA.radMode == 3], df_global, "dvcsGenInb50nA_pPeak")
# df_global = countDF(dvcsSimInb50nA.loc[ dvcsSimInb50nA.radMode == 3], df_global, "dvcsSimInb50nA_pPeak")

# df_global = countGenDF(dvcsGenInb55nA.loc[ dvcsGenInb55nA.radMode == 3], df_global, "dvcsGenInb55nA_pPeak")
# df_global = countDF(dvcsSimInb55nA.loc[ dvcsSimInb55nA.radMode == 3], df_global, "dvcsSimInb55nA_pPeak")

# df_global = countGenDF(dvcsGenInb45nA.loc[ dvcsGenInb45nA.radMode == 3], df_global, "dvcsGenInb45nA_pPeak")
# df_global = countDF(dvcsSimInb45nA.loc[ dvcsSimInb45nA.radMode == 3], df_global, "dvcsSimInb45nA_pPeak")

# df_global = countGenDF(dvcsGenInb0nA.loc[ dvcsGenInb0nA.radMode == 3], df_global, "dvcsGenInb0nA_pPeak")
# df_global = countDF(dvcsSimInb0nA.loc[ dvcsSimInb0nA.radMode == 3], df_global, "dvcsSimInb0nA_pPeak")

# df_global = countGenDF(bhGenInb50nA.loc[ bhGenInb50nA.radMode == 3], df_global, "bhGenInb50nA_pPeak")
# df_global = countDF(bhSimInb50nA.loc[ bhSimInb50nA.radMode == 3], df_global, "bhSimInb50nA_pPeak")

# # df_global = countDF(epgExpInb, df_global, "epgExpInb")
# df_global = countDF(pi0ExpInb, df_global, "pi0ExpInb")

# df_global = countDF(bkgSimInb50nA, df_global, "bkgSimInb50nA")
# df_global = countDF(pi0SimInb50nA, df_global, "pi0SimInb50nA")

# df_global = countDF(bkgSimInb45nA, df_global, "bkgSimInb45nA")
# df_global = countDF(pi0SimInb45nA, df_global, "pi0SimInb45nA")

# df_global = countDF(bkgSimInb55nA, df_global, "bkgSimInb55nA")
# df_global = countDF(pi0SimInb55nA, df_global, "pi0SimInb55nA")

# df_global = countDF(bkgSimInb0nA, df_global, "bkgSimInb0nA")
# df_global = countDF(pi0SimInb0nA, df_global, "pi0SimInb0nA")

# df_global.to_pickle("df_global_inb.pkl")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-if","--ifname", help="a single pickle to be numbered", default="/Users/sangbaek/Dropbox (MIT)/data/project/merged_9628_files.root")
    parser.add_argument("-of","--ofname", help="output with numbered", default="/Users/sangbaek/Dropbox (MIT)/data/project/merged_9628_files.root")
    parser.add_argument("-n","--numbering", help="numbering only", action = "store_true")
    parser.add_argument("-g","--gen", help="gen or rec", action = "store_true")
    parser.add_argument("-i","--inglobal", help="a single pickle file name as an input", default="df_global_Feb.pkl")
    parser.add_argument("-o","--outglobal", help="a single pickle file name as an output", default="df_global_Feb_out.pkl")
    parser.add_argument("-N","--nonrad", help="to make fractional", action="store_true")
    parser.add_argument("-S","--speak", help="to make fractional", action="store_true")
    parser.add_argument("-P","--ppeak", help="to make fractional", action="store_true")
    parser.add_argument("-c","--colName", help="columnName", default="dvcsSimInb50nA")
    
    args = parser.parse_args()

    if args.numbering:
        print("reading..")
        df = pd.read_pickle(args.ifname)
        print("done with reading..")
        if args.gen:
            df.loc[:, "event"] = df.index
            df = df.rename(columns ={"phi2": "phi1", "t2": "t1"})
        print("numbering..")
        df = numberingDF(df)
        print("saving..")
        df.to_pickle(args.ofname)
    else:
        print("reading..")
        df = pd.read_pickle(args.ifname)
        print("reading seed..")
        df_global = pd.read_pickle(args.inglobal)
        if args.gen:
            print("count Gen..")
            if args.nonrad:
                print("select non rad events only")
                df = df.loc[df.radMode == 1, :]
            if args.speak:
                print("select s-peak events only")
                df = df.loc[df.radMode == 2, :]
            if args.ppeak:
                print("select p-peak events only")
                df = df.loc[df.radMode == 3, :]
            df_global = countGenDF(df, df_global, args.colName)
            print("done with counting..")
            df_global.to_pickle(args.outglobal)
        else:
            if args.nonrad:
                print("select non rad events only")
                df = df.loc[df.radMode == 1, :]
            if args.speak:
                print("select s-peak events only")
                df = df.loc[df.radMode == 2, :]
            if args.ppeak:
                print("select p-peak events only")
                df = df.loc[df.radMode == 3, :]
            print("count Rec..")
            df_global = countDF(df, df_global, args.colName)
            print("done with counting..")
        df_global.to_pickle(args.outglobal)