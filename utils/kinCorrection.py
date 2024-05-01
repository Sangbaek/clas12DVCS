'''
independent script to perform the energy loss corection
'''
from utils.const import *
from utils.physics import *
from copy import copy

def electronMomentumCorrection(pol, df_electronRec):
    # https://clasweb.jlab.org/wiki/index.php/CLAS12_Momentum_Corrections#tab=Correction_Code
    df_electronRec = copy(df_electronRec)
    df_electronRec.loc[:, "dp"] = 0

    eleCorr = df_electronRec.loc[:, ["Ep", "Ephi", "Esector"]]
    eleCorr.loc[ ( ((eleCorr.Esector == 4) | (eleCorr.Esector == 3)) & (eleCorr.Ephi < 0)) | ((df_electronRec.Esector > 4) & (eleCorr.Ephi < 90)), "Ephi"] += 360
    eleCorr.loc[:, "Ephi"] = eleCorr.Ephi - (eleCorr.Esector - 1)*60
    if pol == "inbending":
        eleCorr.loc[eleCorr.Esector == 1, "dp"] = ((-4.3303e-06)*eleCorr.loc[eleCorr.Esector==1, "Ephi"]*eleCorr.loc[eleCorr.Esector==1, "Ephi"] + (1.1006e-04)* eleCorr.loc[eleCorr.Esector==1, "Ephi"] + (-5.7235e-04))*eleCorr.loc[eleCorr.Esector==1, "Ep"]*eleCorr.loc[eleCorr.Esector==1, "Ep"] + ((3.2555e-05)* eleCorr.loc[eleCorr.Esector==1, "Ephi"]*eleCorr.loc[eleCorr.Esector==1, "Ephi"] + (-0.0014559)* eleCorr.loc[eleCorr.Esector==1, "Ephi"] + (0.0014878))*  eleCorr.loc[eleCorr.Esector==1, "Ep"] + ((-1.9577e-05)*eleCorr.loc[eleCorr.Esector==1, "Ephi"]*eleCorr.loc[eleCorr.Esector==1, "Ephi"] + (0.0017996)*  eleCorr.loc[eleCorr.Esector==1, "Ephi"] + (0.025963));
        eleCorr.loc[eleCorr.Esector == 2, "dp"] = ((-9.8045e-07)*eleCorr.loc[eleCorr.Esector==2, "Ephi"]*eleCorr.loc[eleCorr.Esector==2, "Ephi"] + (6.7395e-05)* eleCorr.loc[eleCorr.Esector==2, "Ephi"] + (-4.6757e-05))*eleCorr.loc[eleCorr.Esector==2, "Ep"]*eleCorr.loc[eleCorr.Esector==2, "Ep"] + ((-1.4958e-05)*eleCorr.loc[eleCorr.Esector==2, "Ephi"]*eleCorr.loc[eleCorr.Esector==2, "Ephi"] + (-0.0011191)* eleCorr.loc[eleCorr.Esector==2, "Ephi"] + (-0.0025143))* eleCorr.loc[eleCorr.Esector==2, "Ep"] + ((1.2699e-04)* eleCorr.loc[eleCorr.Esector==2, "Ephi"]*eleCorr.loc[eleCorr.Esector==2, "Ephi"] + (0.0033121)*  eleCorr.loc[eleCorr.Esector==2, "Ephi"] + (0.020819));
        eleCorr.loc[eleCorr.Esector == 3, "dp"] = ((-5.9459e-07)*eleCorr.loc[eleCorr.Esector==3, "Ephi"]*eleCorr.loc[eleCorr.Esector==3, "Ephi"] + (-2.8289e-05)*eleCorr.loc[eleCorr.Esector==3, "Ephi"] + (-4.3541e-04))*eleCorr.loc[eleCorr.Esector==3, "Ep"]*eleCorr.loc[eleCorr.Esector==3, "Ep"] + ((-1.5025e-05)*eleCorr.loc[eleCorr.Esector==3, "Ephi"]*eleCorr.loc[eleCorr.Esector==3, "Ephi"] + (5.7730e-04)* eleCorr.loc[eleCorr.Esector==3, "Ephi"] + (-0.0077582))* eleCorr.loc[eleCorr.Esector==3, "Ep"] + ((7.3348e-05)* eleCorr.loc[eleCorr.Esector==3, "Ephi"]*eleCorr.loc[eleCorr.Esector==3, "Ephi"] + (-0.001102)*  eleCorr.loc[eleCorr.Esector==3, "Ephi"] + (0.057052));
        eleCorr.loc[eleCorr.Esector == 4, "dp"] = ((-2.2714e-06)*eleCorr.loc[eleCorr.Esector==4, "Ephi"]*eleCorr.loc[eleCorr.Esector==4, "Ephi"] + (-3.0360e-05)*eleCorr.loc[eleCorr.Esector==4, "Ephi"] + (-8.9322e-04))*eleCorr.loc[eleCorr.Esector==4, "Ep"]*eleCorr.loc[eleCorr.Esector==4, "Ep"] + ((2.9737e-05)* eleCorr.loc[eleCorr.Esector==4, "Ephi"]*eleCorr.loc[eleCorr.Esector==4, "Ephi"] + (5.1142e-04)* eleCorr.loc[eleCorr.Esector==4, "Ephi"] + (0.0045641))*  eleCorr.loc[eleCorr.Esector==4, "Ep"] + ((-1.0582e-04)*eleCorr.loc[eleCorr.Esector==4, "Ephi"]*eleCorr.loc[eleCorr.Esector==4, "Ephi"] + (-5.6852e-04)*eleCorr.loc[eleCorr.Esector==4, "Ephi"] + (0.027506));
        eleCorr.loc[eleCorr.Esector == 5, "dp"] = ((-1.1490e-06)*eleCorr.loc[eleCorr.Esector==5, "Ephi"]*eleCorr.loc[eleCorr.Esector==5, "Ephi"] + (-6.2147e-06)*eleCorr.loc[eleCorr.Esector==5, "Ephi"] + (-4.7235e-04))*eleCorr.loc[eleCorr.Esector==5, "Ep"]*eleCorr.loc[eleCorr.Esector==5, "Ep"] + ((3.7039e-06)* eleCorr.loc[eleCorr.Esector==5, "Ephi"]*eleCorr.loc[eleCorr.Esector==5, "Ephi"] + (-1.5943e-04)*eleCorr.loc[eleCorr.Esector==5, "Ephi"] + (-8.5238e-04))*eleCorr.loc[eleCorr.Esector==5, "Ep"] + ((4.4069e-05)* eleCorr.loc[eleCorr.Esector==5, "Ephi"]*eleCorr.loc[eleCorr.Esector==5, "Ephi"] + (0.0014152)*  eleCorr.loc[eleCorr.Esector==5, "Ephi"] + (0.031933));
        eleCorr.loc[eleCorr.Esector == 6, "dp"] = ((1.1076e-06) *eleCorr.loc[eleCorr.Esector==6, "Ephi"]*eleCorr.loc[eleCorr.Esector==6, "Ephi"] + (4.0156e-05)* eleCorr.loc[eleCorr.Esector==6, "Ephi"] + (-1.6341e-04))*eleCorr.loc[eleCorr.Esector==6, "Ep"]*eleCorr.loc[eleCorr.Esector==6, "Ep"] + ((-2.8613e-05)*eleCorr.loc[eleCorr.Esector==6, "Ephi"]*eleCorr.loc[eleCorr.Esector==6, "Ephi"] + (-5.1861e-04)*eleCorr.loc[eleCorr.Esector==6, "Ephi"] + (-0.0056437))* eleCorr.loc[eleCorr.Esector==6, "Ep"] + ((1.2419e-04)* eleCorr.loc[eleCorr.Esector==6, "Ephi"]*eleCorr.loc[eleCorr.Esector==6, "Ephi"] + (4.9084e-04)* eleCorr.loc[eleCorr.Esector==6, "Ephi"] + (0.049976));
    else:
        eleCorr.loc[eleCorr.Esector==1, "dp"] =     ((1.3189e-06)*eleCorr.loc[eleCorr.Esector==1, "Ephi"]*eleCorr.loc[eleCorr.Esector==1, "Ephi"] +  (4.26057e-05)*eleCorr.loc[eleCorr.Esector==1, "Ephi"] +  (-0.002322628))*eleCorr.loc[eleCorr.Esector==1, "Ep"]*eleCorr.loc[eleCorr.Esector==1, "Ep"] +  ((-1.1409e-05)*eleCorr.loc[eleCorr.Esector==1, "Ephi"]*eleCorr.loc[eleCorr.Esector==1, "Ephi"] +    (2.2188e-05)*eleCorr.loc[eleCorr.Esector==1, "Ephi"] + (0.02878927))*eleCorr.loc[eleCorr.Esector==1, "Ep"] +   ((2.4950e-05)*eleCorr.loc[eleCorr.Esector==1, "Ephi"]*eleCorr.loc[eleCorr.Esector==1, "Ephi"] +   (1.6170e-06)*eleCorr.loc[eleCorr.Esector==1, "Ephi"] + (-0.061816275));
        eleCorr.loc[eleCorr.Esector==2, "dp"] =    ((-2.9240e-07)*eleCorr.loc[eleCorr.Esector==2, "Ephi"]*eleCorr.loc[eleCorr.Esector==2, "Ephi"] +   (3.2448e-07)*eleCorr.loc[eleCorr.Esector==2, "Ephi"] +  (-0.001848308))*eleCorr.loc[eleCorr.Esector==2, "Ep"]*eleCorr.loc[eleCorr.Esector==2, "Ep"] +   ((4.4500e-07)*eleCorr.loc[eleCorr.Esector==2, "Ephi"]*eleCorr.loc[eleCorr.Esector==2, "Ephi"] +   (4.76324e-04)*eleCorr.loc[eleCorr.Esector==2, "Ephi"] + (0.02219469))*eleCorr.loc[eleCorr.Esector==2, "Ep"] +   ((6.9220e-06)*eleCorr.loc[eleCorr.Esector==2, "Ephi"]*eleCorr.loc[eleCorr.Esector==2, "Ephi"] +  (-0.00153517)*eleCorr.loc[eleCorr.Esector==2, "Ephi"] + (-0.0479058));
        eleCorr.loc[eleCorr.Esector==3, "dp"] =    ((2.71911e-06)*eleCorr.loc[eleCorr.Esector==3, "Ephi"]*eleCorr.loc[eleCorr.Esector==3, "Ephi"] + (1.657148e-05)*eleCorr.loc[eleCorr.Esector==3, "Ephi"] +  (-0.001822211))*eleCorr.loc[eleCorr.Esector==3, "Ep"]*eleCorr.loc[eleCorr.Esector==3, "Ep"] + ((-4.96814e-05)*eleCorr.loc[eleCorr.Esector==3, "Ephi"]*eleCorr.loc[eleCorr.Esector==3, "Ephi"] + (-3.761117e-04)*eleCorr.loc[eleCorr.Esector==3, "Ephi"] + (0.02564148))*eleCorr.loc[eleCorr.Esector==3, "Ep"] +  ((1.97748e-04)*eleCorr.loc[eleCorr.Esector==3, "Ephi"]*eleCorr.loc[eleCorr.Esector==3, "Ephi"] +  (9.58259e-04)*eleCorr.loc[eleCorr.Esector==3, "Ephi"] + (-0.05818292));
        eleCorr.loc[eleCorr.Esector==4, "dp"] =    ((1.90966e-06)*eleCorr.loc[eleCorr.Esector==4, "Ephi"]*eleCorr.loc[eleCorr.Esector==4, "Ephi"] +  (-2.4761e-05)*eleCorr.loc[eleCorr.Esector==4, "Ephi"] +   (-0.00231562))*eleCorr.loc[eleCorr.Esector==4, "Ep"]*eleCorr.loc[eleCorr.Esector==4, "Ep"] +  ((-2.3927e-05)*eleCorr.loc[eleCorr.Esector==4, "Ephi"]*eleCorr.loc[eleCorr.Esector==4, "Ephi"] +   (2.25262e-04)*eleCorr.loc[eleCorr.Esector==4, "Ephi"] +  (0.0291831))*eleCorr.loc[eleCorr.Esector==4, "Ep"] +   ((8.0515e-05)*eleCorr.loc[eleCorr.Esector==4, "Ephi"]*eleCorr.loc[eleCorr.Esector==4, "Ephi"] + (-6.42098e-04)*eleCorr.loc[eleCorr.Esector==4, "Ephi"] + (-0.06159197));
        eleCorr.loc[eleCorr.Esector==5, "dp"] = ((-3.6760323e-06)*eleCorr.loc[eleCorr.Esector==5, "Ephi"]*eleCorr.loc[eleCorr.Esector==5, "Ephi"] +  (4.04398e-05)*eleCorr.loc[eleCorr.Esector==5, "Ephi"] + (-0.0021967515))*eleCorr.loc[eleCorr.Esector==5, "Ep"]*eleCorr.loc[eleCorr.Esector==5, "Ep"] +  ((4.90857e-05)*eleCorr.loc[eleCorr.Esector==5, "Ephi"]*eleCorr.loc[eleCorr.Esector==5, "Ephi"] +  (-4.37437e-04)*eleCorr.loc[eleCorr.Esector==5, "Ephi"] + (0.02494339))*eleCorr.loc[eleCorr.Esector==5, "Ep"] + ((-1.08257e-04)*eleCorr.loc[eleCorr.Esector==5, "Ephi"]*eleCorr.loc[eleCorr.Esector==5, "Ephi"] +   (0.00146111)*eleCorr.loc[eleCorr.Esector==5, "Ephi"] + (-0.0648485));
        eleCorr.loc[eleCorr.Esector==6, "dp"] =    ((-6.2488e-08)*eleCorr.loc[eleCorr.Esector==6, "Ephi"]*eleCorr.loc[eleCorr.Esector==6, "Ephi"] +  (2.23173e-05)*eleCorr.loc[eleCorr.Esector==6, "Ephi"] +   (-0.00227522))*eleCorr.loc[eleCorr.Esector==6, "Ep"]*eleCorr.loc[eleCorr.Esector==6, "Ep"] +   ((1.8372e-05)*eleCorr.loc[eleCorr.Esector==6, "Ephi"]*eleCorr.loc[eleCorr.Esector==6, "Ephi"] +   (-7.5227e-05)*eleCorr.loc[eleCorr.Esector==6, "Ephi"] +   (0.032636))*eleCorr.loc[eleCorr.Esector==6, "Ep"] +  ((-6.6566e-05)*eleCorr.loc[eleCorr.Esector==6, "Ephi"]*eleCorr.loc[eleCorr.Esector==6, "Ephi"] +  (-2.4450e-04)*eleCorr.loc[eleCorr.Esector==6, "Ephi"] + (-0.072293));


    df_electronRec.loc[:, "Epx"] = (1 + eleCorr.dp/eleCorr.Ep) * df_electronRec.loc[:, "Epx"]
    df_electronRec.loc[:, "Epy"] = (1 + eleCorr.dp/eleCorr.Ep) * df_electronRec.loc[:, "Epy"]
    df_electronRec.loc[:, "Epz"] = (1 + eleCorr.dp/eleCorr.Ep) * df_electronRec.loc[:, "Epz"]
    df_electronRec.loc[:, "Ep"]  = df_electronRec.Ep + eleCorr.dp
    df_electronRec.loc[:, "Ee"]  = np.sqrt(df_electronRec.Ep**2 + me**2)
    return df_electronRec

def electronMomentumSmearing(df_electronRec):
    #p.49 of inclusive note
    df_electronRec = copy(df_electronRec)
    sigma_theta = 0.00387 - 0.00019 * df_electronRec.Etheta + 1.3021e-5 * df_electronRec.Etheta * df_electronRec.Etheta
    eleCorr = df_electronRec.loc[:, ["Ep", "Ephi", "Esector"]]
    eleCorr.loc[:, "dp"] = df_electronRec.Ep * sigma_theta * np.random.normal(0, 1, len(df_electronRec))

    df_electronRec.loc[:, "Epx"] = (1 + eleCorr.dp/eleCorr.Ep) * df_electronRec.loc[:, "Epx"]
    df_electronRec.loc[:, "Epy"] = (1 + eleCorr.dp/eleCorr.Ep) * df_electronRec.loc[:, "Epy"]
    df_electronRec.loc[:, "Epz"] = (1 + eleCorr.dp/eleCorr.Ep) * df_electronRec.loc[:, "Epz"]
    df_electronRec.loc[:, "Ep"]  = df_electronRec.Ep + eleCorr.dp
    df_electronRec.loc[:, "Ee"]  = np.sqrt(df_electronRec.Ep**2 + me**2)
    return df_electronRec

def protonEnergyLossCorr(pol, df_protonRec):
    '''
    a simple function for the energy loss correction.
    '''
    df_protonRec = copy(df_protonRec)
    pro = [df_protonRec['Ppx'], df_protonRec['Ppy'], df_protonRec['Ppz']]
    df_protonRec.loc[:, "Pp"] = mag(pro)
    df_protonRec.loc[:, "Pe"] = getEnergy(pro, M)
    df_protonRec.loc[:, "Ptheta"] = getTheta(pro)
    df_protonRec.loc[:, "Pphi"] = getPhi(pro)

    df_protonRec.loc[:, "PpOrig"] = mag(pro)
    df_protonRec.loc[:, "PeOrig"] = getEnergy(pro, M)
    df_protonRec.loc[:, "PthetaOrig"] = getTheta(pro)
    df_protonRec.loc[:, "PphiOrig"] = getPhi(pro)

    df_protonRecFD = df_protonRec.loc[df_protonRec.Psector<7, :]
    df_protonRecCD = df_protonRec.loc[(df_protonRec.Psector>7) & (df_protonRec.PthetaOrig<75), :]
    df_protonRecOthers = df_protonRec.loc[ ((df_protonRec.Psector>7) & (df_protonRec.PthetaOrig>=75)) | ((df_protonRec.Psector<7) & (df_protonRec.PpOrig<0.3)), :]

    #two band criterion
    def corr(x, t):
        x0, x1, x2, x3 = x
        return x0 + x1*np.power(t-np.ones(len(t))*0.3, x3)

    df_protonRecFD = df_protonRecFD.loc[df_protonRec.PpOrig > 0.3, :]
    df_protonRecFD.loc[:, "PDc1theta"] = getTheta([df_protonRecFD.PDc1Hitx, df_protonRecFD.PDc1Hity, df_protonRecFD.PDc1Hitz])
    best_params = [-53.14680163254601, 79.61307254040804, 0.3, 0.05739232362022314]
    df_protonRecFD_1 = df_protonRecFD.loc[df_protonRecFD.PDc1theta < corr(best_params, df_protonRecFD.Pp), :]
    df_protonRecFD_2 = df_protonRecFD.loc[df_protonRecFD.PDc1theta >= corr(best_params, df_protonRecFD.Pp), :]

    #inbending proton energy loss correction
    if pol == "inbending":
        #FD part
        const_FD = -0.00051894 - 0.00018104 * df_protonRecFD_1.Ptheta
        coeff_FD = 3.29466917*10**(-3) +  5.73663160*10**(-4) * df_protonRecFD_1.Ptheta - 1.40807209 * 10**(-5) * df_protonRecFD_1.Ptheta * df_protonRecFD_1.Ptheta
        CorrectedPp_FD_1 = np.select([df_protonRecFD_1.Pp<1, df_protonRecFD_1.Pp>=1], [const_FD + coeff_FD/df_protonRecFD_1.loc[:, "Pp"] + df_protonRecFD_1.loc[:, "Pp"], np.exp(-2.739 - 3.932*df_protonRecFD_1.Pp) + 0.002907+df_protonRecFD_1.Pp])

        const_FD = -0.16742969 + 0.00697925 * df_protonRecFD_1.Ptheta
        coeff_FD = 0.23352115 - 0.01338697 * df_protonRecFD_1.Ptheta
        CorrectedPtheta_FD_1 = const_FD + coeff_FD/df_protonRecFD_1.loc[:, "Pp"]/df_protonRecFD_1.loc[:, "Pp"] + df_protonRecFD_1.loc[:, "Ptheta"]

        const_FD = 0.21192125 -0.0115175 * df_protonRecFD_1.Ptheta
        coeff_FD = -8.94307411*0.1 + 1.66349766*0.1 * df_protonRecFD_1.Ptheta -8.90617559*0.001 * df_protonRecFD_1.Ptheta * df_protonRecFD_1.Ptheta + 1.64803754*0.0001 * df_protonRecFD_1.Ptheta * df_protonRecFD_1.Ptheta * df_protonRecFD_1.Ptheta
        CorrectedPphi_FD_1 = const_FD + coeff_FD/df_protonRecFD_1.loc[:, "Pp"]/df_protonRecFD_1.loc[:, "Pp"] + df_protonRecFD_1.loc[:, "Pphi"]

        const_FD = -3.03346359*10**(-1) + 1.83368163*10**(-2)*df_protonRecFD_2.Ptheta - 2.86486404*10**(-4)*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta
        coeff_FD =  2.01023276*10**(-1) - 1.13312215*10**(-2)*df_protonRecFD_2.Ptheta + 1.82487916*10**(-4)*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta
        CorrectedPp_FD_2 = np.select([df_protonRecFD_2.Pp<1, df_protonRecFD_2.Pp>=1], [const_FD + coeff_FD/df_protonRecFD_2.loc[:, "Pp"] + df_protonRecFD_2.loc[:, "Pp"], np.exp(-1.2 - 4.228*df_protonRecFD_2.Pp) + 0.007502+df_protonRecFD_2.Pp])

        const_FD = 2.04334532 * 10 -1.81052405 * df_protonRecFD_2.Ptheta + 5.32556360*0.01 * df_protonRecFD_2.Ptheta * df_protonRecFD_2.Ptheta -5.23157558*0.0001 * df_protonRecFD_2.Ptheta * df_protonRecFD_2.Ptheta * df_protonRecFD_2.Ptheta
        coeff_FD = 8.74233279 -7.63869344 * 0.1 * df_protonRecFD_2.Ptheta + 2.22376362*0.01 * df_protonRecFD_2.Ptheta * df_protonRecFD_2.Ptheta -2.16457260*0.0001 * df_protonRecFD_2.Ptheta * df_protonRecFD_2.Ptheta * df_protonRecFD_2.Ptheta
        CorrectedPtheta_FD_2 = const_FD + coeff_FD/df_protonRecFD_2.loc[:, "Pp"]/df_protonRecFD_2.loc[:, "Pp"] + df_protonRecFD_2.loc[:, "Ptheta"]

        const_FD = 0.54697831 -0.04896981*df_protonRecFD_2.Ptheta +  0.00111376*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta
        coeff_FD = -4.06733541*10**2 + 2.43696202*10*df_protonRecFD_2.Ptheta -3.36144736*10**(-1)*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta
        coeff2_FD = 2.06378660*10 - 1.42866062*df_protonRecFD_2.Ptheta + 2.01085440*10**(-2)*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta
        CorrectedPphi_FD_2 = const_FD + coeff_FD*np.exp(coeff2_FD*df_protonRecFD_2.loc[:, "Pp"]) + df_protonRecFD_2.loc[:, "Pphi"]

        #CD part
        const_CD = 1.93686914 - 0.116288824*np.minimum(df_protonRecCD.Ptheta, 60) + 0.00223685833*np.minimum(df_protonRecCD.Ptheta, 60)**2 - 1.40771969 * 10**(-5)*np.minimum(df_protonRecCD.Ptheta, 60)**3
        coeff_CD = -0.738047800 + 0.0443343685*np.minimum(df_protonRecCD.Ptheta, 60) - 8.50985972*10**(-4)*np.minimum(df_protonRecCD.Ptheta, 60)*np.minimum(df_protonRecCD.Ptheta, 60) + 5.36810280 * 10**(-6) * np.minimum(df_protonRecCD.Ptheta, 60)**3

        CorrectedPp_CD = const_CD + coeff_CD/df_protonRecCD.loc[:, "Pp"] + df_protonRecCD.loc[:, "Pp"]

        const_CD = -1.09849291*100 + 8.86664014 * np.minimum(df_protonRecCD.Ptheta, 60) - 0.26643881 * np.minimum(df_protonRecCD.Ptheta, 60)**2 + 3.53814210 * 10**(-3) * np.minimum(df_protonRecCD.Ptheta, 60)**3 - 1.75297107 * 10**(-5) * np.minimum(df_protonRecCD.Ptheta, 60)**4
        coeff_CD = 9.52034523*100 -5.74808292 * 10 * np.minimum(df_protonRecCD.Ptheta, 60) + 1.15386949 * np.minimum(df_protonRecCD.Ptheta, 60)**2 - 7.57970373 * 0.001 * np.minimum(df_protonRecCD.Ptheta, 60)**3
        coeff2_CD = -2.00387313*100 + 1.18979079 * 10 * np.minimum(df_protonRecCD.Ptheta, 60) - 2.37730217*0.1 * np.minimum(df_protonRecCD.Ptheta, 60)**2 + 1.55153003*0.001*np.minimum(df_protonRecCD.Ptheta, 60)**3

        CorrectedPtheta_CD = const_CD + coeff_CD*np.exp(coeff2_CD*df_protonRecCD.loc[:, "Pp"]) + df_protonRecCD.loc[:, "Ptheta"]

        const_CD = 4.94546178 -3.26662886*0.1 * np.minimum(df_protonRecCD.Ptheta, 60) +  7.39069603 * 0.001 * np.minimum(df_protonRecCD.Ptheta, 60)**2 -6.83599356*10**(-5) * np.minimum(df_protonRecCD.Ptheta, 60)**3 + 2.12303103*10**(-7) * np.minimum(df_protonRecCD.Ptheta, 60)**4
        coeff_CD = 1.72181613*10**(5) -1.36827111*10**(4) * np.minimum(df_protonRecCD.Ptheta, 60) + 4.00923146*10**(2) * np.minimum(df_protonRecCD.Ptheta, 60)**2 - 5.12792347 * np.minimum(df_protonRecCD.Ptheta, 60)**3 + 2.41793167*10**(-2) * np.minimum(df_protonRecCD.Ptheta, 60)**4
        coeff2_CD =  1.20477219*10**(2) -5.86630228 * np.minimum(df_protonRecCD.Ptheta, 60) + 7.44007875*10**(-2) * np.minimum(df_protonRecCD.Ptheta, 60)**2 -2.42652473*10**(-4) * np.minimum(df_protonRecCD.Ptheta, 60)**3
        CorrectedPphi_CD = const_CD + coeff_CD*np.exp(coeff2_CD*df_protonRecCD.loc[:, "Pp"]) + df_protonRecCD.loc[:, "Pphi"]

    #outbending proton energy loss correction
    elif pol == "outbending":
        #FD part
        const_FD = 0.05083242 -0.00469777*df_protonRecFD_1.Ptheta + 0.0001082*df_protonRecFD_1.Ptheta*df_protonRecFD_1.Ptheta
        coeff_FD = -1.47443264*0.01 + 1.58220893*0.001*df_protonRecFD_1.Ptheta -3.19490013*0.00001*df_protonRecFD_1.Ptheta*df_protonRecFD_1.Ptheta
        CorrectedPp_FD_1 = np.select([df_protonRecFD_1.Pp<1, df_protonRecFD_1.Pp>=1], [const_FD + coeff_FD/df_protonRecFD_1.loc[:, "Pp"] + df_protonRecFD_1.loc[:, "Pp"], np.exp(-2.739 - 3.932*df_protonRecFD_1.Pp) + 0.002907 + df_protonRecFD_1.Pp])

        const_FD = -2.56460305*10 + 3.29877542*df_protonRecFD_1.Ptheta -1.43106886*0.1*df_protonRecFD_1.Ptheta*df_protonRecFD_1.Ptheta + 2.08341898*0.001*df_protonRecFD_1.Ptheta*df_protonRecFD_1.Ptheta*df_protonRecFD_1.Ptheta
        coeff_FD =  9.12532740*10 -1.20100762*10*df_protonRecFD_1.Ptheta + 5.27654711*0.1*df_protonRecFD_1.Ptheta*df_protonRecFD_1.Ptheta -7.72656759*0.001*df_protonRecFD_1.Ptheta*df_protonRecFD_1.Ptheta*df_protonRecFD_1.Ptheta
        CorrectedPtheta_FD_1 = const_FD + coeff_FD/df_protonRecFD_1.loc[:, "Pp"]/df_protonRecFD_1.loc[:, "Pp"] + df_protonRecFD_1.loc[:, "Ptheta"]

        const_FD = -20.4780893 + 1.67020488*df_protonRecFD_1.Ptheta - 0.03419348*df_protonRecFD_1.Ptheta*df_protonRecFD_1.Ptheta
        coeff_FD = 35.02807194 - 2.9098043*df_protonRecFD_1.Ptheta +  0.06037906*df_protonRecFD_1.Ptheta*df_protonRecFD_1.Ptheta
        CorrectedPphi_FD_1 = const_FD + coeff_FD/df_protonRecFD_1.loc[:, "Pp"]/df_protonRecFD_1.loc[:, "Pp"] + df_protonRecFD_1.loc[:, "Pphi"]

        const_FD = 0.09832589 -0.0066463*df_protonRecFD_2.Ptheta + 0.00010312*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta
        coeff_FD = -9.61421691*0.01 + 6.85638807*0.001*df_protonRecFD_2.Ptheta -9.75766427*0.00001*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta
        CorrectedPp_FD_2 = np.select([df_protonRecFD_2.Pp<1, df_protonRecFD_2.Pp>=1], [const_FD + coeff_FD/df_protonRecFD_2.loc[:, "Pp"] + df_protonRecFD_2.loc[:, "Pp"], np.exp(-1.871 - 3.063*df_protonRecFD_2.Pp) + 0.007517 + df_protonRecFD_2.Pp])

        const_FD = -1.68873940 + 9.56867163*0.01*df_protonRecFD_2.Ptheta -1.43741464*0.001*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta
        coeff_FD = 1.49978357*10 -1.40137094*df_protonRecFD_2.Ptheta + 4.38501543*0.01*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta -4.57982872*0.0001*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta
        CorrectedPtheta_FD_2 = const_FD + coeff_FD/df_protonRecFD_2.loc[:, "Pp"]/df_protonRecFD_2.loc[:, "Pp"] + df_protonRecFD_2.loc[:, "Ptheta"]

        const_FD = 6.75359137 - 0.43199851*df_protonRecFD_2.Ptheta + 0.0068995*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta
        coeff_FD = -1.68588219 + 1.05609627*0.1*df_protonRecFD_2.Ptheta -1.50452832*0.001*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta
        CorrectedPphi_FD_2 = const_FD + coeff_FD/df_protonRecFD_2.loc[:, "Pp"]/df_protonRecFD_2.loc[:, "Pp"] + df_protonRecFD_2.loc[:, "Pphi"]
        #CD part
        const_CD = 1.92657376 - 0.113836734*np.minimum(df_protonRecCD.Ptheta, 60) + 0.00215038526*np.minimum(df_protonRecCD.Ptheta, 60)**2 - 1.32525053 * 10**(-5)*np.minimum(df_protonRecCD.Ptheta, 60)**3
        coeff_CD = -0.755650043 + 0.0445538936*np.minimum(df_protonRecCD.Ptheta, 60) - 8.38241864*10**(-4)*np.minimum(df_protonRecCD.Ptheta, 60)*np.minimum(df_protonRecCD.Ptheta, 60) + 5.16887255 * 10**(-6) * np.minimum(df_protonRecCD.Ptheta, 60)**3

        CorrectedPp_CD = const_CD + coeff_CD/df_protonRecCD.loc[:, "Pp"] + df_protonRecCD.loc[:, "Pp"]

        const_CD = -5.79024055*10 + 4.67197531 * np.minimum(df_protonRecCD.Ptheta, 60) - 0.140156897 * np.minimum(df_protonRecCD.Ptheta, 60)**2 + 1.85853057 * 10**(-3) * np.minimum(df_protonRecCD.Ptheta, 60)**3 - 9.19989908 * 10**(-6) * np.minimum(df_protonRecCD.Ptheta, 60)**4
        coeff_CD = 2.99700765*1000 - 2.18027982 * 10**2 * np.minimum(df_protonRecCD.Ptheta, 60) + 5.84757503 * np.minimum(df_protonRecCD.Ptheta, 60)**2 - 6.80409195 * 0.01 * np.minimum(df_protonRecCD.Ptheta, 60)**3 + 2.89244618 * 0.0001 * np.minimum(df_protonRecCD.Ptheta, 60)**4
        coeff2_CD = -1.82237904*100 + 1.10153549 * 10 * np.minimum(df_protonRecCD.Ptheta, 60) - 2.24699931*0.1 * np.minimum(df_protonRecCD.Ptheta, 60)**2 + 1.49390960*0.001*np.minimum(df_protonRecCD.Ptheta, 60)**3

        CorrectedPtheta_CD = const_CD + coeff_CD*np.exp(coeff2_CD*df_protonRecCD.loc[:, "Pp"]) + df_protonRecCD.loc[:, "Ptheta"]

        const_CD = 7.58761670 - 5.28224578*0.1 * np.minimum(df_protonRecCD.Ptheta, 60) +  1.31580117 * 0.01 * np.minimum(df_protonRecCD.Ptheta, 60)**2 -1.41738951*10**(-4) * np.minimum(df_protonRecCD.Ptheta, 60)**3 + 5.62884363*10**(-7) * np.minimum(df_protonRecCD.Ptheta, 60)**4
        coeff_CD = 1.07644097*10**(5) - 8.67994639*10**(3) * np.minimum(df_protonRecCD.Ptheta, 60) + 2.57187193*10**(2) * np.minimum(df_protonRecCD.Ptheta, 60)**2 - 3.31379317 * np.minimum(df_protonRecCD.Ptheta, 60)**3 + 1.56896621*10**(-2) * np.minimum(df_protonRecCD.Ptheta, 60)**4
        coeff2_CD =  1.92263184*10**(2) -1.00870704 * 10 * np.minimum(df_protonRecCD.Ptheta, 60) + 1.56575252*10**(-1) * np.minimum(df_protonRecCD.Ptheta, 60)**2 -7.71489734*10**(-4) * np.minimum(df_protonRecCD.Ptheta, 60)**3
        CorrectedPphi_CD = const_CD + coeff_CD*np.exp(coeff2_CD*df_protonRecCD.loc[:, "Pp"]) + df_protonRecCD.loc[:, "Pphi"]

    if len(df_protonRecFD_1):
        df_protonRecFD_1.loc[:, "Pp"] = CorrectedPp_FD_1
        df_protonRecFD_1.loc[:, "Ptheta"] = CorrectedPtheta_FD_1
        df_protonRecFD_1.loc[:, "Pphi"] = CorrectedPphi_FD_1

    if len(df_protonRecFD_2):
        df_protonRecFD_2.loc[:, "Pp"] = CorrectedPp_FD_2
        df_protonRecFD_2.loc[:, "Ptheta"] = CorrectedPtheta_FD_2
        df_protonRecFD_2.loc[:, "Pphi"] = CorrectedPphi_FD_2
    
    if len(df_protonRecCD):
        df_protonRecCD.loc[:, "Pp"] = CorrectedPp_CD
        df_protonRecCD.loc[:, "Ptheta"] = CorrectedPtheta_CD
        df_protonRecCD.loc[:, "Pphi"] = CorrectedPphi_CD
    df_protonRec = pd.concat([df_protonRecFD_1, df_protonRecFD_2, df_protonRecCD, df_protonRecOthers])

    #moduli proton phi
    df_protonRec.loc[:, "Pphi"] = np.where(df_protonRec.loc[:, "Pphi"]%360<180, df_protonRec.loc[:, "Pphi"]%360, df_protonRec.loc[:, "Pphi"]%360-360)

    df_protonRec.loc[:, "Ppx"] = df_protonRec.loc[:, "Pp"]*np.sin(np.radians(df_protonRec.loc[:, "Ptheta"]))*np.cos(np.radians(df_protonRec.loc[:, "Pphi"]))
    df_protonRec.loc[:, "Ppy"] = df_protonRec.loc[:, "Pp"]*np.sin(np.radians(df_protonRec.loc[:, "Ptheta"]))*np.sin(np.radians(df_protonRec.loc[:, "Pphi"]))
    df_protonRec.loc[:, "Ppz"] = df_protonRec.loc[:, "Pp"]*np.cos(np.radians(df_protonRec.loc[:, "Ptheta"]))
    pro = [df_protonRec['Ppx'], df_protonRec['Ppy'], df_protonRec['Ppz']]

    df_protonRec.loc[:, 'Pe'] = getEnergy(pro, M)

    return df_protonRec

def quadratic(x, *par):
    a, b, c = par
    return a*x**2 + b*x + c

def protonMomentumCorrection(pol, df_protonRec):
    df_protonRec = copy(df_protonRec)
    if pol == "inbending":
        dp_sector_1_params =  [4.18310463e-02, -9.75322577e-02, 3.92241528e-02]
        dtheta_sector_1_params =  [1.12117990e+00, -1.45651738e+00, -3.61966664e-01]
        dp_sector_2_params =  [9.72618261e-02, -1.47404991e-01, 4.38454459e-02]
        dtheta_sector_2_params =  [2.23808712e-01, 2.64132351e-01, -1.16260483e+00]
        dp_sector_3_params =  [1.50782345e-01, -2.43829059e-01, 8.30817121e-02]
        dtheta_sector_3_params =  [-5.11801825e-01, 1.24494702e+00, -1.37989466e+00]
        dp_sector_4_params =  [4.24673974e-02, -9.95181630e-02, 3.51056884e-02]
        dtheta_sector_4_params =  [9.62471223e-01, -1.24613297e+00, -3.94245612e-01]
        dp_sector_5_params =  [4.51632734e-02, -1.02533879e-01, 3.94177517e-02]
        dtheta_sector_5_params =  [5.55969390e-01, -6.34977038e-01, -6.31316083e-01]
        dp_sector_6_params =  [-4.44092504e-02, 5.21339250e-02, -2.32922494e-02]
        dtheta_sector_6_params =  [-1.98404520e-02, 1.07806469e+00, -1.78427373e+00]
        dp_sector_CD_params =  [-3.09139398e-02, 4.19387700e-02, 7.96761500e-03]
        dtheta_sector_CD_params =  [2.38931902e+00, -3.86848852e+00, 7.50606916e-01]
    if pol == "outbending":
        dp_sector_1_params =  [7.93399152e-02, -1.63161999e-01, 4.64574242e-02]
        dtheta_sector_1_params =  [1.06523064e+00, -1.76103124e+00, 7.66716947e-01]
        dp_sector_2_params =  [3.78073457e-02, -8.69496631e-02, 6.99047328e-03]
        dtheta_sector_2_params =  [4.38259463e-01, -4.35700753e-01, 1.05241964e-01]
        dp_sector_3_params =  [-2.90113621e-03, -2.94800083e-03, -3.40115366e-02]
        dtheta_sector_3_params =  [2.00545395e-01, 2.00951967e-01, -3.18864848e-01]
        dp_sector_4_params =  [8.73715031e-03, -2.19339609e-02, -2.60223791e-02]
        dtheta_sector_4_params =  [-4.21325885e-01, 1.44347102e+00, -8.03632979e-01]
        dp_sector_5_params =  [5.80017470e-02, -1.14948332e-01, 2.05773011e-02]
        dtheta_sector_5_params =  [-2.39263572e-01, 9.61961521e-01, -5.00992863e-01]
        dp_sector_6_params =  [2.07106073e-02, -4.59478839e-02, -4.07474515e-03]
        dtheta_sector_6_params =  [-3.30239417e-03, 3.56104152e-01, -1.85550026e-01]
        dp_sector_CD_params =  [6.69872161e-02, -1.42457411e-01, 7.15694609e-02]
        dtheta_sector_CD_params =  [-1.33140658e+00, 3.28530583e+00, -2.08330489e+00]

    # dp_sector_CD_params =  [7.80244462e-02, -1.69195711e-01, 8.67323653e-02]
    # dtheta_sector_CD_params =  [7.90127087e-01, -1.03797405e+00, -1.42779248e-01]


    df_protonRec.loc[df_protonRec.Psector == 1, "Pp"] = df_protonRec.loc[df_protonRec.Psector == 1, "Pp"] + quadratic(df_protonRec.loc[df_protonRec.Psector == 1, "Pp"], *dp_sector_1_params)
    df_protonRec.loc[df_protonRec.Psector == 2, "Pp"] = df_protonRec.loc[df_protonRec.Psector == 2, "Pp"] + quadratic(df_protonRec.loc[df_protonRec.Psector == 2, "Pp"], *dp_sector_2_params)
    df_protonRec.loc[df_protonRec.Psector == 3, "Pp"] = df_protonRec.loc[df_protonRec.Psector == 3, "Pp"] + quadratic(df_protonRec.loc[df_protonRec.Psector == 3, "Pp"], *dp_sector_3_params)
    df_protonRec.loc[df_protonRec.Psector == 4, "Pp"] = df_protonRec.loc[df_protonRec.Psector == 4, "Pp"] + quadratic(df_protonRec.loc[df_protonRec.Psector == 4, "Pp"], *dp_sector_4_params)
    df_protonRec.loc[df_protonRec.Psector == 5, "Pp"] = df_protonRec.loc[df_protonRec.Psector == 5, "Pp"] + quadratic(df_protonRec.loc[df_protonRec.Psector == 5, "Pp"], *dp_sector_5_params)
    df_protonRec.loc[df_protonRec.Psector == 6, "Pp"] = df_protonRec.loc[df_protonRec.Psector == 6, "Pp"] + quadratic(df_protonRec.loc[df_protonRec.Psector == 6, "Pp"], *dp_sector_6_params)
    df_protonRec.loc[df_protonRec.Psector >4000, "Pp"] = df_protonRec.loc[df_protonRec.Psector >4000, "Pp"] + quadratic(df_protonRec.loc[df_protonRec.Psector >4000, "Pp"], *dp_sector_CD_params)
    df_protonRec.loc[df_protonRec.Psector == 1, "Ptheta"] = df_protonRec.loc[df_protonRec.Psector == 1, "Ptheta"] + quadratic(df_protonRec.loc[df_protonRec.Psector == 1, "Pp"], *dtheta_sector_1_params)
    df_protonRec.loc[df_protonRec.Psector == 2, "Ptheta"] = df_protonRec.loc[df_protonRec.Psector == 2, "Ptheta"] + quadratic(df_protonRec.loc[df_protonRec.Psector == 2, "Pp"], *dtheta_sector_2_params)
    df_protonRec.loc[df_protonRec.Psector == 3, "Ptheta"] = df_protonRec.loc[df_protonRec.Psector == 3, "Ptheta"] + quadratic(df_protonRec.loc[df_protonRec.Psector == 3, "Pp"], *dtheta_sector_3_params)
    df_protonRec.loc[df_protonRec.Psector == 4, "Ptheta"] = df_protonRec.loc[df_protonRec.Psector == 4, "Ptheta"] + quadratic(df_protonRec.loc[df_protonRec.Psector == 4, "Pp"], *dtheta_sector_4_params)
    df_protonRec.loc[df_protonRec.Psector == 5, "Ptheta"] = df_protonRec.loc[df_protonRec.Psector == 5, "Ptheta"] + quadratic(df_protonRec.loc[df_protonRec.Psector == 5, "Pp"], *dtheta_sector_5_params)
    df_protonRec.loc[df_protonRec.Psector == 6, "Ptheta"] = df_protonRec.loc[df_protonRec.Psector == 6, "Ptheta"] + quadratic(df_protonRec.loc[df_protonRec.Psector == 6, "Pp"], *dtheta_sector_6_params)
    df_protonRec.loc[df_protonRec.Psector >4000, "Ptheta"] = df_protonRec.loc[df_protonRec.Psector >4000, "Ptheta"] + quadratic(df_protonRec.loc[df_protonRec.Psector >4000, "Pp"], *dtheta_sector_CD_params)
    #moduli proton phi
    df_protonRec.loc[:, "Pphi"] = np.where(df_protonRec.loc[:, "Pphi"]%360<180, df_protonRec.loc[:, "Pphi"]%360, df_protonRec.loc[:, "Pphi"]%360-360)

    df_protonRec.loc[:, "Ppx"] = df_protonRec.loc[:, "Pp"]*np.sin(np.radians(df_protonRec.loc[:, "Ptheta"]))*np.cos(np.radians(df_protonRec.loc[:, "Pphi"]))
    df_protonRec.loc[:, "Ppy"] = df_protonRec.loc[:, "Pp"]*np.sin(np.radians(df_protonRec.loc[:, "Ptheta"]))*np.sin(np.radians(df_protonRec.loc[:, "Pphi"]))
    df_protonRec.loc[:, "Ppz"] = df_protonRec.loc[:, "Pp"]*np.cos(np.radians(df_protonRec.loc[:, "Ptheta"]))
    pro = [df_protonRec['Ppx'], df_protonRec['Ppy'], df_protonRec['Ppz']]

    df_protonRec.loc[:, 'Pe'] = getEnergy(pro, M)
    return df_protonRec

def protonMomentumCorrection_legacy(pol, df_protonRec):
    df_protonRec = copy(df_protonRec)
    df_protonRecFD = df_protonRec.loc[df_protonRec.Psector<7, :]
    df_protonRecCD = df_protonRec.loc[(df_protonRec.Psector>7) & (df_protonRec.PthetaOrig<75), :]
    df_protonRecOthers = df_protonRec.loc[ ((df_protonRec.Psector>7) & (df_protonRec.PthetaOrig>=75)) | ((df_protonRec.Psector<7) & (df_protonRec.PpOrig<0.3)), :]

    df_protonRecCD.loc[:, "Pp"] = df_protonRecCD.Pp + 0.01
    df_protonRecCD.loc[:, "Ptheta"] = df_protonRecCD.Ptheta - 0.002129*df_protonRecCD.Ptheta**2 + 0.198*df_protonRecCD.Ptheta - 4.762 -0.2/(1+np.exp((df_protonRecCD.Pp-0.55)/(-0.05)))
    df_protonRecCD.loc[:, "Pphi"] = df_protonRecCD.Pphi
    if pol == "inbending":
        corr = np.poly1d([1.671, -4.918, 5.151, -2.434])(df_protonRecFD.Pp)       
        corr = np.where(corr<0, corr, 0)
        df_protonRecFD.loc[:, "Ptheta"] = df_protonRecFD.Ptheta + corr #corr scale -1 to -0.3 degrees
    if pol == "outbending":
        df_protonRecFD.loc[df_protonRecFD.Psector<7, "Pp"] = df_protonRecFD.loc[df_protonRecFD.Psector<7, "Pp"] - 0.02
        df_protonRecFD.loc[:, "Ptheta"] = df_protonRecFD.Ptheta + 0.05*(np.abs(df_protonRecFD.Ptheta - 27) + (df_protonRecFD.Ptheta - 27))

    df_protonRec = pd.concat([df_protonRecFD, df_protonRecCD, df_protonRecOthers])

    #moduli proton phi
    df_protonRec.loc[:, "Pphi"] = np.where(df_protonRec.loc[:, "Pphi"]%360<180, df_protonRec.loc[:, "Pphi"]%360, df_protonRec.loc[:, "Pphi"]%360-360)

    df_protonRec.loc[:, "Ppx"] = df_protonRec.loc[:, "Pp"]*np.sin(np.radians(df_protonRec.loc[:, "Ptheta"]))*np.cos(np.radians(df_protonRec.loc[:, "Pphi"]))
    df_protonRec.loc[:, "Ppy"] = df_protonRec.loc[:, "Pp"]*np.sin(np.radians(df_protonRec.loc[:, "Ptheta"]))*np.sin(np.radians(df_protonRec.loc[:, "Pphi"]))
    df_protonRec.loc[:, "Ppz"] = df_protonRec.loc[:, "Pp"]*np.cos(np.radians(df_protonRec.loc[:, "Ptheta"]))
    pro = [df_protonRec['Ppx'], df_protonRec['Ppy'], df_protonRec['Ppz']]

    df_protonRec.loc[:, 'Pe'] = getEnergy(pro, M)

    return df_protonRec

def cubic(args, x): #equivalent to poly1d
    a, b, c, d = args
    return a*x**3 +b*x**2 + c*x + d

def quartic(x, sector, pol = "inbending"):
    args_sigmas_FD_inb = [[-0.233, 1.216, -2.279, 1.812, -0.445], [ 0.277, -1.366, 2.318, -1.619,  0.466 ],[ 0.0728, -0.223, 0.0888,  0.225, -0.0889],[-0.204, 0.977, -1.766, 1.411, -0.342], [ 0.277, -1.059, 1.362, -0.641, 0.137], [-0.219, 1.132, -2.153, 1.763, -0.447]]
    args_sigmas_FD_outb = [[0.481,-1.548, 1.524,-0.415, 0.0277], [1.872, -8.054, 12.536, -8.358,  2.083], [-0.0656, 0.480, -1.191, 1.169, -0.315], [-1.559, 7.356, -12.639, 9.312,  -2.405], [ 0.189, -0.344, -0.253,  0.717, -0.238], [0.466, -1.560, 1.622, -0.485, 0.0322]]
    if pol == "inbending":
        a, b, c, d, e = args_sigmas_FD_inb[sector - 1]
        return np.select( [x<0.55, (x>=0.55)& (x < 1.55), x>=1.55], [a*0.55**4+b*0.55**3+c*0.55**2+d*0.55+e, a*x**4 +b*x**3 + c*x**2 + d*x + e, a*1.55**4 + b*1.55**3 + c*1.55**2+d*1.55 +e])
    if pol == "outbending":
        a, b, c, d, e = args_sigmas_FD_outb[sector - 1]
        return np.select( [x<0.65, (x>=0.65)& (x < 1.55), x>=1.55], [a*0.65**4+b*0.65**3+c*0.65**2+d*0.65+e, a*x**4 +b*x**3 + c*x**2 + d*x + e, a*1.55**4 + b*1.55**3 + c*1.55**2+d*1.55 +e])

def sigmaFDOutb(x):
    return np.select([x<.95, (x>=.95) & (x<1.2), (x>=1.2)&(x<1.575), (x>=1.575) & (x<1.9), (x>1.9)], [0.1, -0.045/(1.2-.95)*x+1.2*0.045/(1.2-.95) + 0.055, 0.055, -0.015/(1.9-1.575)*x+1.9*0.015/(1.9-1.575) + 0.04,0.04])


def protonMomentumSmearing(pol, df_protonRec, smearing = 1):
    df_protonRec = copy(df_protonRec)
    regulator = np.abs(2*(1/(1+np.exp(-(df_protonRec.loc[df_protonRec["Psector"]>7, "Pp"]-0.3)/0.01))-0.5))
    sigma1_CD = np.where(df_protonRec.loc[df_protonRec["Psector"]>7, "Pp"]<0.85, cubic([0.0926, 0.137, -0.230, 0.139], df_protonRec.loc[df_protonRec["Psector"]>7, "Pp"]), 0.1)
    sigma2_CD = np.where(df_protonRec.loc[df_protonRec["Psector"]>7, "Pp"]<1.34, cubic([-2.797, 9.351, -9.488, 3.503], df_protonRec.loc[df_protonRec["Psector"]>7, "Pp"]), 0.85)
    sigma3_CD = 0.8 + 2.2/(1+np.exp(5.518*(df_protonRec.loc[df_protonRec.Psector>7, "Pp"]-0.625)))
    #CD proton
    df_protonRec.loc[df_protonRec["Psector"]>7, "Pp"] = df_protonRec.loc[df_protonRec["Psector"]>7, "Pp"]*np.random.normal(1, smearing*regulator*sigma1_CD, len(df_protonRec.loc[df_protonRec.Psector>7]))
    df_protonRec.loc[df_protonRec["Psector"]>7, "Ptheta"] = df_protonRec.loc[df_protonRec["Psector"]>7, "Ptheta"] + np.random.normal(0, smearing*sigma2_CD, len(df_protonRec.loc[df_protonRec.Psector>7]))
    df_protonRec.loc[df_protonRec["Psector"]>7, "Pphi"] = df_protonRec.loc[df_protonRec["Psector"]>7, "Pphi"] + np.random.normal(0, smearing*sigma3_CD, len(df_protonRec.loc[df_protonRec.Psector>7])) 
    #FD proton
    for sector in range(1, 7):
        if pol == "inbending":
            regulator = (1/(1+np.exp(-(df_protonRec.loc[df_protonRec["Psector"]==sector, "Pp"]-0.5)/0.05)))
            sigmas_FD = quartic(df_protonRec.loc[df_protonRec.Psector == sector, "Pp"], sector, pol)
        elif pol == "outbending":
            regulator = (1/(1+np.exp(-(df_protonRec.loc[df_protonRec["Psector"]==sector, "Pp"]-0.6)/0.05)))
            sigmas_FD = sigmaFDOutb(df_protonRec.loc[df_protonRec["Psector"]==sector, "Pp"]) #quartic(df_protonRec.loc[df_protonRec.Psector == sector, "Pp"], sector, pol)
        df_protonRec.loc[df_protonRec["Psector"]==sector, "Pp"] = df_protonRec.loc[df_protonRec["Psector"]==sector, "Pp"]*np.random.normal(1, smearing*regulator*sigmas_FD, len(df_protonRec.loc[df_protonRec["Psector"]==sector, "Pp"]))

    #moduli proton phi
    df_protonRec.loc[:, "Pphi"] = np.where(df_protonRec.loc[:, "Pphi"]%360<180, df_protonRec.loc[:, "Pphi"]%360, df_protonRec.loc[:, "Pphi"]%360-360)

    df_protonRec.loc[:, "Ppx"] = df_protonRec.loc[:, "Pp"]*np.sin(np.radians(df_protonRec.loc[:, "Ptheta"]))*np.cos(np.radians(df_protonRec.loc[:, "Pphi"]))
    df_protonRec.loc[:, "Ppy"] = df_protonRec.loc[:, "Pp"]*np.sin(np.radians(df_protonRec.loc[:, "Ptheta"]))*np.sin(np.radians(df_protonRec.loc[:, "Pphi"]))
    df_protonRec.loc[:, "Ppz"] = df_protonRec.loc[:, "Pp"]*np.cos(np.radians(df_protonRec.loc[:, "Ptheta"]))
    pro = [df_protonRec['Ppx'], df_protonRec['Ppy'], df_protonRec['Ppz']]

    df_protonRec.loc[:, 'Pe'] = getEnergy(pro, M)

    return df_protonRec

def gammaMomentumSmearing(df_gammaRec, smearing = 1):
    df_gammaRec = copy(df_gammaRec)
    gam = [df_gammaRec['Gpx'], df_gammaRec['Gpy'], df_gammaRec['Gpz']]
    df_gammaRec.loc[:, 'Gp'] = mag(gam)
    df_gammaRec.loc[:, 'Gtheta'] = getTheta(gam)
    df_gammaRec.loc[:, 'Gphi'] = getPhi(gam)
    #FT photon
    df_gammaRec.loc[df_gammaRec["Gsector"]>7, "Gp"] = df_gammaRec.loc[df_gammaRec["Gsector"]>7, "Gp"]*np.random.normal(1, smearing*(0.013 + 0.003/(1+np.exp(0.761*(df_gammaRec.loc[df_gammaRec["Gsector"]>7, "Gp"]-6)))), len(df_gammaRec.loc[df_gammaRec.Gsector>7]))
    #FD photon
    df_gammaRec.loc[df_gammaRec["Gsector"]<7, "Gp"] = df_gammaRec.loc[df_gammaRec["Gsector"]<7, "Gp"]*np.random.normal(1, smearing*(0.0395/(1+np.exp(5.308*(df_gammaRec.loc[df_gammaRec["Gsector"]<7, "Gp"]- 8.005)))), len(df_gammaRec.loc[df_gammaRec.Gsector<7]))
    df_gammaRec.loc[:, "Gpx"] = df_gammaRec.loc[:, "Gp"]*np.sin(np.radians(df_gammaRec.loc[:, "Gtheta"]))*np.cos(np.radians(df_gammaRec.loc[:, "Gphi"]))
    df_gammaRec.loc[:, "Gpy"] = df_gammaRec.loc[:, "Gp"]*np.sin(np.radians(df_gammaRec.loc[:, "Gtheta"]))*np.sin(np.radians(df_gammaRec.loc[:, "Gphi"]))
    df_gammaRec.loc[:, "Gpz"] = df_gammaRec.loc[:, "Gp"]*np.cos(np.radians(df_gammaRec.loc[:, "Gtheta"]))

    df_gammaRec.loc[:, "Ge"] = df_gammaRec.loc[:, "Gp"]
    return df_gammaRec

def cubic_without_const(args, x):
    a, b, c = args
    return a*x**3 + b*x**2 + c*x

def quintic_without_const(args, x):
    a, b, c = args
    if b < 0:
        return 0*x
    return a*x*(x-b)**3 * (x-c)

def quartic_without_const(args, x):
    a, b, c, d = args
    x = np.array(x)
    return a*x**4 + b*x**3 + c*x**2 + d*x**1

def gammaMomentumCorrection(pol, df_gg, df_gammaRec):
    print("applying the photon kinematic corrections for " + pol)
    df_gammaRec = copy(df_gammaRec)
    df_gg       = copy(df_gg)

    #photon kinematic correction of df_gg
    gam = [df_gg['Gpx'], df_gg['Gpy'], df_gg['Gpz']]
    df_gg.loc[:, 'Gp'] = mag(gam)
    df_gg.loc[:, 'Gtheta'] = getTheta(gam)
    df_gg.loc[:, 'Gphi'] = getPhi(gam)
    #photon kinematic correction of df_gamma
    gam = [df_gammaRec['Gpx'], df_gammaRec['Gpy'], df_gammaRec['Gpz']]
    df_gammaRec.loc[:, 'Gp'] = mag(gam)
    df_gammaRec.loc[:, 'Gtheta'] = getTheta(gam)
    df_gammaRec.loc[:, 'Gphi'] = getPhi(gam)

    #FT - df_gg: perform correction for only one photon
    FT_phot_corr = 0.02815846*df_gg.loc[df_gg["Gsector"]>7, "Gp"]#(-0.00467*df_gg.loc[df_gg["Gsector"]>7, "Gp"]**2 + 0.0802 *df_gg.loc[df_gg["Gsector"]>7, "Gp"]  -0.352) + 0.25
    df_gg.loc[df_gg["Gsector"]>7, "Gp"] = df_gg.loc[df_gg["Gsector"]>7, "Gp"] + FT_phot_corr
    FT_phot_corr2 = 0.02815846*df_gg.loc[df_gg["Gsector2"]>7, "Gp2"]#(-0.00467*df_gg.loc[df_gg["Gsector"]>7, "Gp"]**2 + 0.0802 *df_gg.loc[df_gg["Gsector"]>7, "Gp"]  -0.352) + 0.25
    df_gg.loc[df_gg["Gsector2"]>7, "Gp2"] = df_gg.loc[df_gg["Gsector2"]>7, "Gp2"] + FT_phot_corr2
    #FT - df_gammaRec: perform every photon
    FT_phot_corr = 0.02815846*df_gammaRec.loc[df_gammaRec["Gsector"]>7, "Gp"]#(-0.00467*df_gammaRec.loc[df_gammaRec["Gsector"]>7, "Gp"]**2 + 0.0802 *df_gammaRec.loc[df_gammaRec["Gsector"]>7, "Gp"]  -0.352) + 0.25
    df_gammaRec.loc[df_gammaRec["Gsector"]>7, "Gp"] = df_gammaRec.loc[df_gammaRec["Gsector"]>7, "Gp"] + FT_phot_corr

    #FD
    if pol == "inbending":
        args = [[-0.0000732, 1.480, 9.344], [-0.000135, 3.070, 9.248], [-0.0000437, 0.719, 9.873], [-0.0000428, 0.00234, 0.0103], [0.000250, -0.00314, 0.0232], [-0.0000454, 0.517, 9.447]]
        funcs = [quintic_without_const, quintic_without_const, quintic_without_const, cubic_without_const, cubic_without_const, quintic_without_const]
        funcs_minor = [quintic_without_const, quintic_without_const, quintic_without_const, cubic_without_const, cubic_without_const, cubic_without_const]
        args_minor = [[-0.0000168, 0.821, 8.894], [-0.0000340, 2.720, 8.419], [-0.0000620, 2.793, 8.865], [ 0.000132, -0.00162,  0.00978], [-0.000135,  0.000282, 0.00650], [ 0.000263,  -0.00293,   0.0139]]
        for sector in range(1, 7):
            cond = df_gg.Gsector == sector
            FD_phot_corr_sector = funcs[sector-1](args[sector-1], df_gg.loc[cond, "Gp"])
            df_gg.loc[cond, "Gp"] = df_gg.loc[cond, "Gp"] + FD_phot_corr_sector
            FD_phot_corr_minor_sector = funcs_minor[sector-1](args_minor[sector-1], df_gg.loc[cond, "Gp"])
            df_gg.loc[cond, "Gp"] = df_gg.loc[cond, "Gp"] + FD_phot_corr_minor_sector

            cond = df_gammaRec.Gsector == sector
            FD_phot_corr_sector = funcs[sector-1](args[sector-1], df_gammaRec.loc[cond, "Gp"])
            df_gammaRec.loc[cond, "Gp"] = df_gammaRec.loc[cond, "Gp"] + FD_phot_corr_sector
            FD_phot_corr_minor_sector = funcs_minor[sector-1](args_minor[sector-1], df_gammaRec.loc[cond, "Gp"])
            df_gammaRec.loc[cond, "Gp"] = df_gammaRec.loc[cond, "Gp"] + FD_phot_corr_minor_sector

    if pol == "outbending":
        args = [[-0.000615,  0.0113, -0.0600,   0.115],[-0.000334,  0.00656, -0.0383,  0.0934],[-0.000911,  0.0157, -0.0806,  0.154],[ 0.000117, -0.000905,  0.00215,  0.0331],[-0.000119,  0.000979, -0.00400,  0.0499],[-0.000893,  0.0131, -0.0580,  0.111]]
        for sector in range(1, 7):
            cond = df_gg.Gsector == sector
            FD_phot_corr_sector = quartic_without_const(args[sector-1], df_gg.loc[cond, "Gp"])/(1+np.exp(-(df_gg.loc[cond, "Gp"]-2.2)/0.15))
            df_gg.loc[cond, "Gp"] = df_gg.loc[cond, "Gp"] + FD_phot_corr_sector

            cond = df_gammaRec.Gsector == sector
            FD_phot_corr_sector = quartic_without_const(args[sector-1], df_gammaRec.loc[cond, "Gp"])/(1+np.exp(-(df_gammaRec.loc[cond, "Gp"]-2.2)/0.15))
            df_gammaRec.loc[cond, "Gp"] = df_gammaRec.loc[cond, "Gp"] + FD_phot_corr_sector
    # return df_gg, df_gammaRec
    df_gg.loc[:, "Gpx"] = df_gg.loc[:, "Gp"]*np.sin(np.radians(df_gg.loc[:, "Gtheta"]))*np.cos(np.radians(df_gg.loc[:, "Gphi"]))
    df_gg.loc[:, "Gpy"] = df_gg.loc[:, "Gp"]*np.sin(np.radians(df_gg.loc[:, "Gtheta"]))*np.sin(np.radians(df_gg.loc[:, "Gphi"]))
    df_gg.loc[:, "Gpz"] = df_gg.loc[:, "Gp"]*np.cos(np.radians(df_gg.loc[:, "Gtheta"]))

    df_gg.loc[:, "Ge"] = df_gg.loc[:, "Gp"]

    df_gg.loc[:, "Gpx2"] = df_gg.loc[:, "Gp2"]*np.sin(np.radians(df_gg.loc[:, "Gtheta2"]))*np.cos(np.radians(df_gg.loc[:, "Gphi2"]))
    df_gg.loc[:, "Gpy2"] = df_gg.loc[:, "Gp2"]*np.sin(np.radians(df_gg.loc[:, "Gtheta2"]))*np.sin(np.radians(df_gg.loc[:, "Gphi2"]))
    df_gg.loc[:, "Gpz2"] = df_gg.loc[:, "Gp2"]*np.cos(np.radians(df_gg.loc[:, "Gtheta2"]))

    df_gg.loc[:, "Ge2"] = df_gg.loc[:, "Gp2"]


    df_gammaRec.loc[:, "Gpx"] = df_gammaRec.loc[:, "Gp"]*np.sin(np.radians(df_gammaRec.loc[:, "Gtheta"]))*np.cos(np.radians(df_gammaRec.loc[:, "Gphi"]))
    df_gammaRec.loc[:, "Gpy"] = df_gammaRec.loc[:, "Gp"]*np.sin(np.radians(df_gammaRec.loc[:, "Gtheta"]))*np.sin(np.radians(df_gammaRec.loc[:, "Gphi"]))
    df_gammaRec.loc[:, "Gpz"] = df_gammaRec.loc[:, "Gp"]*np.cos(np.radians(df_gammaRec.loc[:, "Gtheta"]))

    df_gammaRec.loc[:, "Ge"] = df_gammaRec.loc[:, "Gp"]
    return df_gg, df_gammaRec