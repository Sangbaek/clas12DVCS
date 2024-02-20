'''
independent script to perform the energy loss corection
'''
from utils.const import *

def electronCorrection(pol, df_electronRec):

    df_electronRec.loc[:, "dp"] = 0

    eleCorr = df_electronRec.loc[:, ["Ep", "Ephi", "Esector"]]
    eleCorr.loc[ (df_electronRec.Esector == 4 | df_electronRec.Esector == 3) & eleCorr.Ephi < 0) || (df_electronRec.Esector > 4 & eleCorr.Ephi < 90), "Ephi"] += 360
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
    df_electronRec.loc[:, "Ep"]  = df.electronRec.Ep + eleCorr.dp
    df_electronRec.loc[:, "Ee"]  = np.sqrt(df_electronRec.Ep**2 + me**2)
    return df_electronRec

def protonEnergyLossCorr(pol, topo, df_protonRec):
    '''
    a simple function for the energy loss correction.
    '''
    if pol == "inbending":
        if topo == "FD lower":
            const = -0.00051894 - 0.00018104 * df_protonRec.Ptheta
            coeff = 3.29466917*10**(-3) +  5.73663160*10**(-4) * df_protonRec.Ptheta - 1.40807209 * 10**(-5) * df_protonRec.Ptheta * df_protonRec.Ptheta
            CorrectedPp = np.select([df_protonRec.Pp<1, df_protonRec.Pp>=1], [const + coeff/df_protonRec.loc[:, "Pp"] + df_protonRec.loc[:, "Pp"], np.exp(-2.739 - 3.932*df_protonRec.Pp) + 0.002907+df_protonRec.Pp])

            const = -0.16742969 + 0.00697925 * df_protonRec.Ptheta
            coeff = 0.23352115 - 0.01338697 * df_protonRec.Ptheta
            CorrectedPtheta = const + coeff/df_protonRec.loc[:, "Pp"]/df_protonRec.loc[:, "Pp"] + df_protonRec.loc[:, "Ptheta"]

            const = 0.21192125 -0.0115175 * df_protonRec.Ptheta
            coeff = -8.94307411*0.1 + 1.66349766*0.1 * df_protonRec.Ptheta -8.90617559*0.001 * df_protonRec.Ptheta * df_protonRec.Ptheta + 1.64803754*0.0001 * df_protonRec.Ptheta * df_protonRec.Ptheta * df_protonRec.Ptheta
            CorrectedPphi = const + coeff/df_protonRec.loc[:, "Pp"]/df_protonRec.loc[:, "Pp"] + df_protonRec.loc[:, "Pphi"]

        elif topo == "FD upper":
            const = -3.03346359*10**(-1) + 1.83368163*10**(-2)*df_protonRec.Ptheta - 2.86486404*10**(-4)*df_protonRec.Ptheta*df_protonRec.Ptheta
            coeff =  2.01023276*10**(-1) - 1.13312215*10**(-2)*df_protonRec.Ptheta + 1.82487916*10**(-4)*df_protonRec.Ptheta*df_protonRec.Ptheta
            CorrectedPp = np.select([df_protonRec.Pp<1, df_protonRec.Pp>=1], [const + coeff/df_protonRec.loc[:, "Pp"] + df_protonRec.loc[:, "Pp"], np.exp(-1.2 - 4.228*df_protonRec.Pp) + 0.007502+df_protonRec.Pp])

            const = 2.04334532 * 10 -1.81052405 * df_protonRec.Ptheta + 5.32556360*0.01 * df_protonRec.Ptheta * df_protonRec.Ptheta -5.23157558*0.0001 * df_protonRec.Ptheta * df_protonRec.Ptheta * df_protonRec.Ptheta
            coeff = 8.74233279 -7.63869344 * 0.1 * df_protonRec.Ptheta + 2.22376362*0.01 * df_protonRec.Ptheta * df_protonRec.Ptheta -2.16457260*0.0001 * df_protonRec.Ptheta * df_protonRec.Ptheta * df_protonRec.Ptheta
            CorrectedPtheta = const + coeff/df_protonRec.loc[:, "Pp"]/df_protonRec.loc[:, "Pp"] + df_protonRec.loc[:, "Ptheta"]

            const = 0.54697831 -0.04896981*df_protonRec.Ptheta +  0.00111376*df_protonRec.Ptheta*df_protonRec.Ptheta
            coeff = -4.06733541*10**2 + 2.43696202*10*df_protonRec.Ptheta -3.36144736*10**(-1)*df_protonRec.Ptheta*df_protonRec.Ptheta
            coeff2 = 2.06378660*10 - 1.42866062*df_protonRec.Ptheta + 2.01085440*10**(-2)*df_protonRec.Ptheta*df_protonRec.Ptheta
            CorrectedPphi = const + coeff*np.exp(coeff2*df_protonRec.loc[:, "Pp"]) + df_protonRec.loc[:, "Pphi"]

        elif topo == "CD":
            const = 1.93686914 - 0.116288824*np.minimum(60, df_protonRec.Ptheta) + 0.00223685833*df_protonRec.Ptheta**2 - 1.40771969 * 10**(-5)*df_protonRec.Ptheta**3
            coeff = -0.738047800 + 0.0443343685*np.minimum(60, df_protonRec.Ptheta) - 8.50985972*10**(-4)*df_protonRec.Ptheta*np.minimum(60, df_protonRec.Ptheta) + 5.36810280 * 10**(-6) * df_protonRec.Ptheta**3

            CorrectedPp = const + coeff/df_protonRec.loc[:, "Pp"] + df_protonRec.loc[:, "Pp"]

            const = -1.09849291*100 + 8.86664014 * np.minimum(60, df_protonRec.Ptheta) - 0.26643881 * df_protonRec.Ptheta**2 + 3.53814210 * 10**(-3) * df_protonRec.Ptheta**3 - 1.75297107 * 10**(-5) * df_protonRec.Ptheta**4
            coeff = 9.52034523*100 -5.74808292 * 10 * np.minimum(60, df_protonRec.Ptheta) + 1.15386949 * df_protonRec.Ptheta**2 - 7.57970373 * 0.001 * df_protonRec.Ptheta**3
            coeff2 = -2.00387313*100 + 1.18979079 * 10 * np.minimum(60, df_protonRec.Ptheta) - 2.37730217*0.1 * df_protonRec.Ptheta**2 + 1.55153003*0.001*df_protonRec.Ptheta**3

            CorrectedPtheta = const + coeff*np.exp(coeff2*df_protonRec.loc[:, "Pp"]) + df_protonRec.loc[:, "Ptheta"]

            const = 4.94546178 -3.26662886*0.1 * np.minimum(60, df_protonRec.Ptheta) +  7.39069603 * 0.001 * df_protonRec.Ptheta**2 -6.83599356*10**(-5) * df_protonRec.Ptheta**3 + 2.12303103*10**(-7) * df_protonRec.Ptheta**4
            coeff = 1.72181613*10**(5) -1.36827111*10**(4) * np.minimum(60, df_protonRec.Ptheta) + 4.00923146*10**(2) * df_protonRec.Ptheta**2 - 5.12792347 * df_protonRec.Ptheta**3 + 2.41793167*10**(-2) * df_protonRec.Ptheta**4
            coeff2 =  1.20477219*10**(2) -5.86630228 * np.minimum(60, df_protonRec.Ptheta) + 7.44007875*10**(-2) * df_protonRec.Ptheta**2 -2.42652473*10**(-4) * df_protonRec.Ptheta**3
            CorrectedPphi = const + coeff*np.exp(coeff2*df_protonRec.loc[:, "Pp"]) + df_protonRec.loc[:, "Pphi"]

    #outbending proton energy loss correction
    elif pol == "outbending":
        if topo == "FD lower":
            const = 0.05083242 -0.00469777*df_protonRec.Ptheta + 0.0001082*df_protonRec.Ptheta*df_protonRec.Ptheta
            coeff = -1.47443264*0.01 + 1.58220893*0.001*df_protonRec.Ptheta -3.19490013*0.00001*df_protonRec.Ptheta*df_protonRec.Ptheta
            CorrectedPp = np.select([df_protonRec.Pp<1, df_protonRec.Pp>=1], [const + coeff/df_protonRec.loc[:, "Pp"] + df_protonRec.loc[:, "Pp"], np.exp(-2.739 - 3.932*df_protonRec.Pp) + 0.002907 + df_protonRec.Pp])

            const = -2.56460305*10 + 3.29877542*df_protonRec.Ptheta -1.43106886*0.1*df_protonRec.Ptheta*df_protonRec.Ptheta + 2.08341898*0.001*df_protonRec.Ptheta*df_protonRec.Ptheta*df_protonRec.Ptheta
            coeff =  9.12532740*10 -1.20100762*10*df_protonRec.Ptheta + 5.27654711*0.1*df_protonRec.Ptheta*df_protonRec.Ptheta -7.72656759*0.001*df_protonRec.Ptheta*df_protonRec.Ptheta*df_protonRec.Ptheta
            CorrectedPtheta = const + coeff/df_protonRec.loc[:, "Pp"]/df_protonRec.loc[:, "Pp"] + df_protonRec.loc[:, "Ptheta"]

            const = -20.4780893 + 1.67020488*df_protonRec.Ptheta - 0.03419348*df_protonRec.Ptheta*df_protonRec.Ptheta
            coeff = 35.02807194 - 2.9098043*df_protonRec.Ptheta +  0.06037906*df_protonRec.Ptheta*df_protonRec.Ptheta
            CorrectedPphi = const + coeff/df_protonRec.loc[:, "Pp"]/df_protonRec.loc[:, "Pp"] + df_protonRec.loc[:, "Pphi"]

        elif topo == "FD upper":
            const = 0.09832589 -0.0066463*df_protonRec.Ptheta + 0.00010312*df_protonRec.Ptheta*df_protonRec.Ptheta
            coeff = -9.61421691*0.01 + 6.85638807*0.001*df_protonRec.Ptheta -9.75766427*0.00001*df_protonRec.Ptheta*df_protonRec.Ptheta
            CorrectedPp = np.select([df_protonRec.Pp<1, df_protonRec.Pp>=1], [const + coeff/df_protonRec.loc[:, "Pp"] + df_protonRec.loc[:, "Pp"], np.exp(-1.871 - 3.063*df_protonRec.Pp) + 0.007517 + df_protonRec.Pp])

            const = -1.68873940 + 9.56867163*0.01*df_protonRec.Ptheta -1.43741464*0.001*df_protonRec.Ptheta*df_protonRec.Ptheta
            coeff = 1.49978357*10 -1.40137094*df_protonRec.Ptheta + 4.38501543*0.01*df_protonRec.Ptheta*df_protonRec.Ptheta -4.57982872*0.0001*df_protonRec.Ptheta*df_protonRec.Ptheta*df_protonRec.Ptheta
            CorrectedPtheta = const + coeff/df_protonRec.loc[:, "Pp"]/df_protonRec.loc[:, "Pp"] + df_protonRec.loc[:, "Ptheta"]

            const = 6.75359137 - 0.43199851*df_protonRec.Ptheta + 0.0068995*df_protonRec.Ptheta*df_protonRec.Ptheta
            coeff = -1.68588219 + 1.05609627*0.1*df_protonRec.Ptheta -1.50452832*0.001*df_protonRec.Ptheta*df_protonRec.Ptheta
            CorrectedPphi = const + coeff/df_protonRec.loc[:, "Pp"]/df_protonRec.loc[:, "Pp"] + df_protonRec.loc[:, "Pphi"]

        elif topo == "CD":
            const = 1.92657376 - 0.113836734*np.minimum(60, df_protonRec.Ptheta) + 0.00215038526*np.minimum(60, df_protonRec.Ptheta)**2 - 1.32525053 * 10**(-5)*np.minimum(60, df_protonRec.Ptheta)**3
            coeff = -0.755650043 + 0.0445538936*np.minimum(60, df_protonRec.Ptheta) - 8.38241864*10**(-4)*np.minimum(60, df_protonRec.Ptheta)*np.minimum(60, df_protonRec.Ptheta) + 5.16887255 * 10**(-6) * np.minimum(60, df_protonRec.Ptheta)**3

            CorrectedPp = const + coeff/df_protonRec.loc[:, "Pp"] + df_protonRec.loc[:, "Pp"]

            const = -5.79024055*10 + 4.67197531 * np.minimum(60, df_protonRec.Ptheta) - 0.140156897 * np.minimum(60, df_protonRec.Ptheta)**2 + 1.85853057 * 10**(-3) * np.minimum(60, df_protonRec.Ptheta)**3 - 9.19989908 * 10**(-6) * np.minimum(60, df_protonRec.Ptheta)**4
            coeff = 2.99700765*1000 - 2.18027982 * 10**2 * np.minimum(60, df_protonRec.Ptheta) + 5.84757503 * np.minimum(60, df_protonRec.Ptheta)**2 - 6.80409195 * 0.01 * np.minimum(60, df_protonRec.Ptheta)**3 + 2.89244618 * 0.0001 * np.minimum(60, df_protonRec.Ptheta)**4
            coeff2 = -1.82237904*100 + 1.10153549 * 10 * np.minimum(60, df_protonRec.Ptheta) - 2.24699931*0.1 * np.minimum(60, df_protonRec.Ptheta)**2 + 1.49390960*0.001*np.minimum(60, df_protonRec.Ptheta)**3

            CorrectedPtheta = const + coeff*np.exp(coeff2*df_protonRec.loc[:, "Pp"]) + df_protonRec.loc[:, "Ptheta"]

            const = 7.58761670 - 5.28224578*0.1 * np.minimum(60, df_protonRec.Ptheta) +  1.31580117 * 0.01 * np.minimum(60, df_protonRec.Ptheta)**2 -1.41738951*10**(-4) * np.minimum(60, df_protonRec.Ptheta)**3 + 5.62884363*10**(-7) * np.minimum(60, df_protonRec.Ptheta)**4
            coeff = 1.07644097*10**(5) - 8.67994639*10**(3) * np.minimum(60, df_protonRec.Ptheta) + 2.57187193*10**(2) * np.minimum(60, df_protonRec.Ptheta)**2 - 3.31379317 * np.minimum(60, df_protonRec.Ptheta)**3 + 1.56896621*10**(-2) * np.minimum(60, df_protonRec.Ptheta)**4
            coeff2 =  1.92263184*10**(2) -1.00870704 * 10 * np.minimum(60, df_protonRec.Ptheta) + 1.56575252*10**(-1) * np.minimum(60, df_protonRec.Ptheta)**2 -7.71489734*10**(-4) * np.minimum(60, df_protonRec.Ptheta)**3
            CorrectedPphi = const + coeff*np.exp(coeff2*df_protonRec.loc[:, "Pp"]) + df_protonRec.loc[:, "Pphi"]

    return [CorrectedPp, CorrectedPtheta, CorrectedPphi]