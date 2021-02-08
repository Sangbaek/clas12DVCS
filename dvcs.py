import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#from lmfit import Model, Parameter, report_fit
from scipy.optimize import curve_fit
from copy import copy
from utils.algebra_pd import costheta, dot, cross, angle

file = uproot.open("../IgorConvertsHipotoROOT/root/dvcs_RGA_IN_filtered_with_TOF.root")

M = 0.938272081
me = 0.5109989461 *0.001
ebeam = 10.604
pbeam = np.sqrt(ebeam*ebeam-me*me)

tree = file["T;42"]

#initialize df
df_electron = pd.DataFrame()
df_proton = pd.DataFrame()
df_gamma = pd.DataFrame()

df_electron["Epx"] = tree["Epx"].array(library="pd", entry_stop=100)
df_electron["Epy"] = tree["Epy"].array(library="pd", entry_stop=100)
df_electron["Epz"] = tree["Epz"].array(library="pd", entry_stop=100)

df_proton["Ppx"] = tree["Ppx"].array(library="pd", entry_stop=100)
df_proton["Ppy"] = tree["Ppy"].array(library="pd", entry_stop=100)
df_proton["Ppz"] = tree["Ppz"].array(library="pd", entry_stop=100)
df_proton["Pstat"] = tree["Pstat"].array(library="pd", entry_stop=100)

df_gamma["Gpx"] = tree["Gpx"].array(library="pd", entry_stop=100)
df_gamma["Gpy"] = tree["Gpy"].array(library="pd", entry_stop=100)
df_gamma["Gpz"] = tree["Gpz"].array(library="pd", entry_stop=100)
df_proton["Gstat"] = tree["Gstat"].array(library="pd", entry_stop=100)

df_electron['event']=df_electron.index.get_level_values('entry')
df_proton['event']=df_proton.index.get_level_values('entry')
df_gamma['event']=df_gamma.index.get_level_values('entry')

df_pg=pd.merge(df_proton,df_gamma, how='outer', on='event')
df_epg=pd.merge(df_electron, df_pg, how='outer', on='event')

df_epg=df_epg.astype({"Epx":float, "Epy":float, "Epz":float, "Ppx":float, "Ppy":float, "Ppz":float, "Gpx":float, "Gpy":float, "Gpz":float})

beam = [0, 0, pbeam]
target = [0, 0, 0]
ele = [df_epg['Epx'], df_epg['Epy'], df_epg['Epz']]
pro = [df_epg['Ppx'], df_epg['Ppy'], df_epg['Ppz']]
gam = [df_epg['Gpx'], df_epg['Gpy'], df_epg['Gpz']]

df_epg['Ee'] = np.sqrt(me**2+df_epg["Epx"]**2+df_epg["Epy"]**2+df_epg["Epz"]**2)
df_epg['Pe'] = np.sqrt(M**2+df_epg["Ppx"]**2+df_epg["Ppy"]**2+df_epg["Ppz"]**2)
df_epg['Ge'] = np.sqrt(df_epg["Gpx"]**2+df_epg["Gpy"]**2+df_epg["Gpz"]**2)

#binning kinematics
df_epg['Q2'] = -((ebeam-df_epg['Ee'])**2-df_epg['Epx']**2-df_epg['Epy']**2-(pbeam-df_epg['Epz'])**2)
df_epg['xB'] = df_epg['Q2']/2.0/M/(ebeam-df_epg['Ee'])
df_epg['nu'] = (ebeam-df_epg['Ee'])
VGS = [-df_epg['Epx'], -df_epg['Epy'], pbeam - df_epg['Epz']]
df_epg['t1'] = 2*M*(df_epg['Pe']-M)
costheta = costheta(VGS, gam)
df_epg['t2'] = (M*df_epg['Q2']+2*M*df_epg['nu']*(df_epg['nu']-np.sqrt(df_epg['nu']*df_epg['nu']+df_epg['Q2'])*costheta))/(M+df_epg['nu']-np.sqrt(df_epg['nu']*df_epg['nu']+df_epg['Q2'])*costheta)
df_epg['W'] = np.sqrt((ebeam+M-df_epg['Ee'])**2-df_epg['Epx']**2-df_epg['Epy']**2-(pbeam-df_epg['Epz'])**2)

v3l = cross(beam, ele)
v3h = cross(pro, VGS)
v3g = cross(VGS, gam)

df_epg['phi1'] = angle(v3l, v3h)
df_epg['phi1'] = np.where(dot(v3l, pro)>0, 360.0 - df_epg['phi1'], df_epg['phi1'])
df_epg['phi2'] = angle(v3l, v3g)
df_epg['phi2'] = np.where(dot(VGS, cross(v3l, v3g))<0, 360.0 - df_epg['phi2'], df_epg['phi2'])

#exclusivity variables
df_epg['MM2_epg'] = (-M-ebeam+df_epg["Ee"]+df_epg["Pe"]+df_epg["Ge"])**2 - (df_epg["Epx"]+df_epg["Ppx"]+df_epg["Gpx"])**2 - \
(df_epg["Epy"]+df_epg["Ppy"]+df_epg["Gpy"])**2 -(-pbeam+df_epg["Epz"]+df_epg["Ppz"]+df_epg["Gpz"]) **2
df_epg['ME_epg'] = (M+ebeam-df_epg["Ee"]-df_epg["Pe"]-df_epg["Ge"])
df_epg['MM2_ep'] = (-M-ebeam+df_epg["Ee"]+df_epg["Pe"])**2 - (df_epg["Epx"]+df_epg["Ppx"])**2 - \
(df_epg["Epy"]+df_epg["Ppy"])**2 -(-pbeam+df_epg["Epz"]+df_epg["Ppz"]) **2
df_epg['MM2_eg'] = (-M-ebeam+df_epg["Ee"]+df_epg["Ge"])**2 - (df_epg["Epx"]+df_epg["Gpx"])**2 - \
(df_epg["Epy"]+df_epg["Gpy"])**2 -(-pbeam+df_epg["Epz"]+df_epg["Gpz"]) **2
df_epg['MPt'] = np.sqrt((df_epg["Epx"]+df_epg["Ppx"]+df_epg["Gpx"])**2 + \
(df_epg["Epy"]+df_epg["Ppy"]+df_epg["Gpy"])**2)
df_epg['coneAngle'] = angle(ele, gam)#180.*/np.pi*np.arccos((df_epg["Epx"]*df_epg["Gpx"]+df_epg["Epy"]*df_epg["Gpy"]+df_epg["Epz"]*df_epg["Gpz"])/np.sqrt(df_epg["Epx"]*df_epg["Epx"]+df_epg["Epy"]*df_epg["Epy"]+df_epg["Epz"]*df_epg["Epz"])/np.sqrt(df_epg["Gpx"]*df_epg["Gpx"]+df_epg["Gpy"]*df_epg["Gpy"]+df_epg["Gpz"]*df_epg["Gpz"]))
VmissG = [-df_epg["Epx"]-df_epg["Ppx"], -df_epg["Epy"]-df_epg["Ppy"], pbeam-df_epg["Epz"]-df_epg["Ppz"]]
df_epg['reconGam'] = angle(gam, VmissG)
Vhadr = cross(pro, VGS)
Vhad2 = cross(VGS, gam)
df_epg['coplanarity'] = angle(Vhadr, Vhad2)

#cuts
#Kinematics
#      xB>0 && xB<1 && W>2 && Q2>1 && VE.e()>2 && VG1.e() > 2 && VP.p() > 0.12
cut_forward = df_epg["Pstat"]<4000
cut_central = df_epg["Pstat"]>4000
cut_xBupper = df_epg["xB"]<1
cut_xBlower = df_epg["xB"]>0
cut_Q2 = df_epg["Q2"]>1
cut_W = df_epg["W"]>2
cut_Ee = df_epg["Ee"]>2
cut_Ge = df_epg["Ge"]>2
cut_Pp = np.sqrt(dot([df_epg["Ppx"],df_epg["Ppy"],df_epg["Ppz"]],[df_epg["Ppx"],df_epg["Ppy"],df_epg["Ppz"]]))>0.12

#Exclusive
#mmepg
cut_mmepg = np.abs(df_epg["MM2_epg"])<0.04
# #mmeg
df_epg = df_epg[df_epg["MM2_eg"]>0]
cut_mmegupper = np.sqrt(df_epg["MM2_eg"])<1.7
cut_mmeglower = np.sqrt(df_epg["MM2_eg"])>0.1
#meepg
cut_meepgupper = df_epg["ME_epg"]<1.2
cut_meepglower = df_epg["ME_epg"]>-0.5
#mpt
cut_mpt = df_epg["MPt"]<0.12
#coneangle
cut_cone = df_epg["coneAngle"]>10
#recon gam angle
cut_recon = df_epg["reconGam"]<1.1
#cut on coplanarity
#cut on mmep

df_dvcs = df_epg[cut_xBupper&cut_xBlower&cut_Q2&cut_W&cut_Ee&cut_Ge&cut_Pp &cut_mmepg & cut_mmegupper & cut_mmeglower & cut_meepgupper & cut_meepglower & cut_mpt & cut_cone & cut_recon]
df_dvcsForward = df_epg[cut_forward&cut_xBupper&cut_xBlower&cut_Q2&cut_W&cut_Ee&cut_Ge&cut_Pp &cut_mmepg & cut_mmegupper & cut_mmeglower & cut_meepgupper & cut_meepglower & cut_mpt & cut_cone & cut_recon]
df_dvcsCentral = df_epg[cut_central&cut_xBupper&cut_xBlower&cut_Q2&cut_W&cut_Ee&cut_Ge&cut_Pp &cut_mmepg & cut_mmegupper & cut_mmeglower & cut_meepgupper & cut_meepglower & cut_mpt & cut_cone & cut_recon]

xB_edges = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.85, 1]
Q2_edges = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 12.]
t_edges = [0, 0.15, 0.25, 0.45, 0.80, 1.15,4.0]
phi_edges = np.arange(0,360.1,30)

# fig, axs = plt.subplots(len(Q2_edges)-1, len(xB_edges)-1)
# for xBind in range(0,len(xB_edges)-1):
# 	for Q2ind in range(0, len(Q2_edges)-1):
# 		binQ2lower = df_dvcsCentral["Q2"]>Q2_edges[Q2ind]
# 		binQ2upper = df_dvcsCentral["Q2"]<Q2_edges[Q2ind+1]
# 		binxBlower = df_dvcsCentral["xB"]>xB_edges[xBind]
# 		binxBupper = df_dvcsCentral["xB"]<xB_edges[xBind+1]
# 		plotter = df_dvcsCentral[binQ2lower&binQ2upper&binxBlower&binxBupper]
# 		# plotter.hist(column="Q2", bins=Q2_edges, ax=axs[len(Q2_edges)-2-Q2ind,xBind], histtype='stepfilled', facecolor='none', edgecolor='k')
# 		# axs[len(Q2_edges)-2-Q2ind,xBind].hist2d(plotter["phi2"],plotter["t2"],bins=[phi_edges, t_edges], cmap=plt.cm.Greys)
# 		axs[len(Q2_edges)-2-Q2ind,xBind].hist(plotter["phi2"],bins=phi_edges, histtype='stepfilled', facecolor='none', edgecolor='k')
# 		axs[len(Q2_edges)-2-Q2ind,xBind].set_title("")
# 		axs[len(Q2_edges)-2-Q2ind,xBind].xaxis.set_visible(False)
# 		axs[len(Q2_edges)-2-Q2ind,xBind].yaxis.set_visible(False)
# 		# axs[len(Q2_edges)-2-Q2ind,xBind].set_xticks([])
# 		# axs[len(Q2_edges)-2-Q2ind,xBind].set_yticks([])
# 		# axs[len(Q2_edges)-2-Q2ind,xBind].set_xticklabels([])
# 		# axs[len(Q2_edges)-2-Q2ind,xBind].set_yticklabels([])
# 		axs[len(Q2_edges)-2-Q2ind,xBind].set_ylim([1,10000])
# 		axs[len(Q2_edges)-2-Q2ind,xBind].set_yscale('log')
# 		plt.setp( axs[len(Q2_edges)-2-Q2ind,xBind].get_xticklabels(), visible=False)
# 		plt.setp( axs[len(Q2_edges)-2-Q2ind,xBind].get_yticklabels(), visible=False)
# plt.subplots_adjust(wspace=0, hspace=0)
# plt.tick_params(axis='both', labelsize=0, length = 0)
# plt.savefig("CD_phi.pdf")
# # plt.clf()

# fig, axs = plt.subplots(len(Q2_edges)-1, len(xB_edges)-1)
# for xBind in range(0,len(xB_edges)-1):
# 	for Q2ind in range(0, len(Q2_edges)-1):
# 		binQ2lower = df_dvcsForward["Q2"]>Q2_edges[Q2ind]
# 		binQ2upper = df_dvcsForward["Q2"]<Q2_edges[Q2ind+1]
# 		binxBlower = df_dvcsForward["xB"]>xB_edges[xBind]
# 		binxBupper = df_dvcsForward["xB"]<xB_edges[xBind+1]
# 		plotter = df_dvcsForward[binQ2lower&binQ2upper&binxBlower&binxBupper]
# 		# plotter.hist(column="Q2", bins=Q2_edges, ax=axs[len(Q2_edges)-2-Q2ind,xBind], histtype='stepfilled', facecolor='none', edgecolor='k')
# 		# axs[len(Q2_edges)-2-Q2ind,xBind].hist2d(plotter["phi2"],plotter["t2"],bins=[phi_edges, t_edges], cmap=plt.cm.Greys)
# 		axs[len(Q2_edges)-2-Q2ind,xBind].hist(plotter["phi2"],bins=phi_edges, histtype='stepfilled', facecolor='none', edgecolor='k')
# 		axs[len(Q2_edges)-2-Q2ind,xBind].set_title("")
# 		axs[len(Q2_edges)-2-Q2ind,xBind].xaxis.set_visible(False)
# 		axs[len(Q2_edges)-2-Q2ind,xBind].yaxis.set_visible(False)
# 		axs[len(Q2_edges)-2-Q2ind,xBind].set_ylim([1,10000])
# 		axs[len(Q2_edges)-2-Q2ind,xBind].set_yscale('log')
# 		plt.setp( axs[len(Q2_edges)-2-Q2ind,xBind].get_xticklabels(), visible=False)
# 		plt.setp( axs[len(Q2_edges)-2-Q2ind,xBind].get_yticklabels(), visible=False)
# plt.subplots_adjust(wspace=0, hspace=0)
# plt.savefig("FD_phi.pdf")
# plt.clf()

# cmap = copy(plt.cm.get_cmap("jet"))
# # cmap.set_under('w',0)

# fig, axs = plt.subplots(len(Q2_edges)-1, len(xB_edges)-1)
# for xBind in range(0,len(xB_edges)-1):
# 	for Q2ind in range(0, len(Q2_edges)-1):
# 		binQ2lower = df_dvcsCentral["Q2"]>Q2_edges[Q2ind]
# 		binQ2upper = df_dvcsCentral["Q2"]<Q2_edges[Q2ind+1]
# 		binxBlower = df_dvcsCentral["xB"]>xB_edges[xBind]
# 		binxBupper = df_dvcsCentral["xB"]<xB_edges[xBind+1]
# 		plotter = df_dvcsCentral[binQ2lower&binQ2upper&binxBlower&binxBupper]
# 		axs[len(Q2_edges)-2-Q2ind,xBind].hist2d(plotter["phi2"],plotter["t2"],bins=[phi_edges, t_edges], cmap=cmap, cmin=1)
# 		axs[len(Q2_edges)-2-Q2ind,xBind].set_title("")
# 		axs[len(Q2_edges)-2-Q2ind,xBind].xaxis.set_visible(False)
# 		axs[len(Q2_edges)-2-Q2ind,xBind].yaxis.set_visible(False)
# 		plt.setp( axs[len(Q2_edges)-2-Q2ind,xBind].get_xticklabels(), visible=False)
# 		plt.setp( axs[len(Q2_edges)-2-Q2ind,xBind].get_yticklabels(), visible=False)
# plt.subplots_adjust(wspace=0, hspace=0)
# plt.savefig("CD_tphi.pdf")
# plt.clf()

# fig, axs = plt.subplots(len(Q2_edges)-1, len(xB_edges)-1)
# for xBind in range(0,len(xB_edges)-1):
# 	for Q2ind in range(0, len(Q2_edges)-1):
# 		binQ2lower = df_dvcsForward["Q2"]>Q2_edges[Q2ind]
# 		binQ2upper = df_dvcsForward["Q2"]<Q2_edges[Q2ind+1]
# 		binxBlower = df_dvcsForward["xB"]>xB_edges[xBind]
# 		binxBupper = df_dvcsForward["xB"]<xB_edges[xBind+1]
# 		plotter = df_dvcsForward[binQ2lower&binQ2upper&binxBlower&binxBupper]
# 		axs[len(Q2_edges)-2-Q2ind,xBind].hist2d(plotter["phi2"],plotter["t2"],bins=[phi_edges, t_edges], cmap=cmap, cmin=1)
# 		axs[len(Q2_edges)-2-Q2ind,xBind].set_title("")
# 		axs[len(Q2_edges)-2-Q2ind,xBind].xaxis.set_visible(False)
# 		axs[len(Q2_edges)-2-Q2ind,xBind].yaxis.set_visible(False)
# 		plt.setp( axs[len(Q2_edges)-2-Q2ind,xBind].get_xticklabels(), visible=False)
# 		plt.setp( axs[len(Q2_edges)-2-Q2ind,xBind].get_yticklabels(), visible=False)
# plt.subplots_adjust(wspace=0, hspace=0)
# plt.savefig("FD_tphi.pdf")
# plt.clf()



def dsigma(phi,A,B,C,D,E):
	return A+B*np.cos(phi*np.pi/180.)+C*np.sin(phi*np.pi/180.)+D*np.cos(2*phi*np.pi/180.)+E*np.cos(2*phi*np.pi/180.)

for tind in range(0, len(t_edges)-1):
	fig, axs = plt.subplots(len(Q2_edges)-1, len(xB_edges)-1)
	for xBind in range(0,len(xB_edges)-1):
		for Q2ind in range(0, len(Q2_edges)-1):
			binQ2lower = df_dvcsCentral["Q2"]>Q2_edges[Q2ind]
			binQ2upper = df_dvcsCentral["Q2"]<Q2_edges[Q2ind+1]
			binxBlower = df_dvcsCentral["xB"]>xB_edges[xBind]
			binxBupper = df_dvcsCentral["xB"]<xB_edges[xBind+1]
			bintlower = df_dvcsCentral["t2"]>t_edges[tind]
			bintupper = df_dvcsCentral["t2"]<t_edges[tind+1]
			plotter = df_dvcsCentral[binQ2lower&binQ2upper&binxBlower&binxBupper&bintlower&bintupper]

			# plotter.hist(column="Q2", bins=Q2_edges, ax=axs[len(Q2_edges)-2-Q2ind,xBind], histtype='stepfilled', facecolor='none', edgecolor='k')
			# axs[len(Q2_edges)-2-Q2ind,xBind].hist2d(plotter["phi2"],plotter["t2"],bins=[phi_edges, t_edges], cmap=plt.cm.Greys)
			axs[len(Q2_edges)-2-Q2ind,xBind].hist(plotter["phi2"],bins=phi_edges, linewidth = 0.5, histtype='stepfilled', facecolor='none', edgecolor='k')
			axs[len(Q2_edges)-2-Q2ind,xBind].set_title("")
			axs[len(Q2_edges)-2-Q2ind,xBind].xaxis.set_visible(False)
			axs[len(Q2_edges)-2-Q2ind,xBind].yaxis.set_visible(False)
			
			data_entries, bins = np.histogram(plotter["phi2"], bins=phi_edges)
			binscenters = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])
			popt, pcov = curve_fit(dsigma, xdata=binscenters, ydata=data_entries, p0=[100.0, 1.0, 1.0, 1.0, 1.0])
			print(popt)
			phispace = np.arange(0,360.1,1)
			axs[len(Q2_edges)-2-Q2ind,xBind].plot(phispace, dsigma(phispace, *popt), color='red', linewidth=1.0, label=r'Fitted function')


			# axs[len(Q2_edges)-2-Q2ind,xBind].set_xticks([])
			# axs[len(Q2_edges)-2-Q2ind,xBind].set_yticks([])
			# axs[len(Q2_edges)-2-Q2ind,xBind].set_xticklabels([])
			# axs[len(Q2_edges)-2-Q2ind,xBind].set_yticklabels([])
			axs[len(Q2_edges)-2-Q2ind,xBind].set_ylim([1,5000])
			axs[len(Q2_edges)-2-Q2ind,xBind].set_yscale('log')
			plt.setp( axs[len(Q2_edges)-2-Q2ind,xBind].get_xticklabels(), visible=False)
			plt.setp( axs[len(Q2_edges)-2-Q2ind,xBind].get_yticklabels(), visible=False)
	plt.subplots_adjust(wspace=0, hspace=0)
	plt.tick_params(axis='both', labelsize=0, length = 0)
	plt.savefig("central_{}.pdf".format(tind))
	plt.clf()

for tind in range(0, len(t_edges)-1):
	fig, axs = plt.subplots(len(Q2_edges)-1, len(xB_edges)-1)
	for xBind in range(0,len(xB_edges)-1):
		for Q2ind in range(0, len(Q2_edges)-1):
			binQ2lower = df_dvcsForward["Q2"]>Q2_edges[Q2ind]
			binQ2upper = df_dvcsForward["Q2"]<Q2_edges[Q2ind+1]
			binxBlower = df_dvcsForward["xB"]>xB_edges[xBind]
			binxBupper = df_dvcsForward["xB"]<xB_edges[xBind+1]
			bintlower = df_dvcsForward["t2"]>t_edges[tind]
			bintupper = df_dvcsForward["t2"]<t_edges[tind+1]
			plotter = df_dvcsForward[binQ2lower&binQ2upper&binxBlower&binxBupper&bintlower&bintupper]

			# plotter.hist(column="Q2", bins=Q2_edges, ax=axs[len(Q2_edges)-2-Q2ind,xBind], histtype='stepfilled', facecolor='none', edgecolor='k')
			# axs[len(Q2_edges)-2-Q2ind,xBind].hist2d(plotter["phi2"],plotter["t2"],bins=[phi_edges, t_edges], cmap=plt.cm.Greys)
			axs[len(Q2_edges)-2-Q2ind,xBind].hist(plotter["phi2"],bins=phi_edges, linewidth=0.5, histtype='stepfilled', facecolor='none', edgecolor='k')
			axs[len(Q2_edges)-2-Q2ind,xBind].set_title("")
			axs[len(Q2_edges)-2-Q2ind,xBind].xaxis.set_visible(False)
			axs[len(Q2_edges)-2-Q2ind,xBind].yaxis.set_visible(False)
			
			data_entries, bins = np.histogram(plotter["phi2"], bins=phi_edges)
			binscenters = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])
			popt, pcov = curve_fit(dsigma, xdata=binscenters, ydata=data_entries, p0=[100.0, 1.0, 1.0, 1.0, 1.0])
			print(popt)
			phispace = np.arange(0,360.1,1)
			axs[len(Q2_edges)-2-Q2ind,xBind].plot(phispace, dsigma(phispace, *popt), color='red', linewidth=1.0, label=r'Fitted function')


			# axs[len(Q2_edges)-2-Q2ind,xBind].set_xticks([])
			# axs[len(Q2_edges)-2-Q2ind,xBind].set_yticks([])
			# axs[len(Q2_edges)-2-Q2ind,xBind].set_xticklabels([])
			# axs[len(Q2_edges)-2-Q2ind,xBind].set_yticklabels([])
			axs[len(Q2_edges)-2-Q2ind,xBind].set_ylim([1,5000])
			axs[len(Q2_edges)-2-Q2ind,xBind].set_yscale('log')
			plt.setp( axs[len(Q2_edges)-2-Q2ind,xBind].get_xticklabels(), visible=False)
			plt.setp( axs[len(Q2_edges)-2-Q2ind,xBind].get_yticklabels(), visible=False)
	plt.subplots_adjust(wspace=0, hspace=0)
	plt.tick_params(axis='both', labelsize=0, length = 0)
	plt.savefig("forward_{}.pdf".format(tind))
	plt.clf()

for tind in range(0, len(t_edges)-1):
	fig, axs = plt.subplots(len(Q2_edges)-1, len(xB_edges)-1)
	for xBind in range(0,len(xB_edges)-1):
		for Q2ind in range(0, len(Q2_edges)-1):
			binQ2lower = df_dvcs["Q2"]>Q2_edges[Q2ind]
			binQ2upper = df_dvcs["Q2"]<Q2_edges[Q2ind+1]
			binxBlower = df_dvcs["xB"]>xB_edges[xBind]
			binxBupper = df_dvcs["xB"]<xB_edges[xBind+1]
			bintlower = df_dvcs["t2"]>t_edges[tind]
			bintupper = df_dvcs["t2"]<t_edges[tind+1]
			plotter = df_dvcs[binQ2lower&binQ2upper&binxBlower&binxBupper&bintlower&bintupper]

			# plotter.hist(column="Q2", bins=Q2_edges, ax=axs[len(Q2_edges)-2-Q2ind,xBind], histtype='stepfilled', facecolor='none', edgecolor='k')
			# axs[len(Q2_edges)-2-Q2ind,xBind].hist2d(plotter["phi2"],plotter["t2"],bins=[phi_edges, t_edges], cmap=plt.cm.Greys)
			axs[len(Q2_edges)-2-Q2ind,xBind].hist(plotter["phi2"],bins=phi_edges, linewidth=0.5, histtype='stepfilled', facecolor='none', edgecolor='k')
			axs[len(Q2_edges)-2-Q2ind,xBind].set_title("")
			axs[len(Q2_edges)-2-Q2ind,xBind].xaxis.set_visible(False)
			axs[len(Q2_edges)-2-Q2ind,xBind].yaxis.set_visible(False)
			
			data_entries, bins = np.histogram(plotter["phi2"], bins=phi_edges)
			binscenters = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])
			popt, pcov = curve_fit(dsigma, xdata=binscenters, ydata=data_entries, p0=[100.0, 1.0, 1.0, 1.0, 1.0])
			print(popt)
			phispace = np.arange(0,360.1,1)
			axs[len(Q2_edges)-2-Q2ind,xBind].plot(phispace, dsigma(phispace, *popt), color='red', linewidth=1.0, label=r'Fitted function')


			# axs[len(Q2_edges)-2-Q2ind,xBind].set_xticks([])
			# axs[len(Q2_edges)-2-Q2ind,xBind].set_yticks([])
			# axs[len(Q2_edges)-2-Q2ind,xBind].set_xticklabels([])
			# axs[len(Q2_edges)-2-Q2ind,xBind].set_yticklabels([])
			axs[len(Q2_edges)-2-Q2ind,xBind].set_ylim([1,5000])
			axs[len(Q2_edges)-2-Q2ind,xBind].set_yscale('log')
			plt.setp( axs[len(Q2_edges)-2-Q2ind,xBind].get_xticklabels(), visible=False)
			plt.setp( axs[len(Q2_edges)-2-Q2ind,xBind].get_yticklabels(), visible=False)
	plt.subplots_adjust(wspace=0, hspace=0)
	plt.tick_params(axis='both', labelsize=0, length = 0)
	plt.savefig("dvcs_{}.pdf".format(tind))
	plt.clf()