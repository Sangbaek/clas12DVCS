#!/usr/bin/env python3
"""
Modules help pandas algebra without using ROOT.
"""
import numpy as np
from utils.const import *

def dot(vec1, vec2):
	# dot product of two 3d vectors
    return vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2]

def mag(vec1):
	# L2 norm of vector
	return np.sqrt(dot(vec1, vec1))

def mag2(vec1):
	# square of L2 norm
	return	dot(vec1, vec1)

def cosTheta(vec1, vec2):
	# cosine angle between two 3d vectors
    return dot(vec1,vec2)/np.sqrt(mag2(vec1) * mag2(vec2))

def angle(vec1, vec2):
	# angle between two 3d vectors
	return 180/np.pi*np.arccos(np.minimum(1, cosTheta(vec1, vec2)))

def cross(vec1, vec2):
	# cross product of two 3d vectors
    return [vec1[1]*vec2[2]-vec1[2]*vec2[1], vec1[2]*vec2[0]-vec1[0]*vec2[2], vec1[0]*vec2[1]-vec1[1]*vec2[0]]

def vecAdd(gam1, gam2):
	# add two 3d vectors
	return [gam1[0]+gam2[0], gam1[1]+gam2[1], gam1[2]+gam2[2]]

def pi0Energy(gam1, gam2):
	# reconstructed pi0 energy of two 3d photon momenta
	return mag(gam1)+mag(gam2)

def pi0InvMass(gam1, gam2):
	# pi0 invariant mass of two 3d photon momenta
	pi0mass2 = pi0Energy(gam1, gam2)**2-mag2(vecAdd(gam1, gam2))
	pi0mass2 = np.where(pi0mass2 >= 0, pi0mass2, 10**6)
	pi0mass = np.sqrt(pi0mass2)
	pi0mass = np.where(pi0mass > 100, -1000, pi0mass)
	return pi0mass

def getPhi(vec1):
	# azimuthal angle of one 3d vector
	return 180/np.pi*np.arctan2(vec1[1], vec1[0])

def getTheta(vec1):
	# polar angle of one 3d vector
	return 180/np.pi*np.arctan2(np.sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]), vec1[2])

def getEnergy(vec1, mass):
	# for taken 3d momenta p and mass m, return energy = sqrt(p**2 + m**2)
	return np.sqrt(mag2(vec1)+mass**2)

def nu(xB, Q2, t, phi):
    return Q2/(2*M*xB)

def y(xB, Q2, t, phi):
    return nu(xB, Q2, t, phi)/10.604

def xi(xB, Q2, t, phi):
    return xB*(1-t/2/Q2)/(2-xB+xB*(-t/Q2))

def del2(xB, Q2, t, phi):
    return -t

def del2q2(xB, Q2, t, phi):
    return del2(xB, Q2, t, phi)/Q2

def eps(xB, Q2, t, phi):
    return 2*xB*M/np.sqrt(Q2)

def eps2(xB, Q2, t, phi):
    return eps(xB, Q2, t, phi)**2

def qeps2(xB, Q2, t, phi):
    return 1+eps2(xB, Q2, t, phi)

def sqeps2(xB, Q2, t, phi):
    return np.sqrt(qeps2(xB, Q2, t, phi))

def y1eps(xB, Q2, t, phi):
    return 1 - y(xB, Q2, t, phi) - y(xB, Q2, t, phi)*y(xB, Q2, t, phi)*eps2(xB, Q2, t, phi)/4

def tmin(xB, Q2, t, phi):
    return -Q2*(2*(1-xB)*(1-sqeps2(xB, Q2, t, phi))+eps2(xB, Q2, t, phi))/(4*xB*(1-xB)+eps2(xB, Q2, t, phi))

def tmax(xB, Q2, t, phi):
    return -Q2*(2*(1-xB)*(1+sqeps2(xB, Q2, t, phi))+eps2(xB, Q2, t, phi))/(4*xB*(1-xB)+eps2(xB, Q2, t, phi))

def W2(xB, Q2, t, phi):
    return M*M+2.0*M*nu(xB, Q2, t, phi)-Q2

def W(xB, Q2, t, phi):
    return np.sqrt(W2(xB, Q2, t, phi))

def tmin2(xB, Q2, t, phi):
    return -0.5*((Q2/xB-Q2)*(Q2/xB-np.sqrt((Q2/xB)**2+4*M*M*Q2))+2*M*M*Q2)/W2(xB, Q2, t, phi)

def tmax2(xB, Q2, t, phi):
    return -0.5*((Q2/xB-Q2)*(Q2/xB+np.sqrt((Q2/xB)**2+4*M*M*Q2))+2*M*M*Q2)/W2(xB, Q2, t, phi)

def Kfac2(xB, Q2, t, phi):
    return (-del2q2(xB, Q2, t, phi))*(1 - xB)*y1eps(xB, Q2, t, phi)*(1 - np.abs(tmin(xB, Q2, t, phi))/t)*(np.sqrt(1 + eps2(xB, Q2, t, phi)) + 
            ((4*xB*(1 - xB) + eps2(xB, Q2, t, phi))/(4*(1 - xB)))*(-(t - np.abs(tmin(xB, Q2, t, phi)))/Q2))
def Kfac(xB, Q2, t, phi):
    return np.sqrt(Kfac2(xB, Q2, t, phi))

def Jfac(xB, Q2, t, phi):
    return (1 - y(xB, Q2, t, phi) - y(xB, Q2, t, phi)*eps2(xB, Q2, t, phi)/2)*(1 + del2q2(xB, Q2, t, phi)) - (1 - xB)*(2 - y(xB, Q2, t, phi))*del2q2(xB, Q2, t, phi)

def P1(xB, Q2, t, phi):    
    return -(Jfac(xB, Q2, t, phi) + 2*Kfac(xB, Q2, t, phi)*np.cos(np.pi-np.radians(phi)))/(y(xB, Q2, t, phi)*(1 + eps2(xB, Q2, t, phi)))

def P2(xB, Q2, t, phi):
    return 1 + del2q2(xB, Q2, t, phi) - P1(xB, Q2, t, phi)