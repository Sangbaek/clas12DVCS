import numpy as np

def dot(vec1, vec2):
    return vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2]
def mag(vec1):
	return np.sqrt(dot(vec1, vec1))
def costheta(vec1, vec2):
    return dot(vec1,vec2)/np.sqrt(dot(vec1,vec1))/np.sqrt(dot(vec2,vec2))
#     return (vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2])/np.sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2])/np.sqrt(vec2[0]*vec2[0]+vec2[1]*vec2[1]+vec2[2]*vec2[2])
def cross(vec1, vec2):
    return [vec1[1]*vec2[2]-vec1[2]*vec2[1], vec1[2]*vec2[0]-vec1[0]*vec2[2], vec1[0]*vec2[1]-vec1[1]*vec2[0]]
def angle(vec1, vec2):
    return 180/np.pi*np.arccos(costheta(vec1, vec2))