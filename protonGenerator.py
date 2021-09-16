import numpy as np
import argparse

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument("-o","--out", help="a single pickle file name as an output", default="proton.dat")
	parser.add_argument("-n","--entry", help="entry_stop to stop reading the root file", default = "10000")

	args = parser.parse_args()

	file_out = open(args.out, "w")

	for i in range(int(args.entry)):
		p = 0.5
		costheta = np.random.uniform(np.cos(40.0*np.pi/180.0), np.cos(15.0*np.pi/180.0))
		sintheta = np.sqrt(1-costheta*costheta) 
		phi = np.random.uniform(-np.pi, np.pi)
		px = p * sintheta * np.cos(phi)
		py = p * sintheta * np.sin(phi)
		pz = p * costheta
		file_out.write("2   1   1    0.0   0.0 11   10.604   1       1      0.2373006E-02\n")
		file_out.write("1  -1.  1     11   0    0    1.3862    0.3857    4.2996    4.5334    0.0005    0.0000    0.0000    -3.0000\n")
		file_out.write("2   1.  1   2212   0    0   {0:.4f}   {1:.4f}    {2:.4f}    0.5455    0.9380    0.0000    0.0000    -3.0000\n".format(px, py, pz))