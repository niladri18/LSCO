#!/usr/bin/env python
import os
import sys
import pwd
import glob
import pdb
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import matplotlib.cm as cm
from matplotlib import rc
from matplotlib.colors import LogNorm

def convert(z):
	N = int(len(z)**.5)
	z = np.array(z)
	z =  z.reshape(N,N)
	
	return z

def plot_density(x,y,z,U,s):
	font = {'family':'serif'}
	rc('font',**font)
	dic = {}
	dic["Ehi"] = "E1"
	dic["Elow"] = "E0"
	dic["Ediff"] = r"$\Delta E$"
	plt.title(r"$t_{\sigma}$=1.0"+ "\t"+r"$t_{\pi}$=0.2"+"\t" + "U="+U +"\t"+ dic[s])
	#plt.imshow(z, extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)), norm = LogNorm() )
	plt.xlabel(r"k$_{x}$")
	plt.ylabel(r"k$_{y}$")
	#plt.ylim(-0.2,1.0)
	extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y))
	#levels = -4.0
	plt.imshow(z, extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)) )
	plt.colorbar()
	#plt.contour(z, levels = np.array([-0.01,0.01]), colors = 'white', extent = extent)
	plt.contour(z, levels = np.array([0.0]), colors = 'k', extent = extent)
	plt.show()
	pass


def find_u(fname):
	param = fname.split("_")
	U = param[-1][1:-4]
	return U


def find_gap(fname):
	f0 = open(fname,'r')
	U = find_u(fname)
	x = []
	y = []
	Ediff = []
	while True:
		line = f0.readline().split()
		if len(line) < 1:
			break
		t = [float(i) for i in line]
		x.append(t[0])
		y.append(t[1])
		gap = t[-1]-t[-3]
		Ediff.append(gap)

	ans = [float(U), min(Ediff), max(Ediff)]
	return ans

def sort_file(file_list):
        u = {}
        tmp_list = []
        for i in file_list:
                param = i.split("_")
                x = float( param[-1][1:-4] )
                u[x] = i

        print(sorted(u.keys()))
        for i in sorted(u.keys()):
                tmp_list.append(u[i])
        return tmp_list


def plot_data(gap_data):
        gap_data = np.array(gap_data)
        color_code = ['red', 'green']
	linestyle = ['--',':']
	label = [r"$\pi,\pi$",r"$\pi,0$"]
        plt.xlabel("U", fontsize = 20)
        plt.ylabel(r"$\Delta$E",fontsize = 20)
	plt.ylim(-0.1,4.0)
	#plt.legend(loc = 'upper right')
	#plt.legend()
        for x in range(1,3):
                plt.plot( gap_data[:,0],gap_data[:,x],color_code[x-1], linestyle = linestyle[x-1], linewidth = 3, label = label[x-1])
	plt.legend()
        plt.show()


'''
MAIN PROGRAM
'''

file_list = glob.glob('band_H0.00*')
file_list = sort_file(file_list)
print(file_list)
gap_data = []

for i in file_list:
        ans = find_gap(i)
        gap_data.append(ans)
        #print(ans)



print(gap_data)
plot_data(gap_data)


