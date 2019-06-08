#!/usr/bin/env python
import os
import sys
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
	dic["Eup1"] = r"$E(1)\uparrow $"
        dic["Edn1"] = r"$E(1)\downarrow$"
	dic["Eup0"] = r"$E(0)\uparrow $"
        dic["Edn0"] = r"$E(0)\downarrow $"
	dic["Ediff1"] = r"$\Delta E \uparrow$"
	dic["Ediff0"] = r"$\Delta E\downarrow$"

	plt.title(r"$t_{\sigma}$=1.0"+ "\t"+r"$t_{\pi}$=0.2"+"\t" + "U="+U +"\t"+ dic[s])
	#plt.imshow(z, extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)), norm = LogNorm() )
	plt.xlabel(r"k$_{x}$")
	plt.ylabel(r"k$_{y}$")
	extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y))
	#levels = -4.0
	plt.imshow(z, extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)) )
	plt.colorbar()
	#plt.contour(z, levels = np.array([-0.01,0.01]), colors = 'white', extent = extent)
	cs = plt.contour(z, levels = np.array([0.0]), colors = 'k', extent = extent)
	p = cs.collections[0].get_paths()[0]
	v = p.vertices
	plt.show()
	return v

def plot_gap(x,y,z,U,s):
	#cfname: file name of the contour file
	#cfname = "contourE0.txt"
	cfname = "contourE1.txt"
	f0 = open(cfname,"r")
	kx = []
	ky = []
	while True:
		line = f0.readline().split()
		if len(line) < 1:
			break
		line = [float(i) for i in line]
		kx.append(line[0])
		ky.append(line[1])

	f0.close()

        font = {'family':'serif'}
        rc('font',**font)
	dic = {}
	dic["Eup1"] = "Eup1"
        dic["Edn1"] = "Edn1"
	dic["Eup0"] = "Eup0"
        dic["Edn0"] = "Edn"
	dic["Ediff1"] = r"$\Delta E \uparrow$"
	dic["Ediff0"] = r"$\Delta E \downarrow$"
	s = "Ediff1"
        plt.title(r"$t_{\sigma}$=1.0"+ "\t"+r"$t_{\pi}$=0.2"+"\t" + "U="+U +"\t"+ dic[s])
        #plt.imshow(z, extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)), norm = LogNorm() )
        plt.xlabel(r"k$_{x}$")
        plt.ylabel(r"k$_{y}$")
        extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y))
        #levels = -4.0
        plt.imshow(z, extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)) )
        plt.colorbar()
	kx = np.array(kx)
	ky = np.array(ky)
	plt.scatter(kx,ky,s = 2, linewidth=2)
	#plt.scatter(kx,ky,extent=extent, colors = 'k')
	#plt.show()
	pass
        #plt.contour(z, levels = np.array([-0.01,0.01]), colors = 'white', extent = extent)
        #cs = plt.contour(z, levels = np.array([0.0]), colors = 'k', extent = extent)

def find_u(fname):
	param = fname.split("_")
	U = param[-1][1:-4]
	return U

def save_contour(v):
	kx = v[:,0]
	ky = v[:,1]
	N = len(kx)
	for i in range(N):
		s = str(v[i,0]) + '\t' + str(v[i,1])
		#print(v[i,0],v[i,1])
		print(s)



'''
Main program
'''

fname = sys.argv[1]
U = find_u(fname)
print(U)

f0 = open(fname,'r')

x = []
y = []
Eup1 = [] 
Edn1 = []
Eup0 = []
Edn0 = []
Ediff1 = []
Ediff0 = []

while True:
	line = f0.readline().split()
	if len(line) < 1:
		break
	t = [float(i) for i in line]
	x.append(t[0])
	y.append(t[1])
	Eup1.append(t[-1])
	Edn1.append(t[-2])
	Eup0.append(t[-3])
	Edn0.append(t[-4])
	#Ecdw0.append(t[-4])
	#Ecdw1.append(t[-5])
f0.close()


x = np.array(x)
y = np.array(y)
#N = int(len(Elow)**.5)
#Elow = np.array(Elow)
#Elow =  Elow.reshape(N,N)
Eup1 = convert(Eup1)
Edn1 = convert(Edn1)
Eup0 = convert(Eup0)
Edn0 = convert(Edn0)
#Ecdw0 = convert(Ecdw0)
#Ecdw1 = convert(Ecdw1)
#z = Ediff
#s = "Ediff"

z = Eup1
s = "Eup1"
#z = Edn1
#s = "Edn1"
#z = Eup0
#s = "Eup0"

#z = Elow
#s = "Elow"
#z = Ecdw1
#s = "Ecdw1"
#delta = True

v = plot_density(x,y,z, U,s)
save_contour(v)

#cfname = "contourE1.txt"
#save_contour(v)
#kx = v[:,0]
#ky = v[:,1]
#N = len(kx)
#for i in range(N):
#	print(v[i,0],v[i,1])
