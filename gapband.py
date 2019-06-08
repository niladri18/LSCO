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
	#p = cs.collections[0].get_paths()[0]
	#v = p.vertices
	plt.show()
	return v

def plot_gap(x,y,z,U,s):
	#cfname: file name of the contour file
	#var = "{:.2f}".format(U)
	#var = "%.2f"%U
	dic0 = {}
        dic0["Ediffup"] = "1"
        dic0["Ediffdn"] = "0"
	var = U
	cfname = "contourEup"+dic0[s]+"U"+var+".txt"
	if os.path.isfile("./"+cfname) :
		f0 = open(cfname,"r")
		print("Found file: ",cfname)
	else:
		print("Contour File not found..")
		print("Using default..")
		f0 = open("contourEup1U0.00.txt","r")
	#cfname = "contourEup1.txt"
	#f0 = open(cfname,"r")
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
	dic["Ediffup"] = r"$\Delta E \uparrow$"
	dic["Ediffdn"] = r"$\Delta E \downarrow$"
	#s = "Ediffup"
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
	plt.show()
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

fname = sys.argv[1]
cfname = "contourE1.txt"
#print(fname)
U = find_u(fname)
print(U)

f0 = open(fname,'r')

x = []
y = []
Ediffup = []
Ediffdn = []

while True:
	line = f0.readline().split()
	if len(line) < 1:
		break
	t = [float(i) for i in line]
	x.append(t[0])
	y.append(t[1])
	Ediffup.append(t[-1])
	Ediffdn.append(t[-2])
f0.close()


x = np.array(x)
y = np.array(y)
Ediffup = convert(Ediffup)
Ediffdn = convert(Ediffdn)

z = Ediffup
s = "Ediffup"

#z = Ediffdn
#s = "Ediffdn"


plot_gap(x,y,z, U, s)

