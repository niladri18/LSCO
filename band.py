#!/usr/bin/env python
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy as scp
import pylab
from matplotlib import rc
from matplotlib import cm

def plot_band(x,ans):
	ans = np.array(ans, dtype=float)
	x = np.array(x);
	fermi = np.zeros(len(x))
	xlabel = [r'$\Gamma$',r"($\pi,0$)","($\pi,\pi$)",r'$\Gamma$']
	#print(x)
	#print(ans)
	n = len(ans[0])
	#print(n)
	for i in range( 8, n ):
        	plt.plot(x,ans.transpose()[i],'r')

	plt.plot(x,fermi,'b--' )
	plt.xticks([0,100,200,300],xlabel, fontsize=20)
	plt.ylabel('E-E'+r'$_{F}$', fontsize=20)
	plt.title(r"$t_{\sigma}$=1.0"+ "\t"+r"$t_{\pi}$=0.2"+"\t" + "U="+U)
	#plt.set_xticks((0,100,200,300))
	#fig.show()
	plt.show()

	pass

def plot_gap(x,ans):
        ans = np.array(ans, dtype=float)
        x = np.array(x);
        fermi = np.zeros(len(x))
        xlabel = [r'$\Gamma$',"X","W",r'$\Gamma$']
        #print(x)
        #print(ans)
        n = len(ans[0])
        #print(n)
	y = ans.transpose()[n-1] - ans.transpose()[n-2]
	plt.plot(x,y,'r')
        #for i in range( 10, n ):
        #        plt.plot(x,ans.transpose()[i],'r')

        plt.plot(x,fermi,'b--' )
        plt.xticks([0,100,200,300],xlabel)
        plt.ylabel(r'$E_{1}-E_{0}$')
        plt.title(r"$t_{\sigma}$=1.0"+ "\t"+r"$t_{\pi}$=0.2"+"\t" + "U="+U)
        #plt.set_xticks((0,100,200,300))
        #fig.show()
        plt.show()

        pass

fname = sys.argv[1];
print(fname);

font = {'family':'serif'}
rc('font',**font)

f0 = open(fname,"r")
param = fname.split("_")
U = param[-1][1:-4]
print(U)
ans = []
x = []
while True:
	line = f0.readline().split();
	if len(line) < 1:
		break

	t = [float(i) for i in line]
	x.append(t[0])
	#t = np.asarray(x)
	#print(t)	
	if ans == []:
		ans = [t[1:]]
	else:
		ans.append(t[1:])
		

plot_band(x,ans)
#plot_gap(x,ans)

