#!/usr/bin/env python
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy as scp
import pylab
from matplotlib import rc
from matplotlib import cm


fname = sys.argv[1];
print(fname);

font = {'family':'serif'}
rc('font',**font)

f0 = open(fname,"r")
param = fname.split("_")
U = param[-1][1:-4]
print(U)
ans1 = []
ans2 = []
x = []
while True:
	line = f0.readline().split();
	if len(line) < 1:
		break

	t = [float(i) for i in line]
	x.append(t[0])
	ans2.append(t[-2])
	ans1.append(t[-1])
	#t = np.asarray(x)
	#print(t)	
	#if ans == []:
	#	ans = [t[1:]]
	#else:
	#	ans.append(t[1:])
		

#ans = np.array(ans, dtype=float)
ans1 = np.array(ans1)
ans2 = np.array(ans2)
x = np.array(x);
fermi = np.zeros(len(x))
xlabel = [r'$\Gamma$',"X","W",r'$\Gamma$']
#print(x)
#print(ans1)
plt.plot(x,ans1,'rx')
plt.plot(x,ans2,'bo')
plt.ylabel(r'$\Delta E$', fontsize=20)
plt.xlim(0,90)
plt.xlabel(r'$\theta$',fontsize = 20)
plt.show()
'''
for i in range( 10, n ):
	plt.plot(x,ans.transpose()[i],'r')

plt.plot(x,fermi,'b--' )
plt.xticks([0,100,200,300],xlabel)
plt.ylabel('E-E'+r'$_{F}$')
plt.title(r"$t_{\sigma}$=1.0"+ "\t"+r"$t_{\pi}$=0.2"+"\t" + "U="+U)
#plt.set_xticks((0,100,200,300))
#fig.show()
plt.show()
'''
