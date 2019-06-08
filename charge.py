#!/usr/bin/env python
import os
import sys
import pwd
import glob
import pdb
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


def find_H(file_name):
	param = file_name.split("_")
        H = float(param[1][1:4])
	return H
	

def write_mag(file_name):
	param = file_name.split("_")
	H = param[1][1:4]	
	f0 = open(file_name,"r")
	word = "charge ="
	if float(H) == 0:
		ans = [8]
	else:
		ans = [1.0/float(H)]
	
	while True:
		line = f0.readline()
		if word in line:
			ans.append(float(line.split(word)[-1]) )
		if len(line) < 1:
			break
	return ans 

def sort_file(file_list):
	u = {}
	tmp_list = []
	for i in file_list:
		if "moment" not in i:
			continue
		param = i.split("_")
		x = float( param[1][1:4] )
		u[x] = i
		
	print(sorted(u.keys()))
	for i in sorted(u.keys()):
		tmp_list.append(u[i])
	return tmp_list


def plot_data(magnetization_data):
	magnetization_data = np.array(magnetization_data)
	color_code = ['r--', 'b-.']
	legend = ["O[1]","O[2]"]
	plt.xlabel(r"$1/B$", fontsize = 20)
	plt.ylabel(r"$\langle Q \rangle$", fontsize = 20)
	plt.xlim(0,1)
	for x in range(1,3):
		plt.plot( magnetization_data[:,0],magnetization_data[:,x],color_code[x-1] , linewidth=4.0, label = legend[x-1])
	plt.legend(loc = "upper right")
	plt.show()
	#try:
	#	input("Hit ENTER to close")
	#except:
	#	plt.close()
	#pass

'''
MAIN PROGRAM
'''

file_list = glob.glob('*001_U1.00.txt')
file_list = sort_file(file_list)
print(file_list)
magnetization_data = []

for i in file_list:
	H = find_H(i)
	if H == 0:
		continue
	ans = write_mag(i)
	magnetization_data.append(ans)
	#print(ans)



plot_data(magnetization_data)
