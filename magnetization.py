#!/usr/bin/env python
import os
import sys
import pwd
import glob
import pdb
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc



def write_mag(file_name):
	param = file_name.split("_")
	U = param[-1][1:-4]	
	f0 = open(file_name,"r")
	#word = "total moment = "

	word = "charge ="
	ans = [U]
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
		param = i.split("_")
		x = float( param[-1][1:-4] )
		u[x] = i
		
	print(sorted(u.keys()))
	for i in sorted(u.keys()):
		tmp_list.append(u[i])
	return tmp_list


def plot_data(magnetization_data):
	magnetization_data = np.array(magnetization_data)
	color_code = ['red', 'green']
	plt.xlabel("U")
	plt.ylabel("M")
	plt.xlim(0,3)
	for x in range(1,3):
		plt.plot( magnetization_data[:,0],magnetization_data[:,x],color_code[x-1] )
	plt.show()
	#try:
	#	input("Hit ENTER to close")
	#except:
	#	plt.close()
	#pass

'''
MAIN PROGRAM
'''

file_list = glob.glob('moment_H*')
file_list = sort_file(file_list)
print(file_list)
magnetization_data = []

for i in file_list:
	ans = write_mag(i)
	magnetization_data.append(ans)
	#print(ans)



plot_data(magnetization_data)
