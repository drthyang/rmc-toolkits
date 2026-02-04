#!/usr/bin/env python3
import numpy as np
import sys, glob
import matplotlib.pyplot as plt
from matplotlib import rc
plt.rcParams['font.family'] = 'Dejavu Sans'
plt.rcParams['mathtext.fontset'] = 'dejavusans'
plt.rcParams['lines.linewidth']= 1
plt.rcParams['axes.facecolor'] = 'w'

#fname = glob.glob('./*_FQ1partials.csv')

def read_input() :
	f = open('stog_input.dat','r')
	lines = f.readlines()
	inpname = lines[1].split()[0]
	return inpname

def read_stog(fname) :
	f = open(fname,'r')
	lines = f.readlines()
	#labels = lines[0].split(',')
	#labels[-1] = labels[-1].split('\n')[0]
	data_array = []
	for ii in np.arange(2,len(lines),1) :
		data = lines[ii].split()
		tmp = []
		for jj in np.arange(len(data)) :
			tmp.append(np.float64(data[jj]))
		data_array.append(tmp)
	data_array = np.array(data_array)
	data_array = np.transpose(data_array)
	return data_array

def read_chi(fnames) :
	chi_Q = []
	chi_R = []
	for ii in np.arange(len(fnames)) :
		f = open(fnames[ii],'r')
		lines = f.readlines()
		for jj in np.arange(2,len(lines),1) :
			tmp = lines[jj].split()
			chi_Q.append(np.float64(tmp[-2]))
			chi_R.append(np.float64(tmp[-1]))
	chi_R = np.array(chi_R)
	chi_Q = np.array(chi_Q)
	return chi_Q, chi_R


################
fig = plt.figure(figsize=(3.375*2,3.375*1.4))
ax = fig.add_subplot(111)

fname = glob.glob('./scale_ft.gr')
data = read_stog(fname[0])

#for ii in np.arange(1,len(labels),1) :
ax.plot(data[0],data[1],label='scale_ft.gr',
									lw=1.0,alpha=1.0,color='r')
ax.hlines(1,data[0][0],data[0][-1],ls='--',lw=0.5,alpha=1.0,color='black')
ax.hlines(0,data[0][0],data[0][-1],ls='--',lw=0.5,alpha=1.0,color='black')

ax.set_xlabel(r'r ($\mathrm{\AA}$)',fontsize=11)
ax.set_ylabel(r'G(r)',fontsize=11)
ax.legend(loc=1,fontsize=9,frameon=False)
	
################
fig2 = plt.figure(figsize=(3.375*2,3.375*1.4))
bx = fig2.add_subplot(111)

fname = glob.glob('./scale_ft.sq')
data = read_stog(fname[0])

bx.plot(data[0],data[1],label='scale_ft.sq',
									lw=1.0,alpha=1.0,color='r')
bx.hlines(1,data[0][0],data[0][-1],ls='--',lw=0.5,alpha=1.0,color='black')
#bx.hlines(0,data[0][0],data[0][-1],ls='--',lw=0.5,alpha=1.0,color='black')

bx.set_xlabel(r'Q ($\mathrm{\AA^{-1}}$)',fontsize=11)
bx.set_ylabel(r'S(Q)',fontsize=11)
bx.legend(loc=1,fontsize=9,frameon=False)

################
fig3 = plt.figure(figsize=(3.375*2,3.375*1.4))
cx = fig3.add_subplot(111)

inpname = read_input()
fname = glob.glob('./{}'.format(inpname))
data = read_stog(fname[0])

cx.plot(data[0],data[1],label='{}'.format(fname[0].split('/')[-1]),
									lw=1.0,alpha=1.0,color='r')
cx.hlines(1,data[0][0],data[0][-1],ls='--',lw=0.5,alpha=1.0,color='black')

cx.set_xlabel(r'Q ($\mathrm{\AA^{-1}}$)',fontsize=11)
cx.set_ylabel(r'data',fontsize=11)
cx.legend(loc=1,fontsize=9,frameon=False)

################
fig4 = plt.figure(figsize=(3.375*2,3.375*1.4))
dx = fig4.add_subplot(111)

fname = glob.glob('./scale_ft_rmc.fq')
data = read_stog(fname[0])

#dx.plot(chi_Q,label=r'Q',lw=1.0,alpha=0.5)
#dx.plot(chi_R,label=r'R',lw=1.0,alpha=0.5)
dx.plot(data[0],data[1],label='{}'.format(fname[0].split('/')[-1]),
									lw=1.0,alpha=1.0,color='r')
dx.hlines(0,data[0][0],data[0][-1],ls='--',lw=0.5,alpha=1.0,color='black')

dx.set_xlabel(r'Q ($\mathrm{\AA^{-1}}$)',fontsize=11)
dx.set_ylabel(r'data',fontsize=11)
dx.legend(loc=1,fontsize=9,frameon=False)

################
fig5 = plt.figure(figsize=(3.375*2,3.375*1.4))
ex = fig5.add_subplot(111)

inpname = read_input()
fname = glob.glob('./{}'.format(inpname))
data = read_stog(fname[0])

ex.plot(data[0],data[0] * data[1],label='{}'.format(fname[0].split('/')[-1]),
									lw=1.0,alpha=1.0,color='r')
ex.hlines(1,data[0][0],data[0][-1],ls='--',lw=0.5,alpha=1.0,color='black')

ex.set_xlabel(r'Q ($\mathrm{\AA^{-1}}$)',fontsize=11)
ex.set_ylabel(r'Q*S(Q)',fontsize=11)
ex.legend(loc=1,fontsize=9,frameon=False)

plt.show()
