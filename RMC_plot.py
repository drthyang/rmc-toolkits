#!/usr/bin/env python3
import numpy as np
import sys, glob
import matplotlib.pyplot as plt
from matplotlib import rc
plt.rcParams['font.family'] = 'Dejavu Sans'
plt.rcParams['mathtext.fontset'] = 'dejavusans'
plt.rcParams['lines.linewidth']= 1
plt.rcParams['axes.facecolor'] = 'w'

def read_csv(fname) :
	f = open(fname,'r')
	lines = f.readlines()
	labels = lines[0].split(',')
	labels[-1] = labels[-1].split('\n')[0]
	data_array = []
	for ii in np.arange(1,len(lines),1) :
		data = lines[ii].split(',')
		tmp = []
		for jj in np.arange(len(data)) :
			if data[jj]!=' \n' :
				tmp.append(np.float64(data[jj]))
		data_array.append(tmp)
	data_array = np.array(data_array)
	data_array = np.transpose(data_array)
	return labels, data_array

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

def Rwp(r, gr, grfit, fit_range) :
	idx = [ ii for ii in np.arange(len(r)) if (r[ii]>=fit_range[0]) and (r[ii]<=fit_range[-1]) ]
	gr_sub = np.array(gr[idx])
	grfit_sub = np.array(grfit[idx])
	Rsq = np.sum( (grfit_sub-gr_sub) * (grfit_sub-gr_sub)  ) / np.sum( gr_sub*gr_sub )
	return np.sqrt(Rsq)


################ PDFpartials.csv # Partial PDF
fig = plt.figure(figsize=(3.375*2,3.375*1.2))
ax = fig.add_subplot(111)

fname = glob.glob('./*_PDFpartials.csv')
labels, data = read_csv(fname[0])

for ii in np.arange(1,len(labels),1) :
	ax.plot(data[0],data[ii],label=labels[ii].strip(),lw=1.0,alpha=0.5)

ax.set_xlabel(r'{}'.format(labels[0].strip()),fontsize=11)
ax.set_ylabel(r'data',fontsize=11)
ax.legend(loc=1,fontsize=9,frameon=False)
fig.suptitle('Partial PDF', fontsize=14)
				
################ *FT_XFQ1.csv or *PDF1.csv # Real space G(r)
fig2 = plt.figure(figsize=(3.375*2,3.375*1.2))
bx = fig2.add_subplot(111)

fname_x = glob.glob('./*_FT_XFQ1.csv')
#fname_n = glob.glob('./*_SQ1.csv')
fname_n = glob.glob('./*PDF1.csv')
fname = fname_x + fname_n


for jj in np.arange(len(fname)) :
	labels, data = read_csv(fname[jj])
	colors = ['r','b']

	#fit_range = [2.2,15.0]
	fit_range = [data[0][0],data[0][-1]]
	Rw = Rwp(data[0],data[1],data[2],fit_range)
	print(f"{'Real space:':<20} R = {Rw:.6f}")

	for ii in np.arange(1,len(labels),1) :
		bx.plot(data[0],data[ii],label=labels[ii].strip(),lw=1.0,alpha=0.5)

bx.set_xlabel(r'{}'.format(labels[0].strip()),fontsize=11)
bx.set_ylabel(r'data',fontsize=11)
bx.legend(loc=1,fontsize=9,frameon=False)
fig2.suptitle('PDF', fontsize=14)

################ *_FQ1.csv or *_SQ1.csv # Reciprocal space S(Q)
fig3 = plt.figure(figsize=(3.375*2,3.375*1.2))
cx = fig3.add_subplot(111)

fname_x = glob.glob('./*_FQ1.csv')
fname_n = glob.glob('./*_SQ1.csv')
fname = fname_x + fname_n
labels, data = read_csv(fname[0])

#fit_range = [0.5,26.0]
fit_range = [data[0][0],data[0][-1]]
Rw = Rwp(data[0],data[1],data[2],fit_range)
#print('Reciprocal space: R = {:.6f}'.format(Rw))
print(f"{'Reciprocal space:':<20} R = {Rw:.6f}")

for ii in np.arange(1,len(labels),1) :
	cx.plot(data[0],data[ii],label=labels[ii].strip(),lw=1.0,alpha=0.5)

cx.set_xlabel(r'{}'.format(labels[0].strip()),fontsize=11)
cx.set_ylabel(r'data',fontsize=11)
cx.legend(loc=1,fontsize=9,frameon=False)
fig3.suptitle('S(Q)', fontsize=14)

################ BRAGG
fnames = glob.glob('./*_bragg.csv')
if len(fnames)>0 :
	fig5 = plt.figure(figsize=(3.375*2,3.375*1.2))
	ex = fig5.add_subplot(111)

	labels, data = read_csv(fnames[0])
	fit_range = [data[0][0],data[0][-1]]
	Rw = Rwp(data[0],data[1],data[2],fit_range)
	#print('BRAGG: R = {:.6f}'.format(Rw))
	print(f"{'BRAGG:':<20} R = {Rw:.6f}")
	
	for ii in np.arange(1,len(labels),1) :
		ex.plot(data[0],data[ii],label=labels[ii].strip(),lw=1.0,alpha=0.5)

	ex.set_xlabel(r'Q ($\mathrm{\AA^{-1}}$)',fontsize=11)
	ex.set_ylabel(r'data',fontsize=11)
	ex.legend(loc=1,fontsize=9,frameon=False)
	fig5.suptitle('BRAGG', fontsize=14)

################ chi-values
fig4 = plt.figure(figsize=(3.375*2,3.375*1.2))
dx = fig4.add_subplot(111)

fnames = glob.glob('./*-*.log')
chi_Q, chi_R = read_chi(fnames)

#dx.plot(np.log(chi_Q),label=r'Q',lw=1.0,alpha=0.5)
dx.plot(np.log(chi_R),label=r'R',lw=1.0,alpha=0.5)

dx.set_xlabel(r'Time steps',fontsize=11)
dx.set_ylabel(r'log($\mathrm{\chi}$)',fontsize=11)
dx.legend(loc=1,fontsize=9,frameon=False)
fig4.suptitle('R-value', fontsize=14)


plt.show()
