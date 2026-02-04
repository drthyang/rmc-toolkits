#!/usr/bin/env python3
import numpy as np
from scipy.stats import gaussian_kde
import sys, glob
import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.transforms import Affine2D
import seaborn as sns

plt.rcParams['font.family'] = 'Dejavu Sans'
plt.rcParams['mathtext.fontset'] = 'dejavusans'
plt.rcParams['lines.linewidth']= 1
plt.rcParams['axes.facecolor'] = 'w'

color_map = {
        0 : (0.8, 0.2, 0.2),  # red
        1 :  (0.2, 0.2, 0.8),  # blue
        2 : (0.2, 0.8, 0.2),  # green
        3 : (0.8, 0.8, 0.2),  # yellow
        4 : (0.5, 0.5, 0.5),  # gray
        5 :  (0.3, 0.6, 0.3),  # dark green
        6 : (0.7, 0.4, 0.7),  # purple
        7 : (0.2, 0.6, 0.8),  # teal
		}

#fpath = '/Users/tt9/Research/LacunarSpinels/rmc/server_data/GTS/local_5A/5K_try1/'

fpath = './'
fname = glob.glob(fpath + 'Frac*.txt')
print(fname)



# Read the fractional atomic coordinates from Frac*.txt
def read_frac_atom(fpath,atype='Ta',mode='Frac') :
	# get atom idx
	atom_dic = get_atom_idx(fpath)
	print(atom_dic)
	# get lattice vecs and supercell dims
	v1,v2,v3,dim = read_cell_vec(fpath)
	UB = np.array([v1,v2,v3])
	v1_norm = v1/dim[0]
	v2_norm = v2/dim[1]
	v3_norm = v3/dim[2]
	v1_norm = v1_norm/np.sqrt(np.dot(v1_norm,v1_norm))
	v2_norm = v2_norm/np.sqrt(np.dot(v2_norm,v2_norm))
	v3_norm = v3_norm/np.sqrt(np.dot(v3_norm,v3_norm))
	fname = glob.glob(fpath + 'Frac*.txt')[0]
	f = open(fname,'r')
	lines = f.readlines()
	atmtype = []
	data = []
	for ii in np.arange(5,len(lines),1) :
		ln = lines[ii].split()
		#print(ln[0])
		if atype==0 :
			atmtype.append(ln[0])
			xyz = np.array(np.float64(ln[1:4]))*dim # Fractional coord. for unit cell
			#print(xyz)
			xyz = [x-dim[0] if x > 1 else x for x in xyz]
			if mode=='Frac' :
				data.append(xyz)
			else :
				xyz = cvt_pos(xyz,v1_norm,v2_norm,v3_norm)
				data.append(xyz)
		elif (int(ln[0]) in atom_dic[atype]) :
			atmtype.append(ln[0])
			xyz = np.array(np.float64(ln[1:4]))*dim
			#print(xyz)
			#xyz = [x-dim[0] if x > 1.5 else x for x in xyz]
			xyz = [ -x%1 if x > 1.9 else x for x in xyz]
			if mode=='Frac' :
				data.append(xyz)
			else :
				xyz = cvt_pos(xyz,v1_norm,v2_norm,v3_norm)
				data.append(xyz)
	return atmtype, data

def read_cell_vec(fpath) :
	fname = glob.glob(fpath + '*.rmc6f')[0]
	print(fname)
	f = open(fname,'r')
	lines = f.readlines()
	for ii in np.arange(len(lines)) :
		if (lines[ii].split()[0]=='Lattice') :
			v1 = np.array( np.float64(lines[ii+1].split()) )
			v2 = np.array( np.float64(lines[ii+2].split()) )
			v3 = np.array( np.float64(lines[ii+3].split()) )
			print('Lattice vectors:')
			print(v1)
			print(v2)
			print(v3)
			break
		elif (lines[ii].split()[0]=='Supercell') :
			dim = np.float64(lines[ii].split()[-3:])
			print('Supercell dimensions = {}'.format(dim))
	return v1,v2,v3,np.array(dim)

def get_atom_idx(fpath) :
	fname = glob.glob(fpath + '*.rmc6f')[0]
	print(fname)
	f = open(fname,'r')
	lines = f.readlines()
	atom_dic = {}
	for ii in np.arange(len(lines)) :
		if (lines[ii].split()[0]=='Atoms:') :
			idx_ini = ii
			break
	for ii in np.arange(idx_ini+1,len(lines),1) :
		ln = lines[ii].split()
		#print(ln)
		atom = ln[1]
		atom_idx = int(ln[-4])
		if atom not in atom_dic :
			atom_dic[atom] = [atom_idx]
		else :
			atom_dic[atom].append(atom_idx)
	for key in atom_dic :
		atom_dic[key] = list(set(atom_dic[key]))
	return atom_dic
			
def gen_3d_plot(fpath,atype=0) :
	print('Making folded atomic positions for {} ...'.format(atype))
	# 3D plot
	atmtype, data = read_frac_atom(fpath,atype)
	xyz = np.transpose(data)
	kde = gaussian_kde(xyz)
	density = kde(xyz)

	from mayavi import mlab
	# Plot scatter with mayavi
	figure = mlab.figure('DensityPlot')
	pts = mlab.points3d(xyz[0], xyz[1], xyz[2], density, scale_mode='none', scale_factor=0.01, color=(1,1,0))
	mlab.axes()
	mlab.show()

def gen_3d_plot_all(fpath,atm_dic) :
	print('Making folded atomic positions for all elements ...')
	# 3D plot
	from mayavi import mlab
	# Plot scatter with mayavi
	figure = mlab.figure('DensityPlot')
	icolor = 0
	for atype in atm_dic :
		atmtype, data = read_frac_atom(fpath,atype)
		xyz = np.transpose(data)
		kde = gaussian_kde(xyz)
		density = kde(xyz)
		pts = mlab.points3d(xyz[0], xyz[1], xyz[2], density, scale_mode='none', scale_factor=0.01, color=color_map[icolor])
		icolor += 1
	mlab.axes(extent=[0, 1, 0, 1, 0, 1])
	print('Done. Please check the pop-up windw ...')
	mlab.show()


atm_dic = get_atom_idx('./')

gen_3d_plot_all('./',atm_dic)

#for atype in atm_dic :
#	gen_3d_plot('./',atype)


