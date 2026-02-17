#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RMC_KDE.py
==========

3D visualization tool for Reverse Monte Carlo (RMC) modeling results.

This script provides functions to visualize atomic configurations 
in three dimensions, with options for:
    - Generating atomic density slices (via KDE)
	- Visualizing average atomic displacements (under development)

Dependencies
------------
- Python 3.x
- NumPy
- SciPy
- Matplotlib

Author Information
------------------
Author  : Tsung-Han Yang, Ph.D.  
Affiliation : Oak Ridge National Laboratory (Spallation Neutron Source)  
Contact : tsung-han_yang@alumni.brown.edu  

Version
-------
v1.0.0   (2025-08-27)

License
-------
This code is distributed for academic and research use.  
Please cite the corresponding publication(s) if used in published work.  
See LICENSE file for details.

"""

from typing import Optional, Tuple
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons
from matplotlib.gridspec import GridSpec
from scipy.stats import gaussian_kde
import glob, sys
import argparse

# Read the fractional atomic coordinates from Frac*.txt
def read_frac_atom(fpath,atype='Ta',mode='Latt') :
	# get atom idx
	atom_dic = get_atom_idx(fpath)
	print('>>> Atomic labels: {} '.format(atom_dic))
	# get lattice vecs and supercell dims
	v1,v2,v3,dim = read_cell_vec(fpath)
	UB = np.array([v1,v2,v3])
	v1_norm = v1/dim[0]
	v2_norm = v2/dim[1]
	v3_norm = v3/dim[2]
	fname = glob.glob(fpath + 'Frac*.txt')[0]
	f = open(fname,'r')
	lines = f.readlines()
	atmtype = []
	data = []
	for ii in np.arange(5, len(lines), 1):
		ln = lines[ii].split()
		#print(ln[0])
		if atype == 0:
			atmtype.append(ln[0])
			xyz = np.array(np.float64(ln[1:4])) * dim  # Fractional coord. for unit cell
			#print(xyz)
			xyz = [x - dim[0] if x > 1 else x for x in xyz]
			if mode == 'Frac':
				data.append(xyz)
			else:
				#xyz = cvt_pos(xyz, v1_norm, v2_norm, v3_norm) 
				xyz = xyz[0] * v1_norm + xyz[1] * v2_norm + xyz[2] * v3_norm
				data.append(xyz)
		elif (int(ln[0]) in atom_dic[atype]):
			atmtype.append(ln[0])
			xyz = np.array(np.float64(ln[1:4])) * dim
			#print(xyz)
			#xyz = [x-dim[0] if x > 1.5 else x for x in xyz] 
			# Deal with periodic boundary condition
			xyz = [(x % 1) if x > 1.0 else x for x in xyz]
			if mode == 'Frac':
				data.append(xyz)
			else:
				#xyz = cvt_pos(xyz, v1_norm, v2_norm, v3_norm) 
				xyz = xyz[0] * v1_norm + xyz[1] * v2_norm + xyz[2] * v3_norm
				data.append(xyz)
	#print([v1/dim[0], v2/dim[1], v3/dim[2]])
	return atmtype, data, [v1/dim[0], v2/dim[1], v3/dim[2]]

def read_cell_vec(fpath) :
	fname = glob.glob(fpath + '*.rmc6f')[0] 
	print('>>> Reading lattice vectors from: {}'.format(fname))
	#print(fname)
	f = open(fname,'r')
	lines = f.readlines()
	for ii in np.arange(len(lines)) :
		if (lines[ii].split()[0]=='Lattice') :
			v1 = np.array( np.float64(lines[ii+1].split()) )
			v2 = np.array( np.float64(lines[ii+2].split()) )
			v3 = np.array( np.float64(lines[ii+3].split()) )
			print('>>> Lattice vectors:')
			print('>>> v1 = {} '.format(v1))
			print('>>> v2 = {} '.format(v2))
			print('>>> v3 = {} '.format(v3))
			break
		elif (lines[ii].split()[0]=='Supercell') :
			dim = np.float64(lines[ii].split()[-3:])
			print('>>> Supercell dimensions = {}'.format(dim))
	return v1,v2,v3,np.array(dim)

def get_atom_idx(fpath,verbose=0) :
	fname = glob.glob(fpath + '*.rmc6f')[0] 
	if verbose>0 :
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

def interactive_kde_slice(
    xyz,  # accept list/tuple/ndarray/DataFrame
    slab_thickness: Optional[float] = None,
    gridsize_xy: int = 300,
    gridsize_xz: Tuple[int, int] = (150, 150),
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    zlim: Optional[Tuple[float, float]] = None,
    cmap_xy: str = "seismic",
    cmap_xz: str = "Reds",
    quantile_clip: Tuple[float, float] = (0.01, 0.99),
    z_hist_bins: int = 150,
    normalize_hist: bool = True,
    # contour style
    contour_levels: int = 8,
    contour_color: str = "gray",
    contour_linewidth: float = 1,
) -> None:
    # --- Robust input coercion ---
    if hasattr(xyz, "to_numpy"):  # pandas DataFrame
        xyz = xyz.to_numpy()
    xyz = np.asarray(xyz, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3:
        raise ValueError(f"`xyz` must have shape (N, 3), got {xyz.shape}.")

    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]

    # Axis limits (robust quantiles for x,y)
    if xlim is None:
        xlim = (np.quantile(x, quantile_clip[0]), np.quantile(x, quantile_clip[1]))
    if ylim is None:
        ylim = (np.quantile(y, quantile_clip[0]), np.quantile(y, quantile_clip[1]))
    if zlim is None:
        zlim = (float(np.min(z)), float(np.max(z)))
    zmin, zmax = zlim

    # Default slab thickness
    if slab_thickness is None:
        slab_thickness = 0.05 * (zmax - zmin) if zmax > zmin else 1.0

    # XY grid for KDE
    gx = np.linspace(*xlim, gridsize_xy)
    gy = np.linspace(*ylim, gridsize_xy)
    XX, YY = np.meshgrid(gx, gy)
    grid_xy = np.vstack([XX.ravel(), YY.ravel()])  # (2, M)

    # Figure + axes (3 panels horizontally)
    fig = plt.figure(figsize=(16, 8))
    gs = GridSpec(1, 3, width_ratios=[1,1,0.2], figure=fig, wspace=0.25)
    ax_main = fig.add_subplot(gs[0, 0])  # left: main KDE
    ax_xz   = fig.add_subplot(gs[0, 1])  
    ax_hist = fig.add_subplot(gs[0, 2])  
    plt.subplots_adjust(bottom=0.15)

    z0 = 0.5 * (zmin + zmax)

    def compute_xy_kde(z_center: float, dz: float, bw: float=0.03) -> np.ndarray:
        mask = (z >= z_center - 0.5 * dz) & (z <= z_center + 0.5 * dz)
        pts = np.column_stack([x[mask], y[mask]])
        if pts.shape[0] < 5:
            return np.zeros_like(XX)
        kde2d = gaussian_kde(pts.T, bw_method=bw)
        return kde2d(grid_xy).reshape(XX.shape)

    # --- Main XY heatmap + contours ---
    dens0 = compute_xy_kde(z0, slab_thickness)
    im_xy = ax_main.imshow(
        dens0, origin="lower",
        extent=[xlim[0], xlim[1], ylim[0], ylim[1]],
        aspect="equal", cmap=cmap_xy
    )
    cbar_xy = fig.colorbar(im_xy, ax=ax_main, pad=0.01, aspect=30)
    cbar_xy.set_label("KDE density (slab)", rotation=270, labelpad=14, fontsize=14)
    ax_main.set_xlabel(r"x ($\mathrm{\AA}$)",fontsize=12); ax_main.set_ylabel(r"y ($\mathrm{\AA}$)",fontsize=12)
    title_main = ax_main.set_title(f"XY KDE @ z={z0:.3g}, dz={slab_thickness:.3g}",fontsize=14)

    # contour overlay on XY panel
    contours_xy = ax_main.contour(
        XX, YY, dens0,
        levels=contour_levels,
        colors=contour_color,
        linewidths=contour_linewidth,
    )

    # --- z histogram + slab band ---
    ax_hist.hist(z, bins=z_hist_bins, density=normalize_hist,
                 orientation="horizontal", color="0.6")
    ax_hist.set_ylim(zmin, zmax)
    ax_hist.set_xlabel("count" if not normalize_hist else "density",fontsize=12)
    ax_hist.set_title("z distribution",fontsize=14)
    bandA = ax_hist.axhspan(z0 - 0.5 * slab_thickness, z0 + 0.5 * slab_thickness,
                            color="C3", alpha=0.25)
    ax_hist.set_aspect('equal')
    # --- Global x–z density + band ---
    nx, nz = gridsize_xz
    H, xedges, zedges = np.histogram2d(x, z, bins=[nx, nz], range=[xlim, zlim])
    im_xz = ax_xz.imshow(
        H.T, origin="lower",
        extent=[xedges[0], xedges[-1], zedges[0], zedges[-1]],
        aspect="auto", cmap=cmap_xz
    )
    cbar_xz = fig.colorbar(im_xz, ax=ax_xz, pad=0.01, aspect=30)
    cbar_xz.set_label("counts", rotation=270, labelpad=16,fontsize=12)
    ax_xz.set_xlabel(r"x ($\mathrm{\AA}$)",fontsize=12); ax_xz.set_ylabel(r"z ($\mathrm{\AA}$)",fontsize=12)
    ax_xz.set_title("XZ projection",fontsize=14)
    bandB = ax_xz.axhspan(z0 - 0.5 * slab_thickness, z0 + 0.5 * slab_thickness,
                          color="C3", alpha=0.25)
    ax_xz.set_aspect('equal')

    # --- Sliders ---
    ax_z  = plt.axes([0.25, 0.07, 0.50, 0.06])
    ax_dz = plt.axes([0.25, 0.01, 0.50, 0.06])
    s_z  = Slider(ax_z,  r"z0 ($\mathrm{\AA}$)", zmin, zmax, valinit=z0)
    s_dz = Slider(ax_dz, r"dz ($\mathrm{\AA}$)", max((zmax - zmin) / 1000.0, 1e-9),
                  (zmax - zmin), valinit=slab_thickness)
    # --- Change font size ---
    s_z.label.set_fontsize(14)
    s_z.valtext.set_fontsize(14)
    s_dz.label.set_fontsize(14)
    s_dz.valtext.set_fontsize(14)

    def update(_):
        nonlocal bandA, bandB, contours_xy
        zc = s_z.val
        dz = max(s_dz.val, 1e-12)

        dens = compute_xy_kde(zc, dz)
        im_xy.set_data(dens)
        im_xy.set_clim(vmin=np.min(dens), vmax=np.max(dens) if np.max(dens) > 0 else 1.0)
        title_main.set_text(f"XY KDE @ z={zc:.3g}, dz={dz:.3g}")

        # refresh bands
        try:
            bandA.remove(); bandB.remove()
        except Exception:
            pass
        bandA = ax_hist.axhspan(zc - 0.5 * dz, zc + 0.5 * dz, color="C3", alpha=0.25)
        bandB = ax_xz.axhspan(zc - 0.5 * dz, zc + 0.5 * dz, color="C3", alpha=0.25)

        # refresh contours
        for coll in contours_xy.collections:
            coll.remove()
        contours_xy = ax_main.contour(
            XX, YY, dens,
            levels=contour_levels,
            colors=contour_color,
            linewidths=contour_linewidth,
        )

        fig.canvas.draw_idle()

    s_z.on_changed(update)
    s_dz.on_changed(update)
    plt.show()


def main():
    """
    Entry point for RMC_KDE visualization.

    Usage examples:
      python KDE.py                 # plot all elements
      python KDE.py --el Mn         # plot only Mn
    """
    # Check if Frac*.txt in the current folder
    fpath = './'
    fname = glob.glob(fpath + 'Frac*.txt')
    if len(fname)==0 : 
        print('No Frac*.txt file found in the folder. Please run "thermal_ellipsoid" first ...') 
        sys.exit() 
    else : 
        print('>>> Fetching folded atomic positions from {}'.format(fname))
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--el", "--element", dest="element", default=None,
                        help="Element symbol or index understood by read_frac_atom (e.g. 'Mn'). "
                             "If omitted, use all elements.")
    args = parser.parse_args()

    data_dir = "./" 

    # If no element provided -> use '0' (read_frac_atom convention for 'all')
    element_token = args.element if args.element is not None else 0
    print(">>> Plotting KDE map for >> {} << element(s) ".format("all" if element_token == 0 else element_token))

    # Read data; this already supports strings like 'Mn'
    atmtype, data, lattvecs = read_frac_atom(data_dir, element_token)

    # xyz should be (N, 3).
    xyz = np.asarray(data)
    if xyz.shape[0] == 3 and xyz.shape[1] != 3:
        xyz = xyz.T

    # limits of plots, only supports orthorhombic cell now ...
    xlim = (0, float(lattvecs[0].max()))
    ylim = (0, float(lattvecs[1].max()))
    zlim = (0, float(lattvecs[2].max()))

    # Call your existing interactive plot
    interactive_kde_slice(xyz, xlim=xlim, ylim=ylim, zlim=zlim)

if __name__ == "__main__":
    main()
