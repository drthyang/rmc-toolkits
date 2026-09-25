#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

import numpy as np
import sys, glob, re, os
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import rc

# The R-factor comes from the package (the single source of truth, shared with the
# web dashboard) so this script can never drift from it again.
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))
from rmc_toolkits.parsers import read_chi as _package_read_chi, rwp as _package_rwp, rwp_columns  # noqa: E402

plt.rcParams['font.family'] = 'Dejavu Sans'
plt.rcParams['mathtext.fontset'] = 'dejavusans'
plt.rcParams['lines.linewidth']= 1
plt.rcParams['axes.facecolor'] = 'w'

def parse_arguments():
    parser = argparse.ArgumentParser(description='Plot RMCProfile output files.')
    parser.add_argument('--dir', type=str, default='./', help='Input directory containing RMC output files (default: ./)')
    parser.add_argument('--save', action='store_true', help='Save figures as PNG files')
    parser.add_argument('--no-show', action='store_true', help='Do not display figures (useful for batch processing)')
    return parser.parse_args()

def read_csv(fname) :
    try:
        f = open(fname,'r')
        lines = f.readlines()
        f.close()
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
    except Exception as e:
        print(f"Error reading {fname}: {e}")
        return [], []

def read_chi(fnames) :
    """(second-to-last, last) log columns via rmc_toolkits.parsers.read_chi.

    The package reader checks each row against the header's column count, drops
    a half-written final line and keeps non-finite chi^2 rows as NaN.
    """
    return _package_read_chi(list(fnames))

def Rwp(r, observed, fitted, fit_range=None):
    """R-factor of ``fitted`` against ``observed`` (the EXPERIMENT) over ``fit_range``.

    Delegates to :func:`rmc_toolkits.parsers.rwp`: only rows finite in both columns
    count, and an undefined value (no such row, or an all-zero experiment) is
    ``None`` -- never ``0.0``, which would read as a perfect fit.
    """
    r = np.asarray(r, dtype=float)
    mask = np.ones(r.shape, dtype=bool)
    if fit_range is not None:
        mask = (r >= fit_range[0]) & (r <= fit_range[-1])
    return _package_rwp(
        r[mask],
        np.asarray(observed, dtype=float)[mask],
        np.asarray(fitted, dtype=float)[mask],
    )

def plot_data(fname, title, xlabel, ylabel, args, calc_rwp=False, rwp_label_prefix=""):
    if not fname:
        return
    
    labels, data = read_csv(fname)
    if len(labels) == 0:
        return

    if calc_rwp and len(data) >= 3:
        # RMCProfile writes (x, calculated, experimental); a header naming the
        # roles overrides that order. The experiment is the denominator.
        calculated, experimental = rwp_columns([label.strip() for label in labels], len(data))
        Rw = Rwp(data[0], data[experimental], data[calculated])
        shown = f"{Rw:.6f}" if Rw is not None else "n/a (undefined for this data)"
        print(f"{rwp_label_prefix:<20} R = {shown}")

    fig = plt.figure(figsize=(3.375*2,3.375*1.2))
    ax = fig.add_subplot(111)
    
    for ii in np.arange(1,len(labels),1) :
        ax.plot(data[0],data[ii],label=labels[ii].strip(),lw=1.0,alpha=0.5)

    ax.set_xlabel(xlabel,fontsize=11)
    ax.set_ylabel(ylabel,fontsize=11)
    ax.legend(loc=1,fontsize=9,frameon=False)
    fig.suptitle(title, fontsize=14)
    
    if args.save:
        outname = os.path.splitext(fname)[0] + '.png'
        fig.savefig(outname, dpi=300, bbox_inches='tight')
        print(f"Saved {outname}")

def _pdf_idx(path: str) -> int:
    m = re.search(r'PDF(\d+)\.csv$', path)
    return int(m.group(1)) if m else 0

def main():
    args = parse_arguments()
    input_dir = args.dir
    
    # Real space G(r) - X-ray
    fname_x = glob.glob(os.path.join(input_dir, '*_FT_XFQ1.csv'))
    if fname_x:
        plot_data(fname_x[0], 'xPDF', r'{}'.format('r ($\mathrm{\AA}$)'), 'data', args, calc_rwp=True, rwp_label_prefix="G(r) (x-ray):")

    # Real space G(r) - Neutron (PDF*.csv)
    pdf_files = sorted(glob.glob(os.path.join(input_dir, '*PDF*.csv')), key=_pdf_idx)
    for fpath in pdf_files:
        x = _pdf_idx(fpath)
        plot_title_suffix = fpath.split('.csv')[0].split('_')[-1]
        
        # Determine title and Rwp label
        if 'PDFpartials' in plot_title_suffix:
             # Partials usually don't have Rwp calculated in the same way or it's not requested in original code
             # But we can plot them. Original code commented out partials, but user might want them.
             # Let's stick to original behavior: plot but maybe no Rwp print if it's partials?
             # The original code had partials commented out. The user's code had them commented out.
             # But the user's code had a loop for *PDF*.csv.
             # Let's follow the user's recent logic:
             # "if plot_title!='PDFpartials': ... Rw = Rwp..."
             
             plot_data(fpath, plot_title_suffix, r'r ($\mathrm{\AA}$)', 'data', args, calc_rwp=False)
        else:
             tag = f"G(r) (neutron{'' if x in (0,1) else f' #{x}'})"
             plot_data(fpath, plot_title_suffix, r'r ($\mathrm{\AA}$)', 'data', args, calc_rwp=True, rwp_label_prefix=tag)

    # Reciprocal space fits. RMCProfile writes F(Q) into *_FQ1.csv (header
    # F(Q)_RMC, F(Q)_Expt; -> 0 at high Q), not S(Q): title each file by the
    # function its own header names, else by its name (FQ -> F(Q), SQ -> S(Q)).
    def _function_title(labels, default):
        match = next((re.search(r'([A-Za-z])\(([QqRr])\)', label) for label in labels[1:]
                      if re.search(r'([A-Za-z])\(([QqRr])\)', label)), None)
        return f"{match.group(1)}({match.group(2)})" if match else default

    for pattern, default in (('*_FQ1.csv', 'F(Q)'), ('*_SQ1.csv', 'S(Q)')):
        fnames_q = glob.glob(os.path.join(input_dir, pattern))
        if fnames_q:
            labels, _ = read_csv(fnames_q[0])
            xlabel = labels[0].strip() if labels else r'Q ($\mathrm{\AA^{-1}}$)'
            title = _function_title(labels, default)
            plot_data(fnames_q[0], title, xlabel, title, args, calc_rwp=True, rwp_label_prefix=f"{title}:")

    # BRAGG
    fnames = sorted(set(glob.glob(os.path.join(input_dir, '*_bragg.csv'))
                        + glob.glob(os.path.join(input_dir, '*_bragg_*.csv'))))
    if fnames:
        plot_data(fnames[0], 'BRAGG', r'Q ($\mathrm{\AA^{-1}}$)', 'data', args, calc_rwp=True, rwp_label_prefix="BRAGG:")

    # Chi-values
    def log_sort_key(path):
        name = os.path.basename(path)
        match = re.search(r'^(.+)-(\d+)\.log$', name)
        return (match.group(1).lower(), int(match.group(2))) if match else (name.lower(), -1)

    fnames = sorted(glob.glob(os.path.join(input_dir, '*-*.log')), key=log_sort_key)
    if fnames:
        chi_Q, chi_R = read_chi(fnames)
        if len(chi_R) > 0:
            fig4 = plt.figure(figsize=(3.375*2,3.375*1.2))
            dx = fig4.add_subplot(111)
            dx.plot(np.log(np.maximum(chi_R, 1e-12)),label=r'R',lw=1.0,alpha=0.5)
            dx.set_xlabel(r'Time steps',fontsize=11)
            dx.set_ylabel(r'log($\mathrm{\chi}$)',fontsize=11)
            dx.legend(loc=1,fontsize=9,frameon=False)
            fig4.suptitle('R-value', fontsize=14)
            if args.save:
                fig4.savefig(os.path.join(input_dir, 'R-value.png'), dpi=300, bbox_inches='tight')
                print(f"Saved {os.path.join(input_dir, 'R-value.png')}")

    if not args.no_show:
        plt.show()

if __name__ == "__main__":
    main()
