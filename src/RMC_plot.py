#!/usr/bin/env python3
import numpy as np
import sys, glob, re, os
import argparse
import matplotlib.pyplot as plt
from matplotlib import rc

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
    chi_Q = []
    chi_R = []
    for ii in np.arange(len(fnames)) :
        try:
            f = open(fnames[ii],'r')
            lines = f.readlines()
            f.close()
            for jj in np.arange(2,len(lines),1) :
                tmp = lines[jj].split()
                chi_Q.append(np.float64(tmp[-2]))
                chi_R.append(np.float64(tmp[-1]))
        except Exception as e:
            print(f"Error reading {fnames[ii]}: {e}")
    chi_R = np.array(chi_R)
    chi_Q = np.array(chi_Q)
    return chi_Q, chi_R

def Rwp(r, gr, grfit, fit_range) :
    idx = [ ii for ii in np.arange(len(r)) if (r[ii]>=fit_range[0]) and (r[ii]<=fit_range[-1]) ]
    if len(idx) == 0:
        return 0.0
    gr_sub = np.array(gr[idx])
    grfit_sub = np.array(grfit[idx])
    denom = np.sum( gr_sub*gr_sub )
    if denom == 0:
        return 0.0
    Rsq = np.sum( (grfit_sub-gr_sub) * (grfit_sub-gr_sub)  ) / denom
    return np.sqrt(Rsq)

def plot_data(fname, title, xlabel, ylabel, args, calc_rwp=False, rwp_label_prefix=""):
    if not fname:
        return
    
    labels, data = read_csv(fname)
    if len(labels) == 0:
        return

    if calc_rwp and len(data) >= 3:
        fit_range = [data[0][0],data[0][-1]]
        Rw = Rwp(data[0],data[1],data[2],fit_range)
        print(f"{rwp_label_prefix:<20} R = {Rw:.6f}")

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

    # Reciprocal space S(Q) - X-ray
    fname_x = glob.glob(os.path.join(input_dir, '*_FQ1.csv'))
    if fname_x:
        labels, _ = read_csv(fname_x[0]) # Read just to get label if needed, but we can just use generic
        xlabel = labels[0].strip() if labels else r'Q ($\mathrm{\AA^{-1}}$)'
        plot_data(fname_x[0], 'S(Q) (x-ray)', xlabel, 'data', args, calc_rwp=True, rwp_label_prefix="S(Q) (x-ray):")

    # Reciprocal space S(Q) - Neutron
    fname_n = glob.glob(os.path.join(input_dir, '*_SQ1.csv'))
    if fname_n:
        labels, _ = read_csv(fname_n[0])
        xlabel = labels[0].strip() if labels else r'Q ($\mathrm{\AA^{-1}}$)'
        plot_data(fname_n[0], 'S(Q) (neutron)', xlabel, 'data', args, calc_rwp=True, rwp_label_prefix="S(Q) (neutron):")

    # BRAGG
    fnames = glob.glob(os.path.join(input_dir, '*_bragg.csv'))
    if fnames:
        plot_data(fnames[0], 'BRAGG', r'Q ($\mathrm{\AA^{-1}}$)', 'data', args, calc_rwp=True, rwp_label_prefix="BRAGG:")

    # Chi-values
    fnames = glob.glob(os.path.join(input_dir, '*-*.log'))
    if fnames:
        chi_Q, chi_R = read_chi(fnames)
        if len(chi_R) > 0:
            fig4 = plt.figure(figsize=(3.375*2,3.375*1.2))
            dx = fig4.add_subplot(111)
            dx.plot(np.log(chi_R),label=r'R',lw=1.0,alpha=0.5)
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
