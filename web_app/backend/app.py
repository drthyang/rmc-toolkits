from flask import Flask, jsonify, request, send_file
from flask_cors import CORS
import os
import sys
import glob
import io
import matplotlib
matplotlib.use('Agg') # Use non-interactive backend
import matplotlib.pyplot as plt

# Add project root to sys.path to import RMC_plot
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
import RMC_plot

app = Flask(__name__)
CORS(app) # Enable CORS for frontend communication

@app.route('/api/files', methods=['GET'])
def list_files():
    directory = request.args.get('dir', './')
    try:
        # Expand user path if needed
        directory = os.path.expanduser(directory)
        if not os.path.exists(directory):
             return jsonify({'error': 'Directory not found'}), 404
        
        files = []
        # List relevant files
        patterns = ['*.csv', '*.log', '*.rmc6f']
        for pattern in patterns:
            for fpath in glob.glob(os.path.join(directory, pattern)):
                files.append({
                    'name': os.path.basename(fpath),
                    'path': os.path.abspath(fpath),
                    'type': 'file'
                })
        
        # Also list subdirectories for navigation
        for item in os.listdir(directory):
            item_path = os.path.join(directory, item)
            if os.path.isdir(item_path) and not item.startswith('.'):
                files.append({
                    'name': item,
                    'path': os.path.abspath(item_path),
                    'type': 'directory'
                })
                
        return jsonify({'files': sorted(files, key=lambda x: (x['type'] != 'directory', x['name']))})
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/plot', methods=['GET'])
def plot_file():
    fpath = request.args.get('path')
    if not fpath or not os.path.exists(fpath):
        return jsonify({'error': 'File not found'}), 404

    try:
        # Determine plot type based on filename
        filename = os.path.basename(fpath)
        
        # Clear previous figures
        plt.clf()
        plt.close('all')
        
        # Use RMC_plot logic
        # We need to simulate args object
        class Args:
            save = False
            no_show = True # We handle showing
            
        args = Args()
        
        # Capture the figure
        fig = None
        
        # Logic adapted from RMC_plot.main()
        if '_FT_XFQ1.csv' in filename:
            RMC_plot.plot_data(fpath, 'xPDF', r'r ($\mathrm{\AA}$)', 'data', args, calc_rwp=True, rwp_label_prefix="G(r) (x-ray):")
            fig = plt.gcf()
        elif 'PDF' in filename and '.csv' in filename:
            # Extract index logic
            x = RMC_plot._pdf_idx(fpath)
            plot_title_suffix = filename.split('.csv')[0].split('_')[-1]
            if 'PDFpartials' in plot_title_suffix:
                 RMC_plot.plot_data(fpath, plot_title_suffix, r'r ($\mathrm{\AA}$)', 'data', args, calc_rwp=False)
            else:
                 tag = f"G(r) (neutron{'' if x in (0,1) else f' #{x}'})"
                 RMC_plot.plot_data(fpath, plot_title_suffix, r'r ($\mathrm{\AA}$)', 'data', args, calc_rwp=True, rwp_label_prefix=tag)
            fig = plt.gcf()
        elif '_FQ1.csv' in filename:
            labels, _ = RMC_plot.read_csv(fpath)
            xlabel = labels[0].strip() if labels else r'Q ($\mathrm{\AA^{-1}}$)'
            RMC_plot.plot_data(fpath, 'S(Q) (x-ray)', xlabel, 'data', args, calc_rwp=True, rwp_label_prefix="S(Q) (x-ray):")
            fig = plt.gcf()
        elif '_SQ1.csv' in filename:
            labels, _ = RMC_plot.read_csv(fpath)
            xlabel = labels[0].strip() if labels else r'Q ($\mathrm{\AA^{-1}}$)'
            RMC_plot.plot_data(fpath, 'S(Q) (neutron)', xlabel, 'data', args, calc_rwp=True, rwp_label_prefix="S(Q) (neutron):")
            fig = plt.gcf()
        elif '_bragg.csv' in filename:
            RMC_plot.plot_data(fpath, 'BRAGG', r'Q ($\mathrm{\AA^{-1}}$)', 'data', args, calc_rwp=True, rwp_label_prefix="BRAGG:")
            fig = plt.gcf()
        elif '.log' in filename:
            # Special handling for log files as they are not single plot_data calls in original script
            # But we can reuse the logic
            chi_Q, chi_R = RMC_plot.read_chi([fpath])
            if len(chi_R) > 0:
                fig = plt.figure(figsize=(3.375*2,3.375*1.2))
                dx = fig.add_subplot(111)
                dx.plot(np.log(chi_R),label=r'R',lw=1.0,alpha=0.5)
                dx.set_xlabel(r'Time steps',fontsize=11)
                dx.set_ylabel(r'log($\mathrm{\chi}$)',fontsize=11)
                dx.legend(loc=1,fontsize=9,frameon=False)
                fig.suptitle('R-value', fontsize=14)
        
        if fig:
            img = io.BytesIO()
            fig.savefig(img, format='png', bbox_inches='tight', dpi=150)
            img.seek(0)
            return send_file(img, mimetype='image/png')
        else:
            return jsonify({'error': 'Could not generate plot for this file type'}), 400

    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True, port=5000)
