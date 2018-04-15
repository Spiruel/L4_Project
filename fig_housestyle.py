import numpy as np

def figsize(scale):
	fig_width_pt = 469.75502                          # Get this from LaTeX using \the\textwidth
	inches_per_pt = 1.0/72.27                       # Convert pt to inch
	golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
	fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
	fig_height = fig_width*golden_mean              # height in inches
	fig_size = [fig_width,fig_height]
	return fig_size

import matplotlib as mpl
#http://bkanuka.com/articles/native-latex-plots/
params = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "axes.labelsize": 12,               # LaTeX default is 10pt font.
    "font.size": 10,
    "legend.fontsize": 10,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "figure.figsize": figsize(0.95),
    }
mpl.rcParams.update(params)
font = {'family' : 'monospace'}
print 'set housestyle for all figures!'
