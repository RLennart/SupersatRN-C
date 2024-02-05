import matplotlib as mpl
import  numpy as np;

def figsize(scale=1.3):
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            
    fig_width = 345./72.27*scale   
    fig_height = fig_width*golden_mean
    fig_size = [fig_width,fig_height*2]
    return fig_size

rc = {                      
    "pgf.texsystem": "pdflatex",        
    "text.usetex": False,               
    "font.family": ['Tahoma', 'DejaVu Sans', 'Verdana'],
    "font.serif": [],                   
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 18,  
    "axes.labelpad": 10,               
    "font.size": 15,
    "figure.subplot.left": 0.125,
    "figure.subplot.bottom": 0.175,
    "figure.autolayout" : True,
    "legend.fontsize": 18,   
    "ytick.direction" : "in",
    "xtick.direction" : "in",
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "lines.linewidth" : 3,
    "lines.markersize" : 4,
    "axes.linewidth" : 2,
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.width": 2,
    "ytick.major.width": 2,
    "xtick.minor.width": 2,
    "ytick.minor.width": 2,
    "xtick.major.size": 10,
    "ytick.major.size": 10,
    "xtick.minor.size": 7,
    "ytick.minor.size": 7,
    "figure.figsize": figsize(),  
    "legend.fontsize": 15,
    "legend.framealpha"    : 1.0,
    "legend.shadow"    : False,
    "legend.borderpad" : 0.25,
    "legend.borderaxespad" : 0.75,     
    "axes.grid"     : False, 
    "axes.prop_cycle"  : mpl.cycler('color',['black','firebrick', 'steelblue', 'forestgreen','orange']),

    }
mpl.rcParams.update(rc)

