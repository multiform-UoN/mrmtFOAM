import itertools
import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator
import seaborn as sns
import sys

#sns.set(style='white')
#sns.set()

rcParams.update({'figure.autolayout': True})
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True



# Dims given in [width, height]
fig_dims = [7.33, 5] # height in inches, single column figure
double_dims = [6.66, 2.058] # double width fig
si_dims = [6.66, 4.116] # double size figure

figsize = {'paper':fig_dims,
           'si':si_dims, 
           'double':double_dims}
prefix = {'paper':'',
          'si':'SI_',
          'double':'d_'}
fig_locs = ['paper', 'si']

cols = sns.color_palette('Paired')

beltv_all = ['1mm','3mm','5mm','10mm','25mm','50mm','100mm']
beltv_lim = ['1mm','3mm','5mm','10mm','25mm','50mm']

#colours = {'1mm':cols[1],
#           '3mm':cols[3],
#           '5mm':cols[5],
#           '10mm':cols[7],
#           '25mm':cols[9],
#           '50mm':cols[11],
#           '100mm':cols[0],   
#       }

           
D=['0.125','0.25','0.5']

linestyles = {'1mm':'--', '100mm':'--','0.125':'--','0.25':'--','0.5':'--'}
markers = {'1mm':'o', '3mm':'v','5mm':'s', '10mm':'D','25mm':'>', '50mm':'<','100mm':'^','0.125':'o','0.25':'v','0.5':'s'}

linewidth = 1.1
markersize = 8.0

nullloc = plt.NullLocator()

def smooth(x, y, sigma=None):
    t = np.linspace(0, 1, len(x))
    t2 = np.linspace(0, 1, 100)

    x2 = np.interp(t2, t, x)
    y2 = np.interp(t2, t, y)

    if sigma is None:
        sigma = 1
    x3 = gaussian_filter1d(x2, sigma)
    y3 = gaussian_filter1d(y2, sigma)

    return x3, y3
    

