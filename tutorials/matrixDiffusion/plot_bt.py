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

import scipy.io as spio
plt.style.use('seaborn')




# Import data from breakthrough
codeRes = np.loadtxt(open("./breakthrough.dat", "rb"), delimiter=" ", skiprows=0)


tau = [0.5, 2.0, 3.0, 4.0]
pl = np.power(tau,-3/2)
plsp = np.power(tau,-1/2)

plt.figure(figsize=(250 /25.4, 200 / 25.4))
plt.loglog  (codeRes[:,0],codeRes[:,1],'b-',label='mrmtScalarTransportFoam')
plt.loglog  (tau,pl,'r-',label='power law')
plt.loglog  (tau,plsp,'b-',label='power law 1/2')
plt.ylabel('$c$',fontsize=30)
plt.xlabel('$t$',fontsize=30)
plt.xticks(fontsize=28)
plt.yticks(fontsize=28)
#plt.show()
plt.legend(loc='upper right', frameon=False,fontsize=28)

#    plt.rc('grid', linestyle="-", color='black')
plt.grid(True)
#    plt.tight_layout()
plt.savefig('bt.pdf')
