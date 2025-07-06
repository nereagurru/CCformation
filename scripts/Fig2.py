# to reproduce figure 2
# before running the script, copy "disk_data/*.dat" and "disk_data/*.info" files to ../disk/.

import matplotlib.pyplot as plt
import numpy as np
from Plotting import init_plot
from astropy import units as u
from scipy import interpolate

# where to find disk data
testfile = '../disk/'

# read data from testfile
t = (np.loadtxt(testfile + 'time.dat',dtype=float)*u.s).to(u.Myr) # time in Myr
# read position grid and remove first element
r = np.loadtxt(testfile + 'grid.info',dtype=float)
n = int(r[0])
r = (r[1::]*u.cm).to(u.AU) # position in AU

# read gas surface density over time and across disk
sigma_dat = np.loadtxt(testfile+'sigma.dat', dtype=float,
                   skiprows=2).reshape((n, t.shape[0]), order='F')*u.g/u.cm**2

# create interpolation function for gas surface density
f_sigma = (interpolate.interp2d(t, r, sigma_dat))

# time difference between t0 of the disk and CAI formation
t_CAI = 0.19*u.Myr

# time to evaluate with respect to CAI formation
t_eval = np.array([2.3, 3.5, 3.98, 4.06, 4.5])*u.Myr # time in Myr

# create plot
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 10))
init_plot(ax, None, 'Radial distance (au)', 
          r'Gas surface density (g cm$^{-2}$)', font=30)

# color of lines
color = ['#8BBCCC', '#3E6D9C', '#6D67E4', 'darkred', '#000000'] # '#5C2E7E', 

# order of lines
z = [len(t_eval)-i for i in range(0, len(t_eval))]

# plot r vs. sigma at different times
for i, c in enumerate(color):
    ax.plot(r, f_sigma(t_eval[i]+t_CAI,r), c=color[i],
            label=rf'${t_eval[i].value:.2f}$' + r'$\,\mathrm{Myr}$',
            linewidth=6, zorder=z[i])
   
# set logarithmic scale
ax.set_xscale('log')
ax.set_xlim(1, 120)

# set limits
ax.set_yscale('log')
ax.set_ylim(0.001, 10**3)
ax.set_yticks([0.001, 0.01, 0.1, 1, 10, 100, 1000])
plt.subplots_adjust(left=0.15, bottom=0.2, top=0.8, right=0.73)

# set legend
fig.legend(fontsize=20, ncol=1, bbox_to_anchor=(0.73, 0.46))
plt.show()
