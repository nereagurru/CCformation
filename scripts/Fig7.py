#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 19:50:24 2026

@author: gurrutxaga
"""

import h5py
import os
import glob
import numpy as np
from Plotting import init_plot
import matplotlib.pyplot as plt
from astropy import units as u


output_dir = "../outputs/local_data_reduced"



all_files = sorted(glob.glob(os.path.join(output_dir, "*.h5")))
# colors 
color = {'pebble':'#FCC98A', 'planetesimal': 'blueviolet', 'observation':'grey'}

# read planetesimals
with h5py.File(all_files[0], 'r') as f:

    tfor = f['tfor'][...][:,0] # formation of planetesimals
    rigid_plt = f['fri_plt'][...] # one could also plot f['rigid_plt'][...]

Stmin, Stmax = 0.01, 1

pebbles_fri = np.empty((len(all_files),))

t_arr = np.empty((len(all_files),))
mswarm_arr = np.empty((len(all_files),))
nrl = 100

from_AU_to_cm = (1*u.AU).to(u.cm).value



# read pebble data
for i, file_path in enumerate(all_files[1:]):
    with h5py.File(file_path, 'r') as f:

        # Read compound dataset from swarms/swarmsout
        t_arr[i] = f['snapt'][...]
        rigid = f['fri'][...] # one could also plot f['rigid'][...]
        St = f['Stokes number'][...]
        # Read compound dataset from swarms/swarmsout
        r_arr = f['cylindrical radius [AU]'][...]

        # find dust trap; where dust surface density is maximum
        mswarm = f.attrs['mass_of_swarm']
        sorted_indices = np.argsort(r_arr)
        sorted_positions = r_arr[sorted_indices]

        sigmadl = np.empty((nrl,))
        r_evall = np.empty((nrl,))
        drdis = np.empty((nrl,))
        jump_ = int(len(r_arr)/nrl)
        for j in range(0, nrl):
            r_evall[j] = np.mean(sorted_positions[jump_*j:jump_*(j+1)])
            drdis[j] = sorted_positions[jump_*(j+1)-1]-sorted_positions[jump_*j]


        sigmadl[:] = jump_*mswarm/(2*np.pi*r_evall[:]*drdis[:])/from_AU_to_cm**2
        arg_max = sigmadl.argmax()
        rmin, rmax = r_evall[arg_max]-0.5*drdis[arg_max], r_evall[arg_max]+0.5*drdis[arg_max]
        bool_r = (r_arr>rmin)&(r_arr<rmax)
        bool_St = (St>Stmin)&(St<Stmax)
        pebbles_fri[i] = np.mean(rigid[(bool_r)&(bool_St)])
        mswarm_arr[i] = mswarm
        
        
Mplt_arr = np.zeros((len(all_files),))
fri_arr = np.zeros((len(all_files),))
for it in range(len(t_arr)):
    values = rigid_plt[(tfor==t_arr[it])]
    Mplt_arr[it] = mswarm_arr[it]*len(values)
    fri_arr[it] = np.mean(values)
    
    
    
    

idxs = np.searchsorted(t_arr, tfor)


left  = t_arr[idxs - 1]
right = t_arr[idxs]
idxs -= (tfor - left <= right - tfor)

unique_idxs, inverse = np.unique(idxs, return_inverse=True)

mass_values = np.zeros(len(unique_idxs))

for k, uu in enumerate(unique_idxs):
    with h5py.File(all_files[uu], 'r') as f:
        mass_values[k] = f.attrs['mass_of_swarm']

Mplt = mass_values[inverse] 

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 12), sharey=True)

xlabel = 'Matrix mass fraction'
ylabel = 'Planetesimal mass (Earth mass)'

t_CAI = 0.19 # time difference between disk formation and CAI formation, in Myr
init_plot(ax, title='', xlabel=xlabel, ylabel=ylabel, font = 30)

# fraction ranges


range_ = np.arange(0, 101, 5)*0.01



y_start = 1.35
step = 0.04
fm = np.array([[0.09-0.05, 0.09+0.05], [0.2-0.12, 0.2+0.12], [0.3-0.07, 0.3+0.07],
                [0.41-0.12, 0.41+0.12], [0.64-0.23, 0.64+0.23], [1.-0.02, 1.]]).T

# compare with meteoritic data, from Hellmann et al 2023

CC_label = [r'$\rm{CR}$', r'$\rm{CO}$',
            r'$\rm{CV}$', r'$\rm{CM}$',
            r'$\rm{TL}$', r'$\rm{CI}$']

for i, (fmmin, fmmax) in enumerate(zip(fm[0], fm[1])):
    ax.axvspan(fmmin, fmmax, ymin=y_start-step, ymax=y_start, color='gray', alpha=0.8,
               clip_on=False)
    
    ax.text(s=CC_label[i], x=-0.1, y=y_start-step, color='k',
            ha='center', va='bottom', transform=ax.get_xaxis_transform())
    y_start -=step + 0.01



# add potential of outer disk

# read global data and calculate its dust surface denisty


output_dir = "../outputs/global_data_reduced"

all_files_out = sorted(glob.glob(os.path.join(output_dir, "*.h5")))[1:]



t_arr_out = np.empty((len(all_files_out),))



for i, file_path in enumerate(all_files_out):

    with h5py.File(file_path, 'r') as f:

        # Read compound dataset from swarms/swarmsout
        t_arr_out[i] = f['snapt'][...]


    
arg = np.argmin(np.absolute(t_arr_out-t_CAI*10**6-tfor.max()))


with h5py.File(all_files_out[arg], 'r') as f:

    # Read compound dataset from swarms/swarmsout
    r_arr = f['cylindrical radius [AU]'][...]
    rigid = f['rigid'][...]
    mswarm_out = f.attrs['mass_of_swarm']
    
boolr = (r_arr >6.5)&(rigid==0) # no pure rigids; they do not get trapped in the same loc as varying sizes


ax.hist([1.-fri_arr, 1.-rigid[boolr]], bins=range_,
        weights=[Mplt_arr * u.g.to(u.Mearth), np.full(np.sum(boolr), mswarm_out* u.g.to(u.Mearth))],
        edgecolor='lightgray', lw=2, color=[color['planetesimal'], 'lightskyblue'], stacked=True, alpha=0.8)
Mplt_tot = np.nansum(Mplt_arr)*u.g.to(u.Mearth)
ax.text(x=0, y=1.1*Mplt_tot, s='total planetesimal mass', color=color['planetesimal'])

xmin, xmax = -0.05, 1.05

ax.hlines(xmin=xmin, xmax=xmax, y=Mplt_tot,
          lw=5, color=color['planetesimal'])

# asteroid belt
text = ax.text(x=0, y=1.15*(5e-4), s='asteroid belt', color='k')
text.set_bbox(dict(facecolor='white', alpha=0.5))
ax.hlines(xmin=xmin, xmax=xmax, y=5e-4,
          lw=5, color='k', zorder=10)

# asteroid belt
text = ax.text(x=0, y=0.0225, s='Kuiper belt', color='k')
text.set_bbox(dict(facecolor='white', alpha=0.5))
ax.hlines(xmin=xmin, xmax=xmax, y=0.02,
          lw=5, color='k', zorder=10)

#text.
ax.set_yscale('log')
ax.set_ylim(10**-4, 1.)

ax.set_xlim(xmin,xmax)
fig.savefig('planetesimal_mass.pdf', dpi=300, bbox_inches='tight', format='pdf')
