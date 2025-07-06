# reproduce figure 5

import h5py
import os
import glob
import numpy as np
from Plotting import init_plot
import matplotlib.pyplot as plt


# read data

output_dir = "../outputs/local_data_reduced"

all_files = sorted(glob.glob(os.path.join(output_dir, "*.h5")))



# read planetesimals
with h5py.File(all_files[0], 'r') as f:

    tfor = f['tfor'][...] # formation of planetesimals
    rigid_plt = f['rigid_plt'][...]

Stmin, Stmax = 0.01, 1
rmin, rmax = 6.09, 6.11
pebbles_fri = np.empty((len(all_files),))
t_arr = np.empty((len(all_files),))
# read pebble data
for i, file_path in enumerate(all_files[1:]):
    print(i)
    with h5py.File(file_path, 'r') as f:

        # Read compound dataset from swarms/swarmsout
        t_arr[i] = f['snapt'][...]
        r = f['cylindrical radius [AU]'][...]
        rigid = f['rigid'][...]
        St = f['Stokes number'][...]
        bool_r = (r>rmin)&(r<rmax)
        bool_St = (St>Stmin)&(St<Stmax)
        pebbles_fri[i] = np.mean(rigid[(bool_r)&(bool_St)])
        

# make plot


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 12), sharey=True)

xlabel = 'Matrix mass fraction'
ylabel = 'Time (Myr after CAIs)'

t_CAI = 0.19 # time difference between disk formation and CAI formation, in Myr
init_plot(ax, title='', xlabel=xlabel, ylabel=ylabel, font = 30)


# plot planetesimals

for it in range(len(t_arr)):
    values = rigid_plt[(tfor[:,0]==t_arr[it])]
    if len(values)>100: # to avoid data points with little statistical 
                        # significance (setting it to 0 won't change overall plot)

        ax.scatter(1-np.mean(values), 0.5*(t_arr[it]+t_arr[it-1])/10**6-t_CAI,  # plt formed betweem  t_arr[i-1] and t_arr[I], we take the mean value as the formation time
                   c='#7B5EA4', s=100,
                   alpha=1, zorder=0)

ax.scatter([], [],
           c='#7B5EA4', s=100, label=r'Planetesimals',
           alpha=1, zorder=0)

    
# plot pebbles

ax.scatter(1-pebbles_fri, t_arr/10**6-t_CAI,
            color='#FCC98A', s=100, zorder=-1, marker="s", label='Pebbles')



# compare with meteoritic data, from Hellmann et al 2023

CC_label = [r'$\rm{CR}$', r'$\rm{CO}$',
            r'$\rm{CV}$', r'$\rm{CM}$',
            r'$\rm{TL}$', r'$\rm{CI}$']

tacc = np.array([[3.5, 4.], [2.2, 2.6], [2.3, 3.1],
                 [2.5, 4.2], [3, 4.0], [3.1, 4.1]]).T*10**6

fri = np.array([[0.88, 0.98], [0.68, 0.92], [0.62, 0.77],
                [0.48, 0.7], [0.15, 0.6], [0., 0.]]).T

ax.errorbar(1-fri.mean(axis=0), tacc.mean(axis=0),
               xerr=fri[1,:]-fri.mean(axis=0), 
               yerr=(tacc[1,:]-tacc.mean(axis=0)), 
               ecolor='k', c='k', elinewidth=3,fmt='s')

for i, label in enumerate(CC_label):
    if i ==1 or i ==2 or i ==3 or i==4 or i==5 or i==0:
        ax.text(y=(tacc[:, i].mean()) + 0.05,
                x=(1-fri[:, i].mean(axis=0) + 0.01), s=label)



# set axis limits
ax.set_ylim(2, 4.3)
ax.set_xlim(-0.1, 1.08)

# set legends
ax.legend(ncol=1, fontsize=25)
