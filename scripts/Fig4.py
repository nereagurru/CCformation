# reproduce figure 4
# before running the script, copy "local_data_reduced" directory to "../outputs/."
import h5py
import os
import glob
import numpy as np
from Plotting import init_plot
import matplotlib.pyplot as plt
from astropy import units as u



# read local data and calculate dust surface denisty

output_dir = "../outputs/local_data_reduced"

all_files = sorted(glob.glob(os.path.join(output_dir, "*.h5")))

t_arr = np.empty((len(all_files[1:]),))

# first save time values. Then we will decide the closest ones to our evaluation times
for i, file_path in enumerate(all_files[1:]):

    with h5py.File(file_path, 'r') as f:

        # Read compound dataset from swarms/swarmsout
        t_arr[i] = f['snapt'][...]

# time for evaluating dust surface density
t_CAI = 0.19 # value in Myr
t_eval = np.array([2., 2.3, 3.5])
t_eval += t_CAI

nrl = 100

from_AU_to_cm = (1*u.AU).to(u.cm).value

sigmadl = np.empty((nrl, t_eval.shape[0]))
sigmadl_fri = np.empty((nrl, t_eval.shape[0]))
#drdis = (rl[1:]-rl[:-1]) # size of bin
#r_evall = 0.5*(rl[1:]+rl[:-1]) # radial distance within bin
drdis = np.empty((nrl, t_eval.shape[0]))
r_evall = np.empty((nrl, t_eval.shape[0]))

for i in range(len(t_eval)):
    

    # what file is the nearest to the time we want to evaluate?
    file_id = np.argmin(np.absolute(t_arr-t_eval[i]*10**6))

    with h5py.File(all_files[1:][file_id], 'r') as f:

        # Read compound dataset from swarms/swarmsout
        r_arr = f['cylindrical radius [AU]'][...]

        rigid = f['rigid'][...]
        mswarm = f.attrs['mass_of_swarm']
        sorted_indices = np.argsort(r_arr)
        sorted_positions = r_arr[sorted_indices]
        bins = np.empty((nrl,))
        rcounts_fri = np.empty((nrl,))
        jump_ = int(len(r_arr)/nrl)
        for j in range(0, nrl):
            r_evall[j,i] = np.mean(sorted_positions[jump_*j:jump_*(j+1)])
            drdis[j,i] = sorted_positions[jump_*(j+1)-1]-sorted_positions[jump_*j]
            rcounts_fri[j] = np.sum(rigid[sorted_indices][jump_*j:jump_*(j+1)])
        #rcounts,edges = np.histogram(sorted_positions, bins=bins) # number of swarms per bin
        #rcounts_fri,_ = np.histogram(sorted_positions, bins=bins,
        #                              weights=rigid[sorted_indices])


        

        sigmadl[:, i] = jump_*mswarm/(2*np.pi*r_evall[:,i]*drdis[:,i])/from_AU_to_cm**2
        sigmadl_fri[:, i] = rcounts_fri*mswarm/(2*np.pi*r_evall[:,i]*drdis[:,i])/from_AU_to_cm**2


# read global data and calculate its dust surface denisty


output_dir = "../outputs/global_data_reduced"

all_files = sorted(glob.glob(os.path.join(output_dir, "*.h5")))



t_arr = np.empty((len(all_files),))


for i, file_path in enumerate(all_files):

    with h5py.File(file_path, 'r') as f:

        # Read compound dataset from swarms/swarmsout
        t_arr[i] = f['snapt'][...]



nrg = 30
rg = np.logspace(np.log10(6.5), np.log10(120), nrg)
sigmadg = np.empty((rg.shape[0]-1, t_eval.shape[0]))
sigmadg_fri = np.empty((rg.shape[0]-1, t_eval.shape[0]))
drdis = rg[1:]-rg[:-1]
r_eval = 0.5*(rg[1:]+rg[:-1])

for i in range(len(t_eval)):
    

    # what file is the nearest to the time we want to evaluate?
    file_id = np.argmin(np.absolute(t_arr-t_eval[i]*10**6))

    with h5py.File(all_files[file_id], 'r') as f:

        # Read compound dataset from swarms/swarmsout
        r_arr = f['cylindrical radius [AU]'][...]
        rigid = f['rigid'][...]
        mswarm = f.attrs['mass_of_swarm']
        rcounts,_ = np.histogram(r_arr, bins=rg) # number of swarms per bin
        rcounts_fri,_ = np.histogram(r_arr, bins=rg,
                                      weights=rigid)


        sigmadg[:, i] = rcounts*mswarm/(2*np.pi*r_eval*drdis)/from_AU_to_cm**2
        sigmadg_fri[:, i] = rcounts_fri*mswarm/(2*np.pi*r_eval*drdis)/from_AU_to_cm**2



from scipy import interpolate
# --- indicate minimum sigmad for planetesimal formation
# where to find disk data
testfile = '../disk/PE/disk_midalpha/gap/shallow/'

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
Zcrit = 0.05



# make plot
#%%

fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(25, 18), sharey=True)

# color to indicate fragile and rigid
colors = {'fragile':'#3274B5', 'rigid':'#8B551B', 'planetesimal':'blueviolet'}
ls_arr = ['-', '-', '-']
lw = 8
#Label panel 
labels = ['A', 'B', 'C']
for idx, ls in enumerate(ls_arr):
    title = 'Local simulation' if idx==0 else None
    xlabel = 'Radial distance (AU)' if idx==2 else None
    init_plot(ax[idx, 0], title=title,
              xlabel=xlabel, font=30)
    title = 'Global simulation' if idx==0 else None
    init_plot(ax[idx, 1], title=title,
              xlabel=xlabel, font=30)

    ax[idx, 0].plot(r_evall[:,idx],
               (sigmadl[:, idx]-sigmadl_fri[:, idx]),
               c=colors['fragile'], ls=ls, lw=lw)
    

    
    ax[idx, 1].plot(r_eval,
               (sigmadg[:, idx]-sigmadg_fri[:, idx]),
               colors['fragile'], ls=ls, lw=lw)
    
    

    ax[idx,0].plot(r_evall[:,idx],
               sigmadl_fri[:, idx], 
               ls=ls, c=colors['rigid'], lw=lw)
    
    ax[idx,1].plot(r_eval,
               sigmadg_fri[:, idx],
               ls=ls, c=colors['rigid'], lw=lw)
    #ax[idx,0].plot(r_evall[:,idx],
    #           sigmadl[:, idx], 
    #           ls=ls, c='k', lw=lw)
    ax[idx,0].plot(r_evall[:,idx], Zcrit*f_sigma(t_eval[idx], r_evall[:,idx]), alpha=0.5, c=colors['planetesimal'], 
               lw=lw/2, ls=ls, zorder=-1)

    # indicate time
    ax[idx,1].text(0.75, 0.92, f'{t_eval[idx]-t_CAI:.1f} Myr', transform=ax[idx,1].transAxes,
                   fontsize=30, va='top', ha='left')
    
    

    ax[idx,0].text(0.02, 0.95, labels[idx], transform=ax[idx,0].transAxes,
            fontsize=34, fontweight='bold', va='top', ha='left')
    
# text for planetesimal formation
ax[0,0].text(0.35, 0.87, 'Planetesimal formation threshold', transform=ax[0,0].transAxes,
        fontsize=25, va='top', ha='left', c=colors['planetesimal'])

ymin, ymax = 10**-5, 20

rmax = r_evall[1:-1, 1][sigmadl[1:-1, 1].argmax()]
ax[1,0].axvline(ymin=ymin, ymax=ymax, x=rmax, c=colors['planetesimal'], lw=lw/2, ls='-', alpha=0.5, zorder=-1)

rmax = r_evall[1:-1, 2][sigmadl[1:-1, 2].argmax()]
ax[2,0].axvline(ymin=ymin, ymax=ymax, x=rmax, c=colors['planetesimal'], lw=lw/2, ls='-', alpha=0.5, zorder=-1)



ax[0,1].plot([],[], lw=lw, label='Fragile', c=colors['fragile'])
ax[0,1].plot([],[], lw=lw, label='Rigid', c=colors['rigid'])
ax[0,1].legend(ncol=2, loc='upper left')

for i in range(0, t_eval.shape[0]):
    ax[i,0].set_xlim(r_evall[:,1].min(), r_evall[:,1].max())
    ax[i,1].set_xlim(6.5, 120)
    ax[i,1].set_xscale('log')



ax[0,1].set_yscale('log')

ax[0,1].set_ylim(ymin, ymax)






plt.subplots_adjust(wspace=0.03)

fig.supylabel('Dust surface density (g cm$^{-2}$)', y=0.5, x=0.05, fontsize=30)


plt.show()

fig.savefig(f'local_and_global.pdf', dpi=300, bbox_inches='tight', format='pdf')
