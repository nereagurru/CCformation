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
t_CAI = 0.19
t_eval = np.array([2., 2.3, 3.5])
t_eval += t_CAI

nrl = 16
rl = np.linspace(5.725, 6.5, nrl)
from_AU_to_cm = (1*u.AU).to(u.cm).value

sigmadl = np.empty((rl.shape[0]-1, t_eval.shape[0]))
sigmadl_fri = np.empty((rl.shape[0]-1, t_eval.shape[0]))
drdis = (rl[1:]-rl[:-1]) # size of bin
r_evall = 0.5*(rl[1:]+rl[:-1]) # radial distance within bin

for i in range(len(t_eval)):
    

    # what file is the nearest to the time we want to evaluate?
    file_id = np.argmin(np.absolute(t_arr-t_eval[i]*10**6))

    with h5py.File(all_files[1:][file_id], 'r') as f:

        # Read compound dataset from swarms/swarmsout
        r_arr = f['cylindrical radius [AU]'][...]
        rigid = f['rigid'][...]
        mswarm = f.attrs['mass_of_swarm']
        rcounts,_ = np.histogram(r_arr, bins=rl) # number of swarms per bin
        rcounts_fri,_ = np.histogram(r_arr, bins=rl,
                                      weights=rigid)


        sigmadl[:, i] = rcounts*mswarm/(2*np.pi*r_evall*drdis)/from_AU_to_cm**2
        sigmadl_fri[:, i] = rcounts_fri*mswarm/(2*np.pi*r_evall*drdis)/from_AU_to_cm**2


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



# make plot


fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(20, 10), sharey=True)
init_plot(ax[0], title='Local simulation', ylabel='Dust surface density (g cm$^{-2}$)', font=30)
init_plot(ax[1], title='Global simulation', font=30)

ls_arr = [':', '--', '-']
lw = 5
for idx, ls in enumerate(ls_arr):

    ax[0].plot(r_evall,
               (sigmadl[:, idx]-sigmadl_fri[:, idx]),
               '#3274B5', ls=ls, lw=lw)
    

    
    ax[1].plot(r_eval,
               (sigmadg[:, idx]-sigmadg_fri[:, idx]),
               '#3274B5', ls=ls, lw=lw)
    

    ax[0].plot(r_evall,
               sigmadl_fri[:, idx], 
               ls=ls, c='#8B551B', lw=lw)
    
    ax[1].plot(r_eval,
               sigmadg_fri[:, idx],
               ls=ls, c='#8B551B', lw=lw)

for i in range(0, t_eval.shape[0]):
    ax[1].plot([],[], lw=lw, ls=ls_arr[i],label=f'{t_eval[i]-0.19:.1f} Myr', c='k')

ymin, ymax = 10**-5, 5

rmax = r_evall[1:-1][sigmadl[1:-1, 1].argmax()]
ax[0].axvline(ymin=ymin, ymax=ymax, x=rmax, c='lightgrey', lw=5, ls='-', zorder=-1)


ax[1].plot([],[], lw=lw, label='Fragile', c='#3274B5')
ax[1].plot([],[], lw=lw, label='Rigid', c='#8B551B')

ax[0].set_xlim(5.8, 6.5)
ax[1].set_xlim(6.5, 120)
ax[1].legend(ncol=2)
ymin, ymax = 10**-5, 5

ax[1].set_yscale('log')
ax[1].set_xscale('log')
ax[1].set_ylim(ymin, ymax)


fig.supxlabel('Radial distance (AU)', y=0., fontsize=30)
plt.subplots_adjust(wspace=0.)
plt.show()
