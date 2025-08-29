# reproduce figure 6 (or Supplementary figure 2)
# before running the script, copy "local_data_reduced/" directory to ../outputs/.
# for supplementary plot, copy "local_topheavy_data_reduced/" directory to ../outputs/.

import h5py
import os
import glob
import numpy as np
from Plotting import init_plot
import matplotlib.pyplot as plt
from astropy import units as u


main_manuscript = True  # if True, reproduce figure 6 from main manuscript
                        # if False, reproduce figure 2 from supplementary material
            
# read data

output_dir = "../outputs/local_data_reduced" if main_manuscript else "../outputs/local_topheavy_data_reduced"




all_files = sorted(glob.glob(os.path.join(output_dir, "*.h5")))


# read planetesimals
with h5py.File(all_files[0], 'r') as f:

    tfor = f['tfor'][...] # formation of planetesimals
    rigid_plt = f['fri_plt'][...] # one could also plot f['rigid_plt'][...]

Stmin, Stmax = 0.01, 1

pebbles_fri = np.empty((len(all_files),))

t_arr = np.empty((len(all_files),))

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

        

# make plot


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 12), sharey=True)

xlabel = 'Matrix mass fraction'
ylabel = 'Time (Myr after CAIs)'

t_CAI = 0.19 # time difference between disk formation and CAI formation, in Myr
init_plot(ax, title='', xlabel=xlabel, ylabel=ylabel, font = 30)


# colors 
color = {'pebble':'#FCC98A', 'planetesimal': 'blueviolet', 'observation':'grey'}
# plot planetesimals

for it in range(len(t_arr)):
    values = rigid_plt[(tfor[:,0]==t_arr[it])]
    if len(values)>100: # to avoid data points with little statistical 
                        # significance (setting it to 0 won't change overall plot)

        ax.scatter(1-np.mean(values), 0.5*(t_arr[it]+t_arr[it-1])/10**6-t_CAI,  # plt formed betweem  t_arr[i-1] and t_arr[I], we take the mean value as the formation time
                   c=color['planetesimal'], s=100,
                   alpha=1, zorder=0)
ax.scatter([], [],
           c=color['planetesimal'], s=100, label=r'Planetesimals',
           alpha=1, zorder=0)

    
# plot pebbles

ax.scatter(1-pebbles_fri, t_arr/10**6-t_CAI,
            color=color['pebble'], s=100, zorder=-1, marker="s", label='Pebbles')


# compare with meteoritic data, from Hellmann et al 2023

CC_label = [r'$\rm{CR}$', r'$\rm{CO}$',
            r'$\rm{CV}$', r'$\rm{CM}$',
            r'$\rm{TL}$', r'$\rm{CI}$']

tacc = np.array([[3.5, 4.], [2.2, 2.6], [2.3, 3.1],
                 [2.5, 4.2], [3, 4.2], [3.1, 4.1]]).T


fm = np.array([[0.09-0.05, 0.09+0.05], [0.2-0.12, 0.2+0.12], [0.3-0.07, 0.3+0.07],
                [0.41-0.12, 0.41+0.12], [0.64-0.23, 0.64+0.23], [1., 1.]]).T
ax.errorbar(fm.mean(axis=0), tacc.mean(axis=0),
                xerr=fm[1,:]-fm.mean(axis=0), 
                yerr=(tacc[1,:]-tacc.mean(axis=0)), 
                ecolor=color['observation'], c=color['observation'],
                elinewidth=3,fmt='s', zorder=-1)

for i, label in enumerate(CC_label):
    if i == 1:
        ax.text(y=(tacc[:, i].mean()) + 0.05,
                x=(fm[:, i].mean(axis=0))-0.09, s=label, c=color['observation'])
    if i ==0 or i ==2 or i ==3 or i==4 or i==5:
        ax.text(y=(tacc[:, i].mean()) + 0.05,
                x=(fm[:, i].mean(axis=0))+0.01, s=label, c=color['observation'])


# set axis limits
ax.set_ylim(2, 4.3)
ax.set_xlim(-0.1, 1.08)

# set legends
ax.legend(ncol=1, fontsize=25)

filename = 'plt.pdf' if main_manuscript else 'plt_S.pdf'
fig.savefig(f'{filename}', dpi=300, bbox_inches='tight', format='pdf')
