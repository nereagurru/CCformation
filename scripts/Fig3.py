# to reproduce figure 3 


import h5py
import os
import glob
import numpy as np
from Plotting import init_plot
import matplotlib.pyplot as plt
from astropy import units as u
from scipy import interpolate
from astropy.constants import k_B, N_A




# read disk data and create function to calculate sound speed and gas surface density


     
# location of data where diskevol results are saved
testfile = '../disk/'
# read time grid
t = (np.loadtxt(testfile + 'time.dat',dtype=float)*u.s).to(u.Myr)
# read position grid and remove first element
r = np.loadtxt(testfile + 'grid.info',dtype=float)
n = int(r[0])

r = (r[1::]*u.cm).to(u.AU)

# read temperature across the disk and over time and create its interpolation function
Td_data = np.loadtxt(testfile+ 'temperature.dat', dtype=float, 
                     skiprows=2).reshape((n, t.shape[0]), order='F')*u.K

Td = interpolate.interp2d(t, r, Td_data)


# read gas surface density over time and create its interpolation function
sigma_dat = np.loadtxt(testfile+'sigma.dat', dtype=float,
                   skiprows=2).reshape((n, t.shape[0]), order='F')*u.g/u.cm**2

sigma = (interpolate.interp2d(t, r, sigma_dat))
        

    
# function for estimating the gas surface density. It can include a pressure bump
def f_sigma(t, r):
    if not isinstance(t, u.Quantity) : TypeError('time variable must be a float')
    # calculate surface density without substructure
    sigma_ = np.squeeze(sigma(t, r))*u.g/u.cm**2 if isinstance(r, u.Quantity) \
            else np.squeeze([sigma(t, ri) for ri in r])*u.g/u.cm**2
    return sigma_
        

# midplane temperature 
def f_Td(t, r):
    if not isinstance(t, u.Quantity) : TypeError('time variable must be a float')
    Td_ = np.squeeze(Td(t, r))*u.K if isinstance(r, u.Quantity) \
         else np.squeeze([Td(t, ri) for ri in r])*u.K
    return Td_
    
    
# function for sound speed (t in *u.Myr and r in *u.AU) 
def f_cs(t, r):
    return (np.sqrt(k_B*N_A*f_Td(t, r)/(2.3*u.g/u.mol))).to(u.m/u.s)

# particle size assuming Epstein regime
def a_Epstein(St, rhoi, sigma):
    return 2/np.pi*St*sigma/rhoi
# internal density of particle from two material
def f_densi(CH):
    return 1*u.g/u.cm**3/((CH/3.3) + (1-CH)/1.2)


# read dust data

output_dir = "../outputs/zerod_data_reduced"

all_files = sorted(glob.glob(os.path.join(output_dir, "*.h5")))


t_arr = np.empty((len(all_files),))



# read pebble data
# first save time values. Then we will decide the closest ones to our evaluation times
for i, file_path in enumerate(all_files):

    with h5py.File(file_path, 'r') as f:

        # Read compound dataset from swarms/swarmsout
        t_arr[i] = f['snapt'][...]


time_eval = np.array([2.3, 3.5, 3.98, 4.06])
t_CAI = 0.19

from_AU_to_cm = (1*u.AU).to(u.cm).value


# every how often compute mean value
Nstep = 25



lw = 5 # linewidth
na = 20 # size bins


# array of sizes logarithmically distributed from 0.1 um to 3 cm
a_range = np.logspace(np.log10(0.00001), np.log10(3), num=na, base=10.)


# limits of the y axis
ymin, ymax = 10**-3, 10



# color to indicate fragile and rigid
colors = ['#3274B5', '#8B551B']



fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(20,15),
                       sharex=True, sharey=True)

# set boundaries of the simulation and calculate disk Area
rmin, rmax = 6.09,6.11
rmean, rdif = 0.5*(rmax+rmin), (rmax-rmin)
Area = 2*np.pi*rmean*rdif





for idx_time in range(0,4):
    
    timedif = 1_000 if idx_time > 1 else 25_000
    # time of the simulation with respect to CAI formation
    time = time_eval[idx_time]

    itt = np.absolute((time+t_CAI)*10**6-t_arr).argmin()

    i_init = np.where(t_arr[itt]-timedif < t_arr)[0].min()


    i_end = np.where(t_arr[itt]+timedif > t_arr)[0].max()
    Nstep = i_end - i_init + 1

    r_arr = np.empty((Nstep, 200))
    rigid = np.empty((Nstep, 200))
    grain_size = np.empty((Nstep, 200))
    itt_0 = itt-i_init # from 0
    for i, file_path in enumerate(all_files[i_init:i_end+1]):

        with h5py.File(file_path, 'r') as f:
            r_arr[i,:] = f['cylindrical radius [AU]'][...]
            rigid[i,:] = f['rigid'][...]
            grain_size[i,:] = f['grain size'][...]
            
            if i==itt_0: mswarm = f.attrs['mass_of_swarm']
    
    # set index for ax
    i = 0 if idx_time < 2 else 1
    j = 0 if idx_time==0 or idx_time==2 else 1
    
    # set axes
    xlabel = 'Dust radius (cm)' if i==1 else None
    ylabel= r'$\sigma_{\rm{d}}$(a) (g cm$^{-2}$)' if j==0 else None

    #create subplot
    init_plot(ax[i,j], f'{time:.02f} Myr', xlabel, ylabel, font=30) 
   
    # Stokes number limited by fragmentation
    alphat=10**-4
    vf = 2 # in m/s
    Stfrag = 0.37/3/alphat*(vf/f_cs((t_arr[itt]*u.yr).to(u.Myr), 6.1*u.AU))**2

    # particle size assuming Epstein regime
    CH = 0 if idx_time == 3 else 0.5
    afrag = a_Epstein(Stfrag, f_densi(CH=CH), f_sigma((t_arr[itt]*u.yr).to(u.Myr), 6.1*u.AU))        

    

    
    # initiate empty arrays for dust surface density
    if Nstep==1:
        
        sigma_a = np.zeros((Nstep, a_range.shape[0]+1))
        sigma_a_fri = np.zeros((Nstep, a_range.shape[0]+1))
    
    
    

        # how many particles per size bin 
        counts, _ = np.histogram(grain_size[itt_0,:], bins=a_range)
        # calculate logarithmic surface density per size bin (see Birnstiel 2023)
        sigma_a[0,1:-1] = (counts*mswarm/Area)/from_AU_to_cm**2 \
                            /(a_range[1:]-a_range[:-1])

    
        # how many rigid particles per size bin
        counts_fri, _ = np.histogram(grain_size[itt_0,:], bins=a_range,
                                     weights=rigid[itt_0,:])
        
        # calculate logarithmic surface density of rigids per size bin (see Birnstiel 2023)
        sigma_a_fri[0,1:-1] = (counts_fri*mswarm/Area)/from_AU_to_cm**2 \
                                /(a_range[1:]-a_range[:-1])

        
    else:
        sigma_a = np.zeros((2*Nstep, a_range.shape[0]+1))
        sigma_a_fri = np.zeros((2*Nstep, a_range.shape[0]+1))


        
    
    
        for it in range(0, Nstep):
            
    
            # how many particles per size bin 
            counts, _ = np.histogram(grain_size[it,:], bins=a_range)
            # calculate logarithmic surface density per size bin (see Birnstiel 2023)
            sigma_a[it,1:-1] = (counts*mswarm/Area)/from_AU_to_cm**2 \
                                /(a_range[1:]-a_range[:-1])
    


    

            # how many rigid particles per size bin
            counts_fri, _ = np.histogram(grain_size[it,:], bins=a_range,
                                         weights=rigid[it,:])
            
            # calculate logarithmic surface density of rigids per size bin (see Birnstiel 2023)
            sigma_a_fri[it,1:-1] = (counts_fri*mswarm/Area)/from_AU_to_cm**2 \
                                    /(a_range[1:]-a_range[:-1])
        

    # plot mean value of size bin
    a_plot = (a_range[1:]+a_range[:-1])/2
    a_plot = np.insert(a_plot, 0, a_range[0]) # min size
    a_plot = np.append(a_plot, a_range[-1]) # max size
    
    # plot line of fragile material 
    ax[i,j].plot(a_plot, a_plot*np.mean((sigma_a[:,:]-sigma_a_fri[:,:]), axis=0),
                 lw=lw, c=colors[0])
    
    # plot line of rigid material 
    ax[i,j].plot(a_plot, a_plot*np.mean(sigma_a_fri[:,:], axis=0), lw=lw, c=colors[1])
    


    ax[i,j].vlines(x=afrag.value, ymin=ymin, ymax=ymax, color='k',
                   lw=lw/2, ls='--', zorder=0)
    

        


        

# set legend
ax[0,0].plot([],[], lw=lw, c=colors[0], label='Fragile')
ax[0,0].plot([],[], lw=lw, c=colors[1], label='Rigid')
ax[0,0].legend(ncol=2, loc='upper left')

# set logarithmic scale
ax[0,0].set_yscale('log')
ax[0,0].set_xscale('log')

# set limits
ax[0,0].set_ylim(ymin, ymax)
ax[0,0].set_xlim(0.00005, 5)

# set ticks
ax[0,0].set_xticks([0.0001, 0.001, 0.01, 0.1, 1])
ax[0,0].set_yticks([0.001, 0.01, 0.1, 1, 10])

plt.subplots_adjust(wspace=0.05, hspace=0.15)
plt.show()
