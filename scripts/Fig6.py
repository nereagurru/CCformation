# reproduce figure 6
# before running the script, copy "local_data_reduced/" directory to ../outputs/.
# and the disk data files (*.dat, *.info) to ../disk/.

import h5py
import os
import glob
import numpy as np
from Plotting import init_plot
import matplotlib.pyplot as plt
from scipy import interpolate
from astropy.constants import G, k_B, N_A
from astropy import units as u
import matplotlib.gridspec as gridspec


# read local data and calculate dust surface density

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
t_eval = np.array([3.98, 4.06])
t_eval += t_CAI

nrl = 16
rl = np.linspace(5.725, 6.5, nrl)

from_AU_to_cm = (1*u.AU).to(u.cm).value
from_Myr_to_s = (1*u.Myr).to(u.s).value

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




# read disk data and create function to calculate velocities


class Disk:
    def __init__(self, testfile='../disk/'):
        
        # location of data where diskevol results are saved
        self.testfile = testfile
        # read time grid
        self.t = (np.loadtxt(self.testfile + 'time.dat',dtype=float)*u.s).to(u.Myr)
        # read position grid and remove first element
        self.r = np.loadtxt(self.testfile + 'grid.info',dtype=float)
        self.n = int(self.r[0])

        self.r = (self.r[1::]*u.cm).to(u.AU)
        # read star mass over time and create its interpolation function
        Mstar_dat = (np.loadtxt(self.testfile + 'mstar.dat',dtype=float)*u.g).to(u.Msun)
        self.Mstar = interpolate.interp1d(self.t, Mstar_dat)

        # read temperature across the disk and over time and create its interpolation function
        Td_data = np.loadtxt(self.testfile+ 'temperature.dat', dtype=float, 
                             skiprows=2).reshape((self.n, self.t.shape[0]), order='F')*u.K
        
        self.Td = interpolate.interp2d(self.t, self.r, Td_data)
        
        # read derivative temperature across the disk and over time and create its interpolation function
        dTdr_data = np.loadtxt(self.testfile+ 'dTdr.dat', dtype=float, 
                             skiprows=1).reshape((self.n, self.t.shape[0]), order='F')*(u.K/u.cm).to(u.K/u.AU)
        
        self.dTdr = interpolate.interp2d(self.t, self.r, dTdr_data)
        
        # read derivative temperature across the disk and over time and create its interpolation function
        dsigmadr_data = np.loadtxt(self.testfile+ 'dsigmadr.dat', dtype=float, 
                             skiprows=1).reshape((self.n, self.t.shape[0]), order='F')*(u.g/u.cm**3).to(u.g/u.cm**2/u.AU)
        
        self.dsigmadr = interpolate.interp2d(self.t, self.r, dsigmadr_data)
        

        # read gas velocity across the disk and over time and create its interpolation function
        vg = (np.loadtxt(self.testfile+ 'velo.dat', dtype=float, 
                             skiprows=2).reshape((self.n, self.t.shape[0]), order='F')*u.cm/u.s).to(u.m/u.s)
        self.vg_itp = interpolate.interp2d(self.t, self.r, vg)



        # read gas surface density over time and create its interpolation function
        self.sigma_dat = np.loadtxt(self.testfile+'sigma.dat', dtype=float,
                           skiprows=2).reshape((self.n, self.t.shape[0]), order='F')*u.g/u.cm**2

        self.sigma = (interpolate.interp2d(self.t, self.r, self.sigma_dat))
        

    
    # function for estimating the gas surface density. It can include a pressure bump
    def f_sigma(self, t, r):
        if not isinstance(t, u.Quantity) : TypeError('time variable must be a float')
        # calculate surface density without substructure
        sigma = np.squeeze(self.sigma(t, r))*u.g/u.cm**2 if isinstance(r, u.Quantity) \
                else np.squeeze([self.sigma(t, ri) for ri in r])*u.g/u.cm**2
        return sigma
        

    def f_Mstar(self, t):
        if not isinstance(t, u.Quantity) : TypeError('time variable must be a float')
        return self.Mstar(t)*u.Msun
    
    def f_Td(self, t, r):
        if not isinstance(t, u.Quantity) : TypeError('time variable must be a float')
        Td = np.squeeze(self.Td(t, r))*u.K if isinstance(r, u.Quantity) \
             else np.squeeze([self.Td(t, ri) for ri in r])*u.K
        return Td
    
    def f_dTdr(self, t, r):
        if not isinstance(t, u.Quantity) : TypeError('time variable must be a float')
        dT = np.squeeze(self.dTdr(t, r))*u.K/u.AU if isinstance(r, u.Quantity) \
             else np.squeeze([self.dTdr(t, ri) for ri in r])*u.K/u.AU
        return dT
    
    def f_dsigmadr(self, t, r):
        if not isinstance(t, u.Quantity) : TypeError('time variable must be a float')
        dsigma = np.squeeze(self.dsigmadr(t, r))*u.g/u.cm**2/u.AU if isinstance(r, u.Quantity) \
             else np.squeeze([self.dsigmadr(t, ri) for ri in r])*u.g/u.cm**2/u.AU 
        return dsigma
    
    # function for sound speed (t in *u.Myr and r in *u.AU) 
    def f_cs(self, t, r):
        return (np.sqrt(k_B*N_A*self.f_Td(t, r)/(2.3*u.g/u.mol))).to(u.m/u.s)
    
    # function for keplerian frequency (t in *u.Myr and r in *u.AU) 
    def f_omegak(self, t, r):
        return (np.sqrt(G*self.f_Mstar(t)/r**3)).to(1/u.Myr)

    # function for vertical scale height (t in *u.Myr and r in *u.AU) 
    def f_H(self, t, r):
        return (self.f_cs(t, r)/self.f_omegak(t, r)).to(u.AU)


    # function for inverse of logarithmic pressure gradient (t in *u.Myr and r in *u.AU) 
    def dlnPdlnr(self, t, r):
        #deltar = (5*10**11*u.cm).to(u.AU)
        return r/self.f_sigma(t, r)*self.f_dsigmadr(t, r) + (r/2/self.f_Td(t, r)*self.f_dTdr(t, r) - 3/2) #+ r/self.f_Td(t, r)*self.f_dTdr(t, r)

    # function for Keplerian reduction (t in *u.Myr and r in *u.AU) 
    def delta_v(self, t, r):
        return -(0.5*(self.f_H(t, r)/r)*self.dlnPdlnr(t, r)*self.f_cs(t, r)).to(u.m/u.s)
    
    # radial gas velocity
    def f_vg(self, t, r):
        if not isinstance(t, u.Quantity) : TypeError('time variable must be a float')
        vg = np.squeeze(self.vg_itp(t, r))*u.m/u.s if isinstance(r, u.Quantity) \
             else np.squeeze([self.vg_itp(t, ri) for ri in r])*u.m/u.s

        return vg
    
    # radial solid velocity
    def f_vs(self, t, r, a, rho):
        St =  np.pi/2*a*rho/self.f_sigma(t, r) # assume Epstein regime
        return (self.f_vg(t, r) - 2*self.delta_v(t, r)*St)/(1+St**2)


disk = Disk('../disk/')


# make plot

fig = plt.figure(figsize=(20, 10))
gs = gridspec.GridSpec(2, 2, width_ratios=[2, 1])

ax1 = fig.add_subplot(gs[:, 0])

ax2 = fig.add_subplot(gs[0, 1], sharex=ax1)
ax3 = fig.add_subplot(gs[1, 1], sharex=ax1)

# fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(20, 10), sharey=True)
init_plot(ax1, title='Local model', ylabel='Dust surface density (g cm$^{-2}$)', font=30)
init_plot(ax2, font=25)
init_plot(ax3, font=25)

ls_arr = ['--', '-']
lw = 5
ymin, ymax = 10**-5, 5

for idx, ls in enumerate(ls_arr):

    ax1.plot(r_evall, (sigmadl[:, idx]-sigmadl_fri[:, idx]),
             '#3274B5', ls=ls, lw=lw, zorder=1)
    

    ax1.plot(r_evall, sigmadl_fri[:, idx], 
             ls=ls, c='#8B551B', lw=lw, zorder=2)
    

r_arr = np.arange(r_evall.min(), r_evall.max(), 0.005)
a, rho = 0.001, 3.3 
ax2.plot(r_arr,  disk.f_vs(t_eval[0]*u.Myr, r_arr*u.AU, a*u.cm, rho*u.g/u.cm**3), lw=lw,
         label=f'{a*1000:.0f} $\mu$m', c='#8B551B', ls='--')


# find dust trap for CR
signs = np.sign(  disk.f_vs(t_eval[0]*u.Myr, r_arr*u.AU, a*u.cm, rho*u.g/u.cm**3))
zero_crossings = np.where(np.diff(signs) != 0)[0]
idx2 = zero_crossings[-1]
rdt2 = r_arr[idx2]
ax1.axvline(x=rdt2, ymin=ymin, ymax=ymax, ls='--',  c='#8B551B', zorder=-1)
ax2.axvline(x=rdt2, ymin=-10, ymax=10, ls='--',  c='#8B551B')
ax3.plot(r_arr,  disk.f_vs(t_eval[1]*u.Myr, r_arr*u.AU, a*u.cm, rho*u.g/u.cm**3), lw=lw, c='#8B551B')

a, rho = 10**-4, 1.2
ax2.plot(r_arr,  disk.f_vs(t_eval[0]*u.Myr, r_arr*u.AU, a*u.cm, rho*u.g/u.cm**3), lw=lw,
         label=f'{a*1000:.0f} $\mu$m', c='#3274B5', ls='--')



ax3.plot(r_arr,  disk.f_vs(t_eval[1]*u.Myr, r_arr*u.AU, a*u.cm, rho*u.g/u.cm**3), lw=lw, c='#3274B5')
# find dust trap for CI
signs = np.sign( disk.f_vs(t_eval[1]*u.Myr, r_arr*u.AU, a*u.cm, rho*u.g/u.cm**3))
zero_crossings = np.where(np.diff(signs) != 0)[0]
idx2 = zero_crossings[-1]
rdt2 = r_arr[idx2]
ax1.axvline(x=rdt2, ymin=ymin, ymax=ymax, ls='-',  c='#3274B5', zorder=-1)
ax3.axvline(x=rdt2, ymin=-100, ymax=100, ls='-',  c='#3274B5')

ax2.hlines(y=0, xmin=r_evall.min(), xmax=r_evall.max(), color='lightgrey', zorder=0)
ax3.hlines(y=0, xmin=r_evall.min(), xmax=r_evall.max(), color='lightgrey', zorder=0)
for i in range(0, t_eval.shape[0]):
    ax1.plot([],[], lw=lw, ls=ls_arr[i],label=f'{t_eval[i]-0.19:.2f} Myr', c='k')




ax1.plot([],[], lw=lw, label='Fragile', c='#3274B5')
ax1.plot([],[], lw=lw, label='Rigid', c='#8B551B')

ax1.set_xlim(5.8, 6.5)
ax1.legend(ncol=2)
ax2.set_ylim(-11, 10)
ax3.set_ylim(-55, 70)
ax2.legend()


ax1.set_yscale('log')
ax1.set_ylim(ymin, ymax)

#ax[0].set_yticks([0.00001, 0.0001, 0.001, 0.01, 0.1, 1])
ax1.set_xlabel('Radial distance (AU)', y=0., fontsize=30)
ax3.set_xlabel('Radial distance (AU)', y=0., fontsize=30)
plt.subplots_adjust(wspace=0.3, hspace=0.2)
fig.text(0.6, 0.5, 'Radial velocity (m s$^{-1}$)', va='center', rotation='vertical', fontsize=30)
plt.show()
