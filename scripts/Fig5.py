# reproduce figure 5
# before running the script, copy "local_data_reduced" directory to "../outputs/."
import h5py
import os
import glob
import numpy as np
from Plotting import init_plot
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.constants import G, k_B, N_A
from scipy import interpolate



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
t_CAI = 0.19 # value in Myr
t_eval = np.array([3.98, 4.06])
t_eval += t_CAI
nrl = 30
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



# --- indicate minimum sigmad for planetesimal formation
# where to find disk data
testfile = '../disk/'

Zcrit = 0.05



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


disk = Disk(testfile)


# make plot


fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(25, 12), sharex=True)

# color to indicate fragile and rigid
colors = {'fragile':'#3274B5', 'rigid':'#8B551B', 'planetesimal':'blueviolet'}

lw = 8
ymin, ymax = 10**-5, 20
#Label panel 
labels = ['A', 'B', 'C', 'D']
r_arr = np.arange(r_evall.min(), r_evall.max(), 0.005)

for idx in [0, 1]:
    title = 'Local simulation' if idx==0 else None
    xlabel = 'Radial distance (AU)' if idx==1 else None
    init_plot(ax[idx, 0], title=title,
              xlabel=xlabel, font=30)

    init_plot(ax[idx, 1],
              xlabel=xlabel, font=30)

    ax[idx,0].plot(r_evall,
               (sigmadl[:, idx]-sigmadl_fri[:, idx]),
               c=colors['fragile'],lw=lw)
    ax[idx,0].plot(r_evall,
               sigmadl_fri[:, idx], 
               c=colors['rigid'], lw=lw)
    ax[idx,0].plot(r_evall, Zcrit*disk.f_sigma(t_eval[idx], r_evall),
                   alpha=0.5, c=colors['planetesimal'], 
               lw=lw/2, zorder=-1)

    rmax = r_evall[1:-1][sigmadl[1:-1, idx].argmax()]
    ax[idx,0].axvline(ymin=ymin, ymax=ymax, x=rmax, c=colors['planetesimal'], lw=lw/2, ls='-', alpha=0.5, zorder=-1)

    for a, rho, mat in zip([10**-4, 10**-3], [1.2, 3.3], ['fragile', 'rigid']):
        ax[idx, 1].plot(r_arr,
                        disk.f_vs(t_eval[idx]*u.Myr, r_arr*u.AU, 
                                  a*u.cm, rho*u.g/u.cm**3),
                        c=colors[mat], lw=lw, label=f'{a*10**4:.0f} $\mu$m')
    
    
        signs = np.sign(disk.f_vs(t_eval[idx]*u.Myr, r_arr*u.AU, a*u.cm, rho*u.g/u.cm**3))
        zero_crossings = np.where(np.diff(signs) != 0)[0]
        if len(zero_crossings)>0:
            idx2 = zero_crossings[-1]
            rdt2 = r_arr[idx2]
            ax[idx,1].scatter(rdt2, 0., s=600, edgecolor='k', c=colors[mat], zorder=3)


    # indicate time
    ax[idx,1].text(0.75, 0.92, f'{t_eval[idx]-t_CAI:.2f} Myr', transform=ax[idx,1].transAxes,
                   fontsize=30, va='top', ha='left')
    
    

    ax[idx,0].text(0.02, 0.95, labels[idx*(idx+1)], transform=ax[idx,0].transAxes,
            fontsize=34, fontweight='bold', va='top', ha='left')
    ax[idx,1].text(0.02, 0.95, labels[idx*(idx+1)+1], transform=ax[idx,1].transAxes,
            fontsize=34, fontweight='bold', va='top', ha='left')
    ax[idx,1].hlines(y=0, xmin=r_evall.min(), xmax=r_evall.max(), color='lightgrey', zorder=0)

    




ax[0,0].plot([],[], lw=lw, label='Fragile', c=colors['fragile'])
ax[0,0].plot([],[], lw=lw, label='Rigid', c=colors['rigid'])
ax[0,0].legend(ncol=2, loc='upper right')
ax[0,1].legend(ncol=2, loc='lower center')


ax[0,0].set_xlim(5.830, 6.419) # same boundary as fig. 4

ax[0,1].set_ylim(-11, 10)
ax[1,1].set_ylim(-55, 70)

ax[0,0].set_yscale('log')
ax[1,0].set_yscale('log')

ax[0,0].set_ylim(ymin, ymax)
ax[1,0].set_ylim(ymin, ymax)

plt.subplots_adjust(wspace=0.2)

fig.supylabel('Dust surface density (g cm$^{-2}$)', y=0.5, x=0.05, fontsize=30)
fig.text(s='Radial velocity (m s$^{-1}$)', y=0.5, x=0.49, fontsize=30,
         rotation=90, va='center')




plt.show()

fig.savefig('CRCI.pdf', dpi=300, bbox_inches='tight', format='pdf')
