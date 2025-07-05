! module gathers physical constants values used in the simulations
module constants
   implicit none
   public
   real, parameter        :: pi      = 3.14159265358979323846
   real, parameter        :: kB      = 1.38065e-16               ! Boltzmann constant [erg/K]
   real, parameter        :: third   = 1./3.                     ! http://xkcd.com/1275/
   real, parameter        :: AU      = 1.49597870691e+13         ! AU in cm 
   real, parameter        :: year    = 365.256363051*24.0*3600.0 ! year in sek
   real, parameter        :: NAv     = 6.02214179e+23            ! Avogadro constant
   real, parameter        :: mH2     = 2.3 / NAv                 ! mean particle mass of gas
   real, parameter        :: Ggrav   = 6.67384e-8                ! gravitational constant
   real, parameter        :: Msun    = 1.9891e+33                ! mass of the Sun (used as mass of central star)
   real, parameter        :: AH2     = 2.e-15                    ! geometrical cross-section of H2 molecule
   real, parameter        :: smallv  = 1.e-10                    ! ridiculously small value
end
