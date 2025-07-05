! put your gas disk model here
! x is distance from the star (in cm)
module discstruct
   use constants, only: AU, kB, mH2, pi, third
   use parameters,only: alphat, maxrad0, vfrag, dtg, disk_path
   use utils,     only: interp2d, interp1d
   
   implicit none

   private
   public   :: init_struct, end_struct, alpha, cs, omegaK, sigmag, densg, Pg, vgas
   public   :: ddensgdz, ddensgdr, gasmass, dlogPg, diffcoefgas, deltav, Temp
   public   :: Msolid_PB, particle_formation_PB
   public   :: stokes_frag, stokes_drift, St_limited

   integer                             :: nr, nt
   real, allocatable, dimension(:, :)  :: sigma_itp, Te_itp, vel_itp, alpha_itp
   real, allocatable, dimension(:, :)  :: dsigmadr_itp, dTdr_itp
   real, allocatable, dimension(:)     :: r_itp, t_itp, Mstar_itp
   logical                             :: bounds_error=.False.

   contains

   subroutine init_struct()

      integer                    :: i, j
      character(len = 200)       :: filename, filepath
      filename = "time.info"
      filepath =  trim(disk_path)//trim(filename)
      open(unit=20, action="read", file=filepath, status="old")
      read(20,*) nt
      close(unit=20)
      allocate(t_itp(nt))

      filename = "time.dat"
      filepath =  trim(disk_path)//trim(filename)
      open(unit=21, action="read", file=filepath, status="old")
      do i=1, nt
         read(21, *) t_itp(i)
      end do
      close(unit=21)

      filename = "grid.info"
      filepath =  trim(disk_path)//trim(filename)
      open(unit=22, action="read", file=filepath, status="old")

      read(22, *) nr

      allocate(r_itp(nr))
      do j=1, nr
         read(22, *) r_itp(j)
      end do
      allocate(sigma_itp(nr, nt), Te_itp(nr, nt), vel_itp(nr, nt), Mstar_itp(nt), alpha_itp(nr, nt))
      filename = "sigma.dat"
      filepath =  trim(disk_path)//trim(filename)
      open(unit=23, action="read", file=filepath, status="old")
      read(23,*) ! skip nr line
      read(23,*) ! skip nr line
      do i =1, nt
         do j=1, nr
            read(23,*) sigma_itp(j, i)
         end do
      end do

      close(unit=22)
      close(unit=23)
      
      ! midplane temperature
      filename = "temperature.dat"
      filepath =  trim(disk_path)//trim(filename)
      open(unit=24, action="read", file=filepath, status="old")
      read(24,*) ! skip nr line
      read(24,*) ! skip nr line
      do i =1, nt
         do j=1, nr
            read(24,*) Te_itp(j, i)
         end do
      end do

      close(unit=24)
      
      ! radial velocity of gas
      filename = "velo.dat"
      filepath =  trim(disk_path)//trim(filename)
      open(unit=25, action="read", file=filepath, status="old")
      read(25,*) ! skip nr line
      read(25,*) ! skip nr line
      do i =1, nt
         do j=1, nr
            read(25,*) vel_itp(j, i)
         end do
      end do

      close(unit=25)
      
      ! mass of the central star
      filename = "mstar.dat"
      filepath =  trim(disk_path)//trim(filename)
      open(unit=26, action="read", file=filepath, status="old")
      read(26,*) ! skip nr line
      do i =1, nt
         read(26,*) Mstar_itp(i)
      end do

      close(unit=26)

      ! Temperature derivative
      filename = "dTdr.dat"
      filepath =  trim(disk_path)//trim(filename)
      open(unit=32, action="read", file=filepath, status="old")
      read(32,*) 
      allocate(dTdr_itp(nr, nt))
      do i =1, nt
         do j=1, nr
            read(32,*) dTdr_itp(j, i)
         end do
      end do
      close(unit=32)

      ! surface density derivative
      filename = "dsigmadr.dat"
      filepath =  trim(disk_path)//trim(filename)
      open(unit=33, action="read", file=filepath, status="old")
      read(33,*) 
      allocate(dsigmadr_itp(nr, nt))
      do i =1, nt
         do j=1, nr
            read(33,*) dsigmadr_itp(j, i)
         end do
      end do
      close(unit=33)
      
   end subroutine init_struct

   ! Shakura-Sunyaev's turbulence parameter
   real function alpha(x, time)
      implicit none
      real, intent(in)  :: x, time
      alpha = interp2d(r_itp, t_itp, alpha_itp, x, time, bounds_error, 0.0)
      return
   end function

   ! speed of the sound in gas
   real function cs(x, time)
      implicit none
      real, intent(in)  :: x, time
      cs     = sqrt(kB * Temp(x, time) / mH2)
      return
   end function
   
   ! temperature, isothermal vertically  
   real function Temp(x, time)
      implicit none
      real, intent(in)  :: x, time
      Temp =  interp2d(r_itp, t_itp, Te_itp, x, time, bounds_error, 0.0) 
      return
   end function

   ! temperature derivative, isothermal vertically  
   real function dTdr(x, time)
      implicit none
      real, intent(in)  :: x, time
      dTdr =  interp2d(r_itp, t_itp, dTdr_itp, x, time, bounds_error, 0.0) 
      return
   end function

   ! radial derivative of surface density
   real function dsigmadr(x, time)
      implicit none
      real, intent(in)  :: x, time
      dsigmadr =  interp2d(r_itp, t_itp, dsigmadr_itp, x, time, bounds_error, 0.0) 
      return
   end function
   
   ! the difference between keplerian and gas orbital motion
   real function deltav(x, time)
      implicit none
      real, intent(in)  :: x, time
      deltav = -0.5*dlogPg(x, 0., time)*cs(x, time)**2/x/omegaK(x, time)
      return
   end function

   ! star mass
   real function Mstar(time)
      implicit none
      real, intent(in)  :: time
      Mstar = interp1d(t_itp, Mstar_itp, time, bounds_error, 0.0) 
      return
   end function

   ! Keplerian frequency
   real function omegaK(x, time)
      use constants, only: Ggrav
      implicit none
      real, intent(in) :: x, time
      omegaK = sqrt(Ggrav*Mstar(time) / x**3)
      return
   end function


   ! gas surface density
   real function sigmag(x, time)
      use constants, only: smallv
      implicit none
      real, intent(in) :: x
      real, intent(in) :: time ! in seconds
      sigmag = interp2d(r_itp, t_itp, sigma_itp, x, time, bounds_error, 10.**(-20.))
      sigmag = max(sigmag, smallv)
      return
   end function
   
   ! Stokes number limited by fragmentation; from Birnstiel+2012
   real function stokes_frag(x, time)
      implicit none
      real, intent(in) :: x
      real, intent(in) :: time 
      stokes_frag = 0.37*third*(1./alphat)*(vfrag/cs(x, time))**2.
      return 
   end function

   ! Stokes number limited by rad drift
   real function stokes_drift(x, time)
      implicit none
      real, intent(in) :: x
      real, intent(in) :: time 
      stokes_drift = 0.5*(x*omegaK(x,time))/deltav(x,time)*dtg
      return
   end function 

   ! Stokes number limited by drift or frag
   subroutine St_limited(x, time, Stlim, power)
      implicit none
      real, intent(in) :: x
      real, intent(in) :: time 
      real, intent(out) :: Stlim, power
      real :: Stfrag, Stdrift
      Stfrag = stokes_frag(x, time)
      Stdrift = abs(stokes_drift(x, time))
      if ( Stfrag < Stdrift) then
         Stlim = Stfrag
         power = -3.5 + 4.
      else
         Stlim = Stdrift
         power = -2.5 + 4.
      endif
      return
   end subroutine 
    
   ! gas volume density
   real function densg(x, z, time)
      use constants, only: pi
      implicit none
      real, intent(in)  :: x, z
      real  :: Hg
      real, intent(in)  :: time
      Hg = cs(x, time) / omegaK(x, time) ! gas disk scale-height
      densg = (sigmag(x, time) / (sqrt(2.*pi) * Hg)) * z_exp(z, Hg)
      return
   end function

   ! this function is to make sure we do not get numerical errors when z>3H
   real function z_exp(z, H)
      real, intent(in)  :: z, H
      real, parameter   :: min_exp = 0.01 ! from exp(-3**2/2)
      z_exp = max(exp(-0.5*(z/H)**2), min_exp)
      return
   end function

   ! gas pressure
   real function Pg(x, z, time)
      implicit none
      real, intent(in)  :: x, z, time
      Pg = densg(x,z, time) * cs(x, time)**2.
      return
   end function

   ! radial gas velocity
   real function vgas(x, time)
      implicit none
      real, intent(in)  :: x, time
      vgas = interp2d(r_itp, t_itp, vel_itp, x, time, bounds_error, 0.0)
      return
   end function

   ! turbulent diffusion coefficient
   real function diffcoefgas(x, time)
      implicit none
      real, intent(in)  :: x, time
      diffcoefgas = alpha(x, time) * cs(x, time)**2 / omegaK(x, time)
      return
   end function
   
   ! partial derivative of gas density with respect to vertical height
   real function ddensgdz(x,z, time)
      implicit none
      real, intent(in)  :: x, z, time
      real              :: Hg
      Hg = cs(x, time) / omegaK(x, time)
      ddensgdz = -z * densg(x, z, time) / Hg**2
      return
   end function

   ! partial derivative of gas density with respect to r
   real function ddensgdr(x,z, time)
      implicit none
      real, intent(in)  :: x, z, time
      real  :: Hg
      Hg = cs(x, time) / omegaK(x, time) ! gas disk scale-height
      ddensgdr = densg(x,z,time)*((dsigmadr(x,time)/sigmag(x,time)) - (dTdr(x,time)/2./Temp(x,time) + 3./2./x)*(1-(z/Hg)**2.))
      
      return
   end function

   real function dlogPg(x, z, time)
      implicit none
      real, intent(in)  :: x, time,z
      real  :: Hg
      Hg = cs(x, time) / omegaK(x, time) ! gas disk scaleheight
      dlogPg = x*ddensgdr(x,z,time)/densg(x,z, time) + x*dTdr(x,time)/Temp(x,time)
      return
   end function

   ! mass of the gas included between xmin and xmax
   real function gasmass(xmin,xmax,time)
      use constants, only: pi
      implicit none
      real, intent(in)  :: xmin,xmax, time
      real              :: x
      integer           :: i
      integer, parameter:: N = 10000
      real              :: dx

      dx = (xmax - xmin) / real(N)

      gasmass = 0.0
      x = xmin + 0.5 * dx
      do i = 1, N-1
         gasmass = gasmass + sigmag(x,time) * 2.*pi*x*dx
         x = x + dx
      enddo

      return
   end function

   ! the total mass of solids for r>min_loc at t=time. We use dtg
   real function Msolid_PB(xmin, xmax, time)
      use constants, only: pi
      implicit none
      real, intent(in)  :: xmin, xmax, time
      real              :: x
      integer           :: i
      integer, parameter:: N = 10000
      real              :: dx

      dx = (xmax - xmin) / real(N)

      Msolid_PB = 0.0
      x = xmin + 0.5 * dx
      do i = 1, N-1
         Msolid_PB = Msolid_PB + dtg*sigmag(x,time) * 2.*pi*x*dx
         x = x + dx
      enddo

      return
   end function


   ! furthest distance at which a swarm can be initiated
   function rmax_calculate(r, array, mswarm, nr_new) result (rmax)

      implicit none
      real, dimension(nr_new), intent(in) :: r
      real, dimension(nr_new), intent(in) :: array
      real, intent(in)                    :: mswarm
      integer, intent(in)                 :: nr_new
      real                                :: rmax
      integer                             :: i
      
      rmax = r(1)
      do i=2, nr_new
         if (array(i)>mswarm) then
            rmax = r(i)
         endif
      end do
      return
   end function rmax_calculate

   ! init particles for initproblem_itp for PB
   subroutine particle_formation_PB(Ni, Ntot, mswarm, time, min_loc, rfor_arr, proportion)
      integer, intent(in)                    :: Ni, Ntot ! Ni indicates the index we start saving data
      real, intent(in)                       :: mswarm, time, min_loc, proportion ! proportion = 1 if only one type of particle
      real, dimension(:), intent(out)        :: rfor_arr
      real, dimension(:), allocatable        :: mass_to_fill
      real, dimension(:), allocatable        :: r_mean, r_diff
      integer  :: i, idx, kont, i_min, i_max 
      real     :: ran, b, r0, Mmax, rmin, rmax
   
      rmin = min_loc
      i_min = minloc(abs(r_itp-rmin), dim=1)
      i_max = minloc(abs(r_itp-maxrad0), dim=1)
      write(*,*) "the index is ", i_min, " and min_loc is ", rmin/AU, " while at i_min r_itp=", r_itp(i_min)/AU
      write(*,*) "the index is ", i_max, " and max_loc is ", maxrad0/AU, " while at i_min r_itp=", r_itp(i_max)/AU
      allocate(mass_to_fill(i_max-i_min), r_mean(i_max-i_min), r_diff(i_max-i_min))
      
      do i=Ni, Ntot
         rfor_arr(i) = 0.
      end do

      r_mean(:) = (r_itp(i_min+1:i_max) + r_itp(i_min:i_max-1))/2.
      r_diff(:) = (r_itp(i_min+1:i_max) - r_itp(i_min:i_max-1))

      do i=1, i_max-i_min
         mass_to_fill(i) = proportion*dtg*sigmag(r_mean(i), time)*2*pi*r_mean(i)*r_diff(i)
      enddo
      
      write(*,*) "Total mass is ", sum(mass_to_fill, dim=1)
      write(*,*) "Ni=", Ni, " and Ntot=", Ntot

      do i=Ni, Ntot
         kont = 0
         rmax = 1.01*rmax_calculate(r_mean, mass_to_fill, mswarm, i_max-i_min)
         Mmax = maxval(mass_to_fill, dim=1)
         do while (.True.)
      
            if (kont>10000) then ! in case there is an infinite loop
               !print*, "Issue..."
               !print*, "We are in particle number", i
               !print*, "Total particle number is ", Ntot
               !print*, "rmin is ", rmin/AU, " and rmax is ", rmax/AU

               idx = maxloc(mass_to_fill, dim=1)
               rfor_arr(i) = r_mean(idx)
               mass_to_fill(idx) = mass_to_fill(idx) - mswarm
               exit
            endif
         
            call random_number(ran)
            r0 = (rmax - rmin)*ran + rmin
            call random_number(ran)
            b = ran*Mmax
            idx = minloc(abs(r_mean-r0), dim=1)
            ! acceptance or rejection
            if (b < mass_to_fill(idx)) then
               mass_to_fill(idx) = mass_to_fill(idx) - mswarm
               rfor_arr(i) = r0
               exit
            end if
            
            kont = kont + 1 
         
         end do

      end do
      deallocate(mass_to_fill, r_mean, r_diff)
   end subroutine particle_formation_PB



   ! for deallocating
   subroutine end_struct()
      deallocate(sigma_itp, Te_itp, vel_itp, alpha_itp, r_itp, t_itp, &
                & Mstar_itp, dsigmadr_itp, dTdr_itp)
   end subroutine end_struct



   end module discstruct
