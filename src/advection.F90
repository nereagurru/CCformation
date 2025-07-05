! this module performes advection of every representative particle due to radial drift, vertical settling and turbulent
! diffusion
! in this module all the quantities should be calculated locally (at the particle location) therefore we don't need the
! grid to perform advection

module advection

   use constants,  only: pi, third, AH2, mH2
   use discstruct, only: Pg, dlogPg, densg, omegaK, cs, vgas, ddensgdr, ddensgdz
   use types,      only: swarm
   use parameters, only: smallr, alphat

   implicit none

   private
   public      :: mc_advection, mc_advection_element, stokesnr, vel_rd, vel_vn, update_St

   contains

   ! routine performes radial and vertical advection due to radial drift, vertical settling and turbulent diffusion
   ! for every particle in the simulation
   subroutine mc_advection(swrm, dtime, realtime &
#ifdef GLOBAL
                           &, Nflux &
#endif
                           & )
      implicit none
      type(swarm), dimension(:), allocatable, target :: swrm
      type(swarm), pointer                           :: particle
      real, intent(in)                               :: dtime
      real, intent(in)                               :: realtime
      real                                           :: vs, vr, velr, velv, vn
      integer                                        :: i
#ifdef GLOBAL
      integer :: Nflux
#endif

      ! loop over all the particles: calculating advection velocities
#ifdef GLOBAL
      !$OMP PARALLEL DO SHARED(Nflux) PRIVATE(particle,vs, vr, velr, velv, vn) SCHEDULE(DYNAMIC)
#else
      !$OMP PARALLEL DO PRIVATE(particle,vs, vr, velr, velv, vn) SCHEDULE(DYNAMIC)
#endif
      do i = 1, size(swrm)
         if (swrm(i)%rdis < smallr) cycle
         particle => swrm(i)

         particle%stnr = stokesnr(particle,realtime)

         ! vertical transport
         call vel_vs(particle, vs, realtime)
         call vel_ver(particle, vs, velv, dtime, realtime)

         ! radial transport
         call vel_vn(particle, vn, realtime)
         call vel_rd(particle, vr, vn, realtime)
         call vel_rad(particle, velr, vr, dtime, realtime)
         ! transport particle radially
         particle%rdis = particle%rdis + vr * dtime 

         ! is it still outside evaporation radius?
         if (particle%rdis > smallr) then
            ! transport particle radially
            particle%zdis = particle%zdis + velv * dtime
            particle%stnr = stokesnr(particle,realtime)
         else
#ifdef GLOBAL
            !$OMP critical
            Nflux = Nflux + 1
            !$OMP end critical
            
            ! when did particle cross 6.5AU?
            particle%tfor = realtime + dtime + (smallr-particle%rdis)/vr
            ! z and St of particle when crossing smallr?
            particle%zdis = particle%zdis + velv * (dtime + (smallr-particle%rdis)/vr)
            particle%rdis = smallr
            particle%stnr = stokesnr(particle,realtime)

#endif
            particle%rdis = 0.99 * smallr
            particle%velr = 1.e-20
            particle%velz = 1.e-20
         endif
      enddo
      !$OMP END PARALLEL DO

      return
   end subroutine mc_advection


   ! this function is employed when adding particles in the local simulation
   subroutine mc_advection_element(particle, dtime, realtime)
      implicit none
      type(swarm), target                            :: particle
      real, intent(in)                               :: dtime
      real, intent(in)                               :: realtime
      real                                           :: vs, vr, velr, velv, vn

      if (particle%rdis > smallr) then
         
         call vel_vn(particle, vn, realtime)
         particle%stnr = stokesnr(particle,realtime)

         call vel_vs(particle, vs, realtime)
         call vel_ver(particle, vs, velv, dtime, realtime)
         particle%zdis = particle%zdis + velv * dtime

         call vel_rd(particle, vr, vn, realtime)
         call vel_rad(particle, velr, vr, dtime, realtime)
         particle%rdis = particle%rdis + vr * dtime 

         if (particle%rdis > smallr) then
            particle%stnr = stokesnr(particle,realtime)
         else
            particle%rdis = 0.99 * smallr
            particle%velr = 1.e-20
            particle%velz = 1.e-20
         endif
      endif

      return
   end subroutine mc_advection_element

   ! Keplerian velocity reduction of gas
   subroutine vel_vn(particle, vn, realtime)
      implicit none
      type(swarm)                                    :: particle
      real, intent(in)                               :: realtime
      real                                           :: vn
      vn = 0.5*dlogPg(particle%rdis, particle%zdis, realtime)*Pg(particle%rdis, particle%zdis, realtime)/particle%rdis/ &
           densg(particle%rdis, particle%zdis, realtime)/omegaK(particle%rdis, realtime)
      return
   end subroutine vel_vn

   ! the complete vertical velocity
   subroutine vel_ver(particle, vs, velv, dtime, realtime)
      implicit none
      type(swarm)                                     :: particle
      real, intent(in)                                :: dtime, realtime
      real, intent(in)                                :: vs
      real                                            :: velv
      real                                            :: dz
      real                                            :: Ldiff
      real, dimension(2)                              :: rand
      real                                            :: vD1, vD2, visc

      visc = alphat*(cs(particle%rdis, realtime)**2)/omegaK(particle%rdis, realtime)
      Ldiff = sqrt(2. * dtime * visc  / (1. + particle%stnr**2))
      call random_number(rand)
      dz = Ldiff/(sqrt(2.*log(2.))) * sqrt(-2.*log(rand(1))) * cos(2.*pi*rand(2))
      vD1 = dz / dtime
      vD2 = (visc / (1. + particle%stnr**2) / densg(particle%rdis,particle%zdis, realtime)) * &
            ddensgdz(particle%rdis,particle%zdis, realtime)
      velv = - vs + vD1 + vD2
      particle%velz = vs

      return
   end subroutine vel_ver

   ! the complete radial velocity
   subroutine vel_rad(particle, velr, vr, dtime, realtime)
      implicit none
      type(swarm)                                     :: particle
      real, intent(in)                                :: dtime, realtime
      real, intent(in)                                :: vr
      real                                            :: velr
      real                                            :: Ldiff
      real, dimension(2)                              :: rand
      real                                            :: dr
      real                                            :: vD1, vD2, visc


      visc = alphat*(cs(particle%rdis, realtime)**2)/omegaK(particle%rdis, realtime)
      Ldiff = sqrt(2. * dtime * visc / (1. + particle%stnr**2))
      call random_number(rand)
      dr = Ldiff/(sqrt(2.*log(2.))) * sqrt(-2.*log(rand(1))) * cos(2.*pi*rand(2))
      vD1 = dr / dtime
      vD2 = (visc / (1. + particle%stnr**2) / densg(particle%rdis,particle%zdis, realtime)) * &
            (ddensgdr(particle%rdis,particle%zdis, realtime))
      velr = vD1 + vD2

      particle%velr = vr
      return
   end subroutine vel_rad


   ! calculates the Stokes numbers of one particle locally
   real function stokesnr(particle, realtime)
      implicit none
      type(swarm)                                     :: particle
      real                                            :: lmfp
      real                                            :: rad
      real                                            :: gasdens, css
      real, intent(in)                                :: realtime

      gasdens = densg(particle%rdis,particle%zdis,realtime)
      css = cs(particle%rdis, realtime)
      lmfp = mH2 / ( gasdens * AH2 )

      rad = (0.75 / pi / particle%dens)**third * particle%mass**third
      
      if (rad > 2.25 * lmfp) then ! Stokes regime
         stokesnr = sqrt(2.*pi) * rad**2. *  particle%dens * omegaK(particle%rdis, realtime) * AH2 / (9. * css * mH2)
      else                       ! Epstein regime
         stokesnr = rad *  particle%dens / (sqrt(8./pi) * css * gasdens) * omegaK(particle%rdis, realtime)
      endif

      return
   end function stokesnr

   ! vertical settling velocity
   subroutine vel_vs(particle, vs, realtime)
      implicit none
      type(swarm)                                     :: particle
      real, intent(in)                                :: realtime
      real                                            :: vs
      vs = particle%zdis * omegaK(particle%rdis, realtime) * particle%stnr / (1. + particle%stnr**2.)
      return
   end subroutine vel_vs

   ! velocity due to radial drift and advection
   subroutine vel_rd(particle, vr, vn, realtime)
      implicit none
      type(swarm)                                     :: particle
      real, intent(in)                                :: realtime
      real                                            :: vr
      real, intent(in)                                :: vn
      vr = 2. * vn / (particle%stnr + 1./particle%stnr) + vgas(particle%rdis, realtime) / (1. + particle%stnr* particle%stnr)
      return
   end subroutine vel_rd


   ! update the stokes number of particles
   subroutine update_St(swrm, realtime)
      implicit none
      type(swarm), dimension(:), allocatable, target :: swrm
      type(swarm), pointer                           :: particle
      real, intent(in)                               :: realtime
      integer                                        :: i

      !$OMP PARALLEL DO PRIVATE(particle) SCHEDULE(DYNAMIC)
      do i = 1, size(swrm)
         if (swrm(i)%rdis < smallr) cycle
         particle => swrm(i)
         particle%stnr = stokesnr(particle,realtime)
      enddo
      !$OMP END PARALLEL DO
      return
   end subroutine
end
