! this module performs collisions between representative bodies
! in this module all the quantities should be computed at the center of current cell
! in the arrays: first index -> representative particle, 2nd -> physical


module collisions

   use constants,  only: pi, third, mH2, AH2
   use discstruct, only: cs, omegaK, densg, Pg, dlogPg, vgas
   use grid,       only: g
   use initproblem,only: m0, mswarm
   use parameters, only: matdens, ridens, dmmax, vfrag, alphat

   use types,      only: swarm, list_of_swarms

   implicit none

   private
   public :: mc_collisions

   real, parameter:: eps_fri=1.e-9

   contains

   ! the routine performes collisional evolution on swarms located in the cell nr,ni with the aid of the MC algorithm
   subroutine mc_collisions(nr, ni, bin, swrm, dtime, realtime, kk, index_init)
      implicit none
      type(swarm), dimension(:), allocatable, target            :: swrm
      type(list_of_swarms), dimension(:,:), allocatable, target :: bin
      type(swarm), dimension(:), pointer              :: swarms      ! local rps array
      integer, intent(in)                             :: nr, ni      ! indices of the cell
      real, intent(in)                                :: dtime       ! time step
      real, intent(in)                                :: realtime    ! physical time
      integer, intent(in)                             :: index_init  ! first index of particles in swrm array
      integer                                         :: nsws        ! number of representative particles in given cell
      integer                                         :: nri, nrk    ! indices of physical and representative particles choosen to the next collision
      real, dimension(:,:), allocatable               :: colrates    ! collision rates matrix
      real, dimension(:,:), allocatable               :: relvels     ! relative velocities matrix
      real, dimension(:,:), allocatable               :: accelncol   ! coagulation acceleration matrix (in the case of high mass ratio,
                                                                     ! instead of performing every collision separately, we group the collisions)
      real, dimension(:), allocatable                 :: stokesnr    ! stokes numbers of particles
      real, dimension(:), allocatable                 :: vs, vr      ! vertical settilng and radial drift velocities
      real, dimension(:), allocatable                 :: colri       ! collision rates for rps
      real                                            :: local_time, dt ! physical time spent in the cell, time step between two collisions
      real                                            :: rand        ! random number
      real                                            :: totr        ! total collision rate
      real                                            :: vn          ! maximum radial velocity from the pressure gradient
      real                                            :: Reynolds, v0, Vg2, veta, tL, teta ! relative velocities stuff
      real                                            :: lmfp, gasdens ! mean free path, gas density
      integer                                         :: i, k, l     ! counters
      integer, intent(out)                            :: kk ! collisions counter
      real                                            :: Kepler_freq, cs_speed ! Keplerian frequency and sound speed
      real                                            :: deltar, deltaz ! size of radial and vertical grid; for adaptative dmmax
      
      
      deltar = g%rup(nr) - g%rlo(nr)
      deltaz = g%zup(nr,ni) - g%zlo(nr,ni)
      
      ! calculation of some values needed for calculations
      gasdens = densg(g%rce(nr),g%zce(nr,ni),realtime)         ! gas density in the center of cell
      lmfp = mH2 / ( gasdens * AH2 )                           ! mean free path in gas in the center of cell
      Kepler_freq = omegaK(g%rce(nr), realtime)                ! keplerian frequency at the radial centre of cell
      cs_speed = cs(g%rce(nr), realtime)
      Reynolds = sqrt(0.5 * pi) * alphat * cs_speed / (Kepler_freq * lmfp) ! Reynolds number
      v0 = sqrt(alphat) * cs_speed                                 ! velocity of the largest eddy
      veta = v0 * Reynolds**(-0.25)                                ! velocity of the smallest eddy
      tL = 1./Kepler_freq                                        ! overturn time of the largest eddy
      teta = Reynolds**(-0.5) * tL                                 ! overturn time of the smallest eddy
      vn = 0.5*dlogPg(g%rce(nr), g%zce(nr,ni), realtime)*Pg(g%rce(nr), g%zce(nr,ni), realtime)/g%rce(nr) / &  ! maximum radial velocity from the pressure gradient
              gasdens / Kepler_freq
      Vg2 = 1.5*v0**2.0                                           ! turbulent velocity of gas squared

      !-------------------------------------------------------------------------------------

      swarms => bin(nr,ni)%p       ! points to the swarms array that contains only swarms located in the current cell nr,ni
      swarms(:)%coll_f = 0         ! setting collision flag to zero, will be updated to 1 if collisions happen
      nsws = size(swarms)          ! number of swarms in the current cell

      allocate (colrates(nsws,nsws), accelncol(nsws,nsws), relvels(nsws,nsws), stokesnr(nsws), vs(nsws), vr(nsws), colri(nsws))

      ! calculates initial Stokes number for particles and their settling and radial drift velocities
      do i = 1, nsws
         call stokes_nr_centr(i, swarms, stokesnr, lmfp, gasdens, Kepler_freq, cs_speed)
         call vel_vs_centr(nr, ni, i, stokesnr, Kepler_freq, vs)
         call vel_rd_centr(nr, i, stokesnr, vr, vn, realtime)
      enddo

      ! calculates relative velocities and collision rates matrix
      do i = 1, nsws
         call rel_vels(i, swarms, stokesnr, vr, vs, relvels, vn, Reynolds, veta, Vg2, tL, teta, Kepler_freq, cs_speed)
         call col_rates(i, swarms, relvels, colrates, accelncol, g%vol(nr,ni), deltar, deltaz)
      enddo

      ! collision rate for representative particles
      colri(:) = sum(colrates(:,:),dim=2)

      !------------ MAIN COLLISIONS LOOP ----------------------------------------------------
      local_time = 0.0
      kk = 0
      do while (local_time < dtime)
         totr = sum(colri)                 ! total collision rate
         call random_number(rand)
         dt = -1. * log(rand) / totr       ! time step between this and the next collision
         if (dt > dtime) then ! 0 or 1 collisions, decided by a random number
         call random_number(rand)
         if (rand> dtime/dt) then
            local_time = dtime
            cycle
         endif
         endif
         local_time = local_time + dt      ! update of the local time
         call choose_swarms(nri,nrk,colrates,colri,totr)


         call collision(nri,nrk,swarms,relvels, accelncol)

         call stokes_nr_centr(nri, swarms, stokesnr, lmfp, gasdens, Kepler_freq, cs_speed)
         call vel_vs_centr(nr, ni, nri, stokesnr, Kepler_freq, vs)
         call vel_rd_centr(nr, nri, stokesnr, vr, vn, realtime)

         call rel_vels(nri, swarms, stokesnr, vr, vs, relvels, vn, Reynolds, veta, Vg2, tL, teta, Kepler_freq,cs_speed)
         relvels(:,nri) = relvels(nri,:)
         colri(:) = colri(:) - colrates(:,nri)
         call col_rates(nri, swarms, relvels, colrates, accelncol, g%vol(nr,ni), deltar, deltaz)
         call col_rates_r(nri, swarms, relvels, colrates, accelncol, g%vol(nr,ni), deltar, deltaz)

         colri(:) = colri(:) + colrates(:,nri)
         colri(nri) = sum(colrates(nri,:))
         kk = kk + 1
      enddo
      !-------------------------------------------------------------------------------------

      deallocate (colrates, accelncol, relvels, stokesnr, vs, vr, colri)

      do k = 1, nsws
         if(swarms(k)%coll_f /= 0) then

            l = index_init

            do while (swrm(l)%idnr /= swarms(k)%idnr)
               l = l + 1
            enddo
            swrm(l) = swarms(k)
         else
            cycle
         endif   
      enddo

      nullify(swarms)

      return
   end subroutine mc_collisions

   ! calculating Stokes numbers of particle "i" IN THE CENTER OF CELL nr,ni
   subroutine stokes_nr_centr(i, swarms, stokesnr, lmfp, gasdens, Kepler_freq, cs_speed)
      implicit none
      type(swarm), dimension(:)                       :: swarms
      integer, intent(in)                             :: i
      real, dimension(:), allocatable                 :: stokesnr
      real, intent(in)                                :: lmfp, gasdens, Kepler_freq, cs_speed
      real                                            :: rad


      rad = (0.75 / pi / swarms(i)%dens)**third * swarms(i)%mass**third ! particle radius

      if (rad > 2.25 * lmfp) then ! Stokes regime
         stokesnr(i) = sqrt(2.*pi) * swarms(i)%dens * AH2 * rad**2. * Kepler_freq / (9. * mH2 * cs_speed )
      else                        ! Epstein regime
         stokesnr(i) = rad * swarms(i)%dens / (sqrt(8./pi) * cs_speed  * gasdens) * Kepler_freq
      endif

      return
   end subroutine stokes_nr_centr

   ! calculation of relative velocities of bodies:
   ! we take 5 sources:
   ! Brownian motion vB
   ! turbulence vT (Ormel & Cuzzi 2007), implementation stolen from Til Birnstiel
   ! radial drift vr
   ! vertical settling vs
   ! azimuthal drift vtan
   subroutine rel_vels(nl, swarms, stokesnr, vr, vs, relvels, vn, Reynolds, veta, Vg2, tL, teta, Kepler_freq, cs_speed)
      implicit none
      type(swarm), dimension(:), pointer              :: swarms
      integer, intent(in)                             :: nl
      real, dimension(:,:), allocatable               :: relvels
      real, dimension(:), allocatable                 :: stokesnr, vr, vs
      real, intent(in)                                :: vn
      real, dimension(size(swarms))                   :: vB2, vT2, vtan
      real                                            :: gts, lts
      real, intent(in)                                :: tL, teta
      real, intent(in)                                :: Reynolds
      real, intent(in)                                :: veta, Vg2
      integer                                         :: i
      real                                            :: St1, St2
      real                                            :: y_star, c1, c2, c3, c0, eps, hulp1, hulp2
      real, parameter                                 :: ya = 1.6
      real, intent(in)                                :: Kepler_freq, cs_speed

      ! Brownian motions
      vB2(:) = 8.* mH2 * cs_speed**2 * (swarms(nl)%mass+swarms(:)%mass) / (pi*swarms(nl)%mass*swarms(:)%mass)

      ! turbulence
      do i = 1, size(swarms)
         if (stokesnr(i) > stokesnr(nl)) then
            gts = stokesnr(i) / Kepler_freq
            lts = stokesnr(nl) / Kepler_freq
            St1 = stokesnr(i)
            St2 = stokesnr(nl)
         else
            gts = stokesnr(nl) / Kepler_freq
            lts = stokesnr(i) / Kepler_freq
            St1 = stokesnr(nl)
            St2 = stokesnr(i)
         endif
         if (gts < 0.2*teta) then
             vT2(i) = 1.5 *(veta/teta *(gts - lts))**2.0

         elseif (gts < teta/ya) then
             vT2(i) = Vg2 *(St1-St2)/(St1+St2)*(St1**2.0/(St1+Reynolds**(-0.5)) - St2**2.0/(St2+Reynolds**(-0.5)))
         elseif (gts < 5.0*teta) then
            !Eq. 17 of OC07. The second term with St_i**2.0 is negligible (assuming !Re>>1)
            !hulp1 = Eq. 17; hulp2 = Eq. 18

            hulp1 = ( (St1-St2)/(St1+St2) * (St1**2.0/(St1+ya*St1) - St2**2.0/(St2+ya*St1)) )!note the -sign
            hulp2 = 2.0*(ya*St1-Reynolds**(-0.5)) + St1**2.0/(ya*St1+St1) - St1**2.0/(St1+Reynolds**(-0.5)) +&
                                           St2**2.0/(ya*St1+St2) - St2**2.0/(St2+Reynolds**(-0.5))
            vT2(i) = Vg2 *(hulp1 + hulp2)

         elseif (gts < tL*0.2)  then
            eps=St2/St1!stopping time ratio
            vT2(i) = Vg2 *( St1*(2.0*ya - (1.0+eps) + 2.0/(1.0+eps) *(1.0/(1.0+ya) + eps**3.0/(ya+eps) )) )

         elseif (gts < tL) then
           !now y* lies between 1.6 (St1 << 1) and 1.0 (St1>=1). The fit below fits ystar to less than 1%
           c3 =-0.29847604
           c2 = 0.32938936
           c1 =-0.63119577
           c0 = 1.6015125
           y_star = c0 + c1*St1 + c2*St1**2.0 + c3*St1**3.0
           !we can then employ the same formula as before
           eps=St2/St1
           vT2(i) = Vg2 *( St1*(2.0*y_star - (1.0+eps) + 2.0/(1.0+eps) *(1.0/(1.0+y_star) + eps**3.0/(y_star+eps) )) )

         else
            vT2(i) = Vg2 *( 1.0/(1.0+St1) + 1.0/(1.0+St2) )
         endif

      enddo

      ! tangential
      vtan(:) = vn * ( 1. /(1.+stokesnr(:)**2.) - 1. / (1.+stokesnr(nl)**2.) )

      ! total
      relvels(nl,:) = sqrt(vB2(:)  + vT2(:) + (vs(:) - vs(nl))**2 + (vr(:) - vr(nl))**2 + vtan(:)**2)

      return
   end subroutine rel_vels

   ! calculation of the collision rates between representative particle nl and all physical particles
   subroutine col_rates(nl, swarms, relvels, colrates, accelncol, vol, deltar, deltaz)
      implicit none
      integer, intent(in)                             :: nl       ! number of line of colrates to update
      type(swarm), dimension(:), pointer              :: swarms
      real, dimension(:,:), allocatable               :: colrates
      real, dimension(:,:), allocatable               :: accelncol
      real, dimension(:,:), allocatable               :: relvels
      real, intent(in)                                :: vol      ! volume of the cell
      real, intent(in) :: deltar, deltaz ! width of grid


      ! basic eq. (see e.g. Drazkowska et al 2013, Eqs 19-20)
      colrates(nl,:) = swarms(:)%npar * relvels(nl,:) * pi**third * (0.75)**(2. * third) * &
                      ((swarms(nl)%mass/swarms(nl)%dens)**third + (swarms(:)%mass/swarms(:)%dens)**third)**2./vol

      
      ! instead of perform 1000 identical collisions, we divide the collision rate by 1000 but if the collision happens,
      ! we stick 1000 small particles at once
      where(swarms(:)%mass/swarms(nl)%mass < dmmax)
         accelncol(nl,:) = max(colrates(nl,:)*min(deltar/abs(swarms(nl)%velr),deltaz/abs(swarms(nl)%velz)), 1.)
         accelncol(nl,:) = min(accelncol(nl,:), swarms(nl)%mass*dmmax/swarms(:)%mass) ! we don't want a error larger than dmmax*100%
      elsewhere
         accelncol(nl,:) = 1. 
      endwhere

      colrates(nl,:) = colrates(nl,:) / accelncol(nl,:)
      ! if the representative particle represents less than 1 particle, the method is not valid anymore,
      ! so the collision rate is supressed
      where(swarms(:)%npar <= 1.0) colrates(nl,:) = 0.0

      return
   end subroutine col_rates

   ! for updating the collision rates after collision: calculating the collsion rates between physical particle nl
   ! and all the representative particles
   subroutine col_rates_r(nl, swarms, relvels, colrates, accelncol, vol, deltar, deltaz)
      implicit none
      integer, intent(in)                             :: nl       ! number of column of colrates to update
      type(swarm), dimension(:), pointer              :: swarms
      real, dimension(:,:), allocatable               :: colrates
      real, dimension(:,:), allocatable               :: accelncol
      real, dimension(:,:), allocatable               :: relvels
      real, intent(in)                                :: vol      ! volume of the cell
      real, intent(in) :: deltar, deltaz ! width of grid

      colrates(:,nl) = swarms(nl)%npar * relvels(nl,:) * pi**third * (0.75)**(2. * third) * &
                       ((swarms(nl)%mass/swarms(nl)%dens)**third + (swarms(:)%mass/swarms(:)%dens)**third)**2./vol

      ! adaptive dmmax for grouping collisions
      where(swarms(nl)%mass/swarms(:)%mass < dmmax)
         accelncol(:,nl) = max(colrates(:,nl) * min(deltar/abs(swarms(:)%velr),deltaz/abs(swarms(:)%velz)), 1.)
         accelncol(:,nl) = min(accelncol(:,nl), swarms(:)%mass * dmmax / swarms(nl)%mass)
         
      elsewhere
         accelncol(:,nl) = 1. 
      endwhere

      colrates(:,nl) = colrates(:,nl) / accelncol(:,nl)
      where(swarms(:)%npar <= 1.0) colrates(:,nl) = 0.0

      return
   end subroutine col_rates_r


   ! choosing particles to the next collision
   ! nri -> representative
   ! nrk -> physical
   subroutine choose_swarms(nri, nrk, colrates, ri, totrate)
      implicit none
      integer, intent(out)                            :: nri, nrk
      real, dimension(:,:), allocatable, intent(in)   :: colrates
      integer                                         :: i
      real, dimension(:), allocatable, intent(in)     :: ri
      real, intent(in)                                :: totrate
      real, dimension(2)                              :: rand
      real                                            :: fin

      call random_number(rand) ! drawing 2 random numbers

      ! choosing the representative particle: the higher ri, the higher chance that it's particle "i"
      rand(1) = rand(1) * totrate
      i = 1
      fin = ri(1)
      do while(rand(1) > fin)
         fin = fin + ri(i+1)
         i = i + 1
      enddo
      nri = i

      ! choosing the physical particle
      i = 1
      rand(2) = rand(2) * ri(nri)
      fin = colrates(nri,1)
      do while (rand(2) > fin)
         fin = fin + colrates(nri,i+1)
         i = i + 1
      enddo
      nrk = i

      return
   end subroutine choose_swarms

   ! performing the collision: deciding the collision outcome - put your collision model here
   ! only the representative particle is updated
   subroutine collision(nri,nrk,swarms,relvels, accelncol)
      implicit none
      integer, intent(in)                             :: nri, nrk
      type(swarm), dimension(:), pointer              :: swarms
      real, dimension(:,:), allocatable               :: accelncol
      real, dimension(:,:), allocatable               :: relvels
      real                                            :: rvel
      real                                            :: vmax

      rvel = relvels(nri,nrk)
      vmax = vfrag


      if (rvel < vmax) then
         if ((swarms(nri)%f_matrix<eps_fri) .and. (swarms(nrk)%f_matrix<eps_fri)) then ! if 2 pure rigids collide they bounce
            continue
         else
            call hit_and_stick(nri,nrk,swarms, accelncol)
         endif

      else if ((swarms(nri)%f_matrix > eps_fri)) then ! pure rigid does not fragment

         call fragmentation(nri, swarms)

      endif
      swarms(nri)%npar = mswarm / swarms(nri)%mass
      swarms(nri)%coll_f =  1
      return
   end subroutine collision


   ! sticking collision
   subroutine hit_and_stick(nri,nrk,swarms, accelncol)
      implicit none
      integer, intent(in)                             :: nri, nrk
      type(swarm), dimension(:), pointer              :: swarms
      real, dimension(:,:), allocatable               :: accelncol
      real                                            :: totmass

      totmass = swarms(nri)%mass + accelncol(nri, nrk) * swarms(nrk)%mass

      ! if similar particles (90% of similarity) update rigid mass of fragile. To keep the simmetry of fragmenting population
      if ((swarms(nri)%rigid .eq. 0) .and. (swarms(nrk)%rigid .eq. 1)) then
         if ((1.1*swarms(nri)%mass>swarms(nrk)%mass) .and. (0.9*swarms(nri)%mass<swarms(nrk)%mass)) then ! similar in mass
               swarms(nri)%rigid_mass = swarms(nrk)%rigid_mass
         end if
      endif

      swarms(nri)%f_CH  = (swarms(nri)%f_CH*swarms(nri)%mass + &
                           accelncol(nri, nrk)*swarms(nrk)%f_CH*swarms(nrk)%mass)/totmass
      swarms(nri)%f_matrix = 1.0 - swarms(nri)%f_CH
      
      if (swarms(nri)%f_matrix < 0.) swarms(nri)%f_matrix = 0. ! to avoid numerical problems
      if (swarms(nri)%f_matrix > 1.) swarms(nri)%f_matrix = 1. ! to avoid numerical problems

      ! internal density
      swarms(nri)%dens = 1./(swarms(nri)%f_CH/ridens + swarms(nri)%f_matrix/matdens)

      swarms(nri)%mass = totmass

      return
   end subroutine hit_and_stick

   ! function to calculate the mass of the fragment
   real function f_mf(mmin, mmax, exp)
      implicit none
      real, intent(in) :: mmin, mmax, exp
      real:: ran
      call random_number(ran)
      f_mf = (ran*(mmax**exp - mmin**exp) +  mmin**exp)**(1./exp)
      return 
   end function

   ! fragmentation collision
   ! put your fragment size distribution here:
   ! n(m) ~ m^(kappa - 2)
   subroutine fragmentation(nri,swarms)
      implicit none
      integer, intent(in)                             :: nri
      real                                            :: mf
      type(swarm), dimension(:), pointer              :: swarms
      real, parameter                                 :: kappa = 1./6.  ! n(m) ~ m^(kappa - 2)

      ! calculate mass of fragment
      mf = f_mf(m0, swarms(nri)%mass, kappa)

      if ((swarms(nri)%f_CH>eps_fri)) then ! if False is pure matrix 

         if (mf < swarms(nri)%rigid_mass/swarms(nri)%f_CH) then ! wants to fragment to masses below rigid

            if (swarms(nri)%rigid==1) then
               mf = swarms(nri)%rigid_mass ! always track same size 
               swarms(nri)%f_CH = 1. ! pure rigid
               swarms(nri)%f_matrix = 0.
               swarms(nri)%dens = ridens
            else
               mf = f_mf(m0, min(swarms(nri)%f_matrix*swarms(nri)%mass, swarms(nri)%rigid_mass/swarms(nri)%f_CH), kappa)               
               swarms(nri)%f_CH = 0. ! pure fragile
               swarms(nri)%f_matrix = 1.
               swarms(nri)%dens = matdens
            endif
         endif
      endif

      swarms(nri)%mass = mf
      swarms(nri)%mass = max(swarms(nri)%mass, m0) ! to avoid numerical problems
      
      return
   end subroutine fragmentation

   ! vertical settling velocity at the centre of grid
   subroutine vel_vs_centr(nr, ni, i, stokesnr, Kepler_freq, vs)
      implicit none
      integer, intent(in)                             :: nr, ni
      real, dimension(:), allocatable, intent(in)     :: stokesnr
      real, dimension(:), allocatable                 :: vs
      integer, intent(in)                             :: i       ! index of particle
      real, intent (in)                               :: Kepler_freq
      vs(i) = g%zce(nr,ni) * Kepler_freq * min(stokesnr(i), 0.5)

      return
   end subroutine vel_vs_centr

   ! velocity of radial drift at the centre of grid
   subroutine vel_rd_centr(nr, i, stokesnr, vr, vn, realtime)
      implicit none
      integer, intent(in)                             :: nr
      real, dimension(:), allocatable                 :: stokesnr
      real, dimension(:), allocatable                 :: vr
      real, intent(in)                                :: vn
      integer, intent(in)                             :: i      ! index of particle
      real, intent(in)                                :: realtime

      vr(i) = 2. * vn * stokesnr(i) / (stokesnr(i)**2. + 1.) + &
               vgas(g%rce(nr), realtime) / (1. + stokesnr(i) * stokesnr(i))

      return
   end subroutine vel_rd_centr


end
