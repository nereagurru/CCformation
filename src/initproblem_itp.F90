! this module should contain all the initial conditions for dust
! different funcitons for 0D (init_swarms_OD) and for 2D simulations (init_swarms_PB2)
module initproblem

   use constants,  only: pi, third, AU, year
   use discstruct, only: cs, omegaK, Msolid_PB, particle_formation_PB
   use discstruct, only: Temp, sigmag, dlogPg, stokes_frag, stokes_drift, St_limited
   use parameters, only: minrad0, maxrad0, alphat, r0, ridens, matdens, prop_CH, aCH_min, aCH_max, zeta
   use types,      only: swarm
   use advection,  only: stokesnr
   
   implicit none
   
   private
   public   :: init_random_seed, mswarm, m0, mCH_max, mCH_min
#ifdef ZERO_D
   public   :: init_swarms_OD
#else
   public   ::init_swarms_PB2
#endif

   real     :: mswarm, m0
   real     :: mCH_min, mCH_max

   contains

#ifdef ZERO_D
   subroutine init_swarms_OD(Ntot, swrm, time)  
      implicit none
      type(swarm), dimension(:), allocatable, target   :: swrm        ! list of all swarms in the simulation
      integer, intent(inout)                           :: Ntot        ! number of all swarms
      real, intent(inout)                              :: time        ! initial time
      real                                             :: mdust       ! mass of the dust in the simulatated domain
      real, dimension(2)                               :: rand        ! random numbers
      real                                             :: Hg, Hd      ! pressure height scale of the gas
      integer                                          :: i
      real                                             :: x
      real                                             :: amin ! if a_limited > ari_min amin=ari_min
      real :: Stlim, power
            ! total mass of dust maxrad0>r>min_loc
      mdust = Msolid_PB(minrad0, maxrad0, time) 

      ! mass of one swarm
      mswarm = mdust/ real(Ntot)

      write(*,*) "total solid mass is ", mdust, "and mswarm is ", mswarm
      ! monomer mass
      m0 = 4. * third * pi * r0**3 * matdens

      ! refractory mass minimum and maximum
      mCH_min = 4. * third * pi * aCH_min**3 * ridens
      mCH_max = 4. * third * pi * aCH_max**3 * ridens

      if (.not.allocated(swrm)) allocate(swrm(Ntot))

      ! initializing the particles
      do i = 1, Ntot
         swrm(i)%idnr = i               ! index
         swrm(i)%plt = 0                ! they are not planetesimals
         swrm(i)%tfor = time            ! time when they are initialized

         ! choose initial position randomly (only for 0D, density expected to be near constant in small range)
         call random_number(x)
         swrm(i)%rdis =  minrad0 + (maxrad0 - minrad0)*x 
         call random_number(x)
         amin = (x*(aCH_max**(zeta+4.) - aCH_min**(zeta+4.)) + aCH_min**(zeta+4.))**(1./(zeta+4.))   ! size of the rigid monomer
         swrm(i)%rigid_mass = (4.*third*pi*ridens)*amin**3.                                          ! mass of the rigid monomer

         ! set if static properties are 0 or 1
         call random_number(x)
         if (x<(prop_CH)) then ! rigid
            swrm(i)%rigid = 1
            swrm(i)%f_CH = 1.
            swrm(i)%f_matrix = 0.
            swrm(i)%dens = ridens
            swrm(i)%mass = swrm(i)%rigid_mass

         else ! fragile nature
            swrm(i)%rigid = 0
            swrm(i)%f_CH = 0.
            swrm(i)%f_matrix = 1.
            swrm(i)%dens = matdens
            amin = r0
            swrm(i)%mass = 4. * third * pi * swrm(i)%dens* amin**3.

         endif

         swrm(i)%zdis = 0.0 ! only to compute Hd. Then it will be changed
         swrm(i)%stnr = stokesnr(swrm(i),swrm(i)%tfor)

         swrm(i)%npar = mswarm / swrm(i)%mass
         Hg = cs(swrm(i)%rdis, swrm(i)%tfor) / omegaK(swrm(i)%rdis, swrm(i)%tfor)
         Hd = Hg*sqrt(alphat/(alphat + swrm(i)%stnr))            ! dust vertical scale height
         call random_number(rand)
         swrm(i)%zdis = 0.0005*AU * sqrt(-2.*log(rand(1))) * cos(2.*pi*rand(2))
         swrm(i)%stnr = stokesnr(swrm(i),swrm(i)%tfor) ! recalculate to estimate spreed
         swrm(i)%velr = 1.e-20
         swrm(i)%velz = 1.e-20
         swrm(i)%coll_f = 0

      enddo
   

      return   

   end subroutine init_swarms_OD
#else
      ! init function for executing main_PB
   subroutine init_swarms_PB2(Ntot, swrm, time)
      implicit none
      type(swarm), dimension(:), allocatable, target   :: swrm        ! list of all swarms in the simulation
      integer, intent(inout)                              :: Ntot        ! number of all swarms
      real, intent(inout)                              :: time        ! initial time
      real                                             :: mdust       ! mass of the dust in the simulatated domain
      real, dimension(2)                               :: rand        ! random numbers
      real                                             :: Hg, Hd      ! pressure height scale of the gas
      integer                                          :: i
      real, dimension(:), allocatable                  :: rfor_arr
      real                                             :: x
      real                                             :: amin ! if a_limited > ari_min amin=ari_min
      real :: ari_max, pri, ari_static
      real :: Stlim, power 

      ! total mass of dust maxrad0>r>min_loc
      mdust = Msolid_PB(minrad0, maxrad0, time) 

      ! mass of one swarm
      mswarm = mdust/ real(Ntot)

      ! monomer mass
      m0 = 4. * third * pi * r0**3 * matdens

      ! refractory mass minimum and maximum
      mCH_min = 4. * third * pi * aCH_min**3 * ridens
      mCH_max = 4. * third * pi * aCH_max**3 * ridens

      if (.not.allocated(swrm)) allocate(swrm(Ntot))
      allocate(rfor_arr(Ntot))

      ! calculate positions to initialize particles
      call particle_formation_PB(1, Ntot, mswarm, time, minrad0, rfor_arr, 1.)
    
      ! initializing the particles
      do i = 1, Ntot
         swrm(i)%idnr = i             ! index
         swrm(i)%plt = 0              ! particle is not a planetesimal
         swrm(i)%tfor = time          ! formation time
         swrm(i)%rdis = rfor_arr(i)   ! initial position

         call St_limited(swrm(i)%rdis, swrm(i)%tfor, Stlim, power)
         swrm(i)%stnr = Stlim         ! set initial stokes number. This will be corrected at the end

         ! max size value for rigid size distribution
         ari_max = min(aCH_max, as_limited(swrm(i)%rdis, swrm(i)%tfor, ridens, swrm(i)%stnr))
         ! probability for initiating a rigid particle
         pri = (ari_max**(zeta+4.)-aCH_min**(zeta+4.))/(aCH_max**(zeta+4.) - aCH_min**(zeta+4.))
         swrm(i)%f_CH = prop_CH*pri               ! update probability. Lower in the outer regions
         swrm(i)%f_matrix = 1.- swrm(i)%f_CH         

         swrm(i)%dens = 1./(swrm(i)%f_CH/ridens + swrm(i)%f_matrix/matdens)        ! internal density

         amin = as_limited(swrm(i)%rdis, swrm(i)%tfor, swrm(i)%dens, swrm(i)%stnr) ! we assume Epstein regime

         ! calculate rigid mass of the monomer from size distribution
         call random_number(x)
         ari_static = (x*(ari_max**(zeta+4.) - aCH_min**(zeta+4.)) + aCH_min**(zeta+4.))**(1./(zeta+4.))
         swrm(i)%rigid_mass = (4.*third*pi*ridens)*ari_static**3.


         ! can rigid exist?

         if (St_eps(swrm(i)%rdis, swrm(i)%tfor, ridens, ari_static) .gt. swrm(i)%stnr) then
            ! particle needs to be fragile due to growth barriers
            swrm(i)%rigid = 0
            !swrm(i)%rigid_mass = 0.
            swrm(i)%dens = matdens
            swrm(i)%f_CH = 0.
            swrm(i)%f_matrix = 1.
            call random_number(x)
            amin = (x*(amin**power- r0**power) + r0**power)**(1./power) ! assume size distribution from frag or drift limited
            swrm(i)%mass = 4. * third * pi * swrm(i)%dens* amin**3      ! total mass of aggregate
         else ! rigid can in principle exist

            ! choose particle nature (rigid or fragile)
            call random_number(x)
            if (x<(prop_CH*pri)) then ! rigid
               swrm(i)%rigid = 1
               !choose chondritic aggregate size from frag-sized distribution
               call random_number(x)
               amin = (x*(amin**power- r0**power) + r0**power)**(1./power)

               ! can the particle be part of a chondritic aggregate?
               if (swrm(i)%rigid_mass .gt. 4*pi*third*swrm(i)%dens*swrm(i)%f_CH*amin**3.) then
                  ! particle is a pure rigid
                  swrm(i)%dens = ridens
                  swrm(i)%f_CH = 1.
                  swrm(i)%f_matrix = 0.
                  swrm(i)%mass = swrm(i)%rigid_mass
               else
                  ! particle is a chondritic aggregate
                  swrm(i)%mass = 4. * third * pi * swrm(i)%dens* amin**3 ! total mass of aggregate
               endif

            else ! fragile nature
               swrm(i)%rigid = 0
               !swrm(i)%rigid_mass = 0.
               !choose chondritic aggregate size from frag-sized distribution
               call random_number(x)
               amin = (x*(amin**power- r0**power) + r0**power)**(1./power)

               ! can the particle be part of a chondritic aggregate?
               if (swrm(i)%rigid_mass .gt. 4*pi*third*swrm(i)%dens*swrm(i)%f_CH*amin**3.) then
                  ! particle is pure fragile
                  swrm(i)%dens = matdens
                  swrm(i)%f_CH = 0.
                  swrm(i)%f_matrix = 1.
                  swrm(i)%mass = 4. * third * pi * swrm(i)%dens* amin**3
               else
                  ! particle is a chondritic aggregate
                  swrm(i)%mass = 4. * third * pi * swrm(i)%dens* amin**3
               endif
            endif

         endif

         swrm(i)%zdis = 0.0 ! only to compute Hd. Then it will be changed
         swrm(i)%stnr = stokesnr(swrm(i),swrm(i)%tfor)   

         swrm(i)%npar = mswarm / swrm(i)%mass
         Hg = cs(swrm(i)%rdis, swrm(i)%tfor) / omegaK(swrm(i)%rdis, swrm(i)%tfor)
         Hd = Hg*sqrt(alphat/(alphat + swrm(i)%stnr))
         call random_number(rand)
         swrm(i)%zdis = Hd * sqrt(-2.*log(rand(1))) * cos(2.*pi*rand(2))
         swrm(i)%stnr = stokesnr(swrm(i),swrm(i)%tfor) ! update final stokes number
         swrm(i)%velr = 1.e-20
         swrm(i)%velz = 1.e-20
         swrm(i)%coll_f = 0

      enddo

      deallocate(rfor_arr)
      return   
   end subroutine init_swarms_PB2
#endif
   
   
   ! initialize the random number generator
   subroutine init_random_seed
      implicit none
      integer                            :: i, n, clock
      integer, dimension(:), allocatable :: seed

      call random_seed(size = n)
      allocate(seed(n))

      call system_clock(count=clock)

      seed = 37 * [(i - 1, i = 1, n)]
      seed = seed + clock

      call random_seed(put = seed)
      deallocate(seed)

   end subroutine init_random_seed

   ! Assume Epstein regime
   real function as_limited(x, time, rhoi, St)
      implicit none
      real, intent(in)        :: x, time, rhoi, St

      as_limited = max((2./pi)*sigmag(x, time)*St/rhoi, r0)
      return
   end function 

   ! Assume Epstein regime
   real function St_eps(x, time, rhoi, as)
      implicit none
      real, intent(in)        :: x, time, rhoi, as

      St_eps = (pi/2.)*rhoi*as/sigmag(x, time)
      return
   end function 

   end
