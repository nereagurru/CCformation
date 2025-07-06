! This code is an extended version of:
! Drazkowska, Windmark & Dullemond (2013) A&A 556, A37, and Vaikundaraman et al. (2025, in prep)

! Authors of this version: Nerea Gurrutxaga, Joanna Drazkowska, Vignesh Vaikundaraman
! Max Planck Institute for Solar System Research, Göttingen, Germany
!
! This code performs a 2D simulation of dust evolution in a protoplanetary disk.
! The gas disk is an input and is not evolved by mcdust (discstruct_itp module)
! The dust is treated as representative particles (RPs) undergoing advection (advection module)
! as well as collisions performed with Monte Carlo algorithm (collisions module)
! To perform collisions the RPs are binned using an adaptive grid (grid module)
!
! 

program main
   use constants
   use advection,    only: mc_advection, mc_advection_element, update_St
   use collisions,   only: mc_collisions
   use grid,         only: g, make_grid, deallocate_grid
   use initproblem,  only: init_random_seed, mswarm, m0, mCH_min, mCH_max
#ifdef ZERO_D
   use initproblem,  only: init_swarms_OD
#else
   use initproblem,  only:init_swarms_PB2
#endif
 
   use discstruct,   only: init_struct, end_struct, omegaK, sigmag, &
                           Msolid_PB

   use parameters,   only: read_parameters, Ntot, nz, nr, dtime, fout, tend, smallr, restart, &
                           output_path, db_data, matdens, ridens, r0, stmin, &
                           datadir, alphat, maxrad0, aCH_min, aCH_max
#ifdef ZERO_D
   use parameters,   only: dtg, minrad0
#endif
#ifdef GLOBAL
   use parameters,   only: feedingdir
#endif

#ifdef LOCAL_SIM
   use parameters,   only: feedingdir
   use parameters,   only: Nmax
#endif   

   use timestep,     only: time_step
   use types
   use hdf5
   use hdf5output, only: hdf5_file_write, hdf5_file_t, hdf5_file_read
   use ppplts

   implicit none
   !
   ! Declaration of arrays and simulation variables
   !
   ! Array of individual dust particles (swarms)
   type(swarm), dimension(:), allocatable, target :: temp_swrm     ! Temporary swarms
   type(swarm), dimension(:), allocatable, target :: sim_swrm      ! Main simulation swarms
   type(swarm), dimension(:), allocatable, target :: temp_arr      ! Temporary array for operations
   integer, dimension(:,:), allocatable :: ncolls     ! Collision matrix between particle bins
   
   ! Swarms binned into spatial cells (2D grid)
   type(list_of_swarms), dimension(:,:), allocatable, target :: bin  
   
   ! Swarms binned into radial zones (1D array)
   type(list_of_swarms), dimension(:), allocatable :: rbin           

   ! HDF5 output file handler
   type(hdf5_file_t) :: file                                        
   
   ! General simulation parameters
   real                       :: total               ! Total time of execution
   real(kind=4), dimension(2) :: elapsed             ! Elapsed time (e.g. for profiling)
   integer                   :: i, j, iter, w, nparts ! Loop counters and particle count
   
   ! Physical simulation time parameters
   real :: time = 2.19e6*year        ! Initial simulation time (2 Myr after CAI)
   real :: time_no                   ! Auxiliary time variable
   real :: timeofnextout = 0.0       ! Time scheduled for next output
   real :: resdt = 0.1*year          ! Physical time step
   integer :: nout = 0, kont         ! Output index and control integer
   real :: totmass                   ! Total mass of dust beyond evaporation line
   
   ! Input / Output file controls
   character(len=100) :: ctrl_file     ! Input parameter file name
   character(len=100) :: open_file     ! Possibly the simulation restart file
   
   ! Parameters for local simulation runs
#ifdef LOCAL_SIM
   integer :: imax                   ! Max index for loop/array
   integer :: Nfeed, Nrem, Nin       ! Feeding/removal/input particle counters
   integer :: Nplt_new, Nplt_tot     ! New and total planetesimals
   real :: x                         ! Variable for Random number generation
   real :: mass_add                  ! Added mass from global simulation for each timestep
   
   real :: prob                      ! Probability (e.g., for stochastic process)
   real :: Mfeed, Mfeed_ri           ! Mass entering feeding zones
   real :: Mleak, Mleak_ri           ! Mass leaking out
   real :: Mplt, Mplt_ri             ! Planetesimal mass
   real :: mswarm_old                ! Swarm mass before interaction
   
   integer :: i_plt, jadd            ! Loop variables or particle counters
   integer :: Nleak, Nleak_ri        ! Number of leaking swarms
   integer :: Nfeed_ri               ! Number of feeding swarms in RI region
   real :: mswarm_feeding            ! Feeding swarm mass
   
   ! Arrays for feeding zone particles
   type(swarm), dimension(:), allocatable, target :: feeding_swrm
   type(swarm), dimension(:), allocatable, target :: feeding_swrm_tot
   type(swarm), dimension(:), allocatable, target :: feeding_lagun
   type(swarm), dimension(:), allocatable, target :: temp             ! Temporary swarms

   integer :: Nfluxtot               ! Total number of flux events
   integer :: nout_feed              ! Output index for feeding diagnostics
   integer :: Ncount                 ! General counter
   integer :: Nlagun                 ! Number of particles in laguna
   integer :: idx                    ! General-purpose index
   character(len=15) :: filename     ! Output file name (short)
#endif

   ! related planetesimal formation
   real :: plts_tot = 0.0              ! Total mass in planetesimals
   real :: etamax                      ! Maximum value of metallicity
   real :: Zcrit                    ! Critical metallicity for planetesimal formation
   integer :: stabilize             ! Number of collision loops before acvection starts

! Parameters for global simulation runs
#ifdef GLOBAL
   real :: last_time               ! Time of last update/output
   integer :: Nfluxtot, Nflux      ! Flux-related counters
   
   ! Feeding swarms for global simulation
   type(swarm), dimension(:), allocatable, target :: feeding_swrm
   character(len=15) :: filename    ! Output filename for global mode
#endif

! Parameters for test cases with constant flux
#ifdef TEST_CONST_FLUX
   real :: prop_const               ! Contant rigid mass fraction of thefeeding flux
   integer :: Nfeed_ri_lagun        ! helper variable for calculating feeding rigid mass *"lagun=helper" in basque
#endif

   ! random number generator initialization
   call init_random_seed

   ! reading parameters and initializing the simulation
   call get_command_argument(1, ctrl_file)
   call read_parameters(ctrl_file)
   call init_struct()

   ! restart simulation?
   if (restart) then
      stabilize = 0
      write(*,*) ' Reading restart...'
      open_file = 'restart.h5'
      call hdf5_file_read(open_file, Ntot, temp_swrm, nout, mswarm, time, resdt &
#ifdef LOCAL_SIM
                          &, imax, .False. &
#endif
                          &)
      
      timeofnextout = time + dtime
      m0 = 4. * third * pi * r0**3 * matdens ! CREATE a function in initproblem to calculate this
      mCH_min = 4. * third * pi * (aCH_min)**3 * ridens
      mCH_max = 4. * third * pi * (aCH_max)**3 * ridens

      do i=1, size(temp_swrm)
         temp_swrm(i)%f_matrix = 1. - temp_swrm(i)%f_CH
         temp_swrm(i)%dens = 1./(temp_swrm(i)%f_CH/ridens + temp_swrm(i)%f_matrix/matdens)
         temp_swrm(i)%npar = mswarm / temp_swrm(i)%mass
         temp_swrm(i)%velr = 0.0 
         temp_swrm(i)%velz = 0.0 
      enddo
      write(*,*) 'init time is ', time/year

#ifdef GLOBAL
      ! initiate Nflux, Nfluxtot and last_time for feeding flux
      Nflux = 0
      Nfluxtot = count(temp_swrm(:)%rdis<smallr)
      last_time = time
  
#endif

   else

      stabilize = 100

      write(*,*) 'init time is ', time/year
#ifdef LOCAL_SIM
      Ntot = Nmax


#ifdef ZERO_D
      call init_swarms_OD(Ntot, temp_swrm, time)
#else
      call init_swarms_PB2(Ntot, temp_swrm, time)
#endif


#ifdef GLOBAL
      last_time = time
      Nfluxtot = 0
      Nflux = 0
#endif

   endif

   ! calculate critical dust-to-gas ratio for planetesimal formation
   Zcrit = 10.**(0.15*(log10(alphat))**2. -0.24*log10(stmin)*log10(alphat) -1.48*log10(stmin) &
            + 1.18*log10(alphat)) ! Lim et al. (2024) eq. 19

   


#ifdef LOCAL_SIM
if (restart) then                                  ! restart simulation
   Nplt_tot = count(temp_swrm(:)%plt .eq. 1)       ! How many planetesimals formed so far

   Nrem =  Nmax - (size(temp_swrm) - Nplt_tot)     ! How many need to be added/removed to keep constant number of particles?
     

   allocate(sim_swrm(Nmax))
   if (Nrem>0) then                                ! add Nrem particles by creating copies
      j = 0                                        ! counter for all particles
      jadd = 0                                     ! counter for extra
      i_plt = 0                                    ! counter for particles with plt==1
      allocate(temp_arr(size(temp_swrm)))          ! create a temporary array
      temp_arr = temp_swrm
      deallocate(temp_swrm)
      allocate(temp_swrm(Nplt_tot+Nmax))
      do i=1, size(temp_arr)
         if (temp_arr(i)%rdis > smallr) then       ! if particle outside evaporation radius, save particle
               prob = real(Nrem-jadd)/real(Nmax-Nrem-j)  ! probability of copying this particle depends on how many we copied before (we want an exact number!)
               j = j + 1
               sim_swrm(j+jadd) = temp_arr(i)            ! save particle in sim_swrm
               temp_swrm(Nplt_tot+j+jadd) = temp_arr(i)  ! save particle in temp_swrm (where plt_arr + sim_swrm)
               call random_number(x)
               if (x<prob) then                          ! do we copy this particle?
                  Ntot = Ntot + 1 
                  jadd = jadd + 1
                  sim_swrm(j+jadd) = temp_arr(i)
                  sim_swrm(j+jadd)%idnr = Ntot
                  temp_swrm(Nplt_tot+j+jadd) = temp_arr(i)
                  temp_swrm(Nplt_tot+j+jadd)%idnr = Ntot
                  
               endif
         else if (temp_arr(i)%plt .eq. 1) then    ! if particle is a planetesimal, save it in temp_swrm
            i_plt = i_plt + 1
            temp_swrm(i_plt) = temp_arr(i)
         endif
      enddo
      j = j + jadd
      deallocate(temp_arr)
   else if (Nrem<0) then ! remove extra particles; same as before but choose randomly the one to remove except for copying
      j = 0
      jadd = 0
      i_plt = 0
      allocate(temp_arr(size(temp_swrm)))
      temp_arr = temp_swrm
      deallocate(temp_swrm)
      allocate(temp_swrm(Nplt_tot+Nmax))
      do i=1, size(temp_arr)
         if (temp_arr(i)%rdis > smallr) then
            write(*,*) j, jadd
            prob = real(-Nrem-jadd)/real(Nmax-Nrem-j)
            call random_number(x)
            if (x<prob) then
               jadd = jadd + 1
               j = j + 1
               write(*,*) j, jadd
            else 
               j = j + 1
               sim_swrm(j-jadd) = temp_arr(i)
               temp_swrm(Nplt_tot+j-jadd) = temp_arr(i)
            endif
         else if (temp_arr(i)%plt .eq. 1) then
            i_plt = i_plt + 1
            temp_swrm(i_plt) = temp_arr(i)
         endif
      enddo
      j = j - jadd
      deallocate(temp_arr) 
   else                                      ! don't add/remove particles    
      j = 0
      i_plt = 0
      allocate(temp_arr(size(temp_swrm)))   
      temp_arr = temp_swrm
      do i=1, size(temp_arr)
         if (temp_arr(i)%rdis > smallr) then 
            j = j + 1
            sim_swrm(j) = temp_arr(i)
            temp_swrm(Nplt_tot+j) = temp_arr(i)
         else if (temp_arr(i)%plt .eq. 1) then
            i_plt = i_plt + 1
            temp_swrm(i_plt) = temp_arr(i)
         endif
      enddo
      deallocate(temp_arr) 

   endif

  ! create array of particles that will be adding later from global simulation
  ! find the output that matches with out simulation
  nout_feed = 0
  do while (.True.)
      write(filename,'(a7,i5.5,a3)') 'swarms-',nout_feed,'.h5'
      open_file = trim(feedingdir)//trim(filename)
      write(*,*) open_file
      call hdf5_file_read(open_file, Nfluxtot, feeding_swrm_tot, nout_feed, mswarm_feeding, time_no, time_no, imax, .True.)
      nout_feed = nout_feed + 1
      if (maxval(feeding_swrm_tot(:)%tfor) > time) exit
      deallocate(feeding_swrm_tot)

  end do


#ifdef TEST_CONST_FLUX
   prop_const = 0.5
#endif
   ! correction in case some swarms in the feeding flux were formed before "time"
   Ncount = 0
   allocate(temp(size(feeding_swrm_tot)))
   do i = 1, size(feeding_swrm_tot)
      if ((feeding_swrm_tot(i)%tfor>=time)) then
#ifdef TEST_CONST_FLUX
      feeding_swrm_tot(i)%f_CH = prop_const
      feeding_swrm_tot(i)%f_matrix = 1.-prop_const

      call random_number(x)
      if (x<prop_const) then
         feeding_swrm_tot(i)%rigid = 1
         if (prop_const*feeding_swrm_tot(i)%mass < feeding_swrm_tot(i)%rigid_mass) then
            feeding_swrm_tot(i)%mass = feeding_swrm_tot(i)%rigid_mass
            feeding_swrm_tot(i)%f_matrix = 0.
            feeding_swrm_tot(i)%f_CH = 1.
         endif
      else
         feeding_swrm_tot(i)%rigid = 0
         if (prop_const*feeding_swrm_tot(i)%mass < feeding_swrm_tot(i)%rigid_mass) then
            feeding_swrm_tot(i)%f_matrix = 1.
            feeding_swrm_tot(i)%f_CH = 0.
         endif
      endif
      feeding_swrm_tot(i)%dens = 1./(feeding_swrm_tot(i)%f_CH/ridens + (1.-feeding_swrm_tot(i)%f_matrix)/matdens)

#endif
         Ncount = Ncount + 1
         temp(Ncount) = feeding_swrm_tot(i)
      endif
   enddo
   deallocate(feeding_swrm_tot)
   allocate(feeding_swrm_tot(Ncount))
   feeding_swrm_tot(:) = temp(:Ncount)
   Nfluxtot = Ncount
   deallocate(temp)

else      ! no restart
   Nplt_tot = 0
   allocate(sim_swrm(Nmax))
   sim_swrm = temp_swrm

   ! choose first feeding file
   nout_feed = 0
   write(filename,'(a7,i5.5,a3)') 'swarms-',nout_feed,'.h5'
   open_file = trim(feedingdir)//trim(filename)
   write(*,*) open_file

   call hdf5_file_read(open_file, Nfluxtot, feeding_swrm_tot, nout_feed, mswarm_feeding, time_no, time_no, imax, .True.)
   Nfluxtot = size(feeding_swrm_tot)
   nout_feed = nout_feed + 1


endif

! the following variables save info about cumulative mass that has entered or leaved the simulation
Mfeed = 0.
Mfeed_ri = 0.
Mleak = 0.
Mleak_ri = 0.
Mplt = 0.
Mplt_ri = 0.
imax = size(temp_swrm)

! from LOCAL_SIM ifdef
#else 
   allocate(sim_swrm(size(temp_swrm)))
   sim_swrm = temp_swrm 
   deallocate(temp_swrm)
#endif


   ! initial time, the formation of the first swarms
   timeofnextout = time

   ! Making grid for the first time
   call make_grid(sim_swrm, bin, rbin, nr, nz, smallr,totmass, ncolls)

   ncolls(:,:) = 1
   do i = 1, size(g%rlo)
      kont=0
      do j=1, size(g%zce,dim=2)
         kont = kont + size(bin(i, j)%p)
      enddo
   enddo


   iter = 0
   ! ------------------- MAIN LOOP -------------------------------------------------------------------------------------
   do while (time < tend)

#ifdef ZERO_D
      resdt = 1./omegaK(minval(sim_swrm(:)%rdis), time)
      write(*,*) 'time is ', time/year
      mswarm = dtg*Msolid_PB(minrad0, maxrad0, time)/size(sim_swrm)
      sim_swrm(:)%npar = mswarm/sim_swrm(:)%mass
#else
      ! determining the time step
      if (iter == 0 .or. (stabilize > 0)) then 
         if (restart .eqv. .False.) resdt = 1./omegaK(minval(sim_swrm(:)%rdis), time)
      else 
         call time_step(bin, ncolls, timeofnextout-time, resdt)
      end if
#endif

      if (db_data) then
         open(23,file=trim(output_path)//trim('timestep.dat'),status='unknown',position='append')
         write(23,*) time/year,resdt/year
         close(23)
      endif


      ! producing output
      if ((modulo(iter,fout) == 0 .or. time >=timeofnextout) .and. time >= 2.19e6*year) then
#ifdef LOCAL_SIM
            call update_St(temp_swrm, time)
            call hdf5_file_write(file, datadir, temp_swrm, time, resdt, 'create', nout, mswarm, Ntot, imax)
#else
            call update_St(sim_swrm, time)
            call hdf5_file_write(file, datadir, sim_swrm, time, resdt, 'create', nout, mswarm, Ntot)
#endif



#ifdef GLOBAL
            allocate(feeding_swrm(Nflux))
            Nfluxtot = Nfluxtot + Nflux
            j = 0
            do i = 1, Nfluxtot
               if ((sim_swrm(i)%tfor > last_time) .and. (sim_swrm(i)%tfor <= time)) then
                  j = j + 1
                  feeding_swrm(j) = sim_swrm(i)
               endif
            enddo
            if (j .ne. Nflux) stop

            ! write all of them at once
            if (Nflux > 0) then
               call hdf5_file_write(file, feedingdir, feeding_swrm, time, resdt, 'create', nout-1, mswarm, Nflux)
            endif
            last_time = time
            deallocate(feeding_swrm)
            Nflux = 0
#endif
            write(*,*) 'Time: ', time/year, 'produced output: ',nout
            open(23, file=trim(output_path)//trim('timesout.dat'),status='unknown',position='append')
            write(23,*) 'time: ', time/year, 'produced output: ',nout
            close(23)

            timeofnextout = time+dtime
            nout = nout + 1


            open(25,file=trim(output_path)//trim('collision_info.dat'),status='unknown',position='append')
            
            do i = 1, size(g%rce)
               do j = 1, size(g%zce,dim=2)
                  if (.not.allocated(bin(i, j)%p)) cycle
                  write(25,*) i, j, size(bin(i, j)%p), g%rlo(i)/AU, g%rup(i)/AU, ncolls(i,j)
               enddo
            enddo
            close(25)

#ifdef LOCAL_SIM
            open(26,file=trim(output_path)//trim('feed_leak.dat'),status='unknown',position='append')          
            write(26,*) time/year, Mfeed, Mfeed_ri, Mleak, Mleak_ri, Mplt, Mplt_ri
            close(26)
#endif
      endif

      iter = iter + 1


#ifdef ZERO_D
#else

if (stabilize .le. 1) then
      write(*,*) ' Performing advection: timestep',resdt/year,'yrs'
      ! performing advection
      call mc_advection(sim_swrm, resdt, time &
#ifdef GLOBAL
                        &, Nflux &
#endif
                        & )
      write(*,*) '  advection done'
endif

#endif      

      ! removing old grid and building new one
      call deallocate_grid
      if (allocated(bin))  deallocate(bin)
      if (allocated(rbin)) deallocate(rbin)
      if (allocated(ncolls)) deallocate(ncolls)
      write(*,*) '    Making grid...'
      call make_grid(sim_swrm, bin, rbin, nr, nz, smallr,totmass, ncolls)
      write(*,*) '     grid done'

      ! performing collisions
      write(*,*) '   Performing collisions...'


      !$OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(DYNAMIC)
      do i = 1, size(g%rce)
         do j = 1, size(g%zce,dim=2)
            if (.not.allocated(bin(i, j)%p)) cycle
            call mc_collisions(i, j, bin, sim_swrm, resdt, time, ncolls(i,j), bin(i, 1)%first_idx)
         enddo
      enddo
      !$OMP END PARALLEL DO


    
      write(*,*) '    collisions done!'

      ! this can also be paralelized, but not so relevant right now
      do i = 1, size(g%rce)
         nparts = 0
         do j = 1, size(bin(:,:),dim=2)  
            ! how many particles do I have above the St = stmin limit? and bellow 0.1
            do w = 1, size(bin(i,j)%p(:))
               if ((bin(i,j)%p(w)%stnr>stmin) .and. (bin(i,j)%p(w)%stnr<100.*stmin)) nparts = nparts + 1
            enddo
         enddo
         etamax = nparts*mswarm/pi/(g%rup(i)**2 - g%rlo(i)**2)/sigmag(g%rce(i),time) 

         if (modulo(iter,fout) == 0 .or. time+resdt>=timeofnextout) then
            open(26,file=trim(output_path)//trim('metallicity_integrated.dat'),status='unknown',position='append')
            write(26,*) time/year, g%rce(i)/AU, etamax, g%rlo(i)/AU, g%rup(i)/AU
            close(26)
         endif
         

#ifdef NO_SI
         write(*,*) "planetesimal formation deactivated" 
#else
         if (etamax > Zcrit) then 
            write(*,*) ' Planetesimal formation!'
            call  planetesimal_formation(i, rbin, sim_swrm, resdt, time)
            plts_tot = plts_tot + plts_mass
         endif
#endif
      enddo


#ifdef LOCAL_SIM

   if (stabilize .le. 1) then

      ! how many new plt.
      Nplt_new = count(sim_swrm(:)%plt .eq. 1) ! this ones are always inside smallr
      if (Nplt_new>0) then ! create new array temp_swrm
         allocate(temp_arr(Nplt_tot))
         temp_arr(:) = temp_swrm(:Nplt_tot)
         if (count(temp_arr(:Nplt_tot)%plt .eq. 1) .ne. Nplt_tot) stop
         deallocate(temp_swrm)
         allocate(temp_swrm(Nplt_tot+Nplt_new+Nmax))
         temp_swrm(:Nplt_tot) = temp_arr(:)
         deallocate(temp_arr)
      endif


      ! how many stay
      Nin = count(sim_swrm(:)%rdis > smallr) 

      ! how many leak
      Nleak = Nmax - Nin - Nplt_new

      ! how many feed from feeding file (not to local, that is calculated after mswarm correction)
      Ncount = 0
      allocate(temp(size(feeding_swrm_tot)))
      do i = 1, size(feeding_swrm_tot)
         if ((feeding_swrm_tot(i)%tfor>time) .and. (feeding_swrm_tot(i)%tfor<=time+resdt)) then
            Ncount = Ncount + 1
            temp(Ncount) = feeding_swrm_tot(i)
         endif
      enddo
      Nfluxtot = Nfluxtot - Ncount ! count how many left before reading
      if (Ncount>0) then 
         allocate(feeding_swrm(Ncount))
         feeding_swrm(:) = temp(:Ncount)
         deallocate(temp)
      else 
         deallocate(temp)
      endif

      if (Nfluxtot .eq. 0) then ! read again
         deallocate(feeding_swrm_tot)
         write(filename,'(a7,i5.5,a3)') 'swarms-',nout_feed,'.h5'
         open_file = trim(feedingdir)//trim(filename)
         call hdf5_file_read(open_file, Nfluxtot, feeding_swrm_tot, nout_feed, mswarm_feeding, time_no, time_no, imax, .True.)
         write(*,*) 'Nluxtot ', Nfluxtot
         nout_feed = nout_feed + 1
         Nlagun = 0
         allocate(temp(size(feeding_swrm_tot)))
         do i = 1, size(feeding_swrm_tot)
#ifdef TEST_CONST_FLUX

            feeding_swrm_tot(i)%f_CH = prop_const
            feeding_swrm_tot(i)%f_matrix = 1.-prop_const
   
   
            call random_number(x)
            if (x<prop_const) then
               feeding_swrm_tot(i)%rigid = 1
               if (prop_const*feeding_swrm_tot(i)%mass < feeding_swrm_tot(i)%rigid_mass) then
                  feeding_swrm_tot(i)%mass = feeding_swrm_tot(i)%rigid_mass
                  feeding_swrm_tot(i)%f_matrix = 0.
                  feeding_swrm_tot(i)%f_CH = 1.
               endif
            else
               feeding_swrm_tot(i)%rigid = 0
               if (prop_const*feeding_swrm_tot(i)%mass < feeding_swrm_tot(i)%rigid_mass) then
                  feeding_swrm_tot(i)%f_matrix = 1.
                  feeding_swrm_tot(i)%f_CH = 0.
               endif
            endif
            feeding_swrm_tot(i)%dens = 1./(feeding_swrm_tot(i)%f_CH/ridens + (1.-feeding_swrm_tot(i)%f_matrix)/matdens)

#endif

            if ((feeding_swrm_tot(i)%tfor>time) .and. (feeding_swrm_tot(i)%tfor<=time + resdt)) then
               Nlagun = Nlagun + 1
               temp(Nlagun) = feeding_swrm_tot(i)
            endif
         enddo
         Nfluxtot = Nfluxtot - Nlagun! restart counter
         if (Nlagun>0) then
            Ncount = Ncount + Nlagun ! add to previous Ncount
            allocate(feeding_lagun(size(feeding_swrm)))

            feeding_lagun = feeding_swrm
            deallocate(feeding_swrm)
            allocate(feeding_swrm(Ncount))
            feeding_swrm(:size(feeding_lagun)) = feeding_lagun
            feeding_swrm(size(feeding_lagun)+1:) = temp(:Nlagun)
            deallocate(temp, feeding_lagun)
         else 
            deallocate(temp)
         endif
      endif

      

      mass_add = real(Ncount)*mswarm_feeding

      mswarm_old = mswarm
      mswarm = (real(Nin)*mswarm + mass_add)/real(Nmax)
      Nfeed = int(mass_add/mswarm) 
      call random_number(x)
      if (x < (mass_add/mswarm-real(Nfeed))) Nfeed = Nfeed + 1  !redistribute remanent

      ! Add (Nrem>0) or remove (Nrem<0) particles
      Nrem =  Nmax - Nfeed - Nin
     
      write(*,*) 'Nrem is ', Nrem, ' and Nin is ', Nin, ' and Nplt_new is is ', Nplt_new
      write(*,*) 'Nfeed is ', Nfeed, ' and Nleak', Nleak

      Nleak_ri = 0
      if (Nrem>0) then ! add extra particles
         j = 0 
         jadd = 0
         i_plt = 0
         allocate(temp_arr(Nmax)) ! because we do j+j_add, we dont wanna take any particle that has been modified
         temp_arr = sim_swrm
         do i=1, Nmax
            if (temp_arr(i)%rdis > smallr) then ! bring particles that will still be simulated to front
               temp_arr(i)%npar = mswarm / temp_arr(i)%mass ! this needs to be updated after changing mswarm
               prob = real(Nrem-jadd)/real(Nin-j)
               j = j + 1
               sim_swrm(j+jadd) = temp_arr(i)
               temp_swrm(Nplt_tot+Nplt_new+j+jadd) = temp_arr(i)
               call random_number(x)
               if (x<prob) then
                  Ntot = Ntot + 1 
                  jadd = jadd + 1
                  sim_swrm(j+jadd) = temp_arr(i)
                  sim_swrm(j+jadd)%idnr = Ntot
                  temp_swrm(Nplt_tot+Nplt_new+j+jadd) = temp_arr(i)
                  temp_swrm(Nplt_tot+Nplt_new+j+jadd)%idnr = Ntot
                  
               endif
            else if (temp_arr(i)%plt .eq. 1) then
               i_plt = i_plt + 1
               temp_swrm(Nplt_tot+i_plt) = temp_arr(i)
               Mplt = Mplt + mswarm
               Mplt_ri = Mplt_ri + mswarm*temp_arr(i)%rigid
               
            else
               Nleak_ri = Nleak_ri + temp_arr(i)%rigid
            endif
         enddo
         write(*,*) 'j and jadd ', j, jadd, ' while Nin and Nrem', Nin, Nrem
         j = j + jadd
         deallocate(temp_arr)
      else if (Nrem<0) then ! remove extra particles
         j = 0
         jadd = 0
         i_plt = 0
         allocate(temp_arr(Nmax))
         temp_arr = sim_swrm
         do i=1, Nmax
            if (temp_arr(i)%rdis > smallr) then ! bring particles that will still be simulated to front
               prob = real(-Nrem-jadd)/real(Nin-j)
               call random_number(x)
               if (x<prob) then
                  jadd = jadd + 1
                  j = j + 1
               else 
                  temp_arr(i)%npar = mswarm / temp_arr(i)%mass ! this needs to be updated after changing mswarm
                  j = j + 1
                  sim_swrm(j-jadd) = temp_arr(i)
                  temp_swrm(Nplt_tot+Nplt_new+j-jadd) = temp_arr(i)
               endif
            else if (temp_arr(i)%plt .eq. 1) then
               i_plt = i_plt + 1
               temp_swrm(Nplt_tot+i_plt) = temp_arr(i)
               Mplt = Mplt + mswarm
               Mplt_ri = Mplt_ri + mswarm*temp_arr(i)%rigid
            else
               Nleak_ri = Nleak_ri + temp_arr(i)%rigid
            endif
         enddo
         write(*,*) 'j and jrem', j, jadd, ' while Nin and Nrem', Nin, Nrem
         j = j - jadd
         deallocate(temp_arr)

      else
         j = 0
         i_plt = 0
         allocate(temp_arr(Nmax))
         temp_arr = sim_swrm
         do i=1, Nmax
            if (temp_arr(i)%rdis > smallr) then ! bring particles that will still be simulated to front
               temp_arr(i)%npar = mswarm / temp_arr(i)%mass ! this needs to be updated after changing mswarm
               j = j + 1
               sim_swrm(j) = temp_arr(i)
               temp_swrm(Nplt_tot+Nplt_new+j) = temp_arr(i)
            else if (temp_arr(i)%plt .eq. 1) then
               i_plt = i_plt + 1
               temp_swrm(Nplt_tot+i_plt) = temp_arr(i)
               Mplt = Mplt + mswarm
               Mplt_ri = Mplt_ri + mswarm*temp_arr(i)%rigid
            else
               Nleak_ri = Nleak_ri + temp_arr(i)%rigid
            endif
         enddo
         write(*,*) 'j ', j, ' while Nin and Nrem', Nin, Nrem
         deallocate(temp_arr)
      endif
      
      Nfeed_ri = 0
      ! 2) feed them
      if (Nfeed .ne. Nmax - j) stop
      Nfeed = Nmax - j ! in case not exactly Nrem were added/removed
      write(*,*) 'Nfeed at the end is ', Nfeed

      if (Nfeed>0) then

         do i=1, Nfeed
            
            ! randomly choose with what particle to initialize
            call random_number(x)
            idx = int(x*size(feeding_swrm)) + 1
            sim_swrm(j+i) = feeding_swrm(idx)
            sim_swrm(j+i)%idnr = Ntot+i
            sim_swrm(j+i)%plt = 0
            sim_swrm(j+i)%f_matrix = 1.- sim_swrm(j+i)%f_CH
            sim_swrm(j+i)%dens = 1./(sim_swrm(j+i)%f_CH/ridens + sim_swrm(j+i)%f_matrix/matdens)
            sim_swrm(j+i)%rdis = maxrad0
            sim_swrm(j+i)%npar = mswarm / sim_swrm(j+i)%mass
            sim_swrm(j+i)%velr = 1.e-20
            sim_swrm(j+i)%velz = 1.e-20
            sim_swrm(j+i)%coll_f = 0
            Nfeed_ri = Nfeed_ri + sim_swrm(j+i)%rigid

            call mc_advection_element(sim_swrm(j+i), time+resdt-sim_swrm(j+i)%tfor, sim_swrm(j+i)%tfor)

            temp_swrm(Nplt_tot+Nplt_new+j+i) = sim_swrm(j+i)
         enddo


         Ntot = Ntot + Nfeed
      
      endif


      ! this is not calculated from Nrem correction
      Mfeed = Mfeed + mswarm*Nfeed
      Mfeed_ri = Mfeed_ri + mswarm*Nfeed_ri
      Mleak= Mleak + mswarm_old*Nleak ! this should be with old leak
      Mleak_ri = Mleak_ri + mswarm_old*Nleak_ri

      Nplt_tot = Nplt_tot + Nplt_new
      imax = Nplt_tot + Nmax


      if (Ncount>0) deallocate(feeding_swrm)
   endif
#endif

      if (stabilize > 0) then
         stabilize = stabilize - 1
      endif
      if (stabilize == 0) then 
         time = time + resdt
      endif

      write(*,*) '---------------------------------------'
   enddo
   ! ---- END OF THE MAIN LOOP -----------------------------------------------------------------------------------------

   write(*,*) 'time: ', time/year, 'produced output: ',nout
   open(23,file=trim(output_path)//trim('timesout.dat'),status='unknown',position='append')
   write(23,*) 'time: ', time/year, 'produced output: ',nout
   close(23)
#ifdef LOCAL_SIM
   call hdf5_file_write(file, datadir, temp_swrm, time, resdt, 'create', nout, mswarm, Ntot, imax)
#else
   call hdf5_file_write(file, datadir, sim_swrm, time, resdt, 'create', nout, mswarm, Ntot)
#endif  
   call end_struct()

   deallocate(bin)
   deallocate(rbin)


   ! calclate time of execution
   total = etime(elapsed)
   write(*,*) 'Elapsed time [s]: ', total, ' user:', elapsed(1), ' system:', elapsed(2)

end
