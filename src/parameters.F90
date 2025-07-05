! this module reads parameters from params.par file and will create a file for output
module parameters
    
   use constants,                only: year, AU, pi, third

   implicit none

   private
   public :: read_parameters, Ntot, ncell, nr, nz, fout, dtime, tend, smallr, restart, restime, minrad0, maxrad0, &
             r0, matdens, ridens, dmmax, dtg, vfrag, con1, con2, db_data, &
             alphat, nbins, prop_CH, aCH_min, aCH_max, zeta, &
             stmin, plt_eff, disk_path, output_path, sim_path, datadir, nzone_min
#ifdef GLOBAL
   public :: feedingdir
#endif
#ifdef LOCAL_SIM
   public :: feedingdir, nzone_max, Nmax
#endif

   integer                             :: Ntot     ! total number of representative particles in the simulation
   integer                             :: ncell    ! number of particles per cell
   integer                             :: nr       ! nr of radial zones
   integer                             :: nz       ! nr of zones in z
   real                                :: dtime    ! do output at least every dtime time
   real                                :: tend     ! time to finish the simulation
   real                                :: smallr   ! smallest distance (evaporation radius)
   real                                :: dtg      ! dust to gas ratio
   logical                             :: restart  ! is the simulation run from restart?
   real                                :: restime  ! time of restart
   real                                :: minrad0  ! minimum distance from the star of initial rim
   real                                :: maxrad0  ! maximum distance from the star of initial rim
   real                                :: r0       ! monomer radius
   real, protected                     :: matdens  ! material density of matrix
   real, protected                     :: ridens   ! material density of refractory inclusions (CH, CAI, AOA)
   real                                :: dmmax    ! MC-acceleration parameter
   real                                :: vfrag    ! fragmentation treshold velocity
   integer                             :: nbins    ! number of bins for mass histograms
   integer                             :: fout     ! steps between outputs
   real                                :: alphat   ! gas disk properties
   real                                :: prop_CH  ! how much percentage of the total mass will be rigid
   real                                :: aCH_min  ! minimum size of rigid particles
   real                                :: aCH_max, zeta ! maximum size of rigid particle and distribution index
   real                                :: stmin    ! minimum St to trigger SI
   real                                :: plt_eff  ! planetesimal formation efficiency
   character(len = 200)                :: disk_path = ".disk/"
   character(len = 200)                :: sim_path = "default/", output_path
   logical                             :: db_data  ! write out files for sanity checks?
   character(len=100)                  :: datadir  ! data directory to write the data into
   real, protected                     :: con1  ! optimization
   real, protected                     :: con2  ! optimization
   integer                             :: nzone_min ! for grid
#ifdef GLOBAL
   character(len=100)                  :: feedingdir ! where do we storage data for local simulation later?
#endif
#ifdef LOCAL_SIM
   character(len=100)                  :: feedingdir
   integer :: nzone_max       ! maximum number of particles per grid
   integer :: Nmax            ! constant number of particles throughout the simulation
#endif


   contains

   subroutine read_parameters(ctrl_file)
   
      implicit none
      character(len=100), intent(in)   :: ctrl_file
      character(len=100)               :: buffer, label
      integer                          :: pos
      integer                          :: ios = 0
      integer                          :: line = 0
      integer, parameter               :: fh = 1
      character(len=100)               :: command

      ! setting default values
      minrad0 = 3.0
      maxrad0 = 5.0
      r0 = 1.e-2
      ncell = 200
      nr   = 10
      nz = 1
      fout = 5
      dtime = 100. * year
      tend = 1000. * year
      matdens = 1.2
      ridens = 3.3 ! Helmann+ 2020
      dmmax = 0.001
      smallr = 0.1 !3.0
      dtg = 0.01
      vfrag = 100.
      alphat = 1.e-3
      restart = .false.
      restime = 0.0
      nbins = 200
      output_path = trim('../../outputs/')//trim(sim_path)
      db_data = .true.
      datadir = "data"
#ifdef GLOBAL
      feedingdir = "feeding_data"
#endif
#ifdef LOCAL_SIM
      feedingdir = "feeding_data"
#endif

      prop_CH = 0.
      aCH_min = 0.001
      aCH_max = 0.1

      stmin = 0.01
      plt_eff = 0.001

      nzone_min = 200 ! From Drazkowska et al. 2013 200
#ifdef LOCAL_SIM
      nzone_max = 250
      Nmax = 13000
#endif
      ! reading the parameter file
      open(fh,file=ctrl_file,action='read')

      ! ios is negative if an end of record condition is encountered or if
      ! an endfile condition was detected.  It is positive if an error was
      ! detected.  ios is zero otherwise.

      do while (ios == 0)
         read(fh, '(A)', iostat=ios) buffer
         if (ios == 0) then
               line = line + 1

               pos = scan(buffer, ' 	')
               label = buffer(1:pos)
               buffer = buffer(pos+1:)

               select case (label)
               case ('minimum_radius_[AU]')
                  read(buffer, *, iostat=ios) minrad0
                  print *, 'Read minimum radius: ', minrad0
                  minrad0 = minrad0*AU
               case ('maximum_radius_[AU]')
                  read(buffer, *, iostat=ios) maxrad0
                  print *, 'Read maximum radius: ', maxrad0
                  maxrad0 = maxrad0*AU
               case ('monomer_radius_[cm]')
                  read(buffer, *, iostat=ios) r0
                  print *, 'Read monomer radius: ', r0
               case('number_of_particles_per_cell')
                  read(buffer, *, iostat=ios) ncell
                  print *, 'Read number of particles per cell: ', ncell
               case('nzone_min')
                  read(buffer, *, iostat=ios) nzone_min
                  print *, 'Read min number of particles per cell: ', nzone_min
#ifdef LOCAL_SIM
               case('nzone_max')
                  read(buffer, *, iostat=ios) nzone_max
                  print *, 'Read max number of particles per cell: ', nzone_max
               case('Nmax_local')
                  read(buffer, *, iostat=ios) Nmax
                  print *, 'Maximum numbers of particles in local sim: ', Nmax
#endif
               case('number_of_radial_zones')
                  read(buffer, *, iostat=ios) nr
                  print *, 'Read number of radial zones: ', nr
               case('number_of_vertical_zones')
                  read(buffer, *, iostat=ios) nz
                  print *, 'Read number of vertical zones: ', nz
               case('time_between_outputs_[yrs]')
                  read(buffer, *, iostat=ios) dtime
                  print *, 'Read time between outputs: ', dtime
                  dtime = dtime * year
               case('steps_between_outputs')
                  read(buffer, *, iostat=ios) fout
                  print *, 'Read number of iterations between outputs: ', fout
               case('maximum_time_of_simulation_[yrs]')
                  read(buffer, *, iostat=ios) tend
                  print *, 'Read maximum time of simulation: ', tend
                  tend = tend * year
               case('material_density_[g/cm3]')
                  read(buffer, *, iostat=ios) matdens
                  print *, 'Read material density: ', matdens
               case('refractory_density_[g/cm3]')
                  read(buffer, *, iostat=ios) ridens
                  print *, 'Read refractory density: ', ridens
               case('dmmax')
                  read(buffer, *, iostat=ios) dmmax
                  print *, 'Read dmmax: ', dmmax
               case('evaporation_radius_[AU]')
                  read(buffer, *, iostat=ios) smallr
                  print *, 'Read evaporation radius: ', smallr
                  smallr = smallr * AU
               case('dust_to_gas_ratio')
                  read(buffer, *, iostat=ios) dtg
                  print *, 'Read dust to gas ratio: ', dtg
               case('fragmentation_velocity_[cm/s]')
                  read(buffer, *, iostat=ios) vfrag
                  print *, 'Read fragmentation velocity: ', vfrag
               case('restart')
                  read(buffer, *, iostat=ios) restart
                  print *, 'Restart? ',restart
               case('Restart_time_[yrs]')
                  read(buffer, *, iostat=ios) restime
                  print *, 'Read restart time: ', restime
                  restime = restime * year
               case('alphat')
                  read(buffer, *, iostat=ios) alphat
                  print *, 'alphat ', alphat
               case('plt_eff')
                  read(buffer, *, iostat=ios) plt_eff
                  print *, 'plt_eff', plt_eff             
               case('CH_abundance')
                  read(buffer, *,iostat=ios) prop_CH
                  print *, 'CH?', prop_CH
               case('CH_size_min')
                  read(buffer, *,iostat=ios) aCH_min
                  print *, 'CH min size?', aCH_min
               case('CH_size_max')
                  read(buffer, *,iostat=ios) aCH_max
                  print *, 'CH max size?', aCH_max
               case('zeta')
                  read(buffer, *,iostat=ios) zeta
                  print *, 'zeta power?', zeta
               case('disk_path')
                  read(buffer, *,iostat=ios) disk_path
                  print *, 'disk_path', disk_path
               case('init_path')
                  read(buffer, *,iostat=ios) sim_path
               case('out_path')
                  read(buffer, *,iostat=ios) output_path
                  output_path = trim('../../outputs/')//trim(output_path)
               case('generate_debug_files?')
                  read(buffer, *,iostat=ios) db_data
                  print *, 'generate files for sanity check?', db_data
               case('data_directory')
                  read(buffer, *,iostat=ios) datadir
                  datadir = trim(output_path)//trim(datadir)
                  print *, 'data output directory ', datadir
#ifdef LOCAL_SIM
               case('feeding_directory')
                  read(buffer, *,iostat=ios) feedingdir
                  feedingdir = trim('../../outputs/')//trim(feedingdir)//trim('feeding_data/')
                  print *, 'feeding output directory ', feedingdir
#endif
               ! if you want other parameters add another cases here and in the setup file

               case default
                  print *, 'Skipping invalid label at line', line
               end select
         end if
      end do
      Ntot = ncell * nr * nz
      con1 = pi**third * (0.75 / matdens)**(2. * third) ! due to ridens != matdens, this is not const anymore 
      con2 = (0.75 / pi / matdens)**third

      ! create files for output if they do not exist yet
      close(fh)
      write(command,*) './directory.sh '//trim(output_path)
      CALL SYSTEM(command)
      write(command,*) './directory.sh '//trim(datadir)
      CALL SYSTEM(command)
#ifdef GLOBAL
      feedingdir = trim(output_path)//trim(feedingdir)
      write(command,*) './directory.sh '//trim(feedingdir)
      CALL SYSTEM(command)
#endif

   end subroutine read_parameters

end module parameters
