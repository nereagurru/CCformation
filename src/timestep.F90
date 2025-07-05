! this module estimates maximum time step of the code taking into account Courant-like condition for advection
module timestep

   use grid,       only: g
   use types,      only: list_of_swarms

   private
   public      :: time_step

   contains

   ! gives back timestep as minimal timestep of all processes
   subroutine time_step(bin, ncolls, dtout, final_dtime)
      implicit none
      type(list_of_swarms), dimension(:,:), allocatable, target :: bin
      integer, dimension(:,:), allocatable           :: ncolls
      real, intent(in)                               :: dtout
      real, intent(inout)                            :: final_dtime   ! final timestep
      real                                           :: old_dt
      integer                                        :: i, j, c2=1, c1=10 !c1,c2 - timestep enlargement factors
      real:: eps_time = 1. !in secs

      old_dt = final_dtime
      if (dtout>eps_time) final_dtime = dtout ! avoid dtout == 0 as the timestep drops to 0
      
      do i = 1, size(ncolls(:,:),dim=1)
         do j = 1, size(ncolls(:,:),dim=2)
            final_dtime = min(final_dtime, c1*size(bin(i,j)%p) * old_dt / real(max(1,ncolls(i,j)))) !we don't want more collisions than ~ nr of particles in a zone per timestep
#ifdef ZERO_D
#else
            final_dtime = min(final_dtime, c2*g%dz(i,j)/maxval(abs(bin(i,j)%p(:)%velz))) ! don't let particles jump over cells in z
            final_dtime = min(final_dtime, c2*g%dr(i)/maxval(abs(bin(i,j)%p(:)%velr))) ! don't let particles jump over cells in r
#endif
         enddo
      enddo
   end subroutine time_step


end
