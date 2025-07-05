! this module converts pebbles into planetesimals by setting plt=1 and frezing them
! by changing their radial distance inside the evaporation radius. 

module ppplts
   use constants
   use discstruct
   use types
   use parameters, only: smallr, stmin, plt_eff, output_path
   use initproblem, only: mswarm
   use grid, only: g
   
   implicit none
   private
   public :: planetesimal_formation, plts_mass
   
   real :: plts_mass
   
   contains
   
   subroutine planetesimal_formation(nr, rbin, swrm, dtime, realtime)
      implicit none
      integer, intent(in)                                             :: nr
      type(swarm), dimension(:), allocatable, target                  :: swrm
      type(list_of_swarms), dimension(:),   allocatable               :: rbin
      real, intent(in)                                                :: dtime, realtime
      real                                                            :: Ztot
      integer                                                         :: nparts, nparts_tot, nremv ! number of swarms above the threshold and a total number of swarms
      integer                                                         :: j
      real                                                            :: nremv_float  
      real                                                            :: x ! random variable

      plts_mass = 0.
     
      nparts_tot = size(rbin(nr)%p(:)) 

      ! my total mettalicity; surface dust-to-gas ratio
      Ztot = real(nparts_tot) * mswarm / (pi*g%rup(nr)**2-pi*g%rlo(nr)**2) / sigmag(g%rce(nr),realtime)

      ! how many particles do I have above the St = stmin limit?
      nparts = 0
      do j = 1, size(rbin(nr)%p(:))
         if ((rbin(nr)%p(j)%stnr > stmin) .and. (rbin(nr)%p(j)%stnr < 1000.*stmin)) nparts = nparts + 1
      enddo

      ! estimate number of particles to remove, following Eq. 16 in Drazkowska anbd Dullemond 2016
      nremv_float = dtime*plt_eff/(2.*pi)*omegaK(g%rce(nr), realtime)*real(nparts)

      ! choose how to round the decimal part.
      ! this part of code is specially relevant when dtime or plt_eff are small
      ! if we simply use floor we underestimate planetesimal formation
      nremv = floor(nremv_float)
      call random_number(x)
      if (x<(nremv_float - nremv)) then
         nremv = nremv + 1
      endif

      write(*,*) "We remove ", nremv, " particles"

      ! remove pebbles
      call form_planetesimals(swrm, rbin, nr, nparts, nparts_tot, nremv)
   
      return
   end subroutine planetesimal_formation


   subroutine form_planetesimals(swrm, rbin, nr, nparts, nparts_tot, nremv)
      implicit none
      type(swarm), dimension(:), allocatable                 :: swrm
      type(list_of_swarms), dimension(:),   allocatable      :: rbin
      integer                                                :: nparts, nparts_tot !not sure if it is relevant to update this
      real                                                   :: mclump
      integer                                                :: nremv
      integer                                                :: i, k
      integer, intent(in)                                    :: nr
      real:: rand
      integer:: remove, j
      mclump = 0.0
      
      ! collect particles for 1 clump
      do while (nremv > 0)
         i = 1
         j = 1
         call random_number(rand)
         remove = floor(rand*nparts) + 1
         do while (j .ne. remove)
            i = i + 1
            ! pebbles that can undergo SI have their St between stmin and 100*stmin
            if ((rbin(nr)%p(i)%stnr > stmin) .and. (rbin(nr)%p(i)%stnr < 100.*stmin)) j = j + 1 
         enddo
         mclump = mclump + mswarm
         rbin(nr)%p(i)%rdis = 0.5 * smallr ! remove the particle from simulation (move beyond the evaporation radius)
         rbin(nr)%p(i)%stnr = 0.00
         k = rbin(nr)%first_idx
         do while (swrm(k)%idnr .ne. rbin(nr)%p(i)%idnr)
            k = k + 1
         enddo
         swrm(k)%rdis = 0.5 * smallr 
         swrm(k)%plt = 1
         nparts = nparts - 1
         nparts_tot = nparts_tot-1
         nremv = nremv-1
      enddo
      
      plts_mass = plts_mass + mclump
      
      write(*,*) 'Clump removed', mclump, nremv, real(nremv)*mswarm
      
      return
   end subroutine form_planetesimals

end module
