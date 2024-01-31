! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!PETER
!! This routine evaluates eddy diffusion coefficients,
!! ekz(k) [cm^2 s^-1].
!!
!! Method: First ekz is modeled after the eddy diffusion coefficient profile 
!! from pg. 479 of Hunten, D. M. et al. Venus. 1983. Second ekz is a log-linear
!! extrapolation of the eddy diffusion coefficient profile from Imamura and 
!! Hashimoto 2001. Both are applicable only to the heights between 70 and 100 km.
!! The third ekz is modeled after the eddy diffusion coefficient profile from 
!! Imamura and Hashimoto 2001, and is applicable only to the heights between 40 
!! and 70 km. The fourth ekz is an attempt to add to an extrapolation of the 
!! third ekz a wide gaussian so as to partially-reproduce the shape of the first
!! ekz in the range between 70 and 100km; as such, this ekz is applicable from 
!! 40 to 100km. 
!!
!!
!!
!! @author Peter Gao
!! @version Oct-2011
subroutine setupedif(carma, cstate, rc)

  !     types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_planet_mod
  use carma_condensate_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)         :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations  
  integer                 :: igroup, ibin, iz, igas, k1, k2, nzm1

  !  Define formats
  2 format(/,'Eddy diffusion coefficient (prior to interpolation)')
  3 format(/,'Particle group ',i3,/,' bin   lev  p [dyne/cm2]  T [K]       r [cm]    wet r [cm]    ekz [cm2/s]',/)
  4 format(i3,4x,i3,5(1pe11.3,4x)) 

  ! Loop over all altitudes - now done in the test cases themselves as an input parameter. 
!  do iz = 1, NZ

    ! Vertical eddy diffusion coefficient
!    ekz(iz) = 1000000._f - 990000._f * exp( - ((iz * 100._f - 5000._f) / (12011.224_f)) ** 2) 
!    ekz(iz,ibin,igroup,igas) = 10._f ** (((iz + 300) * 100._f)/38552.9_f)
!    ekz(iz) = 10000._f * 10._f ** ((iz * 100._f)/38552.9_f) + 2500000._f * exp( - ((iz * 100._f - 12500._f) / (1201.1224_f)) ** 2) 
!    ekz(iz) = 10000._f * 10._f ** ((iz * 200._f)/38552.9_f) + 2500000._f * ( exp( - ((iz * 200._f - 12500._f) / (1201.1224_f)) ** 2) + &	!PETER
!              exp( - ((iz * 200._f - 60000._f) / (12011.224_f)) ** 2) )										!PETER

!  enddo

    
  ! Print out diffusivities.
#ifdef DEBUG
  if (do_print_init) then
   
    write(LUNOPRT,2)

    do igroup = 1, NGROUP
        
      write(LUNOPRT,3) igroup
      
      do ibin = 1,NBIN
    
        do iz = NZ, 1, -1
          write(LUNOPRT,4) ibin,iz,p(iz),t(iz),r(ibin,igroup),r_wet(iz,ibin,igroup),ekz(iz)
        end do
      enddo
    enddo
  
    write(LUNOPRT,*) ""
  endif
#endif

  !  Interpolate <ekz> from layer mid-pts to layer boundaries.
  !  <ekz(k)> is the diffusion coefficient at the lower edge of the layer
  nzm1 = max(1, NZ-1)

  ! Set upper boundary before averaging
  ekz(NZP1) = ekz(NZ)

  if (NZ .gt. 1) then
    ekz(NZ) = sqrt(ekz(nzm1) * ekz(NZ))

    if (NZ .gt. 2) then
      do iz = NZ-1, 2, -1
        ekz(iz) = sqrt(ekz(iz-1) * ekz(iz))
      enddo
    endif
  endif
  
  ! Scale cartesian diffusivities to the appropriate vertical coordinate system.
  ! Non--cartesion coordinates are assumed to be positive downward, but
  ! vertical velocities in this model are always assumed to be positive upward. 
  if( igridv .ne. I_CART )then
    ekz(:) = ekz(:) / (zmetl(:)**2)
  endif

  ! Return to caller with fall velocities evaluated.
  return
end
