! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine drives the potentially slower microphysics calculations.
!!
!!  Originally part of microphy.  Now in this separate routine to allow
!!  time splitting of coagulation at a different timestep size from
!!  other microphysical calcs.
!!
!! @author McKie
!! @version Sep-1997
subroutine microslow(carma, cstate, rc)

  ! carma types defs
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

  ! Local Declarations
  integer   :: i		
  integer   :: iz		
  integer   :: ibin
  integer   :: ielem
  integer   :: elemmultiple	

  if (do_printdiag) then										
    write(lundiag,*) 'PART 4: COAGULATION'								
    write(lundiag,*) '****************************************'						 
    elemmultiple = int(NELEM / 3._f)									
    coagprod(:,:,:) = 0._f										
    coagloss(:,:,:) = 0._f										
  end if												



  ! Calculate (implicit) particle loss rates for coagulation.
  call coagl(carma, cstate, rc)
  ! Calculate particle production terms and solve for particle 
  ! concentrations at end of time step.
  !
  ! NOTE: The order of elements required by CARMA to work with the
  ! element loop first is: if you have a group that is both a source
  ! and product of coagulation, then it needs to come after the
  ! other group that participates in that coagulation in the element
  ! table. For example, icoag(2,1) = 1 will not work, but
  ! icoag(2,1) = 2 should work.
  do ielem = 1,NELEM
    do ibin = 1,NBIN
      call coagp(carma, cstate, ibin, ielem, rc)
      call csolve(carma, cstate, ibin, ielem, rc)
    enddo
  enddo

!  write(*,*) 'COAG'
!  write(*,*) coagprod(33,33,2) - coagloss(33,33,2) * pc(33,33,2)

  if (do_printdiag) then
  do i = 1, elemmultiple
    write(lundiag,'(2A5,$)') 'Z', 'BIN'
    do ielem = 3*(i-1)+1,3*i-1
      write(lundiag,'(A12,i2,A12,i2,A12,i2,$)') 'COAGPROD',ielem,'COAGLOSS',ielem,'NETFLUX',ielem
    end do
    write(lundiag,'(A12,i2,A12,i2,A12,i2)') 'COAGPROD',3*i,'COAGLOSS',3*i,'NETFLUX',3*i
    do iz = 1, NZ
      do ibin = 1, NBIN
        write(lundiag,'(2i5,$)') iz, ibin
        do ielem = 3*(i-1)+1,3*i-1
          write(lundiag,'(3e14.3,$)') coagprod(iz,ibin,ielem), pc(iz,ibin,ielem) * coagloss(iz,ibin,igelem(ielem)), &
            coagprod(iz,ibin,ielem)-pc(iz,ibin,ielem) * coagloss(iz,ibin,igelem(ielem))
        end do
        write(lundiag,'(3e14.3)') coagprod(iz,ibin,3*i), pc(iz,ibin,3*i) * coagloss(iz,ibin,igelem(3*i)), &
          coagprod(iz,ibin,3*i)-pc(iz,ibin,3*i) * coagloss(iz,ibin,igelem(3*i))
      end do
    end do
    write(lundiag,*) ' ' 
  end do
  if (3*elemmultiple+1 .le. NELEM) then
    write(lundiag,'(2A5,$)') 'Z', 'BIN'
    do ielem = 3*elemmultiple+1,NELEM-1
      write(lundiag,'(A12,i2,A12,i2,A12,i2,$)') 'COAGPROD',ielem,'COAGLOSS',ielem,'NETFLUX',ielem
    end do
    write(lundiag,'(A12,i2,A12,i2,A12,i2)') 'COAGPROD',NELEM,'COAGLOSS',NELEM,'NETFLUX',NELEM
    do iz = 1, NZ
      do ibin = 1, NBIN
        write(lundiag,'(2i5,$)') iz, ibin
        do ielem = 3*elemmultiple+1,NELEM-1
          write(lundiag,'(3e14.3,$)') coagprod(iz,ibin,ielem), pc(iz,ibin,ielem) * coagloss(iz,ibin,igelem(ielem)), &
            coagprod(iz,ibin,ielem)-pc(iz,ibin,ielem) * coagloss(iz,ibin,igelem(ielem))
        end do
        write(lundiag,'(3e14.3)') coagprod(iz,ibin,NELEM), pc(iz,ibin,NELEM) * coagloss(iz,ibin,igelem(NELEM)), &
          coagprod(iz,ibin,NELEM)-pc(iz,ibin,NELEM) * coagloss(iz,ibin,igelem(NELEM))
      end do
    end do
  end if
  end if

  if (do_printdiag) write(lundiag,*) ' '
  ! Return to caller with new particle concentrations.
  return
end
