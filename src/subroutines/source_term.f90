!===============================================================================
!> \file source_term.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is the subroutine that control the source terms in the equations.
!! \details
!! Contains source_term(), gravity_predictor()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Romain Teyssier, 
!! Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          25-03-2015 
!! \b last \b modified: 25-03-2015
!<
!===============================================================================
!> Compute source term
!===============================================================================
subroutine source_term
  use params
  use variables
  implicit none

  integer :: i, j, k

  !$py begin_statement

  ! Gravity term
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do k = 1, nz
     do j = 1, ny
        do i = 1, nx
           uin(i,j,k,iu) = uin(i,j,k,iu) &
                & + half*(uin(i,j,k,ir) + uin_old(i,j,k,ir))*gravin(i,j,k,1)*dt
           if (ndim > 1) uin(i,j,k,iv) = uin(i,j,k,iv) &
                & + half*(uin(i,j,k,ir) + uin_old(i,j,k,ir))*gravin(i,j,k,2)*dt
           if (ndim > 2) uin(i,j,k,iw) = uin(i,j,k,iw) &
                & + half*(uin(i,j,k,ir) + uin_old(i,j,k,ir))*gravin(i,j,k,3)*dt
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

  return
end subroutine source_term
!===============================================================================
!> Gravity predictor
!===============================================================================
subroutine gravity_predictor(v, igrav, jgrav, kgrav)
  use params
  use variables
  implicit none

  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,3), intent(inout) :: v
  logical, intent(in) :: igrav, jgrav, kgrav
  integer :: ilo, ihi, jlo, jhi, klo, khi
  integer :: i, j, k

  ilo = min(1,iu1+1); ihi = max(1,iu2-1)
  jlo = min(1,ju1+1); jhi = max(1,ju2-1)
  klo = min(1,ku1+1); khi = max(1,ku2-1)

  ! v = v + 1/2*gravin*dt
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           if (igrav) v(i,j,k,1) = v(i,j,k,1) + half*gravin(i,j,k,1)*dt
           if ((jgrav) .and. (ndim > 1)) then
              v(i,j,k,2) = v(i,j,k,2) + half*gravin(i,j,k,2)*dt
           endif
           if ((kgrav) .and. (ndim > 2)) then
              v(i,j,k,3) = v(i,j,k,3) + half*gravin(i,j,k,3)*dt
           endif
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

  return
end subroutine gravity_predictor
