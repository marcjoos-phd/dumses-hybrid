!===============================================================================
!> \file inertial_wave/condinit.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is the initial conditions subroutine for inertial wave test.
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Romain Teyssier, 
!! Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          07-04-2014 
!! \b last \b modified: 07-11-2014
!<
!===============================================================================
subroutine condinit
  use params
  use variables
  use mpi_var
  implicit none

  real(dp) :: d0, delta_vx
  integer  :: i, j, k

  d0 = one
  delta_vx = 1.d-3*ciso

  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2
           uin(i,j,k,ir) = d0
           uin(i,j,k,iu) = d0*delta_vx
           uin(i,j,k,3:7) = zero
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO
  
  return
end subroutine condinit
