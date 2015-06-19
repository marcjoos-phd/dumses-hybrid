!===============================================================================
!> \file incomp_shearwave/condinit.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is the initial conditions subroutine for incompressible shearing wave
!! test.
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Romain Teyssier, 
!! Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          07-07-2014 
!! \b last \b modified: 01-05-2015
!<
!===============================================================================
subroutine condinit
  use params
  use variables
  use mpi_var
  implicit none

  integer  :: i, j, k
  real(dp) :: d0, delta_vx, kx0, ky0, Lx, Ly
  real(dp), dimension(:,:), allocatable :: Az

  d0 = one
  Lx = xmax - xmin
  Ly = ymax - ymin
  delta_vx = 1.d-4*ciso
  kx0 = -8.d0*twopi/Lx
  ky0 = 2.d0*twopi/Ly
  
  allocate(Az(iu1:iu2,ju1:ju2))

  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do j = ju1, ju2
     do i = iu1, iu2
        Az(i,j) = delta_vx/ky0*sin(kx0*x(i) + ky0*y(j))
     enddo
  enddo
  !$OMP END PARALLEL DO

  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do k = 1, nz
     do j = 1, ny
        do i = 1, nx
           uin(i,j,k,ir) = d0
           uin(i,j,k,iu) = d0*(Az(i,j+1) - Az(i,j))/dy
           uin(i,j,k,iv) = -d0*(Az(i+1,j) - Az(i,j))/dx
           uin(i,j,k,iw:iC+3) = zero
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

  deallocate(Az)

  return
end subroutine condinit
