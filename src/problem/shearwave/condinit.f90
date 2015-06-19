!===============================================================================
!> \file shearwave/condinit.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is the initial conditions subroutine for compressible shearing wave
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
!! \b created:          07-08-2014 
!! \b last \b modified: 07-11-2014
!<
!===============================================================================
subroutine condinit
  use params
  use variables
  use mpi_var
  implicit none

  integer  :: i, j, k
  real(dp) :: d0, delta_vx, delta_vy, delta_rho, kx0, ky0, xi0, Lx, Ly

  d0 = one
  Lx = xmax - xmin
  Ly = ymax - ymin
  delta_vx = -4.d-4*ciso
  delta_vy = 1.d-4*ciso
  kx0 = -4.d0*twopi/Lx
  ky0 = 1.d0*twopi/Ly
  delta_rho = kx0*delta_vy - ky0*delta_vx
  xi0 = half*Omega0/d0
  delta_rho = delta_rho/xi0

  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do k = 1, nz
     do j = 1, ny
        do i = 1, nx
           uin(i,j,k,ir) = d0*(one - delta_rho*sin(kx0*x(i) + ky0*y(j)))
           uin(i,j,k,iu) = uin(i,j,k,ir)*delta_vx*cos(kx0*x(i) + ky0*y(j))
           uin(i,j,k,iv) = uin(i,j,k,ir)*delta_vy*cos(kx0*x(i) + ky0*y(j))
           uin(i,j,k,iw:iC+3) = zero
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

  return
end subroutine condinit
