!===============================================================================
!> \file wind_tunnel/condinit.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is the initial conditions subroutine for 2D Sedov blast problem.
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Romain Teyssier, 
!! Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          04-15-2013 
!! \b last \b modified: 07-29-2014
!<
!===============================================================================
subroutine condinit
  use params
  use variables
  implicit none

  integer  :: i, j, k, iseed
  character(len=3) :: direction

  namelist /init_params/ direction
  
  open(unit=1, file='input', status='old')
  read(1, init_params)
  close(1)

  if(verbose) print*, 'Dimension (ndim) of the problem: ', ndim

  ! Wind tunnel with point explosion
  if (direction == "x") then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = ku1, ku2
        do j = ju1, ju2
           do i = iu1, iu2
              uin(i,j,k,ir) = 1.d0
              uin(i,j,k,ip) = 1.d-5
              if((x(i) > (xmin - half*dx)) .and. (x(i) <= (xmin + half*dx)) &
           .and. (y(j) > (ymin - half*dy)) .and. (y(j) <= (ymin + half*dy))) &
              then
                 uin(i,j,k,ip) = 1.d0/dx/dx
              endif
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  else if (direction == "y") then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = ku1, ku2
        do j = ju1, ju2
           do i = iu1, iu2
              uin(i,j,k,ir) = 1.d0
              uin(i,j,k,ip) = 1.d-5
              if((y(j) > (ymin - half*dy)) .and. (y(j) <= (ymin + half*dy)) &
           .and. (z(k) > (zmin - half*dz)) .and. (z(k) <= (zmin + half*dz))) &
              then
                 uin(i,j,k,ip) = 1.d0/dy/dy
              endif
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  else
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = ku1, ku2
        do j = ju1, ju2
           do i = iu1, iu2
              uin(i,j,k,ir) = 1.d0
              uin(i,j,k,ip) = 1.d-5
              if((x(i) > (xmin - half*dx)) .and. (x(i) <= (xmin + half*dx)) &
           .and. (z(k) > (zmin - half*dz)) .and. (z(k) <= (zmin + half*dz))) &
              then
                 uin(i,j,k,ip) = 1.d0/dz/dz
              endif
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif
   
  return
end subroutine condinit
