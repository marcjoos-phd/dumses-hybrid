!===============================================================================
!> \file shock_tube/condinit.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is the initial conditions subroutine for shock tube problem.
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
!! \b last \b modified: 01-27-2014
!<
!===============================================================================
subroutine condinit
  use params
  use variables
  implicit none

  integer  :: i, j, k
  character(len=3) :: direction

  namelist /init_params/ direction
  
  open(unit=1, file='input', status='old')
  read(1, init_params)
  close(1)

  if(verbose) print*, 'Dimension (ndim) of the problem: ', ndim

  ! Sod shock tube initial conditions
  gamma = 5.d0/3.d0
  if (direction == "x") then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = ku1, ku2
        do j = ju1, ju2
           do i = iu1, iu2
              ! Density & energy
              if (x(i) < half) then
                 uin(i,j,k,ir) = 1.d0
                 uin(i,j,k,ip) = 1.d0/(gamma - 1.d0)
              else
                 uin(i,j,k,ir) = 0.125d0
                 uin(i,j,k,ip) = 0.1d0/(gamma - 1.d0)
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
              ! Density & energy
              if (y(j) < half) then
                 uin(i,j,k,ir) = 1.d0
                 uin(i,j,k,ip) = 1.d0/(gamma - 1.d0)
              else
                 uin(i,j,k,ir) = 0.125d0
                 uin(i,j,k,ip) = 0.1d0/(gamma - 1.d0)
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
              ! Density & energy
              if (z(k) < half) then
                 uin(i,j,k,ir) = 1.d0
                 uin(i,j,k,ip) = 1.d0/(gamma - 1.d0)
              else
                 uin(i,j,k,ir) = 0.125d0
                 uin(i,j,k,ip) = 0.1d0/(gamma - 1.d0)
              endif
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif

  return
end subroutine condinit
