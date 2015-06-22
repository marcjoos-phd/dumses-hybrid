!===============================================================================
!> \file fargo/condinit.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is the initial conditions subroutine for FARGO test.
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Romain Teyssier, 
!! Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          05-22-2015
!! \b last \b modified: 06-22-2015
!<
!===============================================================================
subroutine condinit
  use params
  use variables
  implicit none

  real(dp) :: v0, A0, R
  real(dp) :: dist
  real(dp), dimension(:,:), allocatable :: Az
  integer  :: i, j, k

  namelist /init_params/ R, A0, v0
  
  R  = 1.0d0
  A0 = 1.d-3
  v0 = ciso

  open(unit=1, file='input', status='old')
  read(1, init_params)
  close(1)

  if (verbose) print*, 'Dimension (ndim) of the problem: ', ndim

  ! Advection of a magnetic loop
  allocate(Az(iu1:iu2+1,ju1:ju2+1))

  Az = 0.d0

  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dist)
  do j = ju1, ju2
     do i = iu1, iu2
        dist = sqrt(x(i)**2 + y(j)**2)
        if (dist <= R) then
           Az(i,j) = A0*(R - dist)
        else
           Az(i,j) = 0.d0
        endif
     end do
  end do
  !$OMP END PARALLEL DO

  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2
           uin(i,j,k,ir) = 1.d0
           uin(i,j,k,iu) = v0
           uin(i,j,k,iv) = 0.d0
           uin(i,j,k,iw) = 0.d0
           
           uin(i,j,k,iA)   = (Az(i,j+1) - Az(i,j))/dy    
           uin(i,j,k,iB)   = -(Az(i+1,j) - Az(i,j))/dx    
           uin(i,j,k,iC)   = 0.d0    
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

  ! Isothermal EOS
  !$OMP PARALLEL WORKSHARE
  uin(:,:,:,ip) = 0.d0
  !$OMP END PARALLEL WORKSHARE

  deallocate(Az)

  return
end subroutine condinit
