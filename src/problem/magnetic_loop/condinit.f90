!===============================================================================
!> \file magnetic_loop/condinit.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is the initial conditions subroutine for magnetic loop problem.
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
!! \b last \b modified: 07-11-2014
!<
!===============================================================================
subroutine condinit
  use params
  use variables
  implicit none

  real(dp) :: v0, A0, R, cos_theta, sin_theta, y0
  real(dp) :: dist
  real(dp), dimension(:,:), allocatable :: Az
  integer  :: i, j, k
  character(len=3) :: direction

  namelist /init_params/ direction
  
  open(unit=1, file='input', status='old')
  read(1, init_params)
  close(1)

  if (verbose) print*, 'Dimension (ndim) of the problem: ', ndim

  ! Advection of a magnetic loop
  if (direction == "x") then
     allocate(Az(iu1:iu2+1,ju1:ju2+1))
  else if (direction == "y") then
     allocate(Az(ju1:ju2+1,ku1:ku2+1))
  else
     allocate(Az(ku1:ku2+1,iu1:iu2+1))
  endif
  cos_theta = 2.d0/sqrt(5.d0)
  sin_theta = sqrt(1.d0 - cos_theta**2)
  R  = 1.d0 ; y0 = 0.d0
  A0 = 1.d-6; v0 = sqrt(5.d0)
  Az = 0.d0

  if (direction == "x") then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dist)
     do j = ju1, ju2
        do i = iu1, iu2
           dist = sqrt(x(i)**2 + (y(j) - y0)**2)
           if (dist <= R) then
              Az(i,j) = A0*(R - dist)
           else
              Az(i,j) = 0.d0
           endif
        end do
     end do
     !$OMP END PARALLEL DO
  else if (direction == "y") then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dist)
     do k = ku1, ku2
        do j = ju1, ju2
           dist = sqrt(y(j)**2 + (z(k) - y0)**2)
           if (dist <= R) then
              Az(j,k) = A0*(R - dist)
           else
              Az(j,k) = 0.d0
           endif
        end do
     end do
     !$OMP END PARALLEL DO
  else
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dist)
     do i = iu1, iu2
        do k = ku1, ku2
           dist = sqrt(z(k)**2 + (x(i) - y0)**2)
           if (dist <= R) then
              Az(k,i) = A0*(R - dist)
           else
              Az(k,i) = 0.d0
           endif
        end do
     end do
     !$OMP END PARALLEL DO
  endif

  if (direction == "x") then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = ku1, ku2
        do j = ju1, ju2
           do i = iu1, iu2
              uin(i,j,k,ir) = 1.d0
              uin(i,j,k,iu) = uin(i,j,k,1)*v0*cos_theta
              uin(i,j,k,iv) = uin(i,j,k,1)*v0*sin_theta
              uin(i,j,k,iw) = 0.d0
              
              uin(i,j,k,iA)   = (Az(i,j+1) - Az(i,j))/dy    
              uin(i,j,k,iA+3) = (Az(i+1,j+1) - Az(i+1,j))/dy
              uin(i,j,k,iB)   = -(Az(i+1,j) - Az(i,j))/dx    
              uin(i,j,k,iB+3) = -(Az(i+1,j+1) - Az(i,j+1))/dx
              uin(i,j,k,iC)   = 0.d0    
              uin(i,j,k,iC+3) = 0.d0
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
              uin(i,j,k,iv) = uin(i,j,k,ir)*v0*cos_theta
              uin(i,j,k,iw) = uin(i,j,k,ir)*v0*sin_theta
              uin(i,j,k,iu) = 0.d0
              
              uin(i,j,k,iB)   = (Az(j,k+1) - Az(j,k))/dz    
              uin(i,j,k,iB+3) = (Az(j+1,k+1) - Az(j+1,k))/dz
              uin(i,j,k,iC)   = -(Az(j+1,k) - Az(j,k))/dy    
              uin(i,j,k,iC+3) = -(Az(j+1,k+1) - Az(j,k+1))/dy
              uin(i,j,k,iA)   = 0.d0    
              uin(i,j,k,iA+3) = 0.d0
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  else
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = ku1, ku2
        do j = ju1, ju2
           do i = iu1, iu2
              uin(i,j,k,ir) = 1.d0
              uin(i,j,k,iw) = uin(i,j,k,ir)*v0*cos_theta
              uin(i,j,k,iu) = uin(i,j,k,ir)*v0*sin_theta
              uin(i,j,k,iv) = 0.d0
              
              uin(i,j,k,iC)   = (Az(k,i+1) - Az(k,i))/dx    
              uin(i,j,k,iC+3) = (Az(k+1,i+1) - Az(k+1,i))/dx
              uin(i,j,k,iA)   = -(Az(k+1,i) - Az(k,i))/dz    
              uin(i,j,k,iA+3) = -(Az(k+1,i+1) - Az(k,i+1))/dz
              uin(i,j,k,iB)   = 0.d0    
              uin(i,j,k,iB+3) = 0.d0
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  endif

  ! Isothermal EOS
  !$OMP PARALLEL WORKSHARE
  uin(:,:,:,ip) = 0.d0
  !$OMP END PARALLEL WORKSHARE

  deallocate(Az)

  return
end subroutine condinit
