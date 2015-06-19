!===============================================================================
!> \file boundary.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is boundary conditions subroutines.
!! \details
!! Contains boundary(), innerx_boundary(), outerx_boundary(), innery_boundary(),
!! outery_boundary(), innerz_boundary(), outerz_boundary()
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
!! \b last \b modified: 09-22-2014
!<
!===============================================================================
!> Compute boundary conditions for each domain
!===============================================================================
subroutine boundary
  use params
  use variables
  use mpi_var
  implicit none

  !$py begin_statement

  ! Boundary condition: x direction
#if MPI == 1
  call boundary_x
  if (((mype < xleft) .and. (boundary_type(1) /= 'periodic')) &
       .or. (nxslice == 1)) call innerx_boundary
  if (((mype > xright) .and. (boundary_type(1) /= 'periodic')) &
       .or. (nxslice == 1)) call outerx_boundary
#else
  call innerx_boundary
  call outerx_boundary
#endif

  ! Shearing boundary
#if NDIM == 3
  if (boundary_type(1) == "shearingbox") then
#if MPI == 1
     call boundary_y
#else
     call innery_boundary
     call outery_boundary
#endif
     !$py start_timing shearing boundary
     !$acc update host(dt)
     call shearing_boundary(time+dt)
     !$py end_timing shearing boundary
  endif
#endif

  ! Boundary condition: z direction
#if NDIM == 3
#if MPI == 1
  call boundary_z
  if (((mype < zleft) .and. (boundary_type(3) /= 'periodic')) &
       .or. (nzslice == 1)) call innerz_boundary
  if (((mype > zright) .and. (boundary_type(3) /= 'periodic')) &
       .or. (nzslice == 1)) call outerz_boundary
#else
  call innerz_boundary
  call outerz_boundary
#endif
#endif

  ! Boundary condition: y direction
#if NDIM > 1
#if MPI == 1
  call boundary_y
  if (((mype < yleft) .and. (boundary_type(2) /= 'periodic')) &
       .or. (nyslice == 1)) call innery_boundary
  if (((mype > yright) .and. (boundary_type(2) /= 'periodic')) &
       .or. (nyslice == 1)) call outery_boundary
#else
  call innery_boundary
  call outery_boundary
#endif
#endif

  ! Right face B-field
  !$OMP PARALLEL WORKSHARE
  uin(iu1:iu2-1,ju1:ju2,ku1:ku2,iA+3) = uin(iu1+1:iu2,ju1:ju2,ku1:ku2,iA)
  uin(iu1:iu2,ju1:ju2-1,ku1:ku2,iB+3) = uin(iu1:iu2,ju1+1:ju2,ku1:ku2,iB)
  uin(iu1:iu2,ju1:ju2,ku1:ku2-1,iC+3) = uin(iu1:iu2,ju1:ju2,ku1+1:ku2,iC)
  !$OMP END PARALLEL WORKSHARE

  return
end subroutine boundary
!===============================================================================
!> Compute boundary conditions in the x-direction, inner edge
!===============================================================================
subroutine innerx_boundary
  use params
  use variables
  implicit none

  integer :: j, k

  if (boundary_type(1) == 'periodic') then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = ku1,ku2
        do j = ju1,ju2
           uin(iu1  ,j,k,:) = uin(nx-2,j,k,:)
           uin(iu1+1,j,k,:) = uin(nx-1,j,k,:)
           uin(iu1+2,j,k,:) = uin(nx  ,j,k,:)
        end do
     end do
     !$OMP END PARALLEL DO
  else if (boundary_type(1) == 'zero_gradient') then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = ku1,ku2
        do j = ju1,ju2
           uin(iu1  ,j,k,:) = uin(1,j,k,:)
           uin(iu1+1,j,k,:) = uin(1,j,k,:)
           uin(iu1+2,j,k,:) = uin(1,j,k,:)
        end do
     end do
     !$OMP END PARALLEL DO
  else if (boundary_type(1) == 'analytic') then
     call xinner_ana
  endif

  return
end subroutine innerx_boundary
!===============================================================================
!> Compute boundary conditions in the x-direction, outer edge
!===============================================================================
subroutine outerx_boundary
  use params
  use variables
  implicit none

  integer :: j, k
  
  if (boundary_type(1) == 'periodic') then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = ku1,ku2
        do j = ju1,ju2
           uin(iu2-2,j,k,:) = uin(1,j,k,:)
           uin(iu2-1,j,k,:) = uin(2,j,k,:)
           uin(iu2  ,j,k,:) = uin(3,j,k,:)
        end do
     end do
     !$OMP END PARALLEL DO
  else if (boundary_type(1) == 'zero_gradient') then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = ku1,ku2
        do j = ju1,ju2
           uin(iu2-2,j,k,:) = uin(nx,j,k,:)
           uin(iu2-1,j,k,:) = uin(nx,j,k,:)
           uin(iu2  ,j,k,:) = uin(nx,j,k,:)
        end do
     end do
     !$OMP END PARALLEL DO
  else if (boundary_type(1) == 'analytic') then
     call xouter_ana
  endif

  return
end subroutine outerx_boundary
!===============================================================================
#if NDIM > 1
!> Compute boundary conditions in the y-direction, inner edge
!===============================================================================
subroutine innery_boundary
  use params
  use variables
  implicit none

  integer :: i, k
  
  if (boundary_type(2) == 'periodic') then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = ku1,ku2
        do i = iu1,iu2
           uin(i,ju1  ,k,:) = uin(i,ny-2,k,:)
           uin(i,ju1+1,k,:) = uin(i,ny-1,k,:)
           uin(i,ju1+2,k,:) = uin(i,ny  ,k,:)
        end do
     end do
     !$OMP END PARALLEL DO
  else if (boundary_type(2) == 'zero_gradient') then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = ku1,ku2
        do i = iu1,iu2
           uin(i,ju1  ,k,:) = uin(i,1,k,:)
           uin(i,ju1+1,k,:) = uin(i,1,k,:)
           uin(i,ju1+2,k,:) = uin(i,1,k,:)
        end do
     end do
     !$OMP END PARALLEL DO
  else if (boundary_type(2) == 'analytic') then
     call yinner_ana
  endif

  return
end subroutine innery_boundary
!===============================================================================
!> Compute boundary conditions in the y-direction, outer edge
!===============================================================================
subroutine outery_boundary
  use params
  use variables
  implicit none

  integer :: i, k
  
  if (boundary_type(2) == 'periodic') then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = ku1,ku2
        do i = iu1,iu2
           uin(i,ju2-2,k,:) = uin(i,1,k,:)
           uin(i,ju2-1,k,:) = uin(i,2,k,:)
           uin(i,ju2  ,k,:) = uin(i,3,k,:)
        end do
     end do
     !$OMP END PARALLEL DO
  else if (boundary_type(2) == 'zero_gradient') then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = ku1,ku2
        do i = iu1,iu2
           uin(i,ju2-2,k,:) = uin(i,ny,k,:)
           uin(i,ju2-1,k,:) = uin(i,ny,k,:)
           uin(i,ju2  ,k,:) = uin(i,ny,k,:)
        end do
     end do
     !$OMP END PARALLEL DO
  else if (boundary_type(2) == 'analytic') then
     call youter_ana
  endif

  return
end subroutine outery_boundary
#endif
!===============================================================================
#if NDIM == 3
!> Compute boundary conditions in the z-direction, inner edge
!===============================================================================
subroutine innerz_boundary
  use params
  use variables
  implicit none

  integer :: i, j

  if (boundary_type(3) == 'periodic') then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do j = ju1,ju2
        do i = iu1,iu2
           uin(i,j,ku1  ,:) = uin(i,j,nz-2,:)
           uin(i,j,ku1+1,:) = uin(i,j,nz-1,:)
           uin(i,j,ku1+2,:) = uin(i,j,nz  ,:)
        end do
     end do
     !$OMP END PARALLEL DO
  else if (boundary_type(3) == 'zero_gradient') then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do j = ju1,ju2
        do i = iu1,iu2
           uin(i,j,ku1  ,:) = uin(i,j,1,:)
           uin(i,j,ku1+1,:) = uin(i,j,1,:)
           uin(i,j,ku1+2,:) = uin(i,j,1,:)
        end do
     end do
     !$OMP END PARALLEL DO
  else if (boundary_type(3) == 'analytic') then
     call zinner_ana
  endif

  return
end subroutine innerz_boundary
!===============================================================================
!> Compute boundary conditions in the z-direction, outer edge
!===============================================================================
subroutine outerz_boundary
  use params
  use variables
  implicit none

  integer :: i, j

  if (boundary_type(3) == 'periodic') then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do j = ju1,ju2
        do i = iu1,iu2
           uin(i,j,ku2-2,:) = uin(i,j,1,:)
           uin(i,j,ku2-1,:) = uin(i,j,2,:)
           uin(i,j,ku2  ,:) = uin(i,j,3,:)
        end do
     end do
     !$OMP END PARALLEL DO
  else if (boundary_type(3) == 'zero_gradient') then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do j = ju1,ju2
        do i = iu1,iu2
           uin(i,j,ku2-2,:) = uin(i,j,nz,:)
           uin(i,j,ku2-1,:) = uin(i,j,nz,:)
           uin(i,j,ku2  ,:) = uin(i,j,nz,:)
        end do
     end do
     !$OMP END PARALLEL DO
  else if (boundary_type(3) == 'analytic') then
     call zouter_ana
  endif

  return
end subroutine outerz_boundary
#endif
