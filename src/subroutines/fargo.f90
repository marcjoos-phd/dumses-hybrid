!===============================================================================
!> \file fargo.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is FARGO algorithm subroutines.
!! \details
!! Contains fargo_update()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Romain Teyssier
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          01-06-2015
!! \b last \b modified: 05-26-2015
!<
!===============================================================================
!> Update the conservative variables after the MHD step using the FARGO algorithm
!===============================================================================
subroutine fargo_update
  use params
  use variables
  implicit none

  real(dp), dimension(:,:,:,:), allocatable :: duin
  real(dp), dimension(:,:), allocatable     :: demf
  real(dp) :: vK, eps, lambda
  integer  :: ihi, jhi, khi, nshft, jshft, jshftp
  integer  :: i, j, k

  !$py begin_statement

  allocate(duin(iu1:iu2,ju1:ju2,ku1:ku2,nvar))
  emfx = zero; emfz = zero

  ihi = max(1,if2-1); jhi = max(1,jf2-1); khi = max(1, kf2-1)

  ! Boundary conditions along y & z
#if MPI == 1
#if NDIM > 1
  call boundary_y
#endif
#if NDIM == 3
  call boundary_z
#endif
#else
#if NDIM > 1
  call innery_boundary
  call outery_boundary
#endif
#if NDIM == 3
  call innerz_boundary
  call outerz_boundary
#endif
#endif

  call uin_slope(duin)

  uin_old = uin

  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(vK, nshft, eps, jshft, jshftp)
  do i = 1, nx
     vK    = -1.5d0*Omega0*x(i)
     nshft = int(abs(vK*dt)/dy)
     eps   = mod(abs(vK*dt),dy)
     eps   = eps/dy

     do j = 1, ny
        if (vK > zero) then
           jshft = j - nshft - 1
           if (jshft < 1) jshft = jshft + ny
           jshftp = jshft + 1
           if (jshftp > ny) jshftp = jshftp - ny

           uin(i,j,ku1:ku2,ir:ip) = uin_old(i,jshft,ku1:ku2,ir:ip)*eps &
                                & + uin_old(i,jshftp,ku1:ku2,ir:ip)*(one - eps)
           uin(i,j,ku1:ku2,ir:ip) = uin(i,j,ku1:ku2,ir:ip) &
                & + half*duin(i,jshft,ku1:ku2,ir:ip)*(forth - (half - eps)**2)&
                & + half*duin(i,jshftp,ku1:ku2,ir:ip)*((half - eps)**2 - forth)
        else
           jshft = j + nshft
           if (jshft > ny) jshft = jshft - ny
           jshftp = jshft + 1
           if (jshftp > ny) jshftp = jshftp - ny

           uin(i,j,ku1:ku2,ir:ip) = uin_old(i,jshft,ku1:ku2,ir:ip)*(one - eps) &
                                & + uin_old(i,jshftp,ku1:ku2,ir:ip)*eps
           uin(i,j,ku1:ku2,ir:ip) = uin(i,j,ku1:ku2,ir:ip) &
                & + half*duin(i,jshft,ku1:ku2,ir:ip)*(forth - (half - eps)**2) &
                & + half*duin(i,jshftp,ku1:ku2,ir:ip)*((half - eps)**2 - forth)
        endif
     enddo
  enddo
  !$OMP END PARALLEL DO

  ! Compute x & y EMFs
  call get_emfx(duin)
  call get_emfz(duin)

  ! Constrained transport (for 3D & cartesian geometry)
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do k = 1, khi+1
     do j = 1, jhi+1
        do i = 1, ihi+1
           ! Left face B-field
           uin(i,j,k,iA) = uin(i,j,k,iA) + (emfz(i,j+1,k) - emfz(i,j,k))/dy
           uin(i,j,k,iB) = uin(i,j,k,iB) - (emfz(i+1,j,k) - emfz(i,j,k))/dx &
                                       & + (emfx(i,j,k+1) - emfx(i,j,k))/dz
           uin(i,j,k,iC) = uin(i,j,k,iC) - (emfx(i,j+1,k) - emfx(i,j,k))/dy

           ! Right face B-field
           uin(i,j,k,iA+3) = uin(i,j,k,iA) + (emfz(i+1,j+1,k) - emfz(i,j,k))/dy
           uin(i,j,k,iB+3) = uin(i,j,k,iB) - (emfz(i+1,j+1,k) - emfz(i,j,k))/dx&
                                         & + (emfx(i,j+1,k+1) - emfx(i,j,k))/dz
           uin(i,j,k,iC+3) = uin(i,j,k,iC) - (emfx(i,j+1,k+1) - emfx(i,j,k))/dy
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

  deallocate(duin)
  
  return
end subroutine fargo_update
!===============================================================================
!> Compute x EMF
!===============================================================================
subroutine get_emfx(duin)
  use params
  use variables
  implicit none

  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar), intent(in) :: duin
  real(dp), dimension(ku1:ku2) :: demf
  real(dp) :: vK, eps, lambda
  integer  :: nshft, jshft, i, j, jj
  
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(vK, nshft, eps, lambda, jshft) &
  !$OMP PRIVATE(demf)
  do i = 1, nx
     vK = -1.5d0*Omega0*x(i)
     nshft = int(abs(vK*dt)/dy)
     eps   = mod(abs(vK*dt),dy)
     eps   = eps/dy

     ! Compute emfx
     do j = 1, ny+1
        if (vK > zero) then
           lambda = half*eps*(one - eps)
           if (nshft >= 1) then
              do jj = 1, nshft
                 jshft = j - jj
                 if (jshft < 1) jshft = jshft + ny
                 emfx(i,j,:) = emfx(i,j,:) + uin(i,jshft,:,iC)*dy
              enddo
           endif
           jshft = j - nshft - 1
           if (jshft < 1) jshft = jshft + ny
           demf = eps*uin(i,jshft,:,iC) + lambda*duin(i,jshft,:,iC)
           emfx(i,j,:) = emfx(i,j,:) + demf*dy
        else
           lambda = half*eps*(eps - one)
           if (nshft >= 1) then
              do jj = 0, nshft-1
                 jshft = j + jj
                 if (jshft > ny) jshft = jshft - ny
                 emfx(i,j,:) = emfx(i,j,:) - uin(i,jshft,:,iC)*dy
              enddo
           endif
           jshft = j + nshft
           if (jshft > ny) jshft = jshft - ny
           demf = -eps*uin(i,jshft,:,iC) - lambda*duin(i,jshft,:,iC)
           emfx(i,j,:) = emfx(i,j,:) + demf*dy
        endif
     enddo
  enddo
  !$OMP END PARALLEL DO

  return
end subroutine get_emfx
!===============================================================================
!> Compute z EMF
!===============================================================================
subroutine get_emfz(duin)
  use params
  use variables
  implicit none

  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar), intent(in) :: duin
  real(dp), dimension(ku1:ku2) :: demf
  real(dp) :: vK, eps, lambda
  integer  :: nshft, jshft, i, j, jj
  
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(vK, nshft, eps, lambda, jshft) &
  !$OMP PRIVATE(demf)
  do i = 1, nx+1
     vK = -1.5d0*Omega0*half*(x(i) + x(i-1))
     nshft = int(abs(vK*dt)/dy)
     eps   = mod(abs(vK*dt),dy)
     eps   = eps/dy

     ! Compute emfz
     do j = 1, ny+1
        if (vK > zero) then
           lambda = half*eps*(one - eps)
           if (nshft >= 1) then
              do jj = 1, nshft
                 jshft = j - jj
                 if (jshft < 1) jshft = jshft + ny
                 emfz(i,j,:) = emfz(i,j,:) - uin(i,jshft,:,iA)*dy
              enddo
           endif
           jshft = j - nshft - 1
           if (jshft < 1) jshft = jshft + ny
           demf = eps*uin(i,jshft,:,iA) + lambda*duin(i,jshft,:,iA)
           emfz(i,j,:) = emfz(i,j,:) - demf*dy
        else
           lambda = half*eps*(eps - one)
           if (nshft >= 1) then
              do jj = 0, nshft-1
                 jshft = j + jj
                 if (jshft > ny) jshft = jshft - ny
                 emfz(i,j,:) = emfz(i,j,:) + uin(i,jshft,:,iA)*dy
              enddo
           endif
           jshft = j + nshft
           if (jshft > ny) jshft = jshft - ny
           demf = -eps*uin(i,jshft,:,iA) - lambda*duin(i,jshft,:,iA)
           emfz(i,j,:) = emfz(i,j,:) - demf*dy
        endif
     enddo
  enddo
  !$OMP END PARALLEL DO

  return
end subroutine get_emfz
!===============================================================================
!> Compute uin slope
!===============================================================================
subroutine uin_slope(duin)
  use params
  use variables
  implicit none

  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar), intent(inout) :: duin

  real(dp) :: dlft, drgt, dcen, dsgn, slop, dlim
  integer  :: slope_type_loc, i, j, k, ivar

  slope_type_loc = min(slope_type, 2)
  if (slope_type_loc == 0) then
     duin = zero
     return
  endif

  if ((slope_type_loc == 1) .or. (slope_type_loc == 2)) then
     do ivar = 1, nvar
        !$acc kernels loop
        !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dlft, drgt, dcen, dsgn) &
        !$OMP PRIVATE(slop, dlim)
        do k = 1, nz+1
           do j = 1, ny+1
              do i = 1, nx+1
                 dlft = slope_type_loc*(uin(i,j,k,ivar) - uin(i,j-1,k,ivar))
                 drgt = slope_type_loc*(uin(i,j+1,k,ivar) - uin(i,j,k,ivar))
                 dcen = half*(dlft + drgt)/slope_type_loc
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft), abs(drgt))
                 dlim = slop
                 if ((dlft*drgt) <= zero) dlim = zero
                 duin(i,j,k,ivar) = dsgn*min(dlim, abs(dcen))
              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO
     enddo
  else
     print*, "Unknown slope type!"
     stop
  endif

  return
end subroutine uin_slope
