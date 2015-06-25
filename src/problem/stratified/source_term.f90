!===============================================================================
!> \file stratified/source term.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is the source term for the problem of a stratified shearing box
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Romain Teyssier, 
!! Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          04-03-2015 
!! \b last \b modified: 05-18-2015
!<
!===============================================================================
subroutine source_term
  use params
  use mpi_var
  use variables
  use stratified
  implicit none

  integer :: i, j, k, l

  real(dp) :: rhon, rhonp1, vxn, vyn ! Gravity force parameters
  real(dp) :: theta, smallH, mass    ! Addmass parameters

  if (verbose) write (*,*) 'Entering source_term subroutine...'

  ! Gravity term
  !$acc data copyin(uin_old,gravin)
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(rhon, rhonp1, vxn, vyn)
  do k = 1, nz
     !$acc loop
     do j = 1, ny
        do i = 1, nx

           !density
           rhon   = uin_old(i,j,k,1)
           rhonp1 = uin(i,j,k,1)

           !update momentum
           uin(i,j,k,2) = uin(i,j,k,2) + half*(rhon + rhonp1)*gravin(i,j,k,1)*dt
           if (ndim>1) uin(i,j,k,3) = uin(i,j,k,3) + half*(rhon + rhonp1)*gravin(i,j,k,2)*dt
           if (ndim>2) uin(i,j,k,4) = uin(i,j,k,4) + half*(rhon + rhonp1)*gravin(i,j,k,3)*dt

#if ISO == 0
           !update energy
           uin(i,j,k,5) = uin(i,j,k,5) + &
                half*(uin_old(i,j,k,2) + uin(i,j,k,2))*gravin(i,j,k,1)*dt + &
                half*(uin_old(i,j,k,3) + uin(i,j,k,3))*gravin(i,j,k,2)*dt + &
                half*(uin_old(i,j,k,4) + uin(i,j,k,4))*gravin(i,j,k,3)*dt

           !Total energy source term
           vxn = uin(i,j,k,2)/rhon
           vyn = uin(i,j,k,3)/rhon
           uin(i,j,k,5) = uin(i,j,k,5) + 1.5d0*Omega0*dt*rhon*vxn*vyn

#endif
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  !$acc end data

  ! Set density field to floor
  if (floor) then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = 1, nz
        !$acc loop
        do j = 1, ny
           do i = 1, nx
                 
                 !density
#if ISO == 1
                 uin(i,j,k,1) = max(uin(i,j,k,1), exp(-zfloor**2/two))
#else
                 uin(i,j,k,1) = max(uin(i,j,k,1), dfloor)
#endif
                 
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  endif  

  ! Add source term in continuity equation
  if (addmass) then
     smallH = 0.2d0
     mass = zero
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME) REDUCTION(+:mass)
     do k = 1, nz
        !$acc loop
        do j = 1, ny
           do i = 1, nx
              mass = mass + uin(i,j,k,1)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
#if MPI==1
     call sumBcast(mass, 1) !subroutine sumBcast in stratified/condinit.f90...
#endif
     mass  = mass*dx*dy*dz
     theta = (mass0 - mass)/sqrt(two*pi)/smallH/(xmax - xmin)/(ymax - ymin)
     
     ! Restore initial mass
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = 1, nz
        !$acc loop
        do j = 1, ny
           do i = 1, nx
              uin(i,j,k,2:4) = uin(i,j,k,2:4)/uin(i,j,k,1)
              uin(i,j,k,1)   = uin(i,j,k,1) + theta*exp(-z(k)**2/two/smallH**2)
              uin(i,j,k,2:4) = uin(i,j,k,2:4)*uin(i,j,k,1)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  endif

  if (verbose) write (*,*) 'End of source_term subroutine...'

end subroutine source_term
!===============================================================================
!> Check for low beta regions
!===============================================================================
subroutine checkLowBeta(qin, dq, dbf, flagCell)
  use params
  implicit none

  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar       ), intent(in) :: qin 
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim), intent(in) :: dq
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:3   ,1:ndim), intent(in) :: dbf
  logical , dimension(iu1:iu2,ju1:ju2,ku1:ku2), intent(out) :: flagCell

  real(dp) :: rho, bx, by, bz, pmag, pth, beta
  real(dp) :: minBeta
  integer  :: i, j, k, ncell
  integer  :: ilo, ihi, jlo, jhi, klo, khi

  minBeta = 1.d-3; ncell = 0; flagCell = .False.

  ilo = min(1,iu1+3); ihi = max(1,iu2-3)
  jlo = min(1,ju1+3); jhi = max(1,ju2-3)
  klo = min(1,ku1+3); khi = max(1,ku2-3)

  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(rho, bx, by, bz, pth, pmag, beta)
  do k = klo, khi
     !$acc loop
     do j = jlo, jhi
        do i = ilo, ihi
           rho = qin(i,j,k,1)
           bx  = qin(i,j,k,6)
           by  = qin(i,j,k,7)
           bz  = qin(i,j,k,8)
           
           pth  = rho*ciso**2
           pmag = half*(bx*bx + by*by + bz*bz)
           
           beta = pth/pmag
           
           if (beta < minBeta) then
              !dq(i,j,k,:,:)=zero
              !dbf(l,i:i+1,j:j+1,k:k+1,:,:)=zero
              !flagCell(l,i:i+1,j:j+1,k:k+1)=.True.
              !ncell=ncell+1
           endif
           
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  if (ncell > 0) write (*,*) '   Nb of cells affected by low beta:', ncell

  return
end subroutine checkLowBeta
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
