!===============================================================================
!> \file scheme.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is scheme subroutines. This file cannot be compiled as is, it needs to
!! be preprocessed with the DUMSES preprocessor to generate call to Riemann
!! solvers.
!! \details
!! Contains primitive(), godunov(), fc_magnetic(), slope(), trace3d(), trace2d()
!! , trace1d(), compflux(), update(), constraint_transport()
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
!! \b last \b modified: 06-25-2015
!<
!===============================================================================
!> Compute primitive variables
!===============================================================================
subroutine primitive(dt)
  use params
  use variables, only: uin, qin, x
  use oacc_params
  implicit none

  real(dp), intent(in) :: dt
  integer  :: i, j, k
  real(dp) :: eken, emag, eint, dvx, dvy
  real(dp) :: smalle, smallp

  !$py begin_statement

#if ISO == 1
  smalle = one
  smallp = smallr*ciso**2
#else
  smalle = smallc**2/gamma/(gamma-one)
  smallp = smallr*smallc**2/gamma
#endif

  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(eken, emag, eint, dvx, dvy)
  do k = ku1, ku2
     !$acc loop vector(blocky_primitive)
     do j = ju1, ju2
        !$acc loop vector(blockx_primitive)
        do i = iu1, iu2
           qin(i,j,k,ir) = max(uin(i,j,k,ir), smallr)

           qin(i,j,k,iu) = uin(i,j,k,iu)/uin(i,j,k,ir)
           qin(i,j,k,iv) = uin(i,j,k,iv)/uin(i,j,k,ir)
           qin(i,j,k,iw) = uin(i,j,k,iw)/uin(i,j,k,ir)

           qin(i,j,k,iA) = half*(uin(i,j,k,iA) + uin(i,j,k,iA+3))
           qin(i,j,k,iB) = half*(uin(i,j,k,iB) + uin(i,j,k,iB+3))
           qin(i,j,k,iC) = half*(uin(i,j,k,iC) + uin(i,j,k,iC+3))

#if ISO == 1  
           qin(i,j,k,ip) = uin(i,j,k,ir)*ciso**2
#else
           eken = half*(qin(i,j,k,iu)**2 + qin(i,j,k,iv)**2 + qin(i,j,k,iw)**2)
           emag = half*(qin(i,j,k,iA)**2 + qin(i,j,k,iB)**2 + qin(i,j,k,iC)**2)
           eint = (uin(i,j,k,ip) - emag)/uin(i,j,k,ir) - eken
           qin(i,j,k,ip) = max((gamma - one)*qin(i,j,k,ir)*eint, smallp)
#endif

           ! Coriolis force predictor step
#if GEOM == CARTESIAN
           if (Omega0 > zero) then
              dvx = two*Omega0*qin(i,j,k,iv)
              dvy = -half*Omega0*qin(i,j,k,iu)
              qin(i,j,k,iu) = qin(i,j,k,iu) + dvx*dt*half
              qin(i,j,k,iv) = qin(i,j,k,iv) + dvy*dt*half
           endif
#endif

           ! Remove angular velocity of rotating frame in case of cylindrical geom
#if GEOM == CYLINDRICAL
           if (Omega0 > zero) then
              qin(i,j,k,iv) = qin(i,j,k,iv) - x(i)*Omega0
           endif
#endif

        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

  return
end subroutine primitive
!===============================================================================
!> Godunov scheme subroutine; this is the core routine of the numerical scheme.
!!
!! The main steps are:
!!  -# face-centered magnetic field is computed;
!!  -# primitive variables are computed;
!!  -# slopes at interfaces for hydro variables are computed;
!!  -# a half-timestep prediction is performed;
!!  -# Riemann fluxes are estimated using a Riemann solver;
!!  -# variables are updated using the computed fluxes.
!<
!===============================================================================
subroutine godunov
  use params
  use variables
  implicit none

  integer :: i, j, k, idim

  !$py begin_statement

  !$acc data create(qm, qp, qRT, qRB, qLT, qLB, gravin)
  !$acc data create(dq, bfc, dbfc)

  !$py start_timing Face-centered magnetic field
  call fc_magnetic(bfc)
  !$py end_timing Face-centered magnetic field

  !$py start_timing Primitive
  call primitive(dt)
  !$py end_timing Primitive

  !$py start_timing Slope  	 	 
  call slope(bfc, dq, dbfc)  	 	 
  !$py end_timing Slope 

  !$py start_timing Trace
#if NDIM == 3
  call trace3d(bfc, dq, dbfc, qm, qp, qRT, qRB, qLT, qLB)
#endif
#if NDIM == 2
  call trace2d(bfc, dq, dbfc, qm, qp, qRT, qRB, qLT, qLB)
#endif
#if NDIM == 1
  call trace1d(dq, qm, qp)
#endif
  !$py end_timing Trace
  !$acc end data

  !$py start_timing Gravity predictor
  if (rhs) then
     do idim = 1, ndim
        call gravity_predictor(qm(:,:,:,iu:iw,idim), .true., .true., .true.)
        call gravity_predictor(qp(:,:,:,iu:iw,idim), .true., .true., .true.)
        call gravity_predictor(qRT(:,:,:,iu:iw,idim), .true., .true., .true.)
        call gravity_predictor(qRB(:,:,:,iu:iw,idim), .true., .true., .true.)
        call gravity_predictor(qLT(:,:,:,iu:iw,idim), .true., .true., .true.)
        call gravity_predictor(qLB(:,:,:,iu:iw,idim), .true., .true., .true.)
     enddo
  endif
  !$py end_timing Gravity predictor

  !$acc data create(fgodunov, fgodunov_pre, flux)
  !$py start_timing Riemann
  !$py call_solver riemann (qm, qp, fgodunov, fgodunov_pre)
  !$py end_timing Riemann

  !$acc data create(emfx, emfy, emfz)
  !$py start_timing Riemann magnetic
#if NDIM > 1
  !$py call_solver riemann_magnetic (qRT, qRB, qLT, qLB, emfz, 3)
#endif
#if NDIM == 3
  !$py call_solver riemann_magnetic (qRT, qLT, qRB, qLB, emfy, 2)
  !$py call_solver riemann_magnetic (qRT, qRB, qLT, qLB, emfx, 1)
#endif
  !$py end_timing Riemann magnetic

  !$py start_timing Flux
  call compflux(fgodunov)
  !$py end_timing Flux

#if NDIM == 3
  if (boundary_type(1) == 'shearingbox') then
     !$py start_timing Shearing flux
     !$acc update host(dt)
     call shearing_flux(time+half*dt)
     !$py end_timing Shearing flux
  endif
#endif

  !$py start_timing Update
  call update(fgodunov_pre)
  !$py end_timing Update
  !$acc end data
  !$acc end data
  !$acc end data

  return
end subroutine godunov
!===============================================================================
!> Compute face-centered magnetic field
!===============================================================================
subroutine fc_magnetic(bfc)
  use params
  use variables, only: uin
  implicit none

  real(dp), dimension(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3), intent(out) :: bfc
  integer :: i, j, k

  !$py begin_statement

  ! Store face centered magnetic field
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2+1
           if (i <= iu2) then
              bfc(i,j,k,1) = uin(i,j,k,iA)
           else
              bfc(i,j,k,1) = uin(i-1,j,k,iA+3)
           endif
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do k = ku1, ku2
     do j = ju1, ju2+1
        do i = iu1, iu2
           if (j <= ju2) then
              bfc(i,j,k,2) = uin(i,j,k,iB)
           else
              bfc(i,j,k,2) = uin(i,j-1,k,iB+3)
           endif
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do k = ku1, ku2+1
     do j = ju1, ju2
        do i = iu1, iu2
           if (k <= ku2) then
              bfc(i,j,k,3) = uin(i,j,k,iC)
           else
              bfc(i,j,k,3) = uin(i,j,k-1,iC+3)
           endif
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  return
end subroutine fc_magnetic
!===============================================================================
!> Compute slopes at cell interfaces
!===============================================================================
subroutine slope(bfc, dq, dbfc)
  use params
  use variables, only: qin
  use oacc_params
  implicit none

  real(dp), dimension(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,3),intent(in) :: bfc
  real(dp), dimension(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,nvar,ndim),intent(out) :: dq
  real(dp), dimension(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,3,ndim),intent(out) :: dbfc
  integer  :: ilo, ihi, jlo, jhi, klo, khi
  integer  :: slope_type_loc
  real(dp) :: dsgn, dlim, dcen, dlft, drgt, slop
  real(dp) :: dfll, dflm, dflr, dfml, dfmm, dfmr, dfrl, dfrm, dfrr, vmin, vmax
  real(dp) :: dflll, dflml, dflrl, dfmll, dfmml, dfmrl, dfrll, dfrml, dfrrl
  real(dp) :: dfllm, dflmm, dflrm, dfmlm, dfmmm, dfmrm, dfrlm, dfrmm, dfrrm
  real(dp) :: dfllr, dflmr, dflrr, dfmlr, dfmmr, dfmrr, dfrlr, dfrmr, dfrrr
  real(dp) :: dfx, dfy, dfz, dff
  integer  :: i, j, k, ivar

  !$py begin_statement

  ilo=min(1,iu1+1); ihi=max(1,iu2-1)
  jlo=min(1,ju1+1); jhi=max(1,ju2-1)
  klo=min(1,ku1+1); khi=max(1,ku2-1)

  if (slope_type == 0) then
     !$OMP PARALLEL WORKSHARE
     dq   = zero
     dbfc = zero
     !$OMP END PARALLEL WORKSHARE
     return
  endif

  slope_type_loc = min(slope_type, 2)

#if NDIM == 1
  if (slope_type < 3) then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dsgn, dlim, dcen, dlft, drgt) &
     !$OMP PRIVATE(slop)
     do k = klo, khi
        !$acc loop vector(blocky_slope)
        do j = jlo, jhi
           !$acc loop vector(blockx_slope)
           do i = ilo, ihi
              do ivar = 1, nvar
                 dlft = slope_type*(qin(i  ,j,k,ivar) - qin(i-1,j,k,ivar))
                 drgt = slope_type*(qin(i+1,j,k,ivar) - qin(i  ,j,k,ivar))
                 dcen = half*(dlft + drgt)/slope_type
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft), abs(drgt))
                 dlim = slop
                 if ((dlft*drgt) <= zero) dlim = zero
                 dq(i,j,k,ivar,1) = dsgn*min(dlim, abs(dcen))
              end do
           end do
        enddo
     enddo
     !$OMP END PARALLEL DO
  else
     print*, "Unknown slope type!"
     stop
  endif
#endif

#if NDIM == 2
  if (slope_type < 3) then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dsgn, dlim, dcen, dlft, drgt) &
     !$OMP PRIVATE(slop)
     do k = klo, khi
        !$acc loop vector(blocky_slope)
        do j = jlo, jhi
           !$acc loop vector(blockx_slope)
           do i = ilo, ihi
              do ivar = 1, nvar
                 ! slopes in first coordinate direction
                 dlft = slope_type*(qin(i  ,j,k,ivar) - qin(i-1,j,k,ivar))
                 drgt = slope_type*(qin(i+1,j,k,ivar) - qin(i  ,j,k,ivar))
                 dcen = half*(dlft + drgt)/slope_type
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft), abs(drgt))
                 dlim = slop
                 if ((dlft*drgt) <= zero) dlim = zero
                 dq(i,j,k,ivar,1) = dsgn*min(dlim, abs(dcen))
                 ! slopes in second coordinate direction
                 dlft = slope_type*(qin(i,j  ,k,ivar) - qin(i,j-1,k,ivar))
                 drgt = slope_type*(qin(i,j+1,k,ivar) - qin(i,j  ,k,ivar))
                 dcen = half*(dlft + drgt)/slope_type
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft), abs(drgt))
                 dlim = slop
                 if ((dlft*drgt) <= zero) dlim = zero
                 dq(i,j,k,ivar,2) = dsgn*min(dlim, abs(dcen))
              enddo
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  else if (slope_type == 3) then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dfll,dflm,dflr,dfml,dfmm,dfmr)&
     !$OMP PRIVATE(dfrl,dfrm,dfrr,vmin,vmax,dfx,dfy,dff,slop,dlim)
     do k = klo, khi
        !$acc loop vector(blocky_slope)
        do j = jlo, jhi
           !$acc loop vector(blockx_slope)
           do i = ilo, ihi
              do ivar = 1, nvar
                 dfll = qin(i-1,j-1,k,ivar) - qin(i,j,k,ivar)
                 dflm = qin(i-1,j,k,ivar)   - qin(i,j,k,ivar)
                 dflr = qin(i,j-1,k,ivar)   - qin(i,j,k,ivar)
                 dfml = qin(i,j-1,k,ivar) - qin(i,j,k,ivar)
                 dfmm = qin(i,j,k,ivar)   - qin(i,j,k,ivar)
                 dfmr = qin(i,j+1,k,ivar) - qin(i,j,k,ivar)
                 dfrl = qin(i+1,j-1,k,ivar) - qin(i,j,k,ivar)
                 dfrm = qin(i+1,j,k,ivar)   - qin(i,j,k,ivar)
                 dfrr = qin(i+1,j+1,k,ivar) - qin(i,j,k,ivar)

                 vmin = min(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)
                 vmax = max(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)

                 dfx = half*(qin(i+1,j,k,ivar) - qin(i-1,j,k,ivar))
                 dfy = half*(qin(i,j+1,k,ivar) - qin(i,j-1,k,ivar))
                 dff = half*(abs(dfx) + abs(dfy))
                 
                 if (dff > zero) then
                    slop = min(one, min(abs(vmin), abs(vmax))/dff)
                 else
                    slop = one
                 endif
                 
                 dlim = slop
                 dq(i,j,k,ivar,1) = dlim*dfx
                 dq(i,j,k,ivar,2) = dlim*dfy
              enddo
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  else
     print*, "Unknown slope type!"
     stop
  endif

  ! Bx along Y
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dsgn, dlim, dcen, dlft, drgt, slop)
  do k = klo, khi
     !$acc loop vector(blocky_slope)
     do j = jlo, jhi
        !$acc loop vector(blockx_slope)
        do i = ilo, ihi+1
           dlft = slope_type_loc*(bfc(i,j  ,k,1) - bfc(i,j-1,k,1))
           drgt = slope_type_loc*(bfc(i,j+1,k,1) - bfc(i,j  ,k,1))
           dcen = half*(dlft + drgt)/slope_type_loc
           dsgn = sign(one, dcen)
           slop = min(abs(dlft), abs(drgt))
           dlim = slop
           if ((dlft*drgt) <= zero) dlim = zero
           dbfc(i,j,k,1,2) = dsgn*min(dlim, abs(dcen))
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

  ! By along X
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dsgn, dlim, dcen, dlft, drgt, slop)
  do k = klo, khi
     do j = jlo, jhi+1
        do i = ilo, ihi
           dlft = slope_type_loc*(bfc(i  ,j,k,2) - bfc(i-1,j,k,2))
           drgt = slope_type_loc*(bfc(i+1,j,k,2) - bfc(i  ,j,k,2))
           dcen = half*(dlft + drgt)/slope_type_loc
           dsgn = sign(one, dcen)
           slop = min(abs(dlft), abs(drgt))
           dlim = slop
           if ((dlft*drgt) <= zero) dlim = zero
           dbfc(i,j,k,2,1) = dsgn*min(dlim,abs(dcen))
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO
#endif
  
#if NDIM == 3
  if (slope_type < 3) then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dsgn, dlim, dcen, dlft, drgt) &
     !$OMP PRIVATE(slop)
     do k = klo, khi
        !$acc loop vector(blocky_slope)
        do j = jlo, jhi
           !$acc loop vector(blockx_slope)
           do i = ilo, ihi
              do ivar = 1, nvar
                 ! slopes in first coordinate direction
                 dlft = slope_type*(qin(i  ,j,k,ivar) - qin(i-1,j,k,ivar))
                 drgt = slope_type*(qin(i+1,j,k,ivar) - qin(i  ,j,k,ivar))
                 dcen = half*(dlft + drgt)/slope_type
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft), abs(drgt))
                 dlim = slop
                 if ((dlft*drgt) <= zero) dlim = zero
                 dq(i,j,k,ivar,1) = dsgn*min(dlim, abs(dcen))
                 ! slopes in second coordinate direction
                 dlft = slope_type*(qin(i,j  ,k,ivar) - qin(i,j-1,k,ivar))
                 drgt = slope_type*(qin(i,j+1,k,ivar) - qin(i,j  ,k,ivar))
                 dcen = half*(dlft + drgt)/slope_type
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft), abs(drgt))
                 dlim = slop
                 if ((dlft*drgt) <= zero) dlim = zero
                 dq(i,j,k,ivar,2) = dsgn*min(dlim, abs(dcen))
                 ! slopes in third coordinate direction
                 dlft = slope_type*(qin(i,j,k  ,ivar) - qin(i,j,k-1,ivar))
                 drgt = slope_type*(qin(i,j,k+1,ivar) - qin(i,j,k  ,ivar))
                 dcen = half*(dlft + drgt)/slope_type
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft), abs(drgt))
                 dlim = slop
                 if ((dlft*drgt) <= zero) dlim = zero
                 dq(i,j,k,ivar,3) = dsgn*min(dlim, abs(dcen))
              end do
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  else if (slope_type == 3) then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dflll, dflml, dflrl, dfmll) &
     !$OMP PRIVATE(dfmml, dfmrl, dfrll, dfrml, dfrrl, dfllm, dflmm, dflrm) &
     !$OMP PRIVATE(dfmlm, dfmmm, dfmrm, dfrlm, dfrmm, dfrrm, dfllr, dflmr) &
     !$OMP PRIVATE(dflrr, dfmlr, dfmmr, dfmrr, dfrlr, dfrmr, dfrrr, vmin) &
     !$OMP PRIVATE(vmax, dfx, dfy, dfz, dff, slop, dlim)
     do k = klo, khi
        !$acc loop vector(blocky_slope)
        do j = jlo, jhi
           !$acc loop vector(blockx_slope)
           do i = ilo, ihi
              do ivar = 1, nvar
                 dflll = qin(i-1,j-1,k-1,ivar) - qin(i,j,k,ivar)
                 dflml = qin(i-1,j  ,k-1,ivar) - qin(i,j,k,ivar)
                 dflrl = qin(i-1,j+1,k-1,ivar) - qin(i,j,k,ivar)
                 dfmll = qin(i  ,j-1,k-1,ivar) - qin(i,j,k,ivar)
                 dfmml = qin(i  ,j  ,k-1,ivar) - qin(i,j,k,ivar)
                 dfmrl = qin(i  ,j+1,k-1,ivar) - qin(i,j,k,ivar)
                 dfrll = qin(i+1,j-1,k-1,ivar) - qin(i,j,k,ivar)
                 dfrml = qin(i+1,j  ,k-1,ivar) - qin(i,j,k,ivar)
                 dfrrl = qin(i+1,j+1,k-1,ivar) - qin(i,j,k,ivar)

                 dfllm = qin(i-1,j-1,k  ,ivar) - qin(i,j,k,ivar)
                 dflmm = qin(i-1,j  ,k  ,ivar) - qin(i,j,k,ivar)
                 dflrm = qin(i-1,j+1,k  ,ivar) - qin(i,j,k,ivar)
                 dfmlm = qin(i  ,j-1,k  ,ivar) - qin(i,j,k,ivar)
                 dfmmm = qin(i  ,j  ,k  ,ivar) - qin(i,j,k,ivar)
                 dfmrm = qin(i  ,j+1,k  ,ivar) - qin(i,j,k,ivar)
                 dfrlm = qin(i+1,j-1,k  ,ivar) - qin(i,j,k,ivar)
                 dfrmm = qin(i+1,j  ,k  ,ivar) - qin(i,j,k,ivar)
                 dfrrm = qin(i+1,j+1,k  ,ivar) - qin(i,j,k,ivar)

                 dfllr = qin(i-1,j-1,k+1,ivar) - qin(i,j,k,ivar)
                 dflmr = qin(i-1,j  ,k+1,ivar) - qin(i,j,k,ivar)
                 dflrr = qin(i-1,j+1,k+1,ivar) - qin(i,j,k,ivar)
                 dfmlr = qin(i  ,j-1,k+1,ivar) - qin(i,j,k,ivar)
                 dfmmr = qin(i  ,j  ,k+1,ivar) - qin(i,j,k,ivar)
                 dfmrr = qin(i  ,j+1,k+1,ivar) - qin(i,j,k,ivar)
                 dfrlr = qin(i+1,j-1,k+1,ivar) - qin(i,j,k,ivar)
                 dfrmr = qin(i+1,j  ,k+1,ivar) - qin(i,j,k,ivar)
                 dfrrr = qin(i+1,j+1,k+1,ivar) - qin(i,j,k,ivar)

                 vmin = min(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml &
                      & ,dfrrl,dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm &
                      & ,dfrmm,dfrrm,dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr &
                      & ,dfrlr,dfrmr,dfrrr)
                 vmax = max(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml &
                      & ,dfrrl,dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm &
                      & ,dfrmm,dfrrm,dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr &
                      & ,dfrlr,dfrmr,dfrrr)

                 dfx  = half*(qin(i+1,j,k,ivar) - qin(i-1,j,k,ivar))
                 dfy  = half*(qin(i,j+1,k,ivar) - qin(i,j-1,k,ivar))
                 dfz  = half*(qin(i,j,k+1,ivar) - qin(i,j,k-1,ivar))
                 dff  = half*(abs(dfx) + abs(dfy) + abs(dfz))

                 if (dff > zero) then
                    slop = min(one,min(abs(vmin),abs(vmax))/dff)
                 else
                    slop = one
                 endif

                 dlim = slop

                 dq(i,j,k,ivar,1) = dlim*dfx
                 dq(i,j,k,ivar,2) = dlim*dfy
                 dq(i,j,k,ivar,3) = dlim*dfz
              enddo
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  else
     print*, "Unknown slope type!"
     stop
  endif
  
  ! Bx along direction Y and Z
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dsgn, dlim, dcen, dlft, drgt, slop)
  do k = klo, khi
     !$acc loop vector(blocky_slope)
     do j = jlo, jhi
        !$acc loop vector(blockx_slope)
        do i = ilo, ihi+1
           ! slopes along Y
           dlft = slope_type_loc*(bfc(i,j  ,k,1) - bfc(i,j-1,k,1))
           drgt = slope_type_loc*(bfc(i,j+1,k,1) - bfc(i,j  ,k,1))
           dcen = half*(dlft + drgt)/slope_type_loc
           dsgn = sign(one, dcen)
           slop = min(abs(dlft), abs(drgt))
           dlim = slop
           if ((dlft*drgt) <= zero) dlim = zero
           dbfc(i,j,k,1,2) = dsgn*min(dlim, abs(dcen))
           ! slopes along Z
           dlft = slope_type_loc*(bfc(i,j,k  ,1) - bfc(i,j,k-1,1))
           drgt = slope_type_loc*(bfc(i,j,k+1,1) - bfc(i,j,k  ,1))
           dcen = half*(dlft + drgt)/slope_type_loc
           dsgn = sign(one, dcen)
           slop = min(abs(dlft), abs(drgt))
           dlim = slop
           if ((dlft*drgt) <= zero) dlim = zero
           dbfc(i,j,k,1,3) = dsgn*min(dlim, abs(dcen))
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! By along direction X and Z
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dsgn, dlim, dcen, dlft, drgt, slop)
  do k = klo, khi
     !$acc loop vector(blocky_slope)
     do j = jlo, jhi+1
        !$acc loop vector(blockx_slope)
        do i = ilo, ihi
           ! slopes along X
           dlft = slope_type_loc*(bfc(i  ,j,k,2) - bfc(i-1,j,k,2))
           drgt = slope_type_loc*(bfc(i+1,j,k,2) - bfc(i  ,j,k,2))
           dcen = half*(dlft + drgt)/slope_type_loc
           dsgn = sign(one, dcen)
           slop = min(abs(dlft), abs(drgt))
           dlim = slop
           if ((dlft*drgt) <= zero) dlim = zero
           dbfc(i,j,k,2,1) = dsgn*min(dlim,abs(dcen))
           ! slopes along Z
           dlft = slope_type_loc*(bfc(i,j,k  ,2) - bfc(i,j,k-1,2))
           drgt = slope_type_loc*(bfc(i,j,k+1,2) - bfc(i,j,k  ,2))
           dcen = half*(dlft + drgt)/slope_type_loc
           dsgn = sign(one, dcen)
           slop = min(abs(dlft), abs(drgt))
           dlim = slop
           if ((dlft*drgt) <= zero) dlim = zero
           dbfc(i,j,k,2,3) = dsgn*min(dlim, abs(dcen))
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Bz along direction X and Y
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dsgn, dlim, dcen, dlft, drgt, slop)
  do k = klo, khi+1
     !$acc loop vector(blocky_slope)
     do j = jlo, jhi
        !$acc loop vector(blockx_slope)
        do i = ilo, ihi
           ! slopes along X
           dlft = slope_type_loc*(bfc(i  ,j,k,3) - bfc(i-1,j,k,3))
           drgt = slope_type_loc*(bfc(i+1,j,k,3) - bfc(i  ,j,k,3))
           dcen = half*(dlft + drgt)/slope_type_loc
           dsgn = sign(one, dcen)
           slop = min(abs(dlft), abs(drgt))
           dlim = slop
           if ((dlft*drgt) <= zero) dlim = zero
           dbfc(i,j,k,3,1) = dsgn*min(dlim, abs(dcen))
           ! slopes along Y
           dlft = slope_type_loc*(bfc(i,j  ,k,3) - bfc(i,j-1,k,3))
           drgt = slope_type_loc*(bfc(i,j+1,k,3) - bfc(i,j  ,k,3))
           dcen = half*(dlft + drgt)/slope_type_loc
           dsgn = sign(one, dcen)
           slop = min(abs(dlft), abs(drgt))
           dlim = slop
           if ((dlft*drgt) <= zero) dlim = zero
           dbfc(i,j,k,3,2) = dsgn*min(dlim, abs(dcen))
        end do
     end do
  end do
  !$OMP END PARALLEL DO
#endif

  return
end subroutine slope
!===============================================================================
#if NDIM == 1
!> Perform a half-timestep prediction
!===============================================================================
subroutine trace1d(dq, qm, qp)
  use params
  use variables, only: qin, dt, x
  use oacc_params
  implicit none

  real(dp), dimension(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,nvar,ndim), intent(in) :: dq
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar,ndim), intent(out) :: qm 
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar,ndim), intent(out) :: qp 

  integer  :: i, j, k, l, n
  integer  :: ilo, ihi, jlo, jhi, klo, khi
  real(dp) :: dtdx, xc
  real(dp) :: r, u, v, w, p, A, B, C
  real(dp) :: drx, dux, dvx, dwx, dpx, dAx, dBx, dCx
  real(dp) :: sr0, su0, sv0, sw0, sp0, sA0, sB0, sC0

  !$py begin_statement

  ilo  = min(1,iu1+1); ihi = max(1,iu2-1)
  jlo  = min(1,ju1+1); jhi = max(1,ju2-1)
  klo  = min(1,ku1+1); khi = max(1,ku2-1)
  su0  = zero; sv0 = zero; sw0 = zero
  dtdx = dt/dx

  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(xc, r, u, v, w, p, A, B, C, drx) &
  !$OMP PRIVATE(dux, dvx, dwx, dpx, dAx, dBx, dCx, sr0, su0, sv0, sw0, sp0) &
  !$OMP PRIVATE(sA0, sB0, sC0)
  do k = klo, khi
     !$acc loop vector(blocky_trace)
     do j = jlo, jhi
        !$acc loop vector(blockx_trace)
        do i = ilo, ihi
           xc   = x(i)
           
           ! Cell centered values
           r = qin(i,j,k,ir)
           u = qin(i,j,k,iu)
           v = qin(i,j,k,iv)
           w = qin(i,j,k,iw)
           p = qin(i,j,k,ip)
           A = qin(i,j,k,iA)
           B = qin(i,j,k,iB)
           C = qin(i,j,k,iC)
           
           ! TVD slopes in X direction
           drx = half*dq(i,j,k,ir,1)
           dux = half*dq(i,j,k,iu,1)
           dvx = half*dq(i,j,k,iv,1)
           dwx = half*dq(i,j,k,iw,1)
           dpx = half*dq(i,j,k,ip,1)
           dBx = half*dq(i,j,k,iB,1)
           dCx = half*dq(i,j,k,iC,1)
           
           ! Source terms (including transverse derivatives)
#if GEOM == CARTESIAN
           sr0 = -u*drx - r*dux                  ; sr0 = sr0*dtdx
           su0 = -u*dux - (dpx + B*dBx + C*dCx)/r; su0 = su0*dtdx  
           sv0 = -u*dvx + (A*dBx)/r              ; sv0 = sv0*dtdx
           sw0 = -u*dwx + (A*dCx)/r              ; sw0 = sw0*dtdx
           sp0 = -u*dpx - gamma*p*dux            ; sp0 = sp0*dtdx
           sB0 = -u*dBx + A*dvx - B*dux          ; sB0 = sB0*dtdx
           sC0 = -u*dCx + A*dwx - C*dux          ; sC0 = sC0*dtdx
#endif
#if GEOM == CYLINDRICAL 
           ! Note that A=0 by definition in this case
           sr0 = -u*drx - r*dux
           sr0 = sr0*dtdx - r*u/xc*half*dt
           su0 = -u*dux - (dpx + B*dBx + C*dCx)/r
           su0 = su0*dtdx + (v*v - B*B/r)/xc*half*dt
           sv0 = -u*dvx
           sv0 = sv0*dtdx - u*v/xc*half*dt
           sw0 = -u*dwx
           sw0 = sw0*dtdx
           sp0 = -u*dpx - gamma*p*dux     
           sp0 = sp0*dtdx - gamma*p*u/xc*half*dt
           sB0 = -u*dBx - B*dux   
           sB0 = sB0*dtdx
           sC0 = -u*dCx - C*dux   
           sC0 = sC0*dtdx - u*C/xc*half*dt
#endif
#if GEOM == SPHERICAL
           ! Note that theta=pi/2 is assumed
           sr0 = -u*drx - r*dux                  
           sr0 = sr0*dtdx - two*r*u
           su0 = -u*dux - (dpx + B*dBx + C*dCx)/r
           su0 = su0*dtdx + (v*v + w*w - B*B/r - C*C/r)/xc*half*dt
           sv0 = -u*dvx + (A*dBx)/r              
           sv0 = sv0*dtdx - (u*v - A*B/r)/xc*half*dt
           sw0 = -u*dwx + (A*dCx)/r              
           sw0 = sw0*dtdx - (u*w - A*C/r)/xc*half*dt
           sp0 = -u*dpx - gamma*p*dux            
           sp0 = sp0*dtdx - two*gamma*p*u/xc*half*dt
           sB0 = -u*dBx + A*dvx - B*dux          
           sB0 = sB0*dtdx - (u*B - v*A)/xc*half*dt
           sC0 = -u*dCx + A*dwx - C*dux          
           sC0 = sC0*dtdx - (w*A - u*C)/xc*half*dt
#endif

           ! Cell-centered predicted states
           r = r + sr0
           u = u + su0
           v = v + sv0
           w = w + sw0
           p = p + sp0
           B = B + sB0
           C = C + sC0

#if GEOM == CARTESIAN 
          if (Omega0 > zero) then
              sC0 = sc0 - 1.5d0*Omega0*A*dtdx
           endif
#endif
           
           ! Right state at left interface
           qp(i,j,k,ir,1) = r - drx
           qp(i,j,k,iu,1) = u - dux
           qp(i,j,k,iv,1) = v - dvx
           qp(i,j,k,iw,1) = w - dwx
           qp(i,j,k,ip,1) = p - dpx
           qp(i,j,k,iA,1) = A
           qp(i,j,k,iB,1) = B - dBx
           qp(i,j,k,iC,1) = C - dCx
           qp(i,j,k,ir,1) = max(smallr, qp(i,j,k,ir,1))
           
           ! Left state at right interface
           qm(i,j,k,ir,1) = r + drx
           qm(i,j,k,iu,1) = u + dux
           qm(i,j,k,iv,1) = v + dvx
           qm(i,j,k,iw,1) = w + dwx
           qm(i,j,k,ip,1) = p + dpx
           qm(i,j,k,iA,1) = A
           qm(i,j,k,iB,1) = B + dBx
           qm(i,j,k,iC,1) = C + dCx
           qm(i,j,k,ir,1) = max(smallr, qm(i,j,k,ir,1))
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

  return
end subroutine trace1d
#endif
!===============================================================================
#if NDIM == 2
!> Perform a half-timestep prediction
!===============================================================================
subroutine trace2d(bfc, dq, dbfc, qm, qp, qRT, qRB, qLT, qLB)
  use params
  use variables, only: qin, dt, x, y, Ez
  use oacc_params
  implicit none

  real(dp), dimension(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,3), intent(in) :: bfc
  real(dp), dimension(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,nvar,ndim), intent(in) :: dq
  real(dp), dimension(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,3,ndim), intent(in) :: dbfc
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar,ndim), intent(out) :: qm 
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar,ndim), intent(out) :: qp 
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar,3), intent(out) :: qRT
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar,3), intent(out) :: qRB
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar,3), intent(out) :: qLT
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar,3), intent(out) :: qLB

  integer  :: i, j, k, l
  integer  :: ilo, ihi, jlo, jhi, klo, khi
  real(dp) :: dtdx, dtdy, xc, xL, xR, yL, yR, sinxc, sinxL, sinxR, cotanxc
  real(dp) :: smallp, shear
  real(dp) :: r, u, v, w, p, A, B, C
  real(dp) :: ELL, ELR, ERL, ERR
  real(dp) :: drx, dux, dvx, dwx, dpx, dAx, dBx, dCx
  real(dp) :: dry, duy, dvy, dwy, dpy, dAy, dBy, dCy
  real(dp) :: sr0, su0, sv0, sw0, sp0, sA0, sB0, sC0
  real(dp) :: AL, AR, BL, BR
  real(dp) :: dALy, dARy
  real(dp) :: dBLx, dBRx
  real(dp) :: sAL0, sAR0, sBL0, sBR0

  !$py begin_statement

  smallp = smallr*smallc**2/gamma
  ilo  = min(1,iu1+1); ihi = max(1,iu2-1)
  jlo  = min(1,ju1+1); jhi = max(1,ju2-1)
  klo  = min(1,ku1+1); khi = max(1,ku2-1)
  dtdx = dt/dx
  dtdy = dt/dy

  !$acc data pcreate(Ez)
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(u, v, A, B)
  do k = klo, ku2
     !$acc loop vector(blocky_trace)
     do j = jlo, ju2
        !$acc loop vector(blockx_trace)
        do i = ilo, iu2
           u = 0.25*(qin(i-1,j-1,k,iu) + qin(i-1,j,k,iu) + qin(i,j-1,k,iu) &
                & + qin(i,j,k,iu))
           v = 0.25*(qin(i-1,j-1,k,iv) + qin(i-1,j,k,iv) + qin(i,j-1,k,iv) &
                & + qin(i,j,k,iv))
           A = 0.5*(bfc(i,j-1,k,1) + bfc(i,j,k,1))
           B = 0.5*(bfc(i-1,j,k,2) + bfc(i,j,k,2))
           Ez(i,j,k) = u*B - v*A
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(xL, xR, xc, yL, yR, sinxc, sinxL)&
  !$OMP PRIVATE(sinxR, cotanxc, r, u, v, w, p, A, B, C, ELL, ELR, ERL, ERR) &
  !$OMP PRIVATE(drx, dux, dvx, dwx, dpx, dAx, dBx, dCx, dry, duy, dvy, dwy) &
  !$OMP PRIVATE(dpy, dAy, dBy, dCy, sr0, su0, sv0, sw0, sp0, sA0, sB0, sC0, AL)&
  !$OMP PRIVATE(AR, BL, BR, dALy, dARy, dBLx, dBRx, sAL0, sAR0, sBL0, sBR0)
  do k = klo, khi
     !$acc loop vector(blocky_trace)
     do j = jlo, jhi
#if GEOM == SPHERICAL
        yL = half*(y(j-1) + y(j))
        yR = half*(y(j+1) + y(j))
        sinxc = sin(y(j))
        sinxL = sin(yL)
        sinxR = sin(yR)
        cotanxc = cos(y(j))/sinxc
#endif
        !$acc loop vector(blockx_trace)
        do i = ilo, ihi
           xL = half*(x(i-1)+x(i))
           xR = half*(x(i+1)+x(i))
           xc = x(i)

           ! Cell centered values
           r = qin(i,j,k,ir)
           u = qin(i,j,k,iu)
           v = qin(i,j,k,iv)
           w = qin(i,j,k,iw)            
           p = qin(i,j,k,ip)
           A = qin(i,j,k,iA)
           B = qin(i,j,k,iB)
           C = qin(i,j,k,iC)            

           ! Face centered variables
           AL =  bfc(i  ,j  ,k,1)
           AR =  bfc(i+1,j  ,k,1)
           BL =  bfc(i  ,j  ,k,2)
           BR =  bfc(i  ,j+1,k,2)

           ! Electric field
           ELL = Ez(i  ,j  ,k)
           ELR = Ez(i  ,j+1,k)
           ERL = Ez(i+1,j  ,k)
           ERR = Ez(i+1,j+1,k)

           ! Cell centered TVD slopes in X direction
           drx = half*dq(i,j,k,ir,1)
           dux = half*dq(i,j,k,iu,1)
           dvx = half*dq(i,j,k,iv,1)
           dwx = half*dq(i,j,k,iw,1)
           dCx = half*dq(i,j,k,iC,1)
           dpx = half*dq(i,j,k,ip,1)
           dBx = half*dq(i,j,k,iB,1)

           ! Cell centered TVD slopes in Y direction
           dry = half*dq(i,j,k,ir,2)
           duy = half*dq(i,j,k,iu,2)
           dvy = half*dq(i,j,k,iv,2)
           dwy = half*dq(i,j,k,iw,2)
           dCy = half*dq(i,j,k,iC,2)
           dpy = half*dq(i,j,k,ip,2)
           dAy = half*dq(i,j,k,iA,2)

           ! Face centered TVD slopes in transverse direction
           dALy = half*dbfc(i  ,j  ,k,1,2)
           dARy = half*dbfc(i+1,j  ,k,1,2)
           dBLx = half*dbfc(i  ,j  ,k,2,1)
           dBRx = half*dbfc(i  ,j+1,k,2,1)

           ! Cell centered slopes in normal direction
           dAx = half*(AR - AL)
           dBy = half*(BR - BL)

           ! Source terms (including transverse derivatives)
#if GEOM == CARTESIAN
           sr0 = (-u*drx - dux*r)*dtdx + (-v*dry - dvy*r)*dtdy
           su0 = (-u*dux - dpx/r - B*dBx/r - C*dCx/r)*dtdx &
               & + (-v*duy + B*dAy/r)*dtdy
           sv0 = (-u*dvx + A*dBx/r)*dtdx &
               & + (-v*dvy - dpy/r - A*dAy/r - C*dCy/r)*dtdy
           sw0 = (-u*dwx + A*dCx/r)*dtdx + (-v*dwy + B*dCy/r)*dtdy
           sp0 = (-u*dpx - dux*gamma*p)*dtdx + (-v*dpy - dvy*gamma*p)*dtdy
           sA0 = (u*dBy + B*duy - v*dAy - A*dvy)*dtdy
           sB0 = (-u*dBx - B*dux + v*dAx + A*dvx)*dtdx 
           sC0 = (w*dAx + A*dwx - u*dCx - C*dux)*dtdx &
               & + (-v*dCy - C*dvy + w*dBy + B*dwy)*dtdy

           if (Omega0 > zero) then
              shear = -1.5d0*Omega0*x(i)
              sC0 = sC0 + (shear*dAx - 1.5d0*Omega0*A)*dtdx + shear*dBy*dtdy
           endif
           
           ! Face centered B-field
           sAL0 = +(ELR - ELL)*half*dtdy
           sAR0 = +(ERR - ERL)*half*dtdy
           sBL0 = -(ERL - ELL)*half*dtdx
           sBR0 = -(ERR - ELR)*half*dtdx
#endif

#if GEOM == CYLINDRICAL
           sr0 = (-u*drx - dux*r)*dtdx + (-v*dry - dvy*r)/xc*dtdy
           su0 = (-u*dux - dpx/r - B*dBx/r - C*dCx/r)*dtdx &
               & + (-v*duy + B*dAy/r)/xc*dtdy
           sv0 = (-u*dvx + A*dBx/r)*dtdx &
               & + (-v*dvy - dpy/r - A*dAy/r - C*dCy/r)/xc*dtdy
           sw0 = (-u*dwx + A*dCx/r)*dtdx + (-v*dwy + B*dCy/r)/xc*dtdy
           sp0 = (-u*dpx - dux*gamma*p)*dtdx + (-v*dpy - dvy*gamma*p)/xc*dtdy
           sA0 = (u*dBy + B*duy - v*dAy - A*dvy)/xc*dtdy
           sB0 = (-u*dBx - B*dux + v*dAx + A*dvx)*dtdx 
           sC0 = (w*dAx + A*dwx - u*dCx - C*dux)*dtdx &
               & + (-v*dCy - C*dvy + w*dBy + B*dwy)/xc*dtdy

           ! Geometrical terms
           sr0 = sr0 - r*u/xc*half*dt
           su0 = su0 + (v*v - B*B/r)/xc*half*dt
           sv0 = sv0 + (-u*v + A*B/r)/xc*half*dt
           sp0 = sp0 - gamma*p*u/xc*half*dt
           sC0 = sC0 + (w*A - u*C)/xc*half*dt

           ! Face centered B-field
           sAL0 = (ELR - ELL)/xL*half*dtdy
           sAR0 = (ERR - ERL)/xR*half*dtdy
           sBL0 = (ELL - ERL)*half*dtdx
           sBR0 = (ERR - ELR)*half*dtdx
#endif

#if GEOM == SPHERICAL 
           sr0 = (-u*drx - dux*r)*dtdx + (-v*dry - dvy*r)/xc*dtdy
           su0 = (-u*dux - dpx/r - B*dBx/r - C*dCx/r)*dtdx &
               & + (-v*duy + B*dAy/r)/xc*dtdy
           sv0 = (-u*dvx + A*dBx/r)*dtdx &
               & + (-v*dvy - dpy/r - A*dAy/r - C*dCy/r)/xc*dtdy
           sw0 = (-u*dwx + A*dCx/r)*dtdx + (-v*dwy + B*dCy/r)/xc*dtdy
           sp0 = (-u*dpx - dux*gamma*p)*dtdx + (-v*dpy - dvy*gamma*p)/xc*dtdy
           sA0 = (u*dBy + B*duy - v*dAy - A*dvy)/xc*dtdy
           sB0 = (-u*dBx - B*dux + v*dAx + A*dvx)*dtdx 
           sC0 = (w*dAx + A*dwx - u*dCx - C*dux)*dtdx &
               & + (-v*dCy - C*dvy + w*dBy + B*dwy)/xc*dtdy
           
           ! Geometrical terms
           sr0 = sr0 - (two*u + cotanxc*v)*r/xc*half*dt
           su0 = su0 + (v*v + w*w - (B*B + C*C)/r)/xc*half*dt
           sv0 = sv0 + (-u*v + A*B/r + cotanxc*(w*w - C*C/r))/xc*half*dt
           sw0 = sw0 + (-u*w + A*C/r + cotanxc*(v*w - B*C/r))/xc*half*dt
           sp0 = sp0 - gamma*p*(two*u + cotanxc*v)/xc*half*dt
           sA0 = sA0 + cotanxc*(u*B - v*A)/xc*half*dt
           sB0 = sB0 - (u*B - v*A)/xc*half*dt
           sC0 = sC0 + (w*A - u*C)/xc*half*dt
           
           ! Face centered B-field
           sAL0 = (ELR*sinxR - ELL*sinxL)/xL/sinxc*half*dtdy
           sAR0 = (ERR*sinxR - ERL*sinxL)/xR/sinxc*half*dtdy
           sBL0 = -(ERL*xR - ELL*xL)/xc*half*dtdx
           sBR0 = -(ERR*xR - ELR*xL)/xc*half*dtdx
#endif

           ! Update in time the  primitive variables
           r = r + sr0
           u = u + su0
           v = v + sv0
           w = w + sw0
           p = p + sp0
           A = A + sA0
           B = B + sB0
           C = C + sC0

           AL = AL + sAL0
           AR = AR + sAR0
           BL = BL + sBL0
           BR = BR + sBR0

           ! Right state at left interface
           qp(i,j,k,ir,1) = r - drx  
           qp(i,j,k,iu,1) = u - dux
           qp(i,j,k,iv,1) = v - dvx
           qp(i,j,k,iw,1) = w - dwx
           qp(i,j,k,ip,1) = p - dpx
           qp(i,j,k,iA,1) = AL
           qp(i,j,k,iB,1) = B - dBx
           qp(i,j,k,iC,1) = C - dCx
           qp(i,j,k,ir,1) = max(smallr, qp(i,j,k,ir,1))
           qp(i,j,k,ip,1) = max(smallp, qp(i,j,k,ip,1))

           ! Left state at right interface
           qm(i,j,k,ir,1) = r + drx
           qm(i,j,k,iu,1) = u + dux
           qm(i,j,k,iv,1) = v + dvx
           qm(i,j,k,iw,1) = w + dwx
           qm(i,j,k,ip,1) = p + dpx  
           qm(i,j,k,iA,1) = AR
           qm(i,j,k,iB,1) = B + dBx
           qm(i,j,k,iC,1) = C + dCx
           qm(i,j,k,ir,1) = max(smallr, qm(i,j,k,ir,1))
           qm(i,j,k,ip,1) = max(smallp, qm(i,j,k,ip,1))

           ! Top state at bottom interface
           qp(i,j,k,ir,2) = r - dry  
           qp(i,j,k,iu,2) = u - duy 
           qp(i,j,k,iv,2) = v - dvy 
           qp(i,j,k,iw,2) = w - dwy 
           qp(i,j,k,ip,2) = p - dpy  
           qp(i,j,k,iA,2) = A - dAy
           qp(i,j,k,iB,2) = BL
           qp(i,j,k,iC,2) = C - dCy
           qp(i,j,k,ir,2) = max(smallr, qp(i,j,k,ir,2))
           qp(i,j,k,ip,2) = max(smallp, qp(i,j,k,ip,2))

           ! Bottom state at top interface
           qm(i,j,k,ir,2) = r + dry
           qm(i,j,k,iu,2) = u + duy
           qm(i,j,k,iv,2) = v + dvy
           qm(i,j,k,iw,2) = w + dwy
           qm(i,j,k,ip,2) = p + dpy
           qm(i,j,k,iA,2) = A + dAy
           qm(i,j,k,iB,2) = BR
           qm(i,j,k,iC,2) = C + dCy
           qm(i,j,k,ir,2) = max(smallr, qm(i,j,k,ir,2))
           qm(i,j,k,ip,2) = max(smallp, qm(i,j,k,ip,2))

           ! Right-top state (RT->LL)
           qRT(i,j,k,ir,3) = r + (+drx+dry)
           qRT(i,j,k,iu,3) = u + (+dux+duy)
           qRT(i,j,k,iv,3) = v + (+dvx+dvy)
           qRT(i,j,k,iw,3) = w + (+dwx+dwy)
           qRT(i,j,k,ip,3) = p + (+dpx+dpy)
           qRT(i,j,k,iA,3) = AR+ (   +dARy)
           qRT(i,j,k,iB,3) = BR+ (+dBRx   )
           qRT(i,j,k,iC,3) = C + (+dCx+dCy)
           qRT(i,j,k,ir ,3) = max(smallr, qRT(i,j,k,ir,3))
           qRT(i,j,k,ip ,3) = max(smallp, qRT(i,j,k,ip,3))

           ! Right-Bottom state (RB->LR)
           qRB(i,j,k,ir,3) = r + (+drx-dry)
           qRB(i,j,k,iu,3) = u + (+dux-duy)
           qRB(i,j,k,iv,3) = v + (+dvx-dvy)
           qRB(i,j,k,iw,3) = w + (+dwx-dwy)
           qRB(i,j,k,ip,3) = p + (+dpx-dpy)
           qRB(i,j,k,iA,3) = AR+ (   -dARy)
           qRB(i,j,k,iB,3) = BL+ (+dBLx   )
           qRB(i,j,k,iC,3) = C + (+dCx-dCy)
           qRB(i,j,k,ir,3) = max(smallr, qRB(i,j,k,ir,3))
           qRB(i,j,k,ip,3) = max(smallp, qRB(i,j,k,ip,3))

           ! Left-Bottom state (LB->RR)
           qLB(i,j,k,ir,3) = r + (-drx-dry)
           qLB(i,j,k,iu,3) = u + (-dux-duy)
           qLB(i,j,k,iv,3) = v + (-dvx-dvy)
           qLB(i,j,k,iw,3) = w + (-dwx-dwy)
           qLB(i,j,k,ip,3) = p + (-dpx-dpy)
           qLB(i,j,k,iA,3) = AL+ (   -dALy)
           qLB(i,j,k,iB,3) = BL+ (-dBLx   )
           qLB(i,j,k,iC,3) = C + (-dCx-dCy)
           qLB(i,j,k,ir,3) = max(smallr, qLB(i,j,k,ir,3))
           qLB(i,j,k,ip,3) = max(smallp, qLB(i,j,k,ip,3))

           ! Left-Top state (LT->RL)
           qLT(i,j,k,ir,3) = r + (-drx+dry)
           qLT(i,j,k,iu,3) = u + (-dux+duy)
           qLT(i,j,k,iv,3) = v + (-dvx+dvy)
           qLT(i,j,k,iw,3) = w + (-dwx+dwy)
           qLT(i,j,k,ip,3) = p + (-dpx+dpy)
           qLT(i,j,k,iA,3) = AL+ (   +dALy)
           qLT(i,j,k,iB,3) = BR+ (-dBRx   )
           qLT(i,j,k,iC,3) = C + (-dCx+dCy)
           qLT(i,j,k,ir,3) = max(smallr, qLT(i,j,k,ir,3))
           qLT(i,j,k,ip,3) = max(smallp, qLT(i,j,k,ip,3))
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  !$acc end data

  return
end subroutine trace2d
#endif
!===============================================================================
#if NDIM == 3
!> Perform a half-timestep prediction
!===============================================================================
subroutine trace3d(bfc, dq, dbfc, qm, qp, qRT, qRB, qLT, qLB)
  use params
  use variables, only: qin, dt, x, y, Ex, Ey, Ez
  use oacc_params
  implicit none

  real(dp), dimension(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,3), intent(in) :: bfc
  real(dp), dimension(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,nvar,ndim), intent(in) :: dq
  real(dp), dimension(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,3,ndim), intent(in) :: dbfc
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar,ndim), intent(out) :: qm 
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar,ndim), intent(out) :: qp 
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar,3), intent(out) :: qRT
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar,3), intent(out) :: qRB
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar,3), intent(out) :: qLT
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar,3), intent(out) :: qLB

  integer  :: i, j, k
  integer  :: ilo, ihi, jlo, jhi, klo, khi
  real(dp) :: dtdx, dtdy, dtdz, xc, xL, xR, yL, yR, sinxc, sinxL, sinxR, cotanxc
  real(dp) :: smallp, shear
  real(dp) :: r, u, v, w, p, A, B, C
  real(dp) :: ELL, ELR, ERL, ERR
  real(dp) :: FLL, FLR, FRL, FRR
  real(dp) :: GLL, GLR, GRL, GRR
  real(dp) :: drx, dux, dvx, dwx, dpx, dAx, dBx, dCx
  real(dp) :: dry, duy, dvy, dwy, dpy, dAy, dBy, dCy
  real(dp) :: drz, duz, dvz, dwz, dpz, dAz, dBz, dCz
  real(dp) :: sr0, su0, sv0, sw0, sp0, sA0, sB0, sC0
  real(dp) :: AL, AR, BL, BR, CL, CR
  real(dp) :: dALy, dARy, dALz, dARz
  real(dp) :: dBLx, dBRx, dBLz, dBRz
  real(dp) :: dCLx, dCRx, dCLy, dCRy
  real(dp) :: sAL0, sAR0, sBL0, sBR0, sCL0, sCR0

  !$py begin_statement

  smallp = smallr*smallc**2/gamma
  ilo  = min(1,iu1+1); ihi = max(1,iu2-1)
  jlo  = min(1,ju1+1); jhi = max(1,ju2-1)
  klo  = min(1,ku1+1); khi = max(1,ku2-1)
  dtdx = dt/dx
  dtdy = dt/dy
  dtdz = dt/dz

  !$acc data pcreate(Ex, Ey, Ez)
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(u, v, w, A, B, C, shear)
  do k = klo, ku2
     !$acc loop vector(blocky_trace)
     do j = jlo, ju2
        !$acc loop vector(blockx_trace)
        do i = ilo, iu2
           v = forth*(qin(i,j-1,k-1,iv) + qin(i,j-1,k,iv) + qin(i,j,k-1,iv) &
                & + qin(i,j,k,iv))
           w = forth*(qin(i,j-1,k-1,iw) + qin(i,j-1,k,iw) + qin(i,j,k-1,iw) &
                & + qin(i,j,k,iw))
           B = half*(bfc(i,j,k-1,2) + bfc(i,j,k,2))
           C = half*(bfc(i,j-1,k,3) + bfc(i,j,k,3))
           Ex(i,j,k) = v*C - w*B
#if GEOM == CARTESIAN
           if ((Omega0 > zero) .and. (.not. fargo)) then
              shear = -1.5d0*Omega0*x(i)
              Ex(i,j,k) = Ex(i,j,k) + shear*C
           endif
#endif
           u = forth*(qin(i-1,j,k-1,iu) + qin(i-1,j,k,iu) + qin(i,j,k-1,iu) &
                & + qin(i,j,k,iu))
           w = forth*(qin(i-1,j,k-1,iw) + qin(i-1,j,k,iw) + qin(i,j,k-1,iw) &
                & + qin(i,j,k,iw))
           A = half*(bfc(i,j,k-1,1) + bfc(i,j,k,1))
           C = half*(bfc(i-1,j,k,3) + bfc(i,j,k,3))
           Ey(i,j,k) = w*A - u*C
           
           u = forth*(qin(i-1,j-1,k,iu) + qin(i-1,j,k,iu) + qin(i,j-1,k,iu) &
                & + qin(i,j,k,iu))
           v = forth*(qin(i-1,j-1,k,iv) + qin(i-1,j,k,iv) + qin(i,j-1,k,iv) &
                & + qin(i,j,k,iv))
           A = half*(bfc(i,j-1,k,1) + bfc(i,j,k,1))
           B = half*(bfc(i-1,j,k,2) + bfc(i,j,k,2))
           Ez(i,j,k) = u*B - v*A
#if GEOM == CARTESIAN
           if ((Omega0 > zero) .and. (.not. fargo)) then
              shear = -1.5d0*Omega0*half*(x(i) + x(i-1))
              Ez(i,j,k) = Ez(i,j,k) - shear*A
           endif
#endif
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(xL, xR, xc, yL, yR, sinxc, sinxL)&
  !$OMP PRIVATE(sinxR, cotanxc, r, u, v, w, p, A, B, C, ELL, ELR, ERL, ERR) &
  !$OMP PRIVATE(FLL, FLR, FRL, FRR, GLL, GLR, GRL, GRR, drx, dux, dvx, dwx) &
  !$OMP PRIVATE(dpx, dAx, dBx, dCx, dry, duy, dvy, dwy, dpy, dAy, dBy, dCy) &
  !$OMP PRIVATE(drz, duz, dvz, dwz, dpz, dAz, dBz, dCz, sr0, su0, sv0, sw0) &
  !$OMP PRIVATE(sp0, sA0, sB0, sC0, AL, AR, BL, BR, CL, CR, dALy, dARy, dALz) &
  !$OMP PRIVATE(dARz, dBLx, dBRx, dBLz, dBRz, dCLx, dCRx, dCLy, dCRy, sAL0) &
  !$OMP PRIVATE(sAR0, sBL0, sBR0, sCL0, sCR0)
  do k = klo, khi
     !$acc loop vector(blocky_trace)
     do j = jlo, jhi
#if GEOM == SPHERICAL
        yL = half*(y(j-1) + y(j))
        yR = half*(y(j+1) + y(j))
        sinxc = sin(y(j))
        sinxL = sin(yL)
        sinxR = sin(yR)
        cotanxc = cos(y(j))/sinxc
#endif
        !$acc loop vector(blockx_trace)
        do i = ilo, ihi
           xL = half*(x(i-1) + x(i))
           xR = half*(x(i+1) + x(i))
           xc = x(i)
           
           ! Cell centered values
           r = qin(i,j,k,ir)
           u = qin(i,j,k,iu)
           v = qin(i,j,k,iv)
           w = qin(i,j,k,iw)            
           p = qin(i,j,k,ip)
           A = qin(i,j,k,iA)
           B = qin(i,j,k,iB)
           C = qin(i,j,k,iC)            

           ! Cell centered TVD slopes in X direction
           drx = half*dq(i,j,k,ir,1)
           dux = half*dq(i,j,k,iu,1)
           dvx = half*dq(i,j,k,iv,1)
           dwx = half*dq(i,j,k,iw,1)
           dpx = half*dq(i,j,k,ip,1)
           dBx = half*dq(i,j,k,iB,1)
           dCx = half*dq(i,j,k,iC,1)

           ! Cell centered TVD slopes in Y direction
           dry = half*dq(i,j,k,ir,2)
           duy = half*dq(i,j,k,iu,2)
           dvy = half*dq(i,j,k,iv,2)
           dwy = half*dq(i,j,k,iw,2)
           dpy = half*dq(i,j,k,ip,2)
           dAy = half*dq(i,j,k,iA,2)
           dCy = half*dq(i,j,k,iC,2)

           ! Cell centered TVD slopes in Z direction
           drz = half*dq(i,j,k,ir,3)
           duz = half*dq(i,j,k,iu,3)
           dvz = half*dq(i,j,k,iv,3)
           dwz = half*dq(i,j,k,iw,3)
           dpz = half*dq(i,j,k,ip,3)
           dAz = half*dq(i,j,k,iA,3)
           dBz = half*dq(i,j,k,iB,3)

           ! Face centered variables
           AL = bfc(i  ,j  ,k  ,1)
           AR = bfc(i+1,j  ,k  ,1)
           BL = bfc(i  ,j  ,k  ,2)
           BR = bfc(i  ,j+1,k  ,2)
           CL = bfc(i  ,j  ,k  ,3)
           CR = bfc(i  ,j  ,k+1,3)

           ! Cell centered slopes in normal direction
           dAx = half*(AR - AL)
           dBy = half*(BR - BL)
           dCz = half*(CR - CL)

           ! Face centered TVD slopes in transverse direction
           dALy = half*dbfc(i  ,j  ,k  ,1,2)
           dARy = half*dbfc(i+1,j  ,k  ,1,2)
           dALz = half*dbfc(i  ,j  ,k  ,1,3)
           dARz = half*dbfc(i+1,j  ,k  ,1,3)

           dBLx = half*dbfc(i  ,j  ,k  ,2,1)
           dBRx = half*dbfc(i  ,j+1,k  ,2,1)
           dBLz = half*dbfc(i  ,j  ,k  ,2,3)
           dBRz = half*dbfc(i  ,j+1,k  ,2,3)

           dCLx = half*dbfc(i  ,j  ,k  ,3,1)
           dCRx = half*dbfc(i  ,j  ,k+1,3,1)
           dCLy = half*dbfc(i  ,j  ,k  ,3,2)
           dCRy = half*dbfc(i  ,j  ,k+1,3,2)

           ! Edge centered electric field in X, Y and Z directions
           ELL = Ex(i,j  ,k  )
           ELR = Ex(i,j  ,k+1)
           ERL = Ex(i,j+1,k  )
           ERR = Ex(i,j+1,k+1)

           FLL = Ey(i  ,j,k  )
           FLR = Ey(i  ,j,k+1)
           FRL = Ey(i+1,j,k  )
           FRR = Ey(i+1,j,k+1)

           GLL = Ez(i  ,j  ,k)
           GLR = Ez(i  ,j+1,k)
           GRL = Ez(i+1,j  ,k)
           GRR = Ez(i+1,j+1,k)

           ! Source terms (including transverse derivatives)
#if GEOM == CARTESIAN
           sr0 = (-u*drx - dux*r)*dtdx + (-v*dry - dvy*r)*dtdy &
               & + (-w*drz - dwz*r)*dtdz 
           su0 = (-u*dux - (dpx + B*dBx + C*dCx)/r)*dtdx &
               & + (-v*duy + B*dAy/r)*dtdy &
               & + (-w*duz + C*dAz/r)*dtdz 
           sv0 = (-u*dvx + A*dBx/r)*dtdx + (-v*dvy &
               & - (dpy + A*dAy + C*dCy)/r)*dtdy &
               & + (-w*dvz + C*dBz/r)*dtdz
           sw0 = (-u*dwx + A*dCx/r)*dtdx + (-v*dwy + B*dCy/r)*dtdy &
               & + (-w*dwz - (dpz + A*dAz + B*dBz)/r)*dtdz 
           sp0 = (-u*dpx - dux*gamma*p)*dtdx + (-v*dpy - dvy*gamma*p)*dtdy &
               & + (-w*dpz - dwz*gamma*p)*dtdz
           sA0 = (u*dBy + B*duy - v*dAy - A*dvy)*dtdy &
               & + (u*dCz + C*duz - w*dAz - A*dwz)*dtdz
           sB0 = (v*dAx + A*dvx - u*dBx - B*dux)*dtdx &
               & + (v*dCz + C*dvz - w*dBz - B*dwz)*dtdz
           sC0 = (w*dAx + A*dwx - u*dCx - C*dux)*dtdx &
               & + (w*dBy + B*dwy - v*dCy - C*dvy)*dtdy
           
           if ((Omega0 > zero) .and. (.not. fargo)) then
              shear = -1.5d0*Omega0*x(i)
              sr0 = sr0 - shear*dry*dtdy
              su0 = su0 - shear*duy*dtdy
              sv0 = sv0 - shear*dvy*dtdy
              sw0 = sw0 - shear*dwy*dtdy
              sp0 = sp0 - shear*dpy*dtdy
              sA0 = sA0 - shear*dAy*dtdy
              sB0 = sB0 + (shear*dAx - 1.5d0*Omega0*A*dx)*dtdx &
                   & + shear*dBz*dtdz
              sC0 = sC0 -  shear*dCy*dtdy
           endif

           ! Face-centered B-field
           sAL0 = +(GLR - GLL)*dtdy*half - (FLR - FLL)*dtdz*half
           sAR0 = +(GRR - GRL)*dtdy*half - (FRR - FRL)*dtdz*half
           sBL0 = -(GRL - GLL)*dtdx*half + (ELR - ELL)*dtdz*half
           sBR0 = -(GRR - GLR)*dtdx*half + (ERR - ERL)*dtdz*half
           sCL0 = +(FRL - FLL)*dtdx*half - (ERL - ELL)*dtdy*half
           sCR0 = +(FRR - FLR)*dtdx*half - (ERR - ELR)*dtdy*half
#endif

#if GEOM == CYLINDRICAL
           sr0 = (-u*drx - dux*r)*dtdx + (-v*dry - dvy*r)/xc*dtdy &
                & + (-w*drz - dwz*r)*dtdz 
           su0 = (-u*dux - (dpx + B*dBx + C*dCx)/r)*dtdx &
                & + (-v*duy + B*dAy/r)/xc*dtdy + (-w*duz + C*dAz/r)*dtdz 
           sv0 = (-u*dvx + A*dBx/r)*dtdx &
                & + (-v*dvy - (dpy + A*dAy + C*dCy)/r)/xc*dtdy &
                & + (-w*dvz + C*dBz/r)*dtdz
           sw0 = (-u*dwx + A*dCx/r)*dtdx + (-v*dwy + B*dCy/r)/xc*dtdy &
                & + (-w*dwz - (dpz + A*dAz + B*dBz)/r)*dtdz 
           sp0 = (-u*dpx - dux*gamma*p)*dtdx &
                & + (-v*dpy - dvy*gamma*p)/xc*dtdy &
                & + (-w*dpz - dwz*gamma*p)*dtdz
           sA0 = (u*dBy + B*duy - v*dAy - A*dvy)/xc*dtdy &
                & + (u*dCz + C*duz - w*dAz - A*dwz)*dtdz
           sB0 = (v*dAx + A*dvx - u*dBx - B*dux)*dtdx &
                & + (v*dCz + C*dvz - w*dBz - B*dwz)*dtdz 
           sC0 = (w*dAx + A*dwx - u*dCx - C*dux)*dtdx &
                & + (w*dBy + B*dwy - v*dCy - C*dvy)/xc*dtdy

           ! Geometrical terms
           sr0 = sr0 - r*u/xc*half*dt
           su0 = su0 + (v*v - B*B/r)/xc*half*dt
           sv0 = sv0 + (-u*v + A*B/r)/xc*half*dt
           sp0 = sp0 - gamma*p*u/xc*half*dt
           sC0 = sC0 + (w*A - u*C)/xc*half*dt

           ! Face-centered B-field
           sAL0 = (GLR - GLL)/xL*dtdy*half -(FLR - FLL)*dtdz*half
           sAR0 = (GRR - GRL)/xR*dtdy*half -(FRR - FRL)*dtdz*half
           sBL0 = -(GRL - GLL)*dtdx*half +(ELR - ELL)*dtdz*half
           sBR0 = -(GRR - GLR)*dtdx*half +(ERR - ERL)*dtdz*half
           sCL0 = +(FRL*xR - FLL*xL)/xc*dtdx*half - (ERL - ELL)/xc*dtdy*half
           sCR0 = +(FRR*xR - FLR*xL)/xc*dtdx*half - (ERR - ELR)/xc*dtdy*half
#endif

#if GEOM == SPHERICAL
           sr0 = (-u*drx - dux*r)*dtdx + (-v*dry - dvy*r)/xc*dtdy &
                & + (-w*drz - dwz*r)/xc/sinxc*dtdz 
           su0 = (-u*dux - (dpx + B*dBx + C*dCx)/r)*dtdx &
                & + (-v*duy + B*dAy/r)/xc*dtdy &
                & + (-w*duz + C*dAz/r)/xc/sinxc*dtdz 
           sv0 = (-u*dvx + A*dBx/r)*dtdx &
                & + (-v*dvy - (dpy + A*dAy + C*dCy)/r)/xc*dtdy &
                & + (-w*dvz + C*dBz/r)/xc/sinxc*dtdz
           sw0 = (-u*dwx + A*dCx/r)*dtdx + (-v*dwy + B*dCy/r)/xc*dtdy &
                & + (-w*dwz - (dpz + A*dAz + B*dBz)/r)/xc/sinxc*dtdz 
           sp0 = (-u*dpx - dux*gamma*p)*dtdx &
                & + (-v*dpy - dvy*gamma*p)/xc*dtdy &
                & + (-w*dpz - dwz*gamma*p)/xc/sinxc*dtdz
           sA0 = (u*dBy + B*duy - v*dAy - A*dvy)/xc*dtdy &
                & + (u*dCz + C*duz-w*dAz - A*dwz)/xc/sinxc*dtdz
           sB0 = (v*dAx + A*dvx - u*dBx - B*dux)*dtdx &
                & + (v*dCz + C*dvz - w*dBz - B*dwz)/xc/sinxc*dtdz 
           sC0 = (w*dAx + A*dwx - u*dCx - C*dux)*dtdx &
                & + (w*dBy + B*dwy - v*dCy - C*dvy)/xc*dtdy

           ! Geometrical terms
           sr0 = sr0 - (two*u + cotanxc*v)*r/xc*half*dt
           su0 = su0 + (v*v + w*w - (B*B + C*C)/r)/xc*half*dt
           sv0 = sv0 + (-u*v + A*B/r + cotanxc*(w*w - C*C/r))/xc*half*dt
           sw0 = sw0 + (-u*w + A*C/r + cotanxc*(v*w - B*C/r))/xc*half*dt
           sp0 = sp0 - gamma*p*(two*u + cotanxc*v)/xc*half*dt
           sA0 = sA0 + cotanxc*(u*B - v*A)/xc*half*dt
           sB0 = sB0 -(u*B - v*A)/xc*half*dt
           sC0 = sC0 + (w*A - u*C)/xc*half*dt

           ! Face-centered B-field
           sAL0 = (GLR*sinxR - GLL*sinxL)/xL/sinxc*half*dtdy &
                & - (FLR - FLL)/xc/sinxc*half*dtdz
           sAR0 = (GRR*sinxR - GRL*sinxL)/xR/sinxc*half*dtdy &
                & - (FRR - FRL)/xc/sinxc*half*dtdz
           sBL0 = -(GRL*xR - GLL*xL)/xc*half*dtdx &
                & + (ELR - ELL)/xc/sinxc*half*dtdz
           sBR0 = -(GRR*xR - GLR*xL)/xc*half*dtdx &
                & + (ERR - ERL)/xc/sinxc*half*dtdz
           sCL0 = (FRL*xR - FLL*xL)/xc*half*dtdx - (ERL - ELL)/xc*half*dtdy
           sCR0 = (FRR*xR - FLR*xL)/xc*half*dtdx - (ERR - ELR)/xc*half*dtdy
#endif

           ! Update in time the  primitive variables
           r = r + sr0
           u = u + su0
           v = v + sv0
           w = w + sw0
           p = p + sp0
           A = A + sA0
           B = B + sB0
           C = C + sC0

           AL = AL + sAL0
           AR = AR + sAR0
           BL = BL + sBL0
           BR = BR + sBR0
           CL = CL + sCL0
           CR = CR + sCR0

           ! Face averaged right state at left interface
           qp(i,j,k,ir,1) = r - drx
           qp(i,j,k,iu,1) = u - dux
           qp(i,j,k,iv,1) = v - dvx
           qp(i,j,k,iw,1) = w - dwx
           qp(i,j,k,ip,1) = p - dpx
           qp(i,j,k,iA,1) = AL     
           qp(i,j,k,iB,1) = B - dBx
           qp(i,j,k,iC,1) = C - dCx
           qp(i,j,k,ir,1) = max(smallr, qp(i,j,k,ir,1))
           qp(i,j,k,ip,1) = max(smallp, qp(i,j,k,ip,1))

           ! Face averaged left state at right interface
           qm(i,j,k,ir,1) = r + drx
           qm(i,j,k,iu,1) = u + dux
           qm(i,j,k,iv,1) = v + dvx
           qm(i,j,k,iw,1) = w + dwx
           qm(i,j,k,ip,1) = p + dpx
           qm(i,j,k,iA,1) = AR     
           qm(i,j,k,iB,1) = B + dBx
           qm(i,j,k,iC,1) = C + dCx
           qm(i,j,k,ir,1) = max(smallr, qm(i,j,k,ir,1))
           qm(i,j,k,ip,1) = max(smallp, qm(i,j,k,ip,1))

           ! Face averaged top state at bottom interface
           qp(i,j,k,ir,2) = r - dry
           qp(i,j,k,iu,2) = u - duy
           qp(i,j,k,iv,2) = v - dvy
           qp(i,j,k,iw,2) = w - dwy
           qp(i,j,k,ip,2) = p - dpy
           qp(i,j,k,iA,2) = A - dAy
           qp(i,j,k,iB,2) = BL     
           qp(i,j,k,iC,2) = C - dCy
           qp(i,j,k,ir,2) = max(smallr, qp(i,j,k,ir,2))
           qp(i,j,k,ip,2) = max(smallp, qp(i,j,k,ip,2))

           ! Face averaged bottom state at top interface
           qm(i,j,k,ir,2) = r + dry
           qm(i,j,k,iu,2) = u + duy
           qm(i,j,k,iv,2) = v + dvy
           qm(i,j,k,iw,2) = w + dwy
           qm(i,j,k,ip,2) = p + dpy
           qm(i,j,k,iA,2) = A + dAy
           qm(i,j,k,iB,2) = BR     
           qm(i,j,k,iC,2) = C + dCy
           qm(i,j,k,ir,2) = max(smallr, qm(i,j,k,ir,2))
           qm(i,j,k,ip,2) = max(smallp, qm(i,j,k,ip,2))

           ! Face averaged front state at back interface
           qp(i,j,k,ir,3) = r - drz
           qp(i,j,k,iu,3) = u - duz
           qp(i,j,k,iv,3) = v - dvz
           qp(i,j,k,iw,3) = w - dwz
           qp(i,j,k,ip,3) = p - dpz
           qp(i,j,k,iA,3) = A - dAz
           qp(i,j,k,iB,3) = B - dBz
           qp(i,j,k,iC,3) = CL     
           qp(i,j,k,ir,3) = max(smallr, qp(i,j,k,ir,3))
           qp(i,j,k,ip,3) = max(smallp, qp(i,j,k,ip,3))

           ! Face averaged back state at front interface
           qm(i,j,k,ir,3) = r + drz
           qm(i,j,k,iu,3) = u + duz
           qm(i,j,k,iv,3) = v + dvz
           qm(i,j,k,iw,3) = w + dwz
           qm(i,j,k,ip,3) = p + dpz
           qm(i,j,k,iA,3) = A + dAz
           qm(i,j,k,iB,3) = B + dBz
           qm(i,j,k,iC,3) = CR     
           qm(i,j,k,ir,3) = max(smallr, qm(i,j,k,ir,3))
           qm(i,j,k,ip,3) = max(smallp, qm(i,j,k,ip,3))

           ! X-edge averaged right-top corner state (RT->LL)
           qRT(i,j,k,ir,1) = r + (+dry + drz)
           qRT(i,j,k,iu,1) = u + (+duy + duz)
           qRT(i,j,k,iv,1) = v + (+dvy + dvz)
           qRT(i,j,k,iw,1) = w + (+dwy + dwz)
           qRT(i,j,k,ip,1) = p + (+dpy + dpz)
           qRT(i,j,k,iA,1) = A + (+dAy + dAz)
           qRT(i,j,k,iB,1) = BR + dBRz
           qRT(i,j,k,iC,1) = CR + dCRy
           qRT(i,j,k,ir,1) = max(smallr, qRT(i,j,k,ir,1))
           qRT(i,j,k,ip,1) = max(smallp, qRT(i,j,k,ip,1))

           ! X-edge averaged right-bottom corner state (RB->LR)
           qRB(i,j,k,ir,1) = r + (+dry - drz)
           qRB(i,j,k,iu,1) = u + (+duy - duz)
           qRB(i,j,k,iv,1) = v + (+dvy - dvz)
           qRB(i,j,k,iw,1) = w + (+dwy - dwz)
           qRB(i,j,k,ip,1) = p + (+dpy - dpz)
           qRB(i,j,k,iA,1) = A + (+dAy - dAz)
           qRB(i,j,k,iB,1) = BR - dBRz
           qRB(i,j,k,iC,1) = CL + dCLy
           qRB(i,j,k,ir,1) = max(smallr, qRB(i,j,k,ir,1))
           qRB(i,j,k,ip,1) = max(smallp, qRB(i,j,k,ip,1))

           ! X-edge averaged left-top corner state (LT->RL)
           qLT(i,j,k,ir,1) = r + (-dry + drz)
           qLT(i,j,k,iu,1) = u + (-duy + duz)
           qLT(i,j,k,iv,1) = v + (-dvy + dvz)
           qLT(i,j,k,iw,1) = w + (-dwy + dwz)
           qLT(i,j,k,ip,1) = p + (-dpy + dpz)
           qLT(i,j,k,iA,1) = A + (-dAy + dAz)
           qLT(i,j,k,iB,1) = BL + dBLz
           qLT(i,j,k,iC,1) = CR - dCRy
           qLT(i,j,k,ir,1) = max(smallr, qLT(i,j,k,ir,1))
           qLT(i,j,k,ip,1) = max(smallp, qLT(i,j,k,ip,1))

           ! X-edge averaged left-bottom corner state (LB->RR)
           qLB(i,j,k,ir,1) = r + (-dry - drz)
           qLB(i,j,k,iu,1) = u + (-duy - duz)
           qLB(i,j,k,iv,1) = v + (-dvy - dvz)
           qLB(i,j,k,iw,1) = w + (-dwy - dwz)
           qLB(i,j,k,ip,1) = p + (-dpy - dpz)
           qLB(i,j,k,iA,1) = A + (-dAy - dAz)
           qLB(i,j,k,iB,1) = BL - dBLz
           qLB(i,j,k,iC,1) = CL - dCLy
           qLB(i,j,k,ir,1) = max(smallr, qLB(i,j,k,ir,1))
           qLB(i,j,k,ip,1) = max(smallp, qLB(i,j,k,ip,1))

           ! Y-edge averaged right-top corner state (RT->LL)
           qRT(i,j,k,ir,2) = r + (+drx + drz)
           qRT(i,j,k,iu,2) = u + (+dux + duz)
           qRT(i,j,k,iv,2) = v + (+dvx + dvz)
           qRT(i,j,k,iw,2) = w + (+dwx + dwz)
           qRT(i,j,k,ip,2) = p + (+dpx + dpz)
           qRT(i,j,k,iA,2) = AR + dARz
           qRT(i,j,k,iB,2) = B + (+dBx + dBz)
           qRT(i,j,k,iC,2) = CR + dCRx
           qRT(i,j,k,ir,2) = max(smallr, qRT(i,j,k,ir,2))
           qRT(i,j,k,ip,2) = max(smallp, qRT(i,j,k,ip,2))

           ! Y-edge averaged right-bottom corner state (RB->LR)
           qRB(i,j,k,ir,2) = r + (+drx - drz)
           qRB(i,j,k,iu,2) = u + (+dux - duz)
           qRB(i,j,k,iv,2) = v + (+dvx - dvz)
           qRB(i,j,k,iw,2) = w + (+dwx - dwz)
           qRB(i,j,k,ip,2) = p + (+dpx - dpz)
           qRB(i,j,k,iA,2) = AR - dARz
           qRB(i,j,k,iB,2) = B + (+dBx - dBz)
           qRB(i,j,k,iC,2) = CL + dCLx
           qRB(i,j,k,ir,2) = max(smallr, qRB(i,j,k,ir,2))
           qRB(i,j,k,ip,2) = max(smallp, qRB(i,j,k,ip,2))

           ! Y-edge averaged left-top corner state (LT->RL)
           qLT(i,j,k,ir,2) = r + (-drx + drz)
           qLT(i,j,k,iu,2) = u + (-dux + duz)
           qLT(i,j,k,iv,2) = v + (-dvx + dvz)
           qLT(i,j,k,iw,2) = w + (-dwx + dwz)
           qLT(i,j,k,ip,2) = p + (-dpx + dpz)
           qLT(i,j,k,iA,2) = AL + dALz
           qLT(i,j,k,iB,2) = B + (-dBx + dBz)
           qLT(i,j,k,iC,2) = CR - dCRx
           qLT(i,j,k,ir,2) = max(smallr, qLT(i,j,k,ir,2))
           qLT(i,j,k,ip,2) = max(smallp, qLT(i,j,k,ip,2))

           ! Y-edge averaged left-bottom corner state (LB->RR)
           qLB(i,j,k,ir,2) = r + (-drx - drz)
           qLB(i,j,k,iu,2) = u + (-dux - duz)
           qLB(i,j,k,iv,2) = v + (-dvx - dvz)
           qLB(i,j,k,iw,2) = w + (-dwx - dwz)
           qLB(i,j,k,ip,2) = p + (-dpx - dpz)
           qLB(i,j,k,iA,2) = AL - dALz
           qLB(i,j,k,iB,2) = B + (-dBx - dBz)
           qLB(i,j,k,iC,2) = CL - dCLx
           qLB(i,j,k,ir,2) = max(smallr, qLB(i,j,k,ir,2))
           qLB(i,j,k,ip,2) = max(smallp, qLB(i,j,k,ip,2))

           ! Z-edge averaged right-top corner state (RT->LL)
           qRT(i,j,k,ir,3) = r + (+drx + dry)
           qRT(i,j,k,iu,3) = u + (+dux + duy)
           qRT(i,j,k,iv,3) = v + (+dvx + dvy)
           qRT(i,j,k,iw,3) = w + (+dwx + dwy)
           qRT(i,j,k,ip,3) = p + (+dpx + dpy)
           qRT(i,j,k,iA,3) = AR + dARy
           qRT(i,j,k,iB,3) = BR + dBRx
           qRT(i,j,k,iC,3) = C + (+dCx + dCy)
           qRT(i,j,k,ir,3) = max(smallr, qRT(i,j,k,ir,3))
           qRT(i,j,k,ip,3) = max(smallp, qRT(i,j,k,ip,3))

           ! Z-edge averaged right-bottom corner state (RB->LR)
           qRB(i,j,k,ir,3) = r + (+drx - dry)
           qRB(i,j,k,iu,3) = u + (+dux - duy)
           qRB(i,j,k,iv,3) = v + (+dvx - dvy)
           qRB(i,j,k,iw,3) = w + (+dwx - dwy)
           qRB(i,j,k,ip,3) = p + (+dpx - dpy)
           qRB(i,j,k,iA,3) = AR - dARy
           qRB(i,j,k,iB,3) = BL + dBLx
           qRB(i,j,k,iC,3) = C + (+dCx - dCy)
           qRB(i,j,k,ir,3) = max(smallr, qRB(i,j,k,ir,3))
           qRB(i,j,k,ip,3) = max(smallp, qRB(i,j,k,ip,3))

           ! Z-edge averaged left-top corner state (LT->RL)
           qLT(i,j,k,ir,3) = r + (-drx + dry)
           qLT(i,j,k,iu,3) = u + (-dux + duy)
           qLT(i,j,k,iv,3) = v + (-dvx + dvy)
           qLT(i,j,k,iw,3) = w + (-dwx + dwy)
           qLT(i,j,k,ip,3) = p + (-dpx + dpy)
           qLT(i,j,k,iA,3) = AL + dALy
           qLT(i,j,k,iB,3) = BR - dBRx
           qLT(i,j,k,iC,3) = C + (-dCx + dCy)
           qLT(i,j,k,ir,3) = max(smallr, qLT(i,j,k,ir,3))
           qLT(i,j,k,ip,3) = max(smallp, qLT(i,j,k,ip,3))

           ! Z-edge averaged left-bottom corner state (LB->RR)
           qLB(i,j,k,ir,3) = r + (-drx - dry)
           qLB(i,j,k,iu,3) = u + (-dux - duy)
           qLB(i,j,k,iv,3) = v + (-dvx - dvy)
           qLB(i,j,k,iw,3) = w + (-dwx - dwy)
           qLB(i,j,k,ip,3) = p + (-dpx - dpy)
           qLB(i,j,k,iA,3) = AL - dALy
           qLB(i,j,k,iB,3) = BL - dBLx
           qLB(i,j,k,iC,3) = C + (-dCx - dCy)
           qLB(i,j,k,ir,3) = max(smallr, qLB(i,j,k,ir,3))
           qLB(i,j,k,ip,3) = max(smallp, qLB(i,j,k,ip,3))
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  !$acc end data

  return
end subroutine trace3d
#endif
!===============================================================================
!> Compute fluxes from Riemann fluxes ! to merge with update
!===============================================================================
subroutine compflux(fgodunov)
  use params
  use variables, only: flux, dt, ds
  implicit none

  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar,ndim), intent(in) :: fgodunov
  integer :: i, j, k, idim
  integer :: ln, lt1, lt2, bn, bt1, bt2

  !$py begin_statement

  do idim = 1, ndim
     if (idim == 1) then
        ln = 2; lt1 = 3; lt2 = 4
        bn = 6; bt1 = 7; bt2 = 8
     else if (idim == 2) then
        ln = 3; lt1 = 2; lt2 = 4
        bn = 7; bt1 = 6; bt2 = 8
     else
        ln = 4; lt1 = 2; lt2 = 3
        bn = 8; bt1 = 6; bt2 = 7
     endif

     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     !$acc kernels loop
     do k = kf1, kf2
        do j = jf1, jf2
           do i = if1, if2
              flux(i,j,k,ir,idim)  = fgodunov(i,j,k,1,idim)*ds(i,j,k,idim)*dt
              flux(i,j,k,ip,idim)  = fgodunov(i,j,k,2,idim)*ds(i,j,k,idim)*dt
              flux(i,j,k,ln,idim)  = fgodunov(i,j,k,3,idim)*ds(i,j,k,idim)*dt
              flux(i,j,k,bn,idim)  = fgodunov(i,j,k,4,idim)*ds(i,j,k,idim)*dt
              flux(i,j,k,lt1,idim) = fgodunov(i,j,k,5,idim)*ds(i,j,k,idim)*dt
              flux(i,j,k,bt1,idim) = fgodunov(i,j,k,6,idim)*ds(i,j,k,idim)*dt
              flux(i,j,k,lt2,idim) = fgodunov(i,j,k,7,idim)*ds(i,j,k,idim)*dt
              flux(i,j,k,bt2,idim) = fgodunov(i,j,k,8,idim)*ds(i,j,k,idim)*dt
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  enddo
  
  return
end subroutine compflux
!===============================================================================
!> Update variables
!===============================================================================
subroutine update(fgodunov_pre)
  use params
  use variables, only: uin, flux, dt, x, y, dv
  use oacc_params
  implicit none

  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,ndim), intent(in) :: fgodunov_pre
  real(dp), dimension(3,2,3) :: radius
  real(dp) :: xL, xR, xc, sinxL, sinxR, sinxc, yc
  real(dp) :: ratio, lambda
  real(dp) :: dvloc, dsx, dsy, dfn2, dfn3
  integer  :: i, j, k, idim
  integer  :: ihi, jhi, khi
  integer  :: i0, j0, k0

  !$py begin_statement

#if GEOM == CARTESIAN
  if (Omega0 > zero) then
     lambda = Omega0*dt; lambda = 2.5d-1*lambda*lambda
     ratio  = (one - lambda)/(one + lambda)
  endif
#endif

  ihi = max(1, if2 - 1)
  jhi = max(1, jf2 - 1)
  khi = max(1, kf2 - 1)

#if OACC == 1
  if (ndim_act == 2) then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(yc, radius, i0, j0, k0, dsx) &
     !$OMP PRIVATE(dsy, dfn2, dfn3)
     do k = 1, khi
        do j = 1, jhi
           yc = y(j)
           !$acc loop private(radius)
           do i = 1, ihi
              radius = one
              dvloc  = one/dv(i,j,k)
#if GEOM == CYLINDRICAL || GEOM == SPHERICAL
              radius(2,1:2,1:3) = x(i)
              radius(2,1,1) = half*(x(i) + x(i-1))
              radius(2,2,1) = half*(x(i) + x(i+1))
#endif
#if GEOM == SPHERICAL
              radius(3,1,1) = half*(x(i) + x(i-1))*sin(yc)
              radius(3,2,1) = half*(x(i) + x(i+1))*sin(yc)
              radius(3,1,2) = x(i)*sin(half*(y(j) + y(j-1)))
              radius(3,2,2) = x(i)*sin(half*(y(j) + y(j+1)))
              radius(3,1,3) = x(i)*sin(yc)
              radius(3,2,3) = x(i)*sin(yc)
#endif
#if GEOM == CARTESIAN
              if (Omega0 > zero) then
                 ! Source term in case of shearing box
                 dsx = two*Omega0*dt*uin(i,j,k,iv)/(one + lambda)
                 dsy = -half*Omega0*dt*uin(i,j,k,iu)/(one + lambda)
                 
                 uin(i,j,k,iu) = uin(i,j,k,iu)*ratio + dsx 
                 uin(i,j,k,iv) = uin(i,j,k,iv)*ratio + dsy 
              endif
#endif
              
              !$acc loop independent private(i0, j0, k0)
              do idim = 1, ndim
                 i0 = 0; j0 = 0; k0 = 0
                 if (idim == 1) i0 = 1
                 if (idim == 2) j0 = 1
                 if (idim == 3) k0 = 1
                 uin(i,j,k,ir) = uin(i,j,k,ir) + (flux(i,j,k,ir,idim) &
                      & - flux(i+i0,j+j0,k+k0,ir,idim))*dvloc
                 uin(i,j,k,iw) = uin(i,j,k,iw) &
                      & + (flux(i,j,k,iw,idim)*radius(iw-1,1,idim) &
                      & - flux(i+i0,j+j0,k+k0,iw,idim)*radius(iw-1,2,idim)) &
                      & *dvloc/radius(iw-1,2,3)
                 uin(i,j,k,ip) = uin(i,j,k,ip) + (flux(i,j,k,ip,idim) &
                      & - flux(i+i0,j+j0,k+k0,ip,idim))*dvloc
                 
                 if (ndim == 1) then
                    uin(i,j,k,iA:iC) = uin(i,j,k,iA:iC) &
                         & + (flux(i,j,k,iA:iC,idim) &
                         & - flux(i+i0,j+j0,k+k0,iA:iC,idim))*dvloc
                 else if (ndim == 2) then
                    uin(i,j,k,iC) = uin(i,j,k,iC) + (flux(i,j,k,iC,idim) &
                         & - flux(i+i0,j+j0,k+k0,iC,idim))*dvloc
                 endif
                 
#if GEOM == CARTESIAN
                 if (Omega0 > zero) then
                    dfn2 = (flux(i,j,k,iu,idim) - flux(i+i0,j+j0,k+k0,iu,idim)) &
                         & *dvloc
                    if (idim == 1) dfn2 = dfn2 + fgodunov_pre(i,j,k,1)*dvloc
                    dfn3 = (flux(i,j,k,iv,idim) - flux(i+i0,j+j0,k+k0,iv,idim)) &
                         & *dvloc
#if NDIM > 1
                    if (idim == 2) dfn3 = dfn3 + fgodunov_pre(i,j,k,2)*dvloc
#endif
                    
                    uin(i,j,k,iu) = uin(i,j,k,iu) + dfn2/(one + lambda) &
                         & + Omega0*dt/(one + lambda)*dfn3
                    uin(i,j,k,iv) = uin(i,j,k,iv) + dfn3/(one + lambda) &
                         & - 2.5d-1*Omega0*dt/(one + lambda)*dfn2
                 else
#endif
                    uin(i,j,k,iu) = uin(i,j,k,iu) &
                         & + (flux(i,j,k,iu,idim)*radius(iu-1,1,idim) &
                         & - flux(i+i0,j+j0,k+k0,iu,idim)*radius(iu-1,2,idim)) &
                         & *dvloc/radius(iu-1,2,3)
                    uin(i,j,k,iv) = uin(i,j,k,iv)  &
                         & + (flux(i,j,k,iv,idim)*radius(iv-1,1,idim) &
                         & - flux(i+i0,j+j0,k+k0,iv,idim)*radius(iv-1,2,idim)) &
                         & *dvloc/radius(iv-1,2,3)
#if GEOM == CARTESIAN
                 endif
#endif
              enddo
              
#if GEOM == CARTESIAN
              if (Omega0 == zero) then
#endif
                 uin(i,j,k,iu) = uin(i,j,k,iu) + fgodunov_pre(i,j,k,1)*dvloc
#if NDIM > 1
                 uin(i,j,k,iv) = uin(i,j,k,iv) + fgodunov_pre(i,j,k,2)*dvloc
#endif
#if GEOM == CARTESIAN
              endif
#endif
              
#if NDIM > 2 
              uin(i,j,k,iw) = uin(i,j,k,iw) + fgodunov_pre(i,j,k,3)*dvloc
#endif
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  else
#endif
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(yc, radius, i0, j0, k0, dsx) &
     !$OMP PRIVATE(dsy, dfn2, dfn3)
     do k = 1, khi
        !$acc loop vector(blocky_update)
        do j = 1, jhi
           yc = y(j)
           !$acc loop private(radius) vector(blockx_update)
           do i = 1, ihi
              radius = one
              dvloc  = one/dv(i,j,k)
#if GEOM == CYLINDRICAL || GEOM == SPHERICAL
              radius(2,1:2,1:3) = x(i)
              radius(2,1,1) = half*(x(i) + x(i-1))
              radius(2,2,1) = half*(x(i) + x(i+1))
#endif
#if GEOM == SPHERICAL
              radius(3,1,1) = half*(x(i) + x(i-1))*sin(yc)
              radius(3,2,1) = half*(x(i) + x(i+1))*sin(yc)
              radius(3,1,2) = x(i)*sin(half*(y(j) + y(j-1)))
              radius(3,2,2) = x(i)*sin(half*(y(j) + y(j+1)))
              radius(3,1,3) = x(i)*sin(yc)
              radius(3,2,3) = x(i)*sin(yc)
#endif
#if GEOM == CARTESIAN
              if (Omega0 > zero) then
                 ! Source term in case of shearing box
                 dsx = two*Omega0*dt*uin(i,j,k,iv)/(one + lambda)
                 dsy = -half*Omega0*dt*uin(i,j,k,iu)/(one + lambda)
                 
                 uin(i,j,k,iu) = uin(i,j,k,iu)*ratio + dsx 
                 uin(i,j,k,iv) = uin(i,j,k,iv)*ratio + dsy 
              endif
#endif
              
              !$acc loop independent private(i0, j0, k0)
              do idim = 1, ndim
                 i0 = 0; j0 = 0; k0 = 0
                 if (idim == 1) i0 = 1
                 if (idim == 2) j0 = 1
                 if (idim == 3) k0 = 1
                 uin(i,j,k,ir) = uin(i,j,k,ir) + (flux(i,j,k,ir,idim) &
                      & - flux(i+i0,j+j0,k+k0,ir,idim))*dvloc
                 uin(i,j,k,iw) = uin(i,j,k,iw) &
                      & + (flux(i,j,k,iw,idim)*radius(iw-1,1,idim) &
                      & - flux(i+i0,j+j0,k+k0,iw,idim)*radius(iw-1,2,idim)) &
                      & *dvloc/radius(iw-1,2,3)
                 uin(i,j,k,ip) = uin(i,j,k,ip) + (flux(i,j,k,ip,idim) &
                      & - flux(i+i0,j+j0,k+k0,ip,idim))*dvloc
                 
                 if (ndim == 1) then
                    uin(i,j,k,iA:iC) = uin(i,j,k,iA:iC) &
                         & + (flux(i,j,k,iA:iC,idim) &
                         & - flux(i+i0,j+j0,k+k0,iA:iC,idim))*dvloc
                 else if (ndim == 2) then
                    uin(i,j,k,iC) = uin(i,j,k,iC) + (flux(i,j,k,iC,idim) &
                         & - flux(i+i0,j+j0,k+k0,iC,idim))*dvloc
                 endif
                 
#if GEOM == CARTESIAN
                 if (Omega0 > zero) then
                    dfn2 = (flux(i,j,k,iu,idim) - flux(i+i0,j+j0,k+k0,iu,idim)) &
                         & *dvloc
                    if (idim == 1) dfn2 = dfn2 + fgodunov_pre(i,j,k,1)*dvloc
                    dfn3 = (flux(i,j,k,iv,idim) - flux(i+i0,j+j0,k+k0,iv,idim)) &
                         & *dvloc
#if NDIM > 1
                    if (idim == 2) dfn3 = dfn3 + fgodunov_pre(i,j,k,2)*dvloc
#endif
                    
                    uin(i,j,k,iu) = uin(i,j,k,iu) + dfn2/(one + lambda) &
                         & + Omega0*dt/(one + lambda)*dfn3
                    uin(i,j,k,iv) = uin(i,j,k,iv) + dfn3/(one + lambda) &
                         & - 2.5d-1*Omega0*dt/(one + lambda)*dfn2
                 else
#endif
                    uin(i,j,k,iu) = uin(i,j,k,iu) &
                         & + (flux(i,j,k,iu,idim)*radius(iu-1,1,idim) &
                         & - flux(i+i0,j+j0,k+k0,iu,idim)*radius(iu-1,2,idim)) &
                         & *dvloc/radius(iu-1,2,3)
                    uin(i,j,k,iv) = uin(i,j,k,iv)  &
                         & + (flux(i,j,k,iv,idim)*radius(iv-1,1,idim) &
                         & - flux(i+i0,j+j0,k+k0,iv,idim)*radius(iv-1,2,idim)) &
                         & *dvloc/radius(iv-1,2,3)
#if GEOM == CARTESIAN
                 endif
#endif
              enddo
              
#if GEOM == CARTESIAN
              if (Omega0 == 0) then
#endif
                 uin(i,j,k,iu) = uin(i,j,k,iu) + fgodunov_pre(i,j,k,1)*dvloc
#if NDIM > 1
                 uin(i,j,k,iv) = uin(i,j,k,iv) + fgodunov_pre(i,j,k,2)*dvloc
#endif
#if GEOM == CARTESIAN
              endif
#endif

              
#if NDIM > 2 
              uin(i,j,k,iw) = uin(i,j,k,iw) + fgodunov_pre(i,j,k,3)*dvloc
#endif
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
#if OACC == 1
  endif
#endif

  if (ndim_act > 1) then
     call constrained_transport
  endif

  ! Right faces B-field
#if OACC == 1
  !$acc kernels loop async
  do i = iu1, iu2-1
     uin(i,:,:,iA+3) = uin(i+1,:,:,iA)
  enddo
  !$acc kernels loop async
  do j = ju1, ju2-1
     uin(:,j,:,iB+3) = uin(:,j+1,:,iB)
  enddo
  !$acc kernels loop async
  do k = ku1, ku2-1
     uin(:,:,k,iC+3) = uin(:,:,k+1,iC)
  enddo
#else
  !$OMP PARALLEL WORKSHARE
  uin(iu1:iu2-1,ju1:ju2,ku1:ku2,iA+3) = uin(iu1+1:iu2,ju1:ju2,ku1:ku2,iA)
  uin(iu1:iu2,ju1:ju2-1,ku1:ku2,iB+3) = uin(iu1:iu2,ju1+1:ju2,ku1:ku2,iB)
  uin(iu1:iu2,ju1:ju2,ku1:ku2-1,iC+3) = uin(iu1:iu2,ju1:ju2,ku1+1:ku2,iC)
  !$OMP END PARALLEL WORKSHARE
#endif

  return
end subroutine update
!===============================================================================
!> Update magnetic field from constraint transport method
!===============================================================================
subroutine constrained_transport
  use params
  use variables
  implicit none

  real(dp) :: xL, xR, xc, sinxL, sinxR, sinxc
  integer  :: ihi, jhi, khi
  integer  :: i, j, k

  !$py begin_statement

  ihi = max(1, if2-1)
  jhi = max(1, jf2-1)
  khi = max(1, kf2-1)

  ! Constrained transport for face-centered B-field
#if GEOM == CARTESIAN
  if (ndim >= 2) then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = 1, khi
        do j = 1, jhi+1
           do i = 1, ihi+1
              ! Left-face B-field
              uin(i,j,k,iA) = uin(i,j,k,iA) &
                            & + (emfz(i,j+1,k) - emfz(i,j,k))/dy
              uin(i,j,k,iB) = uin(i,j,k,iB) &
                            & - (emfz(i+1,j,k) - emfz(i,j,k))/dx
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  endif
  if (ndim == 3) then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = 1, khi+1
        do j = 1, jhi+1
           do i = 1, ihi+1
              ! Left-face B-field
              uin(i,j,k,iA) = uin(i,j,k,iA) &
                            & - (emfy(i,j,k+1) - emfy(i,j,k))/dz
              uin(i,j,k,iB) = uin(i,j,k,iB) &
                            & + (emfx(i,j,k+1) - emfx(i,j,k))/dz
              uin(i,j,k,iC) = uin(i,j,k,iC) &
                            & + (emfy(i+1,j,k) - emfy(i,j,k))/dx &
                            & - (emfx(i,j+1,k) - emfx(i,j,k))/dy
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  endif
#endif

#if GEOM == CYLINDRICAL
  if (ndim >= 2) then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = 1, khi
        do j = 1, jhi+1
           do i = 1, ihi+1
              xL = half*(x(i) + x(i-1))
              ! Left-face B-field
              uin(i,j,k,iA) = uin(i,j,k,iA) &
                            & + (emfz(i,j+1,k) - emfz(i,j,k))/dy/xL
              uin(i,j,k,iB) = uin(i,j,k,iB) &
                            & - (emfz(i+1,j,k) - emfz(i,j,k))/dx
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  endif
  if (ndim == 3) then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(xL, xR, xc)
     do k = 1, khi+1
        do j = 1, jhi+1
           do i = 1, ihi+1
              xL = half*(x(i) + x(i-1))
              xR = half*(x(i) + x(i+1))
              xc = x(i)
              ! Left-face B-field
              uin(i,j,k,iA) = uin(i,j,k,iA) &
                            & - (emfy(i,j,k+1) - emfy(i,j,k))/dz
              uin(i,j,k,iB) = uin(i,j,k,iB) &
                            & + (emfx(i,j,k+1) - emfx(i,j,k))/dz
              uin(i,j,k,iC) = uin(i,j,k,iC) &
                            & + (emfy(i+1,j,k)*xR - emfy(i,j,k)*xL)/dx/xc &
                            & - (emfx(i,j+1,k) - emfx(i,j,k))/dy/xc
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  endif
#endif

#if GEOM == SPHERICAL
  if (ndim >= 2) then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(sinxL, sinxR, xinxc, xL, xR, xc)
     do k = 1, khi
        do j = 1, jhi+1
           sinxL = sin(half*(y(j) + y(j-1)))
           sinxR = sin(half*(y(j) + y(j+1)))
           sinxc = sin(y(j))
           do i = 1, ihi+1
              xL = half*(x(i) + x(i-1))
              xR = half*(x(i) + x(i+1))
              xc = x(i)
              ! Left-face B-field
              uin(i,j,k,iA) = uin(i,j,k,iA) &
                            & + (emfz(i,j+1,k)*sinxR - emfz(i,j,k)*sinxL) &
                            & /dy/xL/sinxc
              uin(i,j,k,iB) = uin(i,j,k,iB) &
                            & - (emfz(i+1,j,k)*xR - emfz(i,j,k)*xL)/dx/xc
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  endif
  if (ndim == 3) then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(sinxL, xinxc, xL, xR, xc)
     do k = 1, khi+1
        do j = 1, jhi+1
           sinxL = sin(half*(y(j) + y(j-1)))
           sinxc = sin(y(j))
           do i = 1, ihi+1
              xL = half*(x(i) + x(i-1))
              xR = half*(x(i) + x(i+1))
              xc = x(i)
              ! Left-face B-field
              uin(i,j,k,iA) = uin(i,j,k,iA) &
                            & - (emfy(i,j,k+1) - emfy(i,j,k))/dz/xL/sinxc
              uin(i,j,k,iB) = uin(i,j,k,iB) &
                            & + (emfx(i,j,k+1) - emfx(i,j,k))/dz/xc/sinxL
              uin(i,j,k,iC) = uin(i,j,k,iC) &
                            & + (emfy(i+1,j,k)*xR - emfy(i,j,k)*xL)/dx/xc &
                            & - (emfx(i,j+1,k) - emfx(i,j,k))/dy/xc
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  endif
#endif

  return
end subroutine constrained_transport
