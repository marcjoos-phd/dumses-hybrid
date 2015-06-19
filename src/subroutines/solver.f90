!===============================================================================
!> \file solver.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is Riemann solver subroutines. This file cannot be compiled as is, it 
!! needs to be preprocessed with the DUMSES preprocessor. riemann subroutine is 
!! a global_template for solver, and the other subroutines are solver_template.
!! \details
!! Contains riemann(), upwind(), llf(), hll(), hlld(), acoustic(),
!! athena_roe(), mhd_eigenvalues(), roe_eigenvalues()
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
!! \b last \b modified: 05-04-2015
!<
!===============================================================================
!> Riemann solver:
!!  retrieve left and right states from trace, and call the solver specified
!! in \c riemann variable
!<
!===============================================================================
!$py global_template riemann
subroutine riemann(qm, qp, fgodunov, fgodunov_pre)
  use params
  use variables, only: dt, ds, dv, x, y
  use oacc_params
  implicit none

  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar,ndim), intent(in) :: qm 
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar,ndim), intent(in) :: qp 
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar,ndim),intent(out) :: fgodunov
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,ndim), intent(out) :: fgodunov_pre
  real(dp), dimension(:,:,:), allocatable :: rgstar, ugstar, vgstar, wgstar
  real(dp), dimension(:,:,:), allocatable :: bgstar, cgstar, pgstar
  real(dp) :: bn_mean
  real(dp) :: rl, pl, ul, vl, wl, cl, bl, al
  real(dp) :: rr, pr, ur, vr, wr, cr, br, ar
  real(dp) :: ro, uo, vo, wo, bo, co, ptoto, pressure
  real(dp) :: cotanxc, shear, Ekin, Emag, Etot
  integer  :: i, j, k, idim, im, jm, km, ii, ji, ki, ioffset, joffset, koffset
  integer  :: ln, lt1, lt2, bn, bt1, bt2, ihip, jhip, khip
  integer  :: ilo, ihi, jlo, jhi, klo, khi

  !$py begin_statement

  allocate(rgstar(iu1:iu2,ju1:ju2,ku1:ku2))
  allocate(ugstar(iu1:iu2,ju1:ju2,ku1:ku2))
  allocate(vgstar(iu1:iu2,ju1:ju2,ku1:ku2))
  allocate(wgstar(iu1:iu2,ju1:ju2,ku1:ku2))
  allocate(bgstar(iu1:iu2,ju1:ju2,ku1:ku2))
  allocate(cgstar(iu1:iu2,ju1:ju2,ku1:ku2))
  allocate(pgstar(iu1:iu2,ju1:ju2,ku1:ku2))
  
#if OACC == 1
  !$acc kernels loop
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2
           rgstar(i,j,k) = zero
           ugstar(i,j,k) = zero
           vgstar(i,j,k) = zero
           wgstar(i,j,k) = zero
           bgstar(i,j,k) = zero
           cgstar(i,j,k) = zero
           pgstar(i,j,k) = zero
        enddo
     enddo
  enddo
#else
  !$OMP PARALLEL WORKSHARE
  rgstar = zero; ugstar = zero; vgstar = zero; wgstar = zero
  bgstar = zero; cgstar = zero; pgstar = zero
  !$OMP END PARALLEL WORKSHARE
#endif

  !$acc data pcreate(rgstar, ugstar, vgstar, wgstar, bgstar, cgstar, pgstar)

  do idim = 1, ndim
     if(idim == 1) then
        ln = 2; lt1 = 3; lt2 = 4
        bn = 6; bt1 = 7; bt2 = 8
        ilo = if1          ; ihi = if2
        jlo = min(1, ju1+3); jhi = max(1, ju2-3)
        klo = min(1, ku1+3); khi = max(1, ku2-3)
        ioffset = 1; joffset = 0; koffset = 0
     else if(idim == 2) then
        ln = 3; lt1 = 2; lt2 = 4
        bn = 7; bt1 = 6; bt2 = 8
        ilo = min(1, iu1+3); ihi = max(1, iu2-3)
        jlo = jf1          ; jhi = jf2
        klo = min(1, ku1+3); khi = max(1, ku2-3)
        ioffset = 0; joffset = 1; koffset = 0
     else
        ln = 4; lt1 = 2; lt2 = 3
        bn = 8; bt1 = 6; bt2 = 7
        ilo = min(1, iu1+3); ihi = max(1, iu2-3)
        jlo = min(1, ju1+3); jhi = max(1, ju2-3)
        klo = kf1          ; khi = kf2
        ioffset = 0; joffset = 0; koffset = 1
     endif
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(bn_mean, im, jm, km, shear) &
     !$OMP PRIVATE(Ekin, Emag, Etot, rl, pl, ul, vl, wl, cl, bl, al, rr, pr) &
     !$OMP PRIVATE(ur, vr, wr, cr, br, ar, ro, uo, vo, wo, bo, co, ptoto)
     do k = klo, khi
        !$acc loop vector(blocky_solver)
        do j = jlo, jhi
           !$acc loop vector(blockx_solver)
           do i = ilo, ihi
              im = i - ioffset
              jm = j - joffset
              km = k - koffset

              ! Enforce continuity for normal magnetic field
              bn_mean = half*(qm(im,jm,km,bn,idim) + qp(i,j,k,bn,idim))

              ! Left state
              rl = qm(im,jm,km,1  ,idim) ! Mass density
              pl = qm(im,jm,km,5  ,idim) ! Pressure
              ul = qm(im,jm,km,ln ,idim) ! Normal velocity
              al = bn_mean               ! Normal magnetic field
              vl = qm(im,jm,km,lt1,idim) ! Tangential velocity 1
              bl = qm(im,jm,km,bt1,idim) ! Tangential magnetic field 1
              wl = qm(im,jm,km,lt2,idim) ! Tangential velocity 2
              cl = qm(im,jm,km,bt2,idim) ! Tangential magnetic field 2
              
              ! Right state
              rr = qp(i,j,k,1  ,idim) ! Mass density
              pr = qp(i,j,k,5  ,idim) ! Pressure
              ur = qp(i,j,k,ln ,idim) ! Normal velocity
              ar = bn_mean            ! Normal magnetic field
              vr = qp(i,j,k,lt1,idim) ! Tangential velocity 1
              br = qp(i,j,k,bt1,idim) ! Tangential magnetic field 1
              wr = qp(i,j,k,lt2,idim) ! Tangential velocity 2
              cr = qp(i,j,k,bt2,idim) ! Tangential magnetic field 2

              !$py insert_solver

              ! Upwind solver in case of the shearing box (only on the hydro
              ! variables))
#if GEOM == CARTESIAN
#if NDIM == 1
              if (Omega0 > zero) then
                 shear = -1.5d0*Omega0*(x(i) + x(i-1))
                 fgodunov(i,j,k,8,idim) = fgodunov(i,j,k,8,idim) &
                                      & + shear*bn_mean
              endif
#endif
#if NDIM == 2
              if (Omega0 > zero) then
                 shear = -1.5d0*Omega0*x(i)
                 if (idim == 1) shear = -1.5d0*Omega0*(x(i) + x(i-1))
                 if (idim == 2) shear = -1.5d0*Omega0*x(i)
                 fgodunov(i,j,k,8,idim) = fgodunov(i,j,k,8,idim) &
                                      & + shear*bn_mean
              endif
#endif
#if NDIM == 3
              if ((Omega0 > zero) .and. (idim == 2) .and. (.not. fargo)) then
                 shear = -1.5d0*Omega0*x(i)
                 if (shear > zero) then
                    Emag = half*(al**2 + bl**2 + cl**2)
                    Ekin = half*(ul**2 + vl**2 + wl**2)
                    Etot = Ekin + Emag + pl/(gamma - one)
                    fgodunov(i,j,k,1,idim) = fgodunov(i,j,k,1,idim) &
                            & + shear*rl
                    fgodunov(i,j,k,2,idim) = fgodunov(i,j,k,2,idim) &
                            & + shear*(Etot + Emag - bn_mean**2)
                    fgodunov(i,j,k,3,idim) = fgodunov(i,j,k,3,idim) &
                            & + shear*rl*ul
                    fgodunov(i,j,k,5,idim) = fgodunov(i,j,k,5,idim) &
                            & + shear*rl*vl
                    fgodunov(i,j,k,7,idim) = fgodunov(i,j,k,7,idim) &
                            & + shear*rl*wl
                 else
                    Emag = half*(ar**2 + br**2 + cr**2)
                    Ekin = half*(ur**2 + vr**2 + wr**2)
                    Etot = Ekin + Emag + pr/(gamma - one)
                    fgodunov(i,j,k,1,idim) = fgodunov(i,j,k,1,idim) &
                            & + shear*rr
                    fgodunov(i,j,k,2,idim) = fgodunov(i,j,k,2,idim) &
                            & + shear*(Etot + Emag - bn_mean**2)
                    fgodunov(i,j,k,3,idim) = fgodunov(i,j,k,3,idim) &
                            & + shear*rr*ur
                    fgodunov(i,j,k,5,idim) = fgodunov(i,j,k,5,idim) &
                            & + shear*rr*vr
                    fgodunov(i,j,k,7,idim) = fgodunov(i,j,k,7,idim) &
                            & + shear*rr*wr
                 endif
              endif
#endif
#endif
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO

     ! Add geometrical terms
#if GEOM == CARTESIAN
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(ii, ji, ki)
     do k = klo, khi - koffset
        do j = jlo, jhi - joffset
           do i = ilo, ihi - ioffset
              ii = i + ioffset; ji = j + joffset; ki = k + koffset
              fgodunov_pre(i,j,k,idim) = -(pgstar(ii,ji,ki)*ds(ii,ji,ki,idim) &
                                         - pgstar(i,j,k)*ds(i,j,k,idim))*dt
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
#endif

#if GEOM == CYLINDRICAL
     if (idim == 1) then
        !$acc kernels loop
        !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(ro, vo, bo, pressure)
        do k = klo, khi - koffset
           do j = jlo, jhi - joffset
              do i = ilo, ihi - ioffset
                 fgodunov_pre(i,j,k,1) = -dv(i,j,k)*(pgstar(i+1,j,k) &
                                                 & - pgstar(i,j,k))/dx
                 ro = half*(rgstar(i,j,k) + rgstar(i+1,j,k))
                 vo = half*(vgstar(i,j,k) + vgstar(i+1,j,k))
                 bo = half*(bgstar(i,j,k) + bgstar(i+1,j,k))
                 ! pressure = (ro*vo*vo - bo*bo)*(ds(i+1,j,k,1) - ds(i,j,k,1))
                 pressure = -bo*bo*(ds(i+1,j,k,1) - ds(i,j,k,1))
                 fgodunov_pre(i,j,k,1) = fgodunov_pre(i,j,k,1) + pressure
                 fgodunov_pre(i,j,k,1) = fgodunov_pre(i,j,k,1)*dt
              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO
     endif

     if ((idim == 2).or.(idim == 3)) then
        !$acc kernels loop
        !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(ii, ji, ki)
        do k = klo, khi - koffset
           do j = jlo, jhi - joffset
              do i = ilo, ihi - ioffset
                 ii = i + ioffset; ji = j + joffset; ki = k + koffset
                 fgodunov_pre(i,j,k,idim) = -(pgstar(ii,ji,ki)*ds(ii,ji,ki,idim)&
                                            - pgstar(i,j,k)*ds(i,j,k,idim))*dt
              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO
     endif
#endif

#if GEOM == SPHERICAL
     if (idim == 1) then
        !$acc kernels loop
        !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(ro, vo, wo, bo, co, pressure)
        do k = klo, khi - koffset
           do j = jlo, jhi - joffset
              do i = ilo, ihi - ioffset
                 fgodunov_pre(i,j,k,1) = -dv(i,j,k)*(pgstar(i+1,j,k) &
                                                 & - pgstar(i,j,k))/dx
                 ro = half*(rgstar(i,j,k) + rgstar(i+1,j,k))
                 vo = half*(vgstar(i,j,k) + vgstar(i+1,j,k))
                 wo = half*(wgstar(i,j,k) + wgstar(i+1,j,k))
                 bo = half*(bgstar(i,j,k) + bgstar(i+1,j,k))
                 co = half*(cgstar(i,j,k) + cgstar(i+1,j,k))
                 pressure = (ro*(vo*vo + wo*wo) - bo*bo - co*co)*(ds(i+1,j,k,1)&
                                                              & - ds(i,j,k,1))
                 fgodunov_pre(i,j,k,1) = fgodunov_pre(i,j,k,1) + pressure
                 fgodunov_pre(i,j,k,1) = fgodunov_pre(i,j,k,1)*dt
              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO
     endif

     if (idim == 2) then
        !$acc kernels loop
        !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(cotanxc, ro, wo, co, pressure)
        do k = klo, khi - koffset
           do j = jlo, jhi - joffset
              cotanxc = cos(y(j))/sin(y(j))
              do i = ilo, ihi - ioffset
                 fgodunov_pre(i,j,k,idim) = -dv(i,j,k)*(pgstar(i,j+1,k) &
                                                 & - pgstar(i,j,k))/dx
                 ro = half*(rgstar(i,j,k) + rgstar(i,j+1,k))
                 wo = half*(wgstar(i,j,k) + wgstar(i,j+1,k))
                 co = half*(cgstar(i,j,k) + cgstar(i,j+1,k))
                 pressure = (ro*wo*wo - co*co)*cotanxc*half*(ds(i,j+1,k,1) &
                                                         & - ds(i,j,k,1))
                 fgodunov_pre(i,j,k,idim) = fgodunov_pre(i,j,k,idim) + pressure
                 fgodunov_pre(i,j,k,idim) = fgodunov_pre(i,j,k,idim)*dt
              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO
     endif

     if (idim == 3) then
        !$acc kernels loop
        !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(ii, ji, ki)
        do k = klo, khi - koffset
           do j = jlo, jhi - joffset
              do i = ilo, ihi - ioffset
                 ii = i + ioffset; ji = j + joffset; ki = k + koffset
                 fgodunov_pre(i,j,k,idim) = -(pgstar(ii,ji,ki)*ds(ii,ji,ki,idim)&
                                            - pgstar(i,j,k)*ds(i,j,k,idim))*dt
              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO
     endif
#endif
  enddo

  !$acc end data

  deallocate(rgstar, ugstar, vgstar, wgstar, bgstar, cgstar, pgstar)

  return
end subroutine riemann
!$py global_template riemann
!===============================================================================
!> Updwind Riemann solver
!===============================================================================
!$py solver_template upwind
subroutine upwind

  real(dp) :: a, sgnm, entho, ecin, ptot, vleft, ploc, proc
  real(dp) :: fl1, fl2, fl3, fl4, fl5, fl6, fl7, fl8
  real(dp) :: cl1, cl2, cl3, cl4, cl5, cl6, cl7, cl8
  real(dp) :: fr1, fr2, fr3, fr4, fr5, fr6, fr7, fr8
  real(dp) :: cr1, cr2, cr3, cr4, cr5, cr6, cr7, cr8

  entho = one/(gamma - one)
              
  ! enforce continuity of normal component
  a    = half*(al + ar)
  sgnm = sign(one, a)
  al   = a; ar = a
              
#if ISO == 1
  ploc = rl*ciso**2
  proc = rr*ciso**2
#else
  ploc = pl
  proc = pr
#endif  

  ! left variables
  ecin = half*(ul*ul + vl*vl + wl*wl)*rl
  emag = half*(a*a + bl*bl + cl*cl)
  etot = ploc*entho + ecin + emag
  ptot = ploc + emag

  fl1 = rl*ul
  fl2 = (etot + ptot)*ul - al*(ul*al + vl*bl + wl*cl)
  fl3 = rl*ul*ul + ptot - al*al
  fl4 = zero
  fl5 = rl*ul*vl - al*bl
  fl6 = bl*ul - al*vl
  fl7 = rl*ul*wl - al*cl
  fl8 = cl*ul - al*wl

  cl1 = rl
  cl2 = etot
  cl3 = rl*ul
  cl4 = al
  cl5 = rl*vl
  cl6 = bl
  cl7 = rl*wl
  cl8 = cl

  ! right variables
  ecin = half*(ur*ur + vr*vr + wr*wr)*rr
  emag = half*(a*a + br*br + cr*cr)
  etot = proc*entho + ecin + emag
  ptot = proc + emag

  fr1 = rr*ur
  fr2 = (etot + ptot)*ur - ar*(ur*ar + vr*br + wr*cr)
  fr3 = rr*ur*ur + ptot - ar*ar
  fr4 = zero
  fr5 = rr*ur*vr - ar*br
  fr6 = br*ur - ar*vr
  fr7 = rr*ur*wr - ar*cr
  fr8 = cr*ur - ar*wr

  cr1 = rr
  cr2 = etot
  cr3 = rr*ur
  cr4 = ar
  cr5 = rr*vr
  cr6 = br
  cr7 = rr*wr
  cr8 = cr

  vleft = half*(ul + ur)

  ! compute the godunov flux
  fgodunov(i,j,k,1,idim) = half*(fl1 + fr1) - abs(vleft)*half*(cr1 - cl1)
  fgodunov(i,j,k,2,idim) = half*(fl2 + fr2) - abs(vleft)*half*(cr2 - cl2)
  fgodunov(i,j,k,3,idim) = half*(fl3 + fr3) - abs(vleft)*half*(cr3 - cl3)
  fgodunov(i,j,k,4,idim) = half*(fl4 + fr4) - abs(vleft)*half*(cr4 - cl4)
  fgodunov(i,j,k,5,idim) = half*(fl5 + fr5) - abs(vleft)*half*(cr5 - cl5)
  fgodunov(i,j,k,6,idim) = half*(fl6 + fr6) - abs(vleft)*half*(cr6 - cl6)
  fgodunov(i,j,k,7,idim) = half*(fl7 + fr7) - abs(vleft)*half*(cr7 - cl7)
  fgodunov(i,j,k,8,idim) = half*(fl8 + fr8) - abs(vleft)*half*(cr8 - cl8)

end subroutine upwind
!$py solver_template upwind
!===============================================================================
!> Lax-Friedrich Riemann solver
!===============================================================================
!$py solver_template llf
subroutine llf

  real(dp) :: a, sgnm, entho, ecin, ptot, ploc, proc, vleft, vright, vel_info
  real(dp) :: fl1, fl2, fl3, fl4, fl5, fl6, fl7, fl8
  real(dp) :: cl1, cl2, cl3, cl4, cl5, cl6, cl7, cl8
  real(dp) :: fr1, fr2, fr3, fr4, fr5, fr6, fr7, fr8
  real(dp) :: cr1, cr2, cr3, cr4, cr5, cr6, cr7, cr8
  real(dp) :: c2, b2, d2, cf

  entho = one/(gamma - one)
              
  ! enforce continuity of normal component
  a    = half*(al + ar)
  sgnm = sign(one, a)
  al   = a; ar = a
              
#if ISO == 1
  ploc = rl*ciso**2
  proc = rr*ciso**2
#else
  ploc = pl
  proc = pr
#endif  

  ! left variables
  ecin = half*(ul*ul + vl*vl + wl*wl)*rl
  emag = half*(a*a + bl*bl + cl*cl)
  etot = ploc*entho + ecin + emag
  ptot = ploc + emag

  fl1 = rl*ul
  fl2 = (etot + ptot)*ul - al*(ul*al + vl*bl + wl*cl)
  fl3 = rl*ul*ul + ptot - al*al
  fl4 = zero
  fl5 = rl*ul*vl - al*bl
  fl6 = bl*ul - al*vl
  fl7 = rl*ul*wl - al*cl
  fl8 = cl*ul - al*wl

  cl1 = rl
  cl2 = etot
  cl3 = rl*ul
  cl4 = al
  cl5 = rl*vl
  cl6 = bl
  cl7 = rl*wl
  cl8 = cl

  c2 = gamma*ploc/rl
  b2 = al*al + bl*bl + cl*cl
  d2 = half*(c2 + b2/rl)
  cf = sqrt(d2 + sqrt(d2**2 - c2*al*al/rl))
  vleft = cf + abs(ul)

  ! right variables
  ecin = half*(ur*ur + vr*vr + wr*wr)*rr
  emag = half*(a*a + br*br + cr*cr)
  etot = proc*entho + ecin + emag
  ptot = proc + emag

  fr1 = rr*ur
  fr2 = (etot + ptot)*ur - ar*(ur*ar + vr*br + wr*cr)
  fr3 = rr*ur*ur + ptot - ar*ar
  fr4 = zero
  fr5 = rr*ur*vr - ar*br
  fr6 = br*ur - ar*vr
  fr7 = rr*ur*wr - ar*cr
  fr8 = cr*ur - ar*wr

  cr1 = rr
  cr2 = etot
  cr3 = rr*ur
  cr4 = ar
  cr5 = rr*vr
  cr6 = br
  cr7 = rr*wr
  cr8 = cr

  c2 = gamma*proc/rr
  b2 = ar*ar + br*br + cr*cr
  d2 = half*(c2 + b2/rr)
  cf = sqrt(d2 + sqrt(d2**2 - c2*ar*ar/rr))
  vright = cf + abs(ur)

  ! compute the godunov flux
  vel_info = max(vleft, vright)
  fgodunov(i,j,k,1,idim) = half*(fl1 + fr1) - vel_info*half*(cr1 - cl1)
  fgodunov(i,j,k,2,idim) = half*(fl2 + fr2) - vel_info*half*(cr2 - cl2)
  fgodunov(i,j,k,3,idim) = half*(fl3 + fr3) - vel_info*half*(cr3 - cl3)
  fgodunov(i,j,k,4,idim) = half*(fl4 + fr4) - vel_info*half*(cr4 - cl4)
  fgodunov(i,j,k,5,idim) = half*(fl5 + fr5) - vel_info*half*(cr5 - cl5)
  fgodunov(i,j,k,6,idim) = half*(fl6 + fr6) - vel_info*half*(cr6 - cl6)
  fgodunov(i,j,k,7,idim) = half*(fl7 + fr7) - vel_info*half*(cr7 - cl7)
  fgodunov(i,j,k,8,idim) = half*(fl8 + fr8) - vel_info*half*(cr8 - cl8)
  
  return
end subroutine llf
!$py solver_template llf
!===============================================================================
!> HLL Riemann solver
!===============================================================================
!$py solver_template hll
subroutine hll

  real(dp) :: a, sgnm, entho, ecin, ptot, ploc, proc, vleft, vright
  real(dp) :: fl1, fl2, fl3, fl4, fl5, fl6, fl7, fl8
  real(dp) :: cl1, cl2, cl3, cl4, cl5, cl6, cl7, cl8
  real(dp) :: fr1, fr2, fr3, fr4, fr5, fr6, fr7, fr8
  real(dp) :: cr1, cr2, cr3, cr4, cr5, cr6, cr7, cr8
  real(dp) :: c2, b2, d2, cf, sl, sr

  entho = one/(gamma - one)
              
  ! enforce continuity of normal component
  a    = half*(al + ar)
  sgnm = sign(one, a)
  al   = a; ar = a
              
#if ISO == 1
  ploc = rl*ciso**2
  proc = rr*ciso**2
#else
  ploc = pl
  proc = pr
#endif  

  ! left variables
  ecin = half*(ul*ul + vl*vl + wl*wl)*rl
  emag = half*(a*a + bl*bl + cl*cl)
  etot = ploc*entho + ecin + emag
  ptot = ploc + emag

  fl1 = rl*ul
  fl2 = (etot + ptot)*ul - al*(ul*al + vl*bl + wl*cl)
  fl3 = rl*ul*ul + ptot - al*al
  fl4 = zero
  fl5 = rl*ul*vl - al*bl
  fl6 = bl*ul - al*vl
  fl7 = rl*ul*wl - al*cl
  fl8 = cl*ul - al*wl

  cl1 = rl
  cl2 = etot
  cl3 = rl*ul
  cl4 = al
  cl5 = rl*vl
  cl6 = bl
  cl7 = rl*wl
  cl8 = cl

  c2 = gamma*ploc/rl
  b2 = al*al + bl*bl + cl*cl
  d2 = half*(c2 + b2/rl)
  cf = sqrt(d2 + sqrt(d2**2 - c2*al*al/rl))
  vleft = cf

  ! right variables
  ecin = half*(ur*ur + vr*vr + wr*wr)*rr
  emag = half*(a*a + br*br + cr*cr)
  etot = proc*entho + ecin + emag
  ptot = proc + emag

  fr1 = rr*ur
  fr2 = (etot + ptot)*ur - ar*(ur*ar + vr*br + wr*cr)
  fr3 = rr*ur*ur + ptot - ar*ar
  fr4 = zero
  fr5 = rr*ur*vr - ar*br
  fr6 = br*ur - ar*vr
  fr7 = rr*ur*wr - ar*cr
  fr8 = cr*ur - ar*wr

  cr1 = rr
  cr2 = etot
  cr3 = rr*ur
  cr4 = ar
  cr5 = rr*vr
  cr6 = br
  cr7 = rr*wr
  cr8 = cr

  c2 = gamma*proc/rr
  b2 = ar*ar + br*br + cr*cr
  d2 = half*(c2 + b2/rr)
  cf = sqrt(d2 + sqrt(d2**2 - c2*ar*ar/rr))
  vright = cf

  sl = min(min(ul, ur) - max(vleft, vright), zero)
  sr = max(max(ul, ur) + max(vleft, vright), zero)

  ! compute the godunov flux
  fgodunov(i,j,k,1,idim) = (sr*fl1 - sl*fr1 + sr*sl*(cr1 - cl1))/(sr - sl)
  fgodunov(i,j,k,2,idim) = (sr*fl2 - sl*fr2 + sr*sl*(cr2 - cl2))/(sr - sl)
  fgodunov(i,j,k,3,idim) = (sr*fl3 - sl*fr3 + sr*sl*(cr3 - cl3))/(sr - sl)
  fgodunov(i,j,k,4,idim) = (sr*fl4 - sl*fr4 + sr*sl*(cr4 - cl4))/(sr - sl)
  fgodunov(i,j,k,5,idim) = (sr*fl5 - sl*fr5 + sr*sl*(cr5 - cl5))/(sr - sl)
  fgodunov(i,j,k,6,idim) = (sr*fl6 - sl*fr6 + sr*sl*(cr6 - cl6))/(sr - sl)
  fgodunov(i,j,k,7,idim) = (sr*fl7 - sl*fr7 + sr*sl*(cr7 - cl7))/(sr - sl)
  fgodunov(i,j,k,8,idim) = (sr*fl8 - sl*fr8 + sr*sl*(cr8 - cl8))/(sr - sl)
  
  return
end subroutine hll
!$py solver_template hll
!===============================================================================
!> HLLD Riemann solver
!===============================================================================
!$py solver_template hlld
subroutine hlld

  real(dp) :: sl, sr, sal, sar
  real(dp) :: entho, a, sgnm, ploc, proc
  real(dp) :: ecinl, emagl, etotl, ptotl, vdotbl, el, cfastl, calfvenl, rcl
  real(dp) :: ecinr, emagr, etotr, ptotr, vdotbr, er, cfastr, calfvenr, rcr
  real(dp) :: rstarl, vstarl, wstarl, bstarl, cstarl, vdotbstarl
  real(dp) :: rstarr, vstarr, wstarr, bstarr, cstarr, vdotbstarr
  real(dp) :: sqrrstarl, etotstarl, etotstarstarl
  real(dp) :: sqrrstarr, etotstarr, etotstarstarr
  real(dp) :: ustar, ptotstar, estar, vstarstar, wstarstar, bstarstar
  real(dp) :: cstarstar, vdotbstarstar
  real(dp) :: etoto, vdotbo
  real(dp) :: c2, b2, d2, cf

  entho = one/(gamma-one)
  
  ! enforce continuity of normal component
  a    = half*(al + ar)
  sgnm = sign(one, a)
  al   = a; ar = a
  
#if ISO == 1
  ploc = rl*ciso**2
  proc = rr*ciso**2
#else
  ploc = pl
  proc = pr
#endif
  
  ! left variables
  ecinl  = half*(ul*ul + vl*vl + wl*wl)*rl
  emagl  = half*(a*a + bl*bl + cl*cl)
  etotl  = ploc*entho + ecinl + emagl
  ptotl  = ploc + emagl
  vdotbl = ul*a + vl*bl + wl*cl
  
  ! right variables
  ecinr  = half*(ur*ur + vr*vr + wr*wr)*rr
  emagr  = half*(a*a + br*br + cr*cr)
  etotr  = proc*entho + ecinr + emagr
  ptotr  = proc + emagr
  vdotbr = ur*a + vr*br + wr*cr

  c2     = gamma*ploc/rl
  b2     = a*a + bl*bl + cl*cl
  d2     = half*(c2 + b2/rl)
  cfastl = sqrt(d2 + sqrt(d2**2 - c2*a*a/rl))

  c2     = gamma*proc/rr
  b2     = a*a + br*br + cr*cr
  d2     = half*(c2 + b2/rr)
  cfastr = sqrt(d2 + sqrt(d2**2 - c2*a*a/rr))

  ! compute hll wave speed
  sl = min(ul, ur) - max(cfastl, cfastr)
  sr = max(ul, ur) + max(cfastl, cfastr)
  
  ! compute lagrangian sound speed
  rcl = rl*(ul - sl)
  rcr = rr*(sr - ur)
  
  ! compute acoustic star state
  ustar    = (rcr*ur + rcl*ul + (ptotl - ptotr))/(rcr + rcl)
  ptotstar = (rcr*ptotl + rcl*ptotr + rcl*rcr*(ul-ur))/(rcr + rcl)
  
  ! left star region variables
  rstarl = rl*(sl - ul)/(sl - ustar)
  estar  = rl*(sl - ul)*(sl - ustar) - a**2
  el     = rl*(sl - ul)*(sl - ul) - a**2
  if (estar == zero) then
     vstarl = vl
     bstarl = bl
     wstarl = wl
     cstarl = cl
  else
     vstarl = vl - a*bl*(ustar - ul)/estar
     bstarl = bl*el/estar
     wstarl = wl - a*cl*(ustar - ul)/estar
     cstarl = cl*el/estar
  endif
  vdotbstarl = ustar*a + vstarl*bstarl + wstarl*cstarl
  etotstarl  = ((sl - ul)*etotl - ptotl*ul + ptotstar*ustar &
             + a*(vdotbl - vdotbstarl))/(sl - ustar)
  sqrrstarl  = sqrt(rstarl)
  calfvenl   = abs(a)/sqrrstarl
  sal        = ustar - calfvenl
  
  ! right star region variables
  rstarr = rr*(sr - ur)/(sr - ustar)
  estar  = rr*(sr - ur)*(sr - ustar) - a**2
  er     = rr*(sr - ur)*(sr - ur) - a**2
  if (estar == zero) then
     vstarr = vr
     bstarr = br
     wstarr = wr
     cstarr = cr
  else
     vstarr = vr - a*br*(ustar - ur)/estar
     bstarr = br*er/estar
     wstarr = wr - a*cr*(ustar - ur)/estar
     cstarr = cr*er/estar
  endif
  vdotbstarr = ustar*a + vstarr*bstarr + wstarr*cstarr
  etotstarr  = ((sr - ur)*etotr - ptotr*ur + ptotstar*ustar &
             + a*(vdotbr - vdotbstarr))/(sr - ustar)
  sqrrstarr  = sqrt(rstarr)
  calfvenr   = abs(a)/sqrrstarr
  sar        = ustar + calfvenr
  
  ! double star region variables
  vstarstar = (sqrrstarl*vstarl + sqrrstarr*vstarr &
          + sgnm*(bstarr - bstarl))/(sqrrstarl + sqrrstarr)
  wstarstar = (sqrrstarl*wstarl + sqrrstarr*wstarr &
          + sgnm*(cstarr - cstarl))/(sqrrstarl + sqrrstarr)
  bstarstar = (sqrrstarl*bstarr + sqrrstarr*bstarl &
          + sgnm*sqrrstarl*sqrrstarr*(vstarr - vstarl))/(sqrrstarl + sqrrstarr)
  cstarstar = (sqrrstarl*cstarr + sqrrstarr*cstarl &
          + sgnm*sqrrstarl*sqrrstarr*(wstarr - wstarl))/(sqrrstarl + sqrrstarr)
  vdotbstarstar = ustar*a + vstarstar*bstarstar + wstarstar*cstarstar
  etotstarstarl = etotstarl - sgnm*sqrrstarl*(vdotbstarl - vdotbstarstar)
  etotstarstarr = etotstarr + sgnm*sqrrstarr*(vdotbstarr - vdotbstarstar)
  
  ! sample the solution at x/t=0
  if (sl > zero) then
     ro = rl
     uo = ul
     vo = vl
     wo = wl
     bo = bl
     co = cl
     ptoto  = ptotl
     etoto  = etotl
     vdotbo = vdotbl
  else if (sal > zero) then
     ro = rstarl
     uo = ustar
     vo = vstarl
     wo = wstarl
     bo = bstarl
     co = cstarl
     ptoto  = ptotstar
     etoto  = etotstarl
     vdotbo = vdotbstarl
  else if (ustar > zero) then
     ro = rstarl
     uo = ustar
     vo = vstarstar
     wo = wstarstar
     bo = bstarstar
     co = cstarstar
     ptoto  = ptotstar
     etoto  = etotstarstarl
     vdotbo = vdotbstarstar
  else if (sar > zero)then
     ro = rstarr
     uo = ustar
     vo = vstarstar
     wo = wstarstar
     bo = bstarstar
     co = cstarstar
     ptoto  = ptotstar
     etoto  = etotstarstarr
     vdotbo = vdotbstarstar
  else if (sr > zero)then
     ro = rstarr
     uo = ustar
     vo = vstarr
     wo = wstarr
     bo = bstarr
     co = cstarr
     ptoto  = ptotstar
     etoto  = etotstarr
     vdotbo = vdotbstarr
  else
     ro = rr
     uo = ur
     vo = vr
     wo = wr
     bo = br
     co = cr
     ptoto  = ptotr
     etoto  = etotr
     vdotbo = vdotbr
  end if
  
  ! compute the godunov flux
  fgodunov(i,j,k,1,idim) = ro*uo
  fgodunov(i,j,k,2,idim) = (etoto + ptoto)*uo - a*vdotbo
  fgodunov(i,j,k,3,idim) = ro*uo*uo - a*a
  fgodunov(i,j,k,4,idim) = zero
  fgodunov(i,j,k,5,idim) = ro*uo*vo - a*bo
  fgodunov(i,j,k,6,idim) = bo*uo - a*vo
  fgodunov(i,j,k,7,idim) = ro*uo*wo - a*co
  fgodunov(i,j,k,8,idim) = co*uo - a*wo
  
  rgstar(i,j,k) = ro
  ugstar(i,j,k) = uo
  vgstar(i,j,k) = vo
  wgstar(i,j,k) = wo
  bgstar(i,j,k) = bo
  cgstar(i,j,k) = co
  pgstar(i,j,k) = ptoto

  return
end subroutine hlld
!$py solver_template hlld
!===============================================================================
!> Acoustic Riemann solver
!===============================================================================
!$py solver_template acoustic
subroutine acoustic

  real(dp) :: smallp
  real(dp) :: sgnm, rho, p, vn, vt1, vt2, btr1, btr2, rsl, rsr, psl, psr
  real(dp) :: wsl, wsr, csl, csr, po, rstar, pstar, ustar, cstar
  real(dp) :: spout, spin, ushock, frac, entho, ecin, ptot, ploc

  smallp = smallr*smallc**2

  ! Initial states pressure, density and velocity
  rsl = max(rl, smallr)
  rsr = max(rr, smallr)

  psl = max(pl, smallp)
  psr = max(pr, smallp)

  ! Acoustic star state
  csl = sqrt(gamma*psl/rsl)
  csr = sqrt(gamma*psr/rsr)
  wsl = csl*rsl
  wsr = csr*rsr
  pstar = ((wsr*psl + wsl*psr) + wsl*wsr*(ul - ur))/(wsl + wsr)
  ustar = ((wsr*ur + wsl*ul) + (psl - psr))/(wsl + wsr)

  ! Left or right going contact wave
  sgnm = sign(one, ustar)

  ! Left or right unperturbed state
  if (sgnm == one) then
     ro = rsl
     po = psl
     uo = ul
     wo = wsl
     co = csl
  else
     ro = rsr
     po = psr
     uo = ur
     wo = wsr
     co = csr
  endif

  ! Star region density and sound speed
  rstar = ro + (pstar - po)/co**2
  rstar = max(rstar, smallr)
  cstar = sqrt(abs(gamma*pstar/rstar))
  cstar = max(cstar, smallc)

  ! Head and tail speed of rarefaction wave
  spout = co - sgnm*uo
  spin  = cstar - sgnm*ustar

  ! Shock speed
  ushock = half*(spin + spout)
  ushock = max(ushock, -sgnm*ustar)
  if (pstar >= po) then
     spout = ushock
     spin = spout
  endif

  ! Sample the solution at x/t = 0
  if (spout < zero) then
     ! Initial state
     rho = ro
     p   = po
     vn  = uo
  else if (spin >= zero) then
     ! Star region
     rho = rstar
     p   = pstar
     vn  = ustar
  else
     ! Rarefaction
     frac = spout/(spout - spin)
     rho = frac*rstar + (one - frac)*ro
     p   = frac*pstar + (one - frac)*po
     vn  = frac*ustar + (one - frac)*uo
  endif

  ! Passive scalars
  if (sgnm == one) then
     vt1 = vl; btr1 = bl
     vt2 = wl; btr2 = cl
  else
     vt1 = vr; btr1 = br
     vt2 = wr; btr2 = cr
  endif

#if ISO == 1
  ploc = rho*ciso**2
#else
  ploc = p
#endif

  ! Local variables
  entho = one/(gamma - one)
  ecin  = half*(vn*vn + vt1*vt1 + vt2*vt2)*rho
  emag  = half*(bn_mean*bn_mean + btr1*btr1 + btr2*btr2)
  etot  = ploc*entho + ecin + emag
  ptot  = ploc + emag
  
  ! Compute fluxes
  fgodunov(i,j,k,1,idim) = rho*vn
  fgodunov(i,j,k,2,idim) = (etot + ptot)*vn &
                       & - bn_mean*(vn*bn_mean + vt1*btr1 + vt2*btr2)
  fgodunov(i,j,k,3,idim) = rho*vn*vn + ptot - bn_mean*bn_mean
  fgodunov(i,j,k,4,idim) = zero
  fgodunov(i,j,k,5,idim) = rho*vn*vt1 - bn_mean*btr1
  fgodunov(i,j,k,6,idim) = btr1*vn - bn_mean*vt1
  fgodunov(i,j,k,7,idim) = rho*vn*vt2 - bn_mean*btr2
  fgodunov(i,j,k,8,idim) = btr2*vn - bn_mean*vt2

  return
end subroutine acoustic
!$py solver_template acoustic
!===============================================================================
!> Athena Roe Riemann solver
!===============================================================================
!subroutine athena_roe(qleft, qright, flux, zero_flux)
!  use params
!  implicit none
!
!  real(dp), dimension(nvar), intent(out) :: flux
!  real(dp), dimension(nvar), intent(in)  :: qleft, qright
!  real(dp), intent(in) :: zero_flux
!  real(dp), dimension(nvar)          :: cvarleft, cvarright, cvardiff
!  real(dp), dimension(nvar)          :: fleft, fright
!  real(dp), dimension(nvar-1,nvar-1) :: lem, rem
!  real(dp), dimension(nvar-1)        :: lambda, lambdal, lambdar, a
!
!  integer  :: ivar
!  real(dp) :: rl, pl, ul, vl, wl, al, bl, cl, el, mxl, myl, mzl, pbl, velleft
!  real(dp) :: rr, pr, ur, vr, wr, ar, br, cr, er, mxr, myr, mzr, pbr, velright
!  real(dp) :: sqrtrl, sqrtrr, rroe, uroe, vroe, wroe, broe, croe, hroe
!  real(dp) :: xfactor, yfactor, rim, mxm, mym, mzm, eim, bim, cim, etm, l1, l2
!  real(dp) :: fluxr, fluxe, fluxmx, fluxmy, fluxb, fluxmz, fluxc, coeff
!  logical  :: llf
!
!  ! First step:
!  ! Compute fluxes and conserved variables
!  rl = qleft(1); pl = qleft(2)
!  ul = qleft(3); vl = qleft(5)  ; wl = qleft(7)
!  al = half*(qleft(4)+qright(4)); bl = qleft(6); cl = qleft(8)
!
!  call comp_mhd_flux(rl, pl, ul, vl, wl, al, bl, cl, cvarleft, fleft)
!
!  rr = qright(1); pr = qright(2)
!  ur = qright(3); vr = qright(5); wr = qright(7)
!  ar = al       ; br = qright(6); cr = qright(8)
!
!  call comp_mhd_flux(rr, pr, ur, vr, wr, ar, br, cr, cvarright, fright)
!
!  ! Define explicitely conserved quantities
!  el  = cvarleft(2); er  = cvarright(2)
!  mxl = cvarleft(3); mxr = cvarright(3)
!  myl = cvarleft(5); myr = cvarright(5)
!  mzl = cvarleft(7); mzr = cvarright(7)
!
!  ! Magnetic pressure
!  pbl = half*(al*al + bl*bl + cl*cl)
!  pbr = half*(ar*ar + br*br + cr*cr)
!
!  ! Second step:
!  ! Compute Roe-averaged data from left and right states
!  ! These averages will be used as input to the eigen problem
!  sqrtrl = sqrt(rl)
!  sqrtrr = sqrt(rr)
!  rroe   = sqrtrl*sqrtrr
!  uroe   = (sqrtrl*ul + sqrtrr*ur)/(sqrtrl + sqrtrr)
!  vroe   = (sqrtrl*vl + sqrtrr*vr)/(sqrtrl + sqrtrr)
!  wroe   = (sqrtrl*wl + sqrtrr*wr)/(sqrtrl + sqrtrr)
!  broe   = (sqrtrr*bl + sqrtrl*br)/(sqrtrl + sqrtrr)
!  croe   = (sqrtrr*cl + sqrtrl*cr)/(sqrtrl + sqrtrr)
!  hroe   = ((el + pl + pbl)/sqrtrl + (er + pr + pbr)/sqrtrr)/(sqrtrl + sqrtrr)
!
!  xfactor = ((broe*broe - bl*br) + (croe*croe - cl*cr))/(2*rroe)
!  yfactor = (rl + rr)/(2*rroe)
!
!  ! Third step:
!  ! Compute eigenvaules and eigenmatrices from Roe-averaged values
!  call roe_eigenvalues(rroe, uroe, vroe, wroe, hroe, al, broe, croe &
!       & , xfactor, yfactor, lambda, rem, lem)
!
!  ! Compute eigenvalues from left and right states
!  call mhd_eigenvalues(rl, ul, vl, wl, pl, al, bl, cl, lambdal)
!  call mhd_eigenvalues(rr, ur, vr, wr, pr, ar, br, cr, lambdar)
!
!  ! Fourth step:
!  ! Create intermediate states from eigenmatrices
!  ! Remember: convention is (r, u, v,  w, e, b, c)
!  a = 0
!  a = a + (rr - rl)*lem(1,:) 
!  a = a + (mxr - mxl)*lem(2,:) 
!  a = a + (myr - myl)*lem(3,:)
!  a = a + (mzr - mzl)*lem(4,:) 
!  a = a + (er - el)*lem(5,:) 
!  a = a + (br - bl)*lem(6,:)
!  a = a + (cr - cl)*lem(7,:)
!
!  llf = .false.
!  rim = rl; mxm = mxl; mym = myl; mzm = mzl; eim = el; bim = bl; cim = cl
!  do ivar = 1, nvar-1
!     rim = rim + a(ivar)*rem(ivar,1)
!     mxm = mxm + a(ivar)*rem(ivar,2)
!     mym = mym + a(ivar)*rem(ivar,3)
!     mzm = mzm + a(ivar)*rem(ivar,4)
!     eim = eim + a(ivar)*rem(ivar,5)
!     bim = bim + a(ivar)*rem(ivar,6)
!     cim = cim + a(ivar)*rem(ivar,7)
!     etm = eim - half*(mxm*mxm + mym*mym + mzm*mzm)/rim &
!             & - half*(al*al + bim*bim + cim*cim)
!     if ((rim <= zero) .or. (etm <= zero)) llf = .true.
!     if (llf) print "(/ 'Problem in Roe: density or pressure is negative!', / &
!          & , 'd_int: ', E13.5, ', e_int: ', E13.5, ', d_left: ', E13.5 &
!          & , 'p_left: ', E13.5, /)", rim, etm, rl, pl
!  enddo
!
!  if (llf) then
!     flux = half*(fright + fleft)*zero_flux
!     call comp_speed(rl, pl, ul, al, bl, cl, velleft)
!     call comp_speed(rr, pr, ur, ar, br, cr, velright)
!     cvardiff = half*(cvarright - cvarleft)
!     flux    = flux - max(velleft, velright)*cvardiff
!     return
!  endif
!
!  ! Fifth step:
!  ! Entropy fix for genuinely non-linear waves
!  do ivar = 1, nvar-1, 2
!     l1 = min(lambdal(ivar), lambda(ivar))
!     l2 = max(lambdar(ivar), lambda(ivar))
!     if ((l1 < zero) .and. (l2 > zero)) then
!        lambda(ivar) = (lambda(ivar)*(l2 + l1) - two*l2*l1)/(l2 - l1)
!     endif
!  enddo
!
!  ! Sixth step:
!  ! Compute fluxes at interface using Roe solver
!
!  ! Add left and right fluxes
!  ! Remember convention: (r, e, u, a, v, b, w, c)
!  fleft  = fleft*zero_flux
!  fright = fright*zero_flux
!  fluxr  = fleft(1) + fright(1)
!  fluxe  = fleft(2) + fright(2)
!  fluxmx = fleft(3) + fright(3)
!  fluxmy = fleft(5) + fright(5)
!  fluxb  = fleft(6) + fright(6)
!  fluxmz = fleft(7) + fright(7)
!  fluxc  = fleft(8) + fright(8)
!  
!  ! Compute Roe fluxes
!  do ivar = 1, nvar-1
!     coeff  = abs(lambda(ivar))*a(ivar)
!     fluxr  = fluxr  - coeff*rem(ivar,1)
!     fluxe  = fluxe  - coeff*rem(ivar,5)
!     fluxmx = fluxmx - coeff*rem(ivar,2)
!     fluxmy = fluxmy - coeff*rem(ivar,3)
!     fluxb  = fluxb  - coeff*rem(ivar,6)
!     fluxmz = fluxmz - coeff*rem(ivar,4)
!     fluxc  = fluxc  - coeff*rem(ivar,7)
!  enddo
!
!  ! Take half of it and store it in flux
!  flux(1) = half*fluxr  
!  flux(2) = half*fluxe  
!  flux(3) = half*fluxmx 
!  flux(4) = zero
!  flux(5) = half*fluxmy 
!  flux(6) = half*fluxb  
!  flux(7) = half*fluxmz 
!  flux(8) = half*fluxc  
!
!  return
!end subroutine athena_roe
!!===============================================================================
!!> Compute MHD adiabatic eigenvalues
!!===============================================================================
!subroutine mhd_eigenvalues(r, u, v, w, p, a, b, c, lambda)
!  use params
!  implicit none
!
!  real(dp), intent(in) :: r, u, v, w, p, a, b, c
!  real(dp), dimension(nvar-1), intent(out) :: lambda
!
!  real(dp) :: vsq, btsq, bt, vaxsq, vax, cssq, astarsq
!  real(dp) :: cfastsq, cfast, cslowsq, cslow
!
!  vsq   = u**2 + v**2 + w**2
!  btsq  = b**2 + c**2
!  bt    = sqrt(btsq)
!  vaxsq = a**2/r
!  vax   = sqrt(vaxsq)
!  cssq  = gamma*p/r
!  cssq  = max(cssq, smallc**2)
!  astarsq = cssq + vaxsq + btsq/r
!
!  cfastsq = half*(astarsq + sqrt(astarsq**2 - four*cssq*vaxsq))
!  cfast   = sqrt(cfastsq)
!
!  cslowsq = half*(astarsq - sqrt(astarsq**2 - four*cssq*vaxsq))
!  if (cslowsq <= zero) cslowsq = zero
!  cslow   = sqrt(cslowsq)
!
!  lambda(1) = u - cfast
!  lambda(2) = u - vax
!  lambda(3) = u - cslow
!  lambda(4) = u
!  lambda(5) = u + cslow
!  lambda(6) = u + vax
!  lambda(7) = u + cfast
!
!  return
!end subroutine mhd_eigenvalues
!!===============================================================================
!!> Compute Roe adiabatic eigenvalues
!!===============================================================================
!subroutine roe_eigenvalues(r, u, v, w, h, a, b, c , xfactor, yfactor, lambda &
!     & , rem, lem)
!  use params
!  implicit none
!
!  real(dp), intent(in) :: r, u, v, w, h, a, b, c, xfactor, yfactor
!  real(dp), dimension(nvar-1), intent(out)        :: lambda
!  real(dp), dimension(nvar-1,nvar-1), intent(out) :: rem, lem
!
!  real(dp) :: vsq, btsq, bt, btstarsq, btstar, vaxsq, vax, hp, twidasq, qstarsq
!  real(dp) :: cfastsq, cfast, cslowsq, cslow, betay, betaz, betaystar, betazstar
!  real(dp) :: betastarsq, vbeta, alphaf, alphas, rroot, s, twida, qfast, qslow
!  real(dp) :: afprime, asprime, afpbb, aspbb, na, cff, css, af, as, afpb, aspb
!  real(dp) :: gammana, qystar, qzstar, vqstar, norm
!
!  vsq      = u**2 + v**2 + w**2
!  btsq     = b**2 + c**2
!  bt       = sqrt(btsq)
!  btstarsq = (gamma - one - (gamma - two)*yfactor)*btsq
!  btstar   = sqrt(btstarsq)
!  vaxsq    = a**2/r
!  vax      = sqrt(vaxsq)
!
!  hp = h - (vaxsq + btsq/r)
!  twidasq = ((gamma - one)*(hp - half*vsq) - (gamma - two)*xfactor)
!  twidasq = max(twidasq, smallc**2)
!  qstarsq = twidasq + (vaxsq + btstarsq/r)
!
!  cfastsq = half*(qstarsq + sqrt(qstarsq**2 - four*twidasq*vaxsq))
!  cfast   = sqrt(cfastsq)
!
!  cslowsq = half*(qstarsq - sqrt(qstarsq**2 - four*twidasq*vaxsq))
!  if (cslowsq <= zero) cslowsq = zero
!  cslow   = sqrt(cslowsq)
!
!  if (bt == zero) then
!     betay = half*sqrt(two)
!     betaz = half*sqrt(two)
!     betaystar = betay
!     betazstar = betaz
!  else
!     betay = b/bt
!     betaz = c/bt
!     betaystar = b/btstar
!     betazstar = c/btstar
!  endif
!  betastarsq = betaystar**2 + betazstar**2
!  vbeta = v*betaystar + w*betazstar
!
!  if ((cfastsq - cslowsq) == zero) then
!     alphaf = one
!     alphas = zero
!  else if ((twidasq - cslowsq) <= zero) then
!     alphaf = zero
!     alphas = one
!  else if ((cfastsq - twidasq) <= zero) then
!     alphaf = one
!     alphas = zero
!  else
!     alphaf = sqrt((twidasq - cslowsq)/(cfastsq - cslowsq))
!     alphas = sqrt((cfastsq - twidasq)/(cfastsq - cslowsq))
!  endif
!
!  ! Compute qs and as for eigenmatrices
!  rroot   = sqrt(r)
!  s       = sign(one, a)
!  twida   = sqrt(twidasq)
!  qfast   = s*cfast*alphaf
!  qslow   = s*cslow*alphas
!  afprime = twida*alphaf/rroot
!  asprime = twida*alphas/rroot
!  afpbb   = afprime*btstar*betastarsq
!  aspbb   = asprime*btstar*betastarsq
!
!  ! Eigenvalues
!  lambda(1) = u - cfast
!  lambda(2) = u - vax
!  lambda(3) = u - cslow
!  lambda(4) = u
!  lambda(5) = u + cslow
!  lambda(6) = u + vax
!  lambda(7) = u + cfast
!
!  ! Right eigenmatrix
!  rem(1,1) = alphaf
!  rem(1,2) = alphaf*(u - cfast)
!  rem(1,3) = alphaf*v + qslow*betaystar
!  rem(1,4) = alphaf*w + qslow*betazstar
!  rem(1,5) = alphaf*(hp - u*cfast) + qslow*vbeta + aspbb
!  rem(1,6) = asprime*betaystar
!  rem(1,7) = asprime*betazstar
!
!  rem(2,1) = zero
!  rem(2,2) = zero
!  rem(2,3) = -betaz
!  rem(2,4) = betay
!  rem(2,5) = -(v*betaz - w*betay)
!  rem(2,6) = -s*betaz/rroot
!  rem(2,7) = s*betay/rroot
!
!  rem(3,1) = alphas
!  rem(3,2) = alphas*(u - cslow)
!  rem(3,3) = alphas*v - qfast*betaystar
!  rem(3,4) = alphas*w - qfast*betazstar
!  rem(3,5) = alphas*(hp - u*cslow) - qfast*vbeta - afpbb
!  rem(3,6) = -afprime*betaystar
!  rem(3,7) = -afprime*betazstar
!
!  rem(4,1) = one
!  rem(4,2) = u
!  rem(4,3) = v
!  rem(4,4) = w
!  rem(4,5) = half*vsq + (gamma - two)*xfactor/(gamma - one)
!  rem(4,6) = zero
!  rem(4,7) = zero
!
!  rem(5,1) = alphas
!  rem(5,2) = alphas*(u + cslow)
!  rem(5,3) = alphas*v + qfast*betaystar
!  rem(5,4) = alphas*w + qfast*betazstar
!  rem(5,5) = alphas*(hp + u*cslow) + qfast*vbeta - afpbb
!  rem(5,6) = -afprime*betaystar
!  rem(5,7) = -afprime*betazstar
!
!  rem(6,1) = zero
!  rem(6,2) = zero
!  rem(6,3) = betaz
!  rem(6,4) = -betay
!  rem(6,5) = (v*betaz - w*betay)
!  rem(6,6) = -s*betaz/rroot
!  rem(6,7) = s*betay/rroot
!
!  rem(7,1) = alphaf
!  rem(7,2) = alphaf*(u + cfast)
!  rem(7,3) = alphaf*v - qslow*betaystar
!  rem(7,4) = alphaf*w - qslow*betazstar
!  rem(7,5) = alphaf*(hp + u*cfast) - qslow*vbeta + aspbb
!  rem(7,6) = asprime*betaystar
!  rem(7,7) = asprime*betazstar
!
!  ! Left eigenmatrix
!  ! Start by normalizing some quantities by 1/(2a**2) or (gamma-1)/(2a**2)
!  na    = half/twidasq
!  cff   = na*alphaf*cfast
!  css   = na*alphas*cslow
!  qfast = qfast*na
!  qslow = qslow*na
!  af    = na*afprime*r
!  as    = na*asprime*r
!  afpb  = na*afprime*btstar
!  aspb  = na*asprime*btstar
!
!  gammana = (gamma - one)*na
!  alphaf  = gammana*alphaf
!  alphas  = gammana*alphas
!  qystar  = betaystar/betastarsq
!  qzstar  = betazstar/betastarsq
!  vqstar  = (v*qystar + w*qzstar)
!  norm    = two*gammana
!
!  lem(1,1) = alphaf*(vsq - hp) + cff*(cfast + u) - qslow*vqstar - aspb
!  lem(2,1) = -alphaf*u - cff
!  lem(3,1) = -alphaf*v + qslow*qystar
!  lem(4,1) = -alphaf*w + qslow*qzstar
!  lem(5,1) = alphaf
!  lem(6,1) = as*qystar - alphaf*b
!  lem(7,1) = as*qzstar - alphaf*c
! 
!  lem(1,2) = half*(v*betaz - w*betay)
!  lem(2,2) = zero
!  lem(3,2) = -half*betaz
!  lem(4,2) = half*betay
!  lem(5,2) = zero
!  lem(6,2) = -half*rroot*betaz*s
!  lem(7,2) = half*rroot*betay*s
!
!  lem(1,3) = alphas*(vsq - hp) + css*(cslow + u) + qfast*vqstar + afpb
!  lem(2,3) = -alphas*u - css
!  lem(3,3) = -alphas*v - qfast*qystar
!  lem(4,3) = -alphas*w - qfast*qzstar
!  lem(5,3) = alphas
!  lem(6,3) = -af*qystar - alphas*b
!  lem(7,3) = -af*qzstar - alphas*c
!
!  ! CAUTION! There is a difference in sign compared to RAMSES
!  lem(1,4) = one - norm*(half*vsq - (gamma - two)*xfactor/(gamma - one))
!  ! Old version:
!  ! lem(1,4) = one - norm*(half*vsq + (gamma - two)*xfactor/(gamma - one))
!  lem(2,4) = norm*u
!  lem(3,4) = norm*v
!  lem(4,4) = norm*w
!  lem(5,4) = -norm
!  lem(6,4) = norm*b
!  lem(7,4) = norm*c
!
!  lem(1,5) = alphas*(vsq - hp) + css*(cslow - u) - qfast*vqstar + afpb
!  lem(2,5) = -alphas*u + css
!  lem(3,5) = -alphas*v + qfast*qystar
!  lem(4,5) = -alphas*w + qfast*qzstar
!  lem(5,5) = alphas
!  lem(6,5) = -af*qystar - alphas*b
!  lem(7,5) = -af*qzstar - alphas*c
!
!  lem(1,6) = -half*(v*betaz - w*betay)
!  lem(2,6) = zero
!  lem(3,6) = half*betaz
!  lem(4,6) = -half*betay
!  lem(5,6) = zero
!  lem(6,6) = -half*rroot*betaz*s
!  lem(7,6) = half*rroot*betay*s
!  
!  lem(1,7) = alphaf*(vsq - hp) + cff*(cfast - u) + qslow*vqstar - aspb
!  lem(2,7) = -alphaf*u + cff
!  lem(3,7) = -alphaf*v - qslow*qystar
!  lem(4,7) = -alphaf*w - qslow*qzstar
!  lem(5,7) = alphaf
!  lem(6,7) = as*qystar - alphaf*b
!  lem(7,7) = as*qzstar - alphaf*c
!
!  return
!end subroutine roe_eigenvalues
