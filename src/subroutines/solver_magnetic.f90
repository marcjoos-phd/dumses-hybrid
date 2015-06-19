!===============================================================================
!> \file solver_magnetic.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is Riemann magnetic solver subroutines. This file cannot be compiled as 
!! is, it needs to be preprocessed with the DUMSES preprocessor. 
!! riemann_magnetic subroutine is a global_template for solver, and the other 
!! subroutines are solver_template.
!! \details
!! Contains riemann_magnetic(), upwind(), llf(), hll(), hlld(), acoustic(),
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
!! \b created:          01-12-2014 
!! \b last \b modified: 04-28-2015
!<
!===============================================================================
!> 2D Riemann solver to compute EMF at cell edges
!===============================================================================
!$py global_template riemann_magnetic
subroutine riemann_magnetic(qRT, qRB, qLT, qLB, emf, idim)
  use params
  use variables, only: x, dt
  use oacc_params
  implicit none

  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar,3), intent(in) :: qRT, qRB
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nvar,3), intent(in) :: qLT, qLB
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2), intent(out)       :: emf
  integer, intent(in) :: idim

  ! local variables
  real(dp) :: rLL, rLR, rRL, rRR, pLL, pLR, pRL, pRR, uLL, uLR, uRL, uRR
  real(dp) :: vLL, vLR, vRL, vRR, ALL, ALR, ARL, ARR, BLL, BLR, BRL, BRR
  real(dp) :: CLL, CLR, CRL, CRR, ELL, ELR, ERL, ERR
  real(dp) :: emf_tmp, shear
  integer  :: i, j, k
  integer  :: ilo, ihi, jlo, jhi, klo, khi
  integer  :: ln, lt1, lt2, bn, bt1, bt2
  integer  :: irtoffset, irboffset, iltoffset, ilboffset
  integer  :: jrtoffset, jrboffset, jltoffset, jlboffset
  integer  :: krtoffset, krboffset, kltoffset, klboffset
  integer  :: irt, irb, ilt, ilb, jrt, jrb, jlt, jlb, krt, krb, klt, klb

  !$py begin_statement

  if(idim == 1) then
     ilo = min(1, iu1+3); ihi = max(1, iu2-3)
     jlo = jf1          ; jhi = jf2
     klo = kf1          ; khi = kf2
     ln = 2; lt1 = 3; lt2 = 4
     bn = 6; bt1 = 7; bt2 = 8
     irtoffset = 0; irboffset = 0; iltoffset = 0; ilboffset = 0
     jrtoffset = -1; jrboffset = -1; jltoffset = 0; jlboffset = 0
     krtoffset = -1; krboffset = 0; kltoffset = -1; klboffset = 0
  else if(idim == 2) then
     ilo = if1          ; ihi = if2
     jlo = min(1, ju1+3); jhi = max(1, ju2-3)
     klo = kf1          ; khi = kf2
     ln = 3; lt1 = 4; lt2 = 2
     bn = 7; bt1 = 8; bt2 = 6
     irtoffset = -1; irboffset = 0; iltoffset = -1; ilboffset = 0
     jrtoffset = 0; jrboffset = 0; jltoffset = 0; jlboffset = 0
     krtoffset = -1; krboffset = -1; kltoffset = 0; klboffset = 0
  else
     ilo = if1          ; ihi = if2
     jlo = jf1          ; jhi = jf2
     klo = min(1, ku1+3); khi = max(1, ku2-3)
     ln = 4; lt1 = 2; lt2 = 3
     bn = 8; bt1 = 6; bt2 = 7
     irtoffset = -1; irboffset = -1; iltoffset = 0; ilboffset = 0
     jrtoffset = -1; jrboffset = 0; jltoffset = -1; jlboffset = 0
     krtoffset = 0; krboffset = 0; kltoffset = 0; klboffset = 0
  endif

  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(irt, ilt, jrt, jlt, krt, klt) &
  !$OMP PRIVATE(rLL, rLR, rRL, rRR, pLL, pLR, pRL, pRR, uLL, uLR, uRL, uRR) &
  !$OMP PRIVATE(vLL, vLR, vRL, vRR, ALL, ALR, ARL, ARR, BLL, BLR, BRL, BRR) &
  !$OMP PRIVATE(CLL, CLR, CRL, CRR, ELL, ELR, ERL, ERR, emf_tmp, shear)
  do k = klo, khi
     !$acc loop vector(blocky_solver_mag)
     do j = jlo, jhi
        !$acc loop vector(blockx_solver_mag)
        do i = ilo, ihi
           irt = i + irtoffset; irb = i + irboffset
           ilt = i + iltoffset; ilb = i + ilboffset
           jrt = j + jrtoffset; jrb = j + jrboffset
           jlt = j + jltoffset; jlb = j + jlboffset
           krt = k + krtoffset; krb = k + krboffset
           klt = k + kltoffset; klb = k + klboffset

           ! State vectors; 1=LL, 2=RL, 3=LR, 4=RR
           rLL = qRT(irt,jrt,krt,1,idim)
           rRL = qLT(ilt,jlt,klt,1,idim)
           rLR = qRB(irb,jrb,krb,1,idim)
           rRR = qLB(ilb,jlb,klb,1,idim)
           
#if ISO == 1
           pLL = rLL*ciso**2
           pRL = rRL*ciso**2
           pLR = rLR*ciso**2
           pRR = rRR*ciso**2
#else
           pLL = qRT(irt,jrt,krt,5,idim)
           pRL = qLT(ilt,jlt,klt,5,idim)
           pLR = qRB(irb,jrb,krb,5,idim)
           pRR = qLB(ilb,jlb,klb,5,idim)
#endif
           
           ! First parallel velocity 
           uLL = qRT(irt,jrt,krt,lt1,idim)
           uRL = qLT(ilt,jlt,klt,lt1,idim)
           uLR = qRB(irb,jrb,krb,lt1,idim)
           uRR = qLB(ilb,jlb,klb,lt1,idim)
           
           ! Second parallel velocity 
           vLL = qRT(irt,jrt,krt,lt2,idim)
           vRL = qLT(ilt,jlt,klt,lt2,idim)
           vLR = qRB(irb,jrb,krb,lt2,idim)
           vRR = qLB(ilb,jlb,klb,lt2,idim)
   
           ! First parallel magnetic field (enforce continuity)
           ALL = half*(qRT(irt,jrt,krt,bt1,idim) + qLT(ilt,jlt,klt,bt1,idim))
           ARL = half*(qRT(irt,jrt,krt,bt1,idim) + qLT(ilt,jlt,klt,bt1,idim))
           ALR = half*(qRB(irb,jrb,krb,bt1,idim) + qLB(ilb,jlb,klb,bt1,idim))
           ARR = half*(qRB(irb,jrb,krb,bt1,idim) + qLB(ilb,jlb,klb,bt1,idim))
   
           ! Second parallel magnetic field (enforce continuity)
           BLL = half*(qRT(irt,jrt,krt,bt2,idim) + qRB(irb,jrb,krb,bt2,idim))
           BRL = half*(qLT(ilt,jlt,klt,bt2,idim) + qLB(ilb,jlb,klb,bt2,idim))
           BLR = half*(qRT(irt,jrt,krt,bt2,idim) + qRB(irb,jrb,krb,bt2,idim))
           BRR = half*(qLT(ilt,jlt,klt,bt2,idim) + qLB(ilb,jlb,klb,bt2,idim))
   
           ! Orthogonal magnetic Field
           CLL = qRT(irt,jrt,krt,bn,idim)
           CRL = qLT(ilt,jlt,klt,bn,idim)
           CLR = qRB(irb,jrb,krb,bn,idim)
           CRR = qLB(ilb,jlb,klb,bn,idim)
   
           ! Compute final fluxes
           ! vx*by - vy*bx at the four edge centers
           ELL = uLL*BLL - vLL*ALL
           ERL = uRL*BRL - vRL*ARL
           ELR = uLR*BLR - vLR*ALR
           ERR = uRR*BRR - vRR*ARR

           !$py insert_solver

           ! Upwind solver in case of the shearing box
#if GEOM == CARTESIAN
           if ((Omega0 > zero) .and. (.not. fargo)) then
              if (idim == 1) then
                 shear = -1.5d0*Omega0*x(i)
                 if (shear > zero) then
                    emf_tmp = emf_tmp + shear*BLL*dt
                 else
                    emf_tmp = emf_tmp + shear*BRR*dt
                 endif
              endif
              if (idim == 3) then
                 shear = -1.5d0*Omega0*half*(x(i) + x(i-1))
                 if (shear > zero) then
                    emf_tmp = emf_tmp - shear*ALL*dt
                 else
                    emf_tmp = emf_tmp - shear*ARR*dt
                 endif
              endif
           endif
#endif
           emf(i,j,k) = emf_tmp
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  return
end subroutine riemann_magnetic
!$py global_template riemann_magnetic
!===============================================================================
!> 2D HLLD Riemann solver
!===============================================================================
!$py solver_template hlld
subroutine hlld

  real(dp) :: SL, SR, SB, ST, SAL, SAR, SAT, SAB, E, c2, b2, d2
  real(dp) :: cLLx, cRLx, cLRx, cRRx, cLLy, cRLy, cLRy, cRRy
  real(dp) :: cfastLLx, cfastRLx, cfastLRx, cfastRRx, cfastLLy, cfastRLy
  real(dp) :: cfastLRy, cfastRRy
  real(dp) :: calfvenR, calfvenL, calfvenT, calfvenB
  real(dp) :: vLLx, vRLx, vLRx, vRRx, vLLy, vRLy, vLRy, vRRy
  real(dp) :: PtotLL, PtotLR, PtotRL, PtotRR, rcLLx, rcLRx, rcRLx, rcRRx
  real(dp) :: rcLLy, rcLRy, rcRLy, rcRRy
  real(dp) :: ustar, vstar, rstarLLx, rstarLRx, rstarRLx, rstarRRx, rstarLLy
  real(dp) :: rstarLRy, rstarRLy, rstarRRy
  real(dp) :: rstarLL, rstarLR, rstarRL, rstarRR, AstarLL, AstarLR, AstarRL
  real(dp) :: AstarRR, BstarLL, BstarLR, BstarRL, BstarRR
  real(dp) :: EstarLLx, EstarLRx, EstarRLx, EstarRRx, EstarLLy, EstarLRy
  real(dp) :: EstarRLy, EstarRRy, EstarLL, EstarLR, EstarRL, EstarRR
  real(dp) :: AstarT, AstarB, BstarR, BstarL
  real(dp) :: SLustar, SRustar, STvstar, SBvstar
  real(dp) :: SARSAL, SATSAB, sqrtrstarLL, sqrtrstarRR, sqrtrstarRL, sqrtrstarLR

  ! Compute 4 fast magnetosonic velocity relative to x direction
  c2       = gamma*pLL/rLL
  b2       = ALL*ALL + BLL*BLL + CLL*CLL
  d2       = half*(c2 + b2/rLL)
  cfastLLx = sqrt(d2 + sqrt(d2**2 - c2*ALL*ALL/rLL))

  c2       = gamma*pLR/rLR
  b2       = ALR*ALR + BLR*BLR + CLR*CLR
  d2       = half*(c2 + b2/rLR)
  cfastLRx = sqrt(d2 + sqrt(d2**2 - c2*ALR*ALR/rLR))

  c2       = gamma*pRL/rRL
  b2       = ARL*ARL + BRL*BRL + CRL*CRL
  d2       = half*(c2 + b2/rRL)
  cfastRLx = sqrt(d2 + sqrt(d2**2 - c2*ARL*ARL/rRL))

  c2       = gamma*pRR/rRR
  b2       = ARR*ARR + BRR*BRR + CRR*CRR
  d2       = half*(c2 + b2/rRR)
  cfastRRx = sqrt(d2 + sqrt(d2**2 - c2*ARR*ARR/rRR))

  ! Compute 4 fast magnetosonic velocity relative to y direction
  c2       = gamma*pLL/rLL
  b2       = BLL*BLL + ALL*ALL + CLL*CLL
  d2       = half*(c2 + b2/rLL)
  cfastLLy = sqrt(d2 + sqrt(d2**2 - c2*BLL*BLL/rLL))

  c2       = gamma*pLR/rLR
  b2       = BLR*BLR + ALR*ALR + CLR*CLR
  d2       = half*(c2 + b2/rLR)
  cfastLRy = sqrt(d2 + sqrt(d2**2 - c2*BLR*BLR/rLR))

  c2       = gamma*pRL/rRL
  b2       = BRL*BRL + ARL*ARL + CRL*CRL
  d2       = half*(c2 + b2/rRL)
  cfastRLy = sqrt(d2 + sqrt(d2**2 - c2*BRL*BRL/rRL))

  c2       = gamma*pRR/rRR
  b2       = BRR*BRR + ARR*ARR + CRR*CRR
  d2       = half*(c2 + b2/rRR)
  cfastRRy = sqrt(d2 + sqrt(d2**2 - c2*BRR*BRR/rRR))
  
  SL = min(uLL,uLR,uRL,uRR) - max(cfastLLx,cfastLRx,cfastRLx,cfastRRx)
  SR = max(uLL,uLR,uRL,uRR) + max(cfastLLx,cfastLRx,cfastRLx,cfastRRx)
  SB = min(vLL,vLR,vRL,vRR) - max(cfastLLy,cfastLRy,cfastRLy,cfastRRy)
  ST = max(vLL,vLR,vRL,vRR) + max(cfastLLy,cfastLRy,cfastRLy,cfastRRy)
  
  ELL = uLL*BLL - vLL*ALL
  ELR = uLR*BLR - vLR*ALR
  ERL = uRL*BRL - vRL*ARL
  ERR = uRR*BRR - vRR*ARR
  
  PtotLL = pLL + half*(ALL*ALL + BLL*BLL + CLL*CLL)
  PtotLR = pLR + half*(ALR*ALR + BLR*BLR + CLR*CLR)
  PtotRL = pRL + half*(ARL*ARL + BRL*BRL + CRL*CRL)
  PtotRR = pRR + half*(ARR*ARR + BRR*BRR + CRR*CRR)
  
  rcLLx = rLL*(uLL - SL); rcRLx = rRL*(SR - uRL) 
  rcLRx = rLR*(uLR - SL); rcRRx = rRR*(SR - uRR)
  rcLLy = rLL*(vLL - SB); rcLRy = rLR*(ST - vLR) 
  rcRLy = rRL*(vRL - SB); rcRRy = rRR*(ST - vRR)
  
  ustar = (rcLLx*uLL + rcLRx*uLR + rcRLx*uRL + rcRRx*uRR &
       & + (PtotLL - PtotRL + PtotLR - PtotRR))/(rcLLx + rcLRx + rcRLx + rcRRx)
  vstar = (rcLLy*vLL + rcLRy*vLR + rcRLy*vRL + rcRRy*vRR &
       & + (PtotLL - PtotLR + PtotRL - PtotRR))/(rcLLy + rcLRy + rcRLy + rcRRy)
              
  SLustar = one/(SL - ustar)
  SRustar = one/(SR - ustar)
  SBvstar = one/(SB - vstar)
  STvstar = one/(ST - vstar)

  rstarLLx = rLL*(SL - uLL)*SLustar
  BstarLL  = BLL*(SL - uLL)*SLustar
  rstarLLy = rLL*(SB - vLL)*SBvstar
  AstarLL  = ALL*(SB - vLL)*SBvstar
  rstarLL  = rLL*(SL - uLL)*SLustar*(SB - vLL)*SBvstar
  EstarLLx = ustar*BstarLL - vLL  *ALL
  EstarLLy = uLL*BLL - vstar*AstarLL
  EstarLL  = ustar*BstarLL - vstar*AstarLL
  
  rstarLRx = rLR*(SL - uLR)*SLustar
  BstarLR  = BLR*(SL - uLR)*SLustar
  rstarLRy = rLR*(ST - vLR)*STvstar
  AstarLR  = ALR*(ST - vLR)*STvstar
  rstarLR  = rLR*(SL - uLR)*SLustar*(ST - vLR)*STvstar
  EstarLRx = ustar*BstarLR - vLR  *ALR
  EstarLRy = uLR*BLR - vstar*AstarLR
  EstarLR  = ustar*BstarLR - vstar*AstarLR
  
  rstarRLx = rRL*(SR - uRL)*SRustar
  BstarRL  = BRL*(SR - uRL)*SRustar
  rstarRLy = rRL*(SB - vRL)*SBvstar
  AstarRL  = ARL*(SB - vRL)*SBvstar
  rstarRL  = rRL*(SR - uRL)*SRustar*(SB - vRL)*SBvstar
  EstarRLx = ustar*BstarRL - vRL  *ARL
  EstarRLy = uRL*BRL - vstar*AstarRL
  EstarRL  = ustar*BstarRL - vstar*AstarRL
  
  rstarRRx = rRR*(SR - uRR)*SRustar
  BstarRR  = BRR*(SR - uRR)*SRustar
  rstarRRy = rRR*(ST - vRR)*STvstar
  AstarRR  = ARR*(ST - vRR)*STvstar
  rstarRR  = rRR*(SR - uRR)*SRustar*(ST - vRR)*STvstar
  EstarRRx = ustar*BstarRR - vRR  *ARR
  EstarRRy = uRR*BRR - vstar*AstarRR
  EstarRR  = ustar*BstarRR - vstar*AstarRR
  
  sqrtrstarLR = one/sqrt(rstarLR)
  sqrtrstarRL = one/sqrt(rstarRL)
  sqrtrstarRR = one/sqrt(rstarRR)
  sqrtrstarLL = one/sqrt(rstarLL)

  calfvenL = max(abs(ALR)/sqrt(rstarLRx),abs(AstarLR)*sqrtrstarLR, &
             abs(ALL)/sqrt(rstarLLx),abs(AstarLL)*sqrtrstarLL, smallc)
  calfvenR = max(abs(ARR)/sqrt(rstarRRx),abs(AstarRR)*sqrtrstarRR, &
             abs(ARL)/sqrt(rstarRLx),abs(AstarRL)*sqrtrstarRL, smallc)
  calfvenB = max(abs(BLL)/sqrt(rstarLLy),abs(BstarLL)*sqrtrstarLL, &
             abs(BRL)/sqrt(rstarRLy),abs(BstarRL)*sqrtrstarRL, smallc)
  calfvenT = max(abs(BLR)/sqrt(rstarLRy),abs(BstarLR)*sqrtrstarLR, &
             abs(BRR)/sqrt(rstarRRy),abs(BstarRR)*sqrtrstarRR, smallc)
  SAL = min(ustar - calfvenL, zero); SAR = max(ustar + calfvenR, zero)
  SAB = min(vstar - calfvenB, zero); SAT = max(vstar + calfvenT, zero)

  SARSAL = one/(SAR - SAL)
  SATSAB = one/(SAT - SAB)
  
  AstarT = (SAR*AstarRR - SAL*AstarLR)*SARSAL
  AstarB = (SAR*AstarRL - SAL*AstarLL)*SARSAL
  BstarR = (SAT*BstarRR - SAB*BstarRL)*SATSAB
  BstarL = (SAT*BstarLR - SAB*BstarLL)*SATSAB
  
  if (SB > zero) then
     if (SL > zero) then
        E = ELL
     else if (SR < zero) then
        E = ERL
     else
        E = (SAR*EstarLLx - SAL*EstarRLx + SAR*SAL*(BRL - BLL))*SARSAL
     endif 
  else if (ST < zero) then
     if (SL > zero) then
        E = ELR
     else if (SR < zero) then
        E = ERR
     else
        E = (SAR*EstarLRx - SAL*EstarRRx + SAR*SAL*(BRR - BLR))*SARSAL
     endif 
  else if (SL > zero) then
     E = (SAT*EstarLLy - SAB*EstarLRy - SAT*SAB*(ALR - ALL))*SATSAB
  else if (SR < zero) then
     E = (SAT*EstarRLy - SAB*EstarRRy - SAT*SAB*(ARR - ARL))*SATSAB
  else
     E = (SAL*SAB*EstarRR - SAL*SAT*EstarRL &
     & - SAR*SAB*EstarLR + SAR*SAT*EstarLL)*SARSAL*SATSAB &
     & - SAT*SAB*SATSAB*(AstarT - AstarB) &
     & + SAR*SAL*SARSAL*(BstarR - BstarL)
  endif

  emf_tmp = E*dt

  return
end subroutine hlld
!$py solver_template hlld
!===============================================================================
!> 2D HLLF Riemann solver
!===============================================================================
!$py solver_template hllf
subroutine hllf

  real(dp) :: SL, SR, SB, ST, SAL, SAR, SAT, SAB , c2, b2, d2
  real(dp) :: cfastLLx, cfastRLx, cfastLRx, cfastRRx, cfastLLy, cfastRLy
  real(dp) :: cfastLRy, cfastRRy

  ! Compute 4 fast magnetosonic velocity relative to x direction
  c2       = gamma*pLL/rLL
  b2       = ALL*ALL + BLL*BLL + CLL*CLL
  d2       = half*(c2 + b2/rLL)
  cfastLLx = sqrt(d2 + sqrt(d2**2 - c2*ALL*ALL/rLL))

  c2       = gamma*pLR/rLR
  b2       = ALR*ALR + BLR*BLR + CLR*CLR
  d2       = half*(c2 + b2/rLR)
  cfastLRx = sqrt(d2 + sqrt(d2**2 - c2*ALR*ALR/rLR))

  c2       = gamma*pRL/rRL
  b2       = ARL*ARL + BRL*BRL + CRL*CRL
  d2       = half*(c2 + b2/rRL)
  cfastRLx = sqrt(d2 + sqrt(d2**2 - c2*ARL*ARL/rRL))

  c2       = gamma*pRR/rRR
  b2       = ARR*ARR + BRR*BRR + CRR*CRR
  d2       = half*(c2 + b2/rRR)
  cfastRRx = sqrt(d2 + sqrt(d2**2 - c2*ARR*ARR/rRR))

  ! Compute 4 fast magnetosonic velocity relative to y direction
  c2       = gamma*pLL/rLL
  b2       = BLL*BLL + ALL*ALL + CLL*CLL
  d2       = half*(c2 + b2/rLL)
  cfastLLy = sqrt(d2 + sqrt(d2**2 - c2*BLL*BLL/rLL))

  c2       = gamma*pLR/rLR
  b2       = BLR*BLR + ALR*ALR + CLR*CLR
  d2       = half*(c2 + b2/rLR)
  cfastLRy = sqrt(d2 + sqrt(d2**2 - c2*BLR*BLR/rLR))

  c2       = gamma*pRL/rRL
  b2       = BRL*BRL + ARL*ARL + CRL*CRL
  d2       = half*(c2 + b2/rRL)
  cfastRLy = sqrt(d2 + sqrt(d2**2 - c2*BRL*BRL/rRL))

  c2       = gamma*pRR/rRR
  b2       = BRR*BRR + ARR*ARR + CRR*CRR
  d2       = half*(c2 + b2/rRR)
  cfastRRy = sqrt(d2 + sqrt(d2**2 - c2*BRR*BRR/rRR))
  
  SL = min(uLL,uLR,uRL,uRR) - max(cfastLLx,cfastLRx,cfastRLx,cfastRRx,smallc)
  SR = max(uLL,uLR,uRL,uRR) + max(cfastLLx,cfastLRx,cfastRLx,cfastRRx,smallc)
  SB = min(vLL,vLR,vRL,vRR) - max(cfastLLy,cfastLRy,cfastRLy,cfastRRy,smallc)
  ST = max(vLL,vLR,vRL,vRR) + max(cfastLLy,cfastLRy,cfastRLy,cfastRRy,smallc)

  emf_tmp = (SL*SB*ERR - SL*ST*ERL - SR*SB*ELR + SR*ST*ELL) &
             &  /(SR - SL)/(ST - SB) &
             & - ST*SB/(ST - SB)*(ARR - ALL) &
             & + SR*SL/(SR - SL)*(BRR - BLL)
  emf_tmp = emf_tmp*dt

  return
end subroutine hllf
!$py solver_template hllf
!===============================================================================
!> 2D HLLA Riemann solver
!===============================================================================
!$py solver_template hlla
subroutine hlla

  real(dp) :: SL, SR, SB, ST, SAL, SAR, SAT, SAB 
  real(dp) :: cfastLLx, cfastRLx, cfastLRx, cfastRRx, cfastLLy, cfastRLy
  real(dp) :: cfastLRy, cfastRRy

  ! Compute 4 Alfven velocity relative to x direction
  cfastLLx = sqrt(ALL*ALL/rLL)
  cfastLRx = sqrt(ALR*ALR/rLR)
  cfastRLx = sqrt(ARL*ARL/rRL)
  cfastRRx = sqrt(ARR*ARR/rRR)

  ! Compute 4 Alfven velocity relative to y direction
  cfastLLy = sqrt(BLL*BLL/rLL)
  cfastLRy = sqrt(BLR*BLR/rLR)
  cfastRLy = sqrt(BRL*BRL/rRL)
  cfastRRy = sqrt(BRR*BRR/rRR)
  
  SL = min(min(uLL,uLR,uRL,uRR) - max(cfastLLx,cfastLRx,cfastRLx,cfastRRx,smallc), zero)
  SR = max(max(uLL,uLR,uRL,uRR) + max(cfastLLx,cfastLRx,cfastRLx,cfastRRx,smallc), zero)
  SB = min(min(vLL,vLR,vRL,vRR) - max(cfastLLy,cfastLRy,cfastRLy,cfastRRy,smallc), zero)
  ST = max(max(vLL,vLR,vRL,vRR) + max(cfastLLy,cfastLRy,cfastRLy,cfastRRy,smallc), zero)

  emf_tmp = (SL*SB*ERR - SL*ST*ERL - SR*SB*ELR + SR*ST*ELL) &
      &  /(SR - SL)/(ST - SB) &
      & - ST*SB/(ST - SB)*(ARR - ALL) &
      & + SR*SL/(SR - SL)*(BRR - BLL)
  emf_tmp = emf_tmp*dt

  return
end subroutine hlla
!$py solver_template hlla
!===============================================================================
!> 2D Lax-Friedrich Riemann solver
!===============================================================================
!$py solver_template llf
subroutine llf

  real(dp) :: E, rl, pl, ul, al, bl, cl, rr, pr, ur, ar, br, cr, ploc, proc
  real(dp) :: c2, b2, d2, cf, vleft, vright, vel_info
  real(dp) :: flux_x, flux_y

  E = forth*(ELL + ERL + ELR + ERR)

  ! Call the solver in the first direction
  ! Left state
  rl = half*(rLL + rLR)
  pl = half*(pLL + pLR)
  ul = half*(uLL + uLR)
  al = half*(ALL + ALR)
  bl = half*(BLL + BLR)
  cl = half*(CLL + CLR)

  ! Right state
  rr = half*(rRR + rRL)
  pr = half*(pRR + pRL)
  ur = half*(uRR + uRL)
  ar = half*(ARR + ARL)
  br = half*(BRR + BRL)
  cr = half*(CRR + CRL)

#if ISO == 1
  ploc = rl*ciso**2
  proc = rr*ciso**2
#else
  ploc = pl
  proc = pr
#endif  
              
  ! enforce continuity of normal component
  al = half*(al + ar)
  ar = al
              
  ! left variables
  c2 = gamma*ploc/rl
  b2 = al*al + bl*bl + cl*cl
  d2 = half*(c2 + b2/rl)
  cf = sqrt(d2 + sqrt(d2**2 - c2*al*al/rl))
  vleft = cf + abs(ul)

  ! right variables
  c2 = gamma*proc/rr
  b2 = ar*ar + br*br + cr*cr
  d2 = half*(c2 + b2/rr)
  cf = sqrt(d2 + sqrt(d2**2 - c2*ar*ar/rr))
  vright = cf + abs(ur)

  ! compute the godunov flux
  vel_info = max(vleft, vright)
  flux_x   = -vel_info*half*(br - bl)

  ! Call the solver in the second direction
  ! Left state
  rl = half*(rLL + rRL)
  pl = half*(pLL + pRL)
  ul = half*(vLL + vRL)
  al = half*(BLL + BRL)
  bl = half*(ALL + ARL)
  cl = half*(CLL + CRL)

  ! Right state
  rr = half*(rRR + rLR)
  pr = half*(pRR + pLR)
  ur = half*(vRR + vLR)
  ar = half*(BRR + BLR)
  br = half*(ARR + ALR)
  cr = half*(CRR + CLR)

#if ISO == 1
  ploc = rl*ciso**2
  proc = rr*ciso**2
#else
  ploc = pl
  proc = pr
#endif  
              
  ! enforce continuity of normal component
  al = half*(al + ar)
  ar = al
              
  ! left variables
  c2 = gamma*ploc/rl
  b2 = al*al + bl*bl + cl*cl
  d2 = half*(c2 + b2/rl)
  cf = sqrt(d2 + sqrt(d2**2 - c2*al*al/rl))
  vleft = cf + abs(ul)

  ! right variables
  c2 = gamma*proc/rr
  b2 = ar*ar + br*br + cr*cr
  d2 = half*(c2 + b2/rr)
  cf = sqrt(d2 + sqrt(d2**2 - c2*ar*ar/rr))
  vright = cf + abs(ur)

  ! compute the godunov flux
  vel_info = max(vleft, vright)
  flux_y   = - vel_info*half*(br - bl)

  emf_tmp = E + (flux_x - flux_y)
  emf_tmp = emf_tmp*dt

  return
end subroutine llf
!$py solver_template llf
!===============================================================================
!> 2D Upwind Riemann solver
!===============================================================================
!$py solver_template upwind
subroutine upwind

  real(dp) :: E, ul, ur, bl, br, vel_info, flux_x, flux_y

  E = forth*(ELL + ERL + ELR + ERR)

  ! Call the solver in the first direction
  ! Left state
  ul = half*(uLL + uLR)
  bl = half*(BLL + BLR)

  ! Right state
  ur = half*(uRR + uRL)
  br = half*(BRR + BRL)

  ! compute the godunov flux
  vel_info = half*(ul + ur)
  flux_x   = -abs(vel_info)*half*(br - bl)

  ! Call the solver in the second direction
  ! Left state
  ul = half*(vLL + vRL)
  bl = half*(ALL + ARL)

  ! Right state
  ur = half*(vRR + vLR)
  br = half*(ARR + ALR)

  ! compute the godunov flux
  vel_info = half*(ul + ur)
  flux_y   = - abs(vel_info)*half*(br - bl)

  emf_tmp = E + (flux_x - flux_y)
  emf_tmp = emf_tmp*dt

  return
end subroutine upwind
!$py solver_template upwind
!!===============================================================================
!!> 2D Athena-Roe Riemann solver
!!===============================================================================
!subroutine athena_roe2d(qc, emf)
!  use params
!  use variables, only: dt
!  implicit none
!
!  real(dp), dimension(4,nvar+1), intent(in) :: qc
!  real(dp), intent(out) :: emf
!
!  real(dp), dimension(nvar) :: qleft, qright, flux_x, flux_y
!  real(dp) :: rLL, rLR, rRL, rRR, pLL, pLR, pRL, pRR, uLL, uLR, uRL, uRR, vLL
!  real(dp) :: vLR, vRL, vRR, wLL, wLR, wRL, wRR, ALL, ALR, ARL, ARR, BLL, BLR
!  real(dp) :: BRL, BRR, CLL, CLR, CRL, CRR, ELL, ELR, ERL, ERR, E, zero_flux
!
!  zero_flux = zero
!
!  rLL = qc(1,ir); pLL = qc(1,ip); uLL = qc(1,iu); vLL = qc(1,iv)
!  rRL = qc(2,ir); pRL = qc(2,ip); uRL = qc(2,iu); vRL = qc(2,iv)
!  rLR = qc(3,ir); pLR = qc(3,ip); uLR = qc(3,iu); vLR = qc(3,iv)
!  rRR = qc(4,ir); pRR = qc(4,ip); uRR = qc(4,iu); vRR = qc(4,iv)
!          
!  wLL = qc(1,iw); ALL = qc(1,iA); BLL = qc(1,iB); CLL = qc(1,iC)
!  wRL = qc(2,iw); ARL = qc(2,iA); BRL = qc(2,iB); CRL = qc(2,iC)
!  wLR = qc(3,iw); ALR = qc(3,iA); BLR = qc(3,iB); CLR = qc(3,iC)
!  wRR = qc(4,iw); ARR = qc(4,iA); BRR = qc(4,iB); CRR = qc(4,iC)
!
!  ELL = qc(1,iC+1)
!  ERL = qc(2,iC+1)
!  ELR = qc(3,iC+1)
!  ERR = qc(4,iC+1)
!
!  E = forth*(ELL + ERL + ELR + ERR)
!
!  ! Call the solver in the first direction
!  ! Left state
!  qleft(1) = half*(rLL + rLR)
!  qleft(2) = half*(pLL + pLR)
!  qleft(3) = half*(uLL + uLR)
!  qleft(4) = half*(ALL + ALR)
!  qleft(5) = half*(vLL + vLR)
!  qleft(6) = half*(BLL + BLR)
!  qleft(7) = half*(wLL + wLR)
!  qleft(8) = half*(CLL + CLR)
!
!  ! Right state
!  qright(1) = half*(rRR + rRL)
!  qright(2) = half*(pRR + pRL)
!  qright(3) = half*(uRR + uRL)
!  qright(4) = half*(ARR + ARL)
!  qright(5) = half*(vRR + vRL)
!  qright(6) = half*(BRR + BRL)
!  qright(7) = half*(wRR + wRL)
!  qright(8) = half*(CRR + CRL)
!
!  call athena_roe(qleft, qright, flux_x, zero_flux)
!
!  ! Call the solver in the second direction
!  ! Left state
!  qleft(1) = half*(rLL + rRL)
!  qleft(2) = half*(pLL + pRL)
!  qleft(3) = half*(vLL + vRL)
!  qleft(4) = half*(BLL + BRL)
!  qleft(5) = half*(uLL + uRL)
!  qleft(6) = half*(ALL + ARL)
!  qleft(7) = half*(wLL + wRL)
!  qleft(8) = half*(CLL + CRL)
!
!  ! Right state
!  qright(1) = half*(rRR + rLR)
!  qright(2) = half*(pRR + pLR)
!  qright(3) = half*(vRR + vLR)
!  qright(4) = half*(BRR + BLR)
!  qright(5) = half*(uRR + uLR)
!  qright(6) = half*(ARR + ALR)
!  qright(7) = half*(wRR + wLR)
!  qright(8) = half*(CRR + CLR)
!
!  call athena_roe(qleft, qright, flux_y, zero_flux)
!
!  emf = E + (flux_x(6) - flux_y(6))
!  emf = emf*dt
!
!  return
!end subroutine athena_roe2d
!!===============================================================================
!!> Updwind Riemann solver
!!===============================================================================
!subroutine upwind(qleft, qright, flux, zero_flux)
!  use params
!  implicit none
!
!  real(dp), dimension(nvar) :: flux
!  real(dp), dimension(nvar) :: qleft, qright
!  real(dp), dimension(nvar) :: fleft, fright, cvarleft, cvarright, cvardiff
!  real(dp) :: zero_flux
!  real(dp) :: rho, p, vn, vt1, vt2, bn, bt1, bt2, vleft
!
!  ! Find the largest eigenvalue in the normal direction to the interface
!  ! and compute MHD flux
!  rho = qleft(1); p   = qleft(2)
!  vn  = qleft(3); vt1 = qleft(5) ; vt2 = qleft(7)
!  bn  = half*(qleft(4)+qright(4)); bt1 = qleft(6); bt2 = qleft(8)
!
!  call comp_mhd_flux(rho, p, vn, vt1, vt2, bn, bt1, bt2, cvarleft, fleft)
!  
!  rho = qright(1); p   = qright(2)
!  vn  = qright(3); vt1 = qright(5); vt2 = qright(7)
!  bn  = bn       ; bt1 = qright(6); bt2 = qright(8)
!
!  call comp_mhd_flux(rho, p, vn, vt1, vt2, bn, bt1, bt2, cvarright, fright)
!
!  ! Mean flux
!  flux     = half*(fright + fleft)*zero_flux
!
!  ! Mean normal velocity
!  vleft    = half*(qleft(3) + qright(3))
!
!  ! Difference between the 2 states 
!  cvardiff = half*(cvarright - cvarleft)
!
! ! The upwind flux
!  flux     = flux - abs(vleft)*cvardiff
!
!  return
!end subroutine upwind
!!===============================================================================
!!> Lax-Friedrich Riemann solver
!!===============================================================================
!subroutine lax_friedrich(qleft, qright, flux, zero_flux)
!  use params
!  implicit none
!
!  real(dp), dimension(nvar), intent(out) :: flux
!  real(dp), dimension(nvar), intent(in)  :: qleft, qright
!  real(dp), dimension(nvar) :: fleft, fright, cvarleft, cvarright, cvardiff
!  real(dp), intent(in) :: zero_flux
!  real(dp) :: rho, p, vn, vt1, vt2, bn, bt1, bt2
!  real(dp) :: vel_left, vel_right
!  real(dp) :: vel_info
!
!  ! Find the largest eigenvalue in the normal direction to the interface
!  ! and compute MHD flux
!  rho = qleft(1); p   = qleft(2)
!  vn  = qleft(3); vt1 = qleft(5) ; vt2 = qleft(7)
!  bn  = half*(qleft(4)+qright(4)); bt1 = qleft(6); bt2 = qleft(8)
!
!  call comp_mhd_flux(rho, p, vn, vt1, vt2, bn, bt1, bt2, cvarleft, fleft)
!  call comp_speed(rho, p, vn, bn, bt1, bt2, vel_left)
!
!  rho = qright(1); p   = qright(2)
!  vn  = qright(3); vt1 = qright(5); vt2 = qright(7)
!  bn  = bn       ; bt1 = qright(6); bt2 = qright(8)
!
!  call comp_mhd_flux(rho, p, vn, vt1, vt2, bn, bt1, bt2, cvarright, fright)
!  call comp_speed(rho, p, vn, bn, bt1, bt2, vel_right)
!  vel_info = max(vel_left, vel_right)
!
!  ! Find the mean flux
!  flux = half*(fright + fleft)*zero_flux
!
!  ! Difference between the 2 states
!  cvardiff = half*(cvarright - cvarleft)
!
!  ! Local Lax-Friedrich flux
!  flux = flux - vel_info*cvardiff
!  
!  return
!end subroutine lax_friedrich
!!===============================================================================
!!> Athena Roe Riemann solver
!!===============================================================================
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
