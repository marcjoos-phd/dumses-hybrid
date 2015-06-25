!===============================================================================
!> \file stratified/condinit.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is the initial conditions subroutine for MRI problem.
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
!! \b last \b modified: 05-27-2015
!<
!===============================================================================
!> setup the initial conditions for a stratified shearing box
!===============================================================================
subroutine condinit
  use params
  use variables
  use mpi_var
#if MPI == 1
  use mpi
#endif
  use stratified
  implicit none

  integer  :: i,j,k,l,iseed
  real(dp) :: H,beta,B0,p0,amp,rhofloor
  real(dp) :: rho,rhovx,rhovy,rhovz,Bxc,Byc,Bzc,E
  real :: rvalue
  !
  ! Eigenmode properties
  !
  ! 'l': Legendre poly, 'r': random, 'h': henrik
  ! 'z': zero net flux random (with Bz)
  ! 'y': zero net flux random (with By)
  character(len=10) :: type
  character(len=10) :: nchar,Rm,Re !mode number (character)
  character(len=256) :: filename,dirname
  integer :: n,nz_eigen
  real(dp) :: kn,alpha,va,xi,s,u0,Bnorm
  real, dimension(:), allocatable :: F,dF,G
  real(dp), dimension(:), allocatable :: z_eigen,F_eigen,G_eigen,ux_eigen,uy_eigen,bx_eigen,by_eigen
  real(dp) :: Fvalue,Gvalue,uxvalue,uyvalue,bxvalue,byvalue
  !
#if MPI==1
  real(dp) :: Bnorm_tmp
  integer :: ierr
#endif
  namelist /init_params/d0,beta,amp,n,type,dirname,floor,zfloor,dfloor,smooth,addmass,massDiff

  if(verbose) print*, 'Dimension (ndim) of the problem: ', ndim

  ! MRI initial conditions
  if (verbose) print*, "> Define initial conditions"
  d0=1.d0
  beta=100.d0
  amp=1.d-1
  n=10
  type='l'
  dirname=''
  floor=.false. ; smooth=.false. ; addmass=.false. ; massDiff=.false.
  zfloor=5.d0

  open(unit=1,file='input' ,status='old')
  read(1,init_params)
  close(1)

  H=ciso/Omega0
#if ISO==0
  H=H/sqrt(gamma)
#endif
  B0=sqrt(2.d0*d0*ciso**2/beta)

  !Read eigenmodes from data files
  if (type=='h') then
     write (nchar,'(I2)') n
     filename = trim(dirname) // 'FG' // trim(adjustl(nchar)) //'.txt'
     open(unit=1,file=filename,status='old')
     read(1,*) kn
     read(1,*) nz_eigen
     allocate(z_eigen(nz_eigen),F_eigen(nz_eigen),G_eigen(nz_eigen))
     do k=1,nz_eigen
        read(1,*) z_eigen(k),F_eigen(k),G_eigen(k)
     end do
     close(1)
  endif
  !Read eigenmodes from data files (case with nonzero resistivity)
  if (type=='resist') then
     write (nchar,'(I5)') int(beta)
     write (Rm   ,'(I10)') int(ciso*H/eta)
     if (nu==0) then
        filename = trim(dirname) // 'beta' // trim(adjustl(nchar)) // &
             & 'Rm' // trim(adjustl(Rm)) //'.txt'
     else
        write (Re   ,'(I10)') int(ciso*H/nu)
        filename = trim(dirname) // 'beta' // trim(adjustl(nchar)) // &
             & 'Rm' // trim(adjustl(Rm)) // &
             & 'Re' // trim(adjustl(Re)) //'.txt'
     endif
     open(unit=1,file=filename,status='old')
     read(1,*) nz_eigen
     allocate(z_eigen(nz_eigen))
     allocate(ux_eigen(nz_eigen),uy_eigen(nz_eigen))
     allocate(bx_eigen(nz_eigen),by_eigen(nz_eigen))
     do k=1,nz_eigen
        read(1,*) z_eigen(k),ux_eigen(k),uy_eigen(k),bx_eigen(k),by_eigen(k)
     end do
     close(1)
  endif


  if (type=='l') kn=1.d0/H*sqrt(real(n)+real(n)**2)
  va=B0/sqrt(d0)
  xi=va/H/Omega0
  s=Omega0*sqrt((-(1.d0+2.d0*xi**2*kn**2*H**2)+&
       & sqrt(1.d0+16.d0*xi**2*kn**2*H**2))/2.d0)
  alpha=asin(sqrt(1.d0/6.d0*(sqrt(1.d0+32.d0*kn**2*H**2/beta)-1.d0)))
  u0=s/kn*tan(alpha)

  if (type=='l') then
     allocate(F(0:n)) ; allocate(dF(0:n)) ; allocate(G(0:n))
  endif

  if (type=='z') B0=2.0*sqrt(d0*ciso**2/beta)
  if (type=='hydro') B0=0.d0

  uin=0.d0
  iseed=mype

  !Set vertical profile
  do k=ku1,ku2

     !density
     if (floor) then
#if ISO==1
        uin(:,:,k,1)=d0*max(exp(-z(k)**2/2.D0/H**2),exp(-zfloor**2/2.d0/H**2))
#else
        uin(:,:,k,1)=max(d0*exp(-z(k)**2/2.D0/H**2),dfloor)
#endif
     else
        uin(:,:,k,1)=d0*exp(-z(k)**2/2.D0/H**2)
     endif

  end do

  do k=1,nz+1

     if ((type=='r').or.(type=='z').or.(type=='y')) then
        !Momentum (Random condition)
        do j=1,ny
           do i=1,nx
              call ran2(iseed,rvalue)
              uin(i,j,k,2) = uin(i,j,k,1)*amp*(rvalue-0.5d0)*sqrt(ciso**2)
              call ran2(iseed,rvalue)
              uin(i,j,k,4) = uin(i,j,k,1)*amp*(rvalue-0.5d0)*sqrt(ciso**2)
              call ran2(iseed,rvalue)
              uin(i,j,k,3) = uin(i,j,k,1)*amp*(rvalue-0.5d0)*sqrt(ciso**2)
           end do
        end do
        !Magnetic field (pure vertical or zero net flux)
        if (type=='r') uin(:,:,k,8:nvar+3:3)= B0
        if (type=='z') then
           do j=ju1,ju2
              uin(iu1:iu2,j,k,8     ) = B0*sin(2.d0*pi*x(iu1:iu2))
              uin(iu1:iu2,j,k,nvar+3) = B0*sin(2.d0*pi*x(iu1:iu2))
           end do
        endif
        if (type=='y') then
           if (abs(z(k))<2.d0) uin(:,:,k,7:nvar+3:3)=B0
        endif
     endif

     if (type=='l') then
        !Momentum (Legendre polynomial)
        F(0:n)=0.d0 ; dF(0:n)=0.d0
        call legendre(n,real(tanh(z(k)/H)),F(0:n),dF(0:n))
        G=-dF/kn*(1.-tanh(z(k)/H)**2)
        uin(:,:,k,2) = uin(:,:,k,1)*u0*cos(alpha)*F(n)
        uin(:,:,k,3) = uin(:,:,k,1)*u0*sin(alpha)*F(n)
        uin(:,:,k,4) = 0.d0
        !Magnetic field
        uin(:,:,k,6)=-B0*sin(alpha)*G(n)
        uin(:,:,k,7)= B0*cos(alpha)*G(n)
        uin(:,:,k,8)= B0
     endif

     if (type=='h') then
        call GetInterpolated(z(k),Fvalue,z_eigen,F_eigen,nz_eigen)
        call GetInterpolated(z(k),Gvalue,z_eigen,G_eigen,nz_eigen)
        !Momentum (Legendre polynomial)
        uin(:,:,k,2) = uin(:,:,k,1)*u0*cos(alpha)*Fvalue
        uin(:,:,k,3) = uin(:,:,k,1)*u0*sin(alpha)*Fvalue
        uin(:,:,k,4) = 0.d0
        !Magnetic field
        uin(:,:,k,6)=+B0*sin(alpha)*Gvalue
        uin(:,:,k,7)=-B0*cos(alpha)*Gvalue
        uin(:,:,k,8)= B0
     endif

     if (type=='resist') then
        call GetInterpolated(z(k),uxvalue,z_eigen,ux_eigen,nz_eigen)
        call GetInterpolated(z(k),uyvalue,z_eigen,uy_eigen,nz_eigen)
        call GetInterpolated(z(k),bxvalue,z_eigen,bx_eigen,nz_eigen)
        call GetInterpolated(z(k),byvalue,z_eigen,by_eigen,nz_eigen)
        !Momentum (Legendre polynomial)
        uin(:,:,k,2) = uin(:,:,k,1)*uxvalue
        uin(:,:,k,3) = uin(:,:,k,1)*uyvalue
        uin(:,:,k,4) = 0.d0
        !Magnetic field
        uin(:,:,k,6)= bxvalue
        uin(:,:,k,7)= byvalue
        uin(:,:,k,8)= B0
     endif

  enddo

  !Normalise eigenvectors
  if ((type=='l').or.(type=='h').or.(type=='resist')) then
     Bnorm=maxval(sqrt(uin(1,1,1:nz,6)**2+uin(1,1,1:nz,7)**2))
#if MPI==1
     CALL MPI_Allreduce(Bnorm,Bnorm_tmp,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
          &                MPI_COMM_WORLD,ierr)
     Bnorm=Bnorm_tmp
#endif
     amp=amp*B0/Bnorm
     if (type=='l') amp=(-1)**(n+1)*amp
     !Momentum
     uin(:,:,:,2)=amp*uin(:,:,:,2)
     uin(:,:,:,3)=amp*uin(:,:,:,3)
     !Magnetic field
     uin(:,:,:,6)=amp*uin(:,:,:,6)
     uin(:,:,:,7)=amp*uin(:,:,:,7)
  endif

  if (type=='l') deallocate(F,dF,G)
  if (type=='h') deallocate(z_eigen,F_eigen,G_eigen)

  !Set total energy density
#if ISO==0
  do k=ku1,ku2-1
     do j=ju1,ju2-1
        do i=iu1,iu2-1

           rho=uin(i,j,k,1)
           rhovx=uin(i,j,k,2) ; Bxc=half*(uin(i,j,k,6)+uin(i+1,j  ,k  ,6))
           rhovy=uin(i,j,k,3) ; Byc=half*(uin(i,j,k,7)+uin(i  ,j+1,k  ,7))
           rhovz=uin(i,j,k,4) ; Bzc=half*(uin(i,j,k,8)+uin(i  ,j  ,k+1,8))

           call get_E(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,ciso*ciso)
           uin(i,j,k,5)=E

        end do
     end do
  end do
#endif


  !Calculate initial mass in computational box
  mass0=0.d0
  do k=1,nz
     do j=1,ny
        do i=1,nx
           mass0 = mass0 + uin(i,j,k,1)
        end do
     end do
  end do
#if MPI==1
  call sumBcast(mass0,1)
#endif
  mass0=mass0*dx*dy*dz

  return
end subroutine condinit
!===============================================================================
!> This routine computes and broadcast the sum of an array 
!===============================================================================
subroutine sumBcast(quant,nvar)
#if MPI == 1
  use precision
  use mpi
  implicit none

  integer :: nvar,ierr
  real(dp), dimension(nvar) :: quant
  real(dp), dimension(nvar) :: vol_quant

  call MPI_Allreduce(quant,vol_quant,nvar,MPI_DOUBLE_PRECISION,MPI_SUM, &
   &                 MPI_COMM_WORLD, ierr)
  quant=vol_quant

#endif
  return
end subroutine sumBcast
!====================================================================
!> Numerical recipes random number generator ran2
!>    requires input seed value=iseed
!>    returns real random number=rvalue
!>    Also updates iseed for next call 
!====================================================================
      subroutine ran2(iseed,rvalue)
     
      integer iseed
      real rvalue
      integer idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real AM,EPS,RNMX
      parameter (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,&
               & IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,&
               & IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      integer idum2,jj,kk,iv(NTAB),iy
      data idum2/123456789/, iv/NTAB*0/, iy/0/
!
      idum=iseed
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 jj=NTAB+8,1,-1
          kk=idum/IQ1
          idum=IA1*(idum-kk*IQ1)-kk*IR1
          if (idum.lt.0) idum=idum+IM1
          if (jj.le.NTAB) iv(jj)=idum
11      continue
        iy=iv(1)
      endif
      kk=idum/IQ1
      idum=IA1*(idum-kk*IQ1)-kk*IR1
      if (idum.lt.0) idum=idum+IM1
      kk=idum2/IQ2
      idum2=IA2*(idum2-kk*IQ2)-kk*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      jj=1+iy/NDIV
      iy=iv(jj)-idum2
      iv(jj)=idum
      if(iy.lt.1)iy=iy+IMM1
      rvalue=min(AM*iy,RNMX)
      iseed=idum
      return
   end if

 end subroutine ran2
!====================================================================
!> Subroutine to generate Legendre polynomials P_L(X)
!> for L = 0,1,...,LMAX with given X.
!====================================================================
subroutine legendre(LMAX, X, P_L, dP_L)
  use precision
  implicit none

  integer, intent (IN) :: LMAX
  integer :: L
  real(dp), intent (IN) :: X
  real(dp), intent (OUT), dimension (0:LMAX) :: P_L,dP_L
!
  P_L(0) = 1.0d0 ; dP_L(0) = 0.0d0
  P_L(1) = X     ; dP_L(1) = 1.0d0
  do L = 1, LMAX-1
     P_L(L+1) = ((2.0*L+1)*X* P_L(L)         -L* P_L(L-1))/(L+1)
    dP_L(L+1) = ((2.0*L+1)  *(P_L(L)+X*dP_L(L))-L*dP_L(L-1))/(L+1)
  end do
end subroutine legendre
!====================================================================
!> Subroutine to interpolate linearly a 1D function at location zloc.
!====================================================================
subroutine GetInterpolated(zloc,floc,z,func,n)
  use precision
  implicit none
  integer :: n
  real(dp) :: zloc,floc
  real(dp), dimension(n) :: z,func

  integer :: k
  real(dp) :: z0,f0,df,dz

  k=1
  do while (z(k)<zloc)
     k=k+1
  end do
  f0=func(k) ; z0=z(k)
  df=func(k+1)-func(k)
  dz=z(k+1)-z(k)
  floc=f0+(zloc-z0)*df/dz

  return
end subroutine GetInterpolated
