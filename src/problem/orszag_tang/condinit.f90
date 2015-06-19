!===============================================================================
!> \file orszag_tang/condinit.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is the initial conditions subroutine for Orszag Tang problem.
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

  real(dp) :: d0, beta, p0, B0
  real(dp) :: Ek, Em
  integer  :: i, j, k
  character(len=3) :: direction

  namelist /init_params/ direction
  
  open(unit=1, file='input', status='old')
  read(1, init_params)
  close(1)

  if(verbose) print*, 'Dimension (ndim) of the problem: ', ndim

  ! Orszag-Tang problem
  gamma = 5.d0/3.d0
  B0    = (2.d0*twopi)**(-half)
  beta  = 2.d0*gamma
  p0    = gamma/(2.d0*twopi)
  
  if(direction .eq. 'x') then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = ku1, ku2
        do j = ju1, ju2
           do i = iu1, iu2
              uin(i,j,k,ir) = gamma*p0
              uin(i,j,k,iu) = -1.d0*uin(i,j,k,ir)*sin(y(j)*twopi)
              uin(i,j,k,iv) = uin(i,j,k,ir)*sin(x(i)*twopi)
              
              uin(i,j,k,iA)   = -B0*sin(y(j)*twopi)
              uin(i,j,k,iA+3) = uin(i,j,k,iA)
              
              uin(i,j,k,iB)   = B0*sin(2.d0*x(i)*twopi)
              uin(i,j,k,iB+3) = uin(i,j,k,iB)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO

  else if (direction .eq. 'y') then
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = ku1, ku2
        do j = ju1, ju2
           do i = iu1, iu2
              uin(i,j,k,ir) = gamma*p0
              uin(i,j,k,iv) = -1.d0*uin(i,j,k,ir)*sin(z(k)*twopi)
              uin(i,j,k,iw)  = uin(i,j,k,ir)*sin(y(j)*twopi)
              
              uin(i,j,k,iB)   = -B0*sin(z(k)*twopi)
              uin(i,j,k,iB+3) = uin(i,j,k,iB)
              
              uin(i,j,k,iC)   = B0*sin(2.d0*y(j)*twopi)
              uin(i,j,k,iC+3) = uin(i,j,k,iC)
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
              uin(i,j,k,ir) = gamma*p0
              uin(i,j,k,iw)  = -1.d0*uin(i,j,k,ir)*sin(x(i)*twopi)
              uin(i,j,k,iu) = uin(i,j,k,ir)*sin(z(k)*twopi)
              
              uin(i,j,k,iC)   = -B0*sin(x(i)*twopi)
              uin(i,j,k,iC+3) = uin(i,j,k,iC)
              
              uin(i,j,k,iA)   = B0*sin(2.d0*z(k)*twopi)
              uin(i,j,k,iA+3) = uin(i,j,k,iA)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif
  
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(Ek, Em)
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2
           Ek = uin(i,j,k,iu)**2/uin(i,j,k,ir) + uin(i,j,k,iv)**2/uin(i,j,k,ir)&
              + uin(i,j,k,iw)**2/uin(i,j,k,ir)
           if((i.lt.iu2).and.(j.lt.ju2).and.(k.lt.ku2)) then
              Em = forth*((uin(i,j,k,iA) + uin(i+1,j,k,iA))**2 &
                        + (uin(i,j,k,iB) + uin(i,j+1,k,iB))**2 &
                        + (uin(i,j,k,iC) + uin(i,j,k+1,iC))**2)
           else if((i.lt.iu2).and.(j.eq.ju2).and.(k.lt.ku2)) then
              Em = forth*((uin(i,j,k,iA) + uin(i+1,j,k,iA))**2 &
                        + (uin(i,j,k,iB) + uin(i,4,k,iB))**2 &
                        + (uin(i,j,k,iC) + uin(i,j,k+1,iC))**2)
           else if((i.eq.iu2).and.(j.lt.ju2).and.(k.lt.ku2)) then
              Em = forth*((uin(i,j,k,iA) + uin(4,j,k,iA))**2 &
                        + (uin(i,j,k,iB) + uin(i,j+1,k,iB))**2 &
                        + (uin(i,j,k,iC) + uin(i,j,k+1,iC))**2)
           else if((i.eq.iu2).and.(j.eq.ju2).and.(k.lt.ku2)) then
              Em = forth*((uin(i,j,k,iA) + uin(4,j,k,iA))**2 &
                        + (uin(i,j,k,iB) + uin(i,4,k,iB))**2 &
                        + (uin(i,j,k,iC) + uin(i,j,k+1,iC))**2)
           else if((i.lt.iu2).and.(j.lt.ju2).and.(k.eq.ku2)) then
              Em = forth*((uin(i,j,k,iA) + uin(i+1,j,k,iA))**2 &
                        + (uin(i,j,k,iB) + uin(i,j+1,k,iB))**2 &
                        + (uin(i,j,k,iC) + uin(i,j,4,iC))**2)
           else if((i.lt.iu2).and.(j.eq.ju2).and.(k.eq.ku2)) then
              Em = forth*((uin(i,j,k,iA) + uin(i+1,j,k,iA))**2 &
                        + (uin(i,j,k,iB) + uin(i,4,k,iB))**2 &
                        + (uin(i,j,k,iC) + uin(i,j,4,iC))**2)
           else if((i.eq.iu2).and.(j.lt.ju2).and.(k.eq.ku2)) then
              Em = forth*((uin(i,j,k,iA) + uin(4,j,k,iA))**2 &
                        + (uin(i,j,k,iB) + uin(i,j+1,k,iB))**2 &
                        + (uin(i,j,k,iC) + uin(i,j,4,iC))**2) 
           else if((i.eq.iu2).and.(j.eq.ju2).and.(k.eq.ku2)) then
              Em = forth*((uin(i,j,k,iA) + uin(4,j,k,iA))**2 &
                        + (uin(i,j,k,iB) + uin(i,4,k,iB))**2 &
                        + (uin(i,j,k,iC) + uin(i,j,4,iC))**2)
           endif
           uin(i,j,k,ip) = p0/(gamma - 1.d0) + half*(Ek + Em)
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

  return
end subroutine condinit
