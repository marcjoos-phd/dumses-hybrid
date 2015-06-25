!===============================================================================
!> \file stratified/user_init.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is user's initialization subroutines for the problem of a stratified
!! shearing box.
!! \details
!! Contains user_init(), get_eta()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Romain Teyssier, 
!! Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          03-04-2015
!! \b last \b modified: 03-04-2015
!<
!===============================================================================
!> User initialization
!===============================================================================
subroutine user_init
  use params
  use variables
  use stratified
  implicit none

  integer :: i,j,k,l,ipos,idim
  integer :: ilo,ihi,jlo,jhi,klo,khi
  real(dp) :: xc,yc,phic,dist
  real(dp), dimension(2,ndim) :: phi

  ! set rhs variables (==.true. when there are source terms)
  rhs=.true.

  !Compute gravin array
  ilo = min(1,iu1+1) ; ihi = max(1,iu2-1)
  jlo = min(1,ju1+1) ; jhi = max(1,ju2-1)
  klo = min(1,ku1+1) ; khi = max(1,ku2-1)

  gravin=0.d0 ; phi=0.d0
  do k=klo,khi
     do j=jlo,jhi
        do i=ilo,ihi

           phi(1,3)=half*Omega0**2*z(k-1)**2
           phi(2,3)=half*Omega0**2*z(k+1)**2
           if (smooth) then
              if (abs(z(k-1))>zfloor) phi(1,3)=half*Omega0**2*zfloor**2
              if (abs(z(k+1))>zfloor) phi(2,3)=half*Omega0**2*zfloor**2
           endif

           !First coordinate gravity force
           gravin(i,j,k,1)=-half*(phi(2,1)-phi(1,1))/dx
           !Second coordinate gravity force
           if (ndim>1) gravin(i,j,k,2)=-half*(phi(2,2)-phi(1,2))/dy
           !Third coordinate gravity force
           if (ndim>2) gravin(i,j,k,3)=-half*(phi(2,3)-phi(1,3))/dz

        end do
     end do
  end do

  return
end subroutine user_init
!===============================================================================
!> Get eta parameter
!===============================================================================
subroutine get_eta(etaval,r)
  use precision
  implicit none
  real(dp) :: etaval,r
  return
end subroutine get_eta
