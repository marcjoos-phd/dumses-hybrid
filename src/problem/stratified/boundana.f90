!===============================================================================
!> \file boundana.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This subroutine contains the analytical boundary conditions.
!! \details
!! Contains xinner_ana(), xouter_ana(), yinner_ana(), youter_ana(),
!! zinner_ana(), zouter_ana()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Romain Teyssier, 
!! Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          05-07-2015 
!! \b last modified:    05-18-2015 
!! \b last \b modified: 
!<
!===============================================================================
!> Compute boundary conditions in the x-direction, inner edge
!===============================================================================
subroutine xinner_ana
  use params
  use variables
  implicit none

  integer :: j, k

  ! Periodic BC
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do k = ku1, ku2
     do j = ju1, ju2      
        uin(iu1  ,j,k,:) = uin(nx-2,j,k,:)
        uin(iu1+1,j,k,:) = uin(nx-1,j,k,:)
        uin(iu1+2,j,k,:) = uin(nx  ,j,k,:)
     end do
  end do
  !$OMP END PARALLEL DO

  return
end subroutine xinner_ana
!===============================================================================
!> Compute boundary conditions in the x-direction, outer edge
!===============================================================================
subroutine xouter_ana
  use params
  use variables
  implicit none

  integer :: j, k

  ! Periodic BC
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do k = ku1, ku2
     do j = ju1, ju2
        uin(iu2-2,j,k,:) = uin(1 ,j,k,:)
        uin(iu2-1,j,k,:) = uin(2 ,j,k,:)
        uin(iu2  ,j,k,:) = uin(3 ,j,k,:)
     end do
  end do
  !$OMP END PARALLEL DO

  return
end subroutine xouter_ana
!===============================================================================
!> Compute boundary conditions in the y-direction, inner edge
!===============================================================================
subroutine yinner_ana
#if NDIM > 1
  use params
  use variables
  implicit none

  integer :: k, i

  ! Periodic BC
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do k = ku1, ku2
     do i = iu1, iu2
        uin(i,ju1  ,k,:) = uin(i,ny-2,k,:)
        uin(i,ju1+1,k,:) = uin(i,ny-1,k,:)
        uin(i,ju1+2,k,:) = uin(i,ny  ,k,:)
     end do
  end do
  !$OMP END PARALLEL DO

#endif
  return
end subroutine yinner_ana
!===============================================================================
!> Compute boundary conditions in the y-direction, outer edge
!===============================================================================
subroutine youter_ana
#if NDIM > 1
  use params
  use variables
  implicit none

  integer :: k, i

  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do k = ku1, ku2
     do i = iu1, iu2
        uin(i,ju2-2,k,:) = uin(i,1 ,k,:)
        uin(i,ju2-1,k,:) = uin(i,2 ,k,:)
        uin(i,ju2  ,k,:) = uin(i,3 ,k,:)
     end do
  end do
  !$OMP END PARALLEL DO

#endif
  return
end subroutine youter_ana
!===============================================================================
!> Compute boundary conditions in the z-direction, inner edge
!===============================================================================
subroutine zinner_ana
#if NDIM > 1
  use params
  use variables
  use stratified
  implicit none

  integer  :: i, j, k
  real(dp) :: H, factor, ratio_nyp1, ratio_nyp2, ratio_nyp3, dbxdx, dbydy
  real(dp) :: rho, rhovx, rhovy, rhovz, E, Bxc, Byc, Bzc, csq

  H      = ciso/Omega0
  factor = -dz/two/H/H
  if (floor) then
     ratio_nyp1 = one
     ratio_nyp2 = one
     ratio_nyp3 = one
  else
     ratio_nyp1 = exp(factor*(-two*(zmin + half*dz) +      dz))
     ratio_nyp2 = exp(factor*(-two*(zmin + half*dz) + 3.d0*dz))
     ratio_nyp3 = exp(factor*(-two*(zmin + half*dz) + 5.d0*dz))
  endif
  
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(rho, E, rhovx, rhovy, rhovz) &
  !$OMP PRIVATE(Bxc, Byc, Bzc) &
  !$OMP FIRSTPRIVATE(H, factor, ratio_nyp1, ratio_nyp2, ratio_nyp3)
  do j = ju1, ju2-1
     do i = iu1, iu2-1
#if ISO == 0
        rho = uin(i,j,ku1+3,1)
        E = uin(i,j,ku1+3,5)
        rhovx = uin(i,j,ku1+3,2)
        rhovy = uin(i,j,ku1+3,3)
        rhovz = uin(i,j,ku1+3,4)
        Bxc = half*(uin(i,j,ku1+3,6) + uin(i+1,j  ,ku1+3,6))
        Byc = half*(uin(i,j,ku1+3,7) + uin(i  ,j+1,ku1+3,7))
        Bzc = half*(uin(i,j,ku1+3,8) + uin(i  ,j  ,ku1+4,8))
        !write (*,*) i,j,Bxc,Byc,Bzc
        call get_csq(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,csq)
        H = sqrt(csq)/Omega0/sqrt(gamma)
        factor = -dz/two/H/H
        ratio_nyp1 = one !exp(factor*(-2*(zmin+half*dz)+     dz))
        ratio_nyp2 = one !exp(factor*(-2*(zmin+half*dz)+3.d0*dz))
        ratio_nyp3 = one !exp(factor*(-2*(zmin+half*dz)+5.d0*dz))
#endif
        !extrapolating density
        !factor=(log(uin(i,j,ku1+4,1))-log(uin(i,j,ku1+3,1)))/log(uin(i,j,ku1+3,1))
        !uin(i,j,ku1+2,1) = uin(i,j,ku1+3,1)**(1.d0-factor)
        !factor=(log(uin(i,j,ku1+3,1))-log(uin(i,j,ku1+2,1)))/log(uin(i,j,ku1+2,1))
        !uin(i,j,ku1+1,1) = uin(i,j,ku1+2,1)**(1.d0-factor)
        !factor=(log(uin(i,j,ku1+2,1))-log(uin(i,j,ku1+1,1)))/log(uin(i,j,ku1+1,1))
        !uin(i,j,ku1  ,1) = uin(i,j,ku1+1,1)**(1.d0-factor)
        uin(i,j,ku1+2,1) = max(uin(i,j,ku1+3,1)*ratio_nyp1,smallr)
        uin(i,j,ku1+1,1) = max(uin(i,j,ku1+2,1)*ratio_nyp2,smallr)
        uin(i,j,ku1  ,1) = max(uin(i,j,ku1+1,1)*ratio_nyp3,smallr)
     end do
  end do
  !$OMP END PARALLEL DO
  
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do j = ju1, ju2
     do i = iu1, iu2
        !zero gradient BC for velocities
        uin(i,j,ku1  ,2:3) = uin(i,j,ku1+3,2:3)/uin(i,j,ku1+3,1)
        uin(i,j,ku1+1,2:3) = uin(i,j,ku1+3,2:3)/uin(i,j,ku1+3,1)
        uin(i,j,ku1+2,2:3) = uin(i,j,ku1+3,2:3)/uin(i,j,ku1+3,1)
        uin(i,j,ku1  ,2:3) = uin(i,j,ku1  ,2:3)*uin(i,j,ku1  ,1)
        uin(i,j,ku1+1,2:3) = uin(i,j,ku1+1,2:3)*uin(i,j,ku1+1,1)
        uin(i,j,ku1+2,2:3) = uin(i,j,ku1+2,2:3)*uin(i,j,ku1+2,1)
        uin(i,j,ku1  ,4) = min(uin(i,j,ku1+3,4), zero)
        uin(i,j,ku1+1,4) = min(uin(i,j,ku1+3,4), zero)
        uin(i,j,ku1+2,4) = min(uin(i,j,ku1+3,4), zero)
        !Tangential magnetic field
        uin(i,j,ku1  ,6:7) = zero !Vanishing tangential field
        uin(i,j,ku1+1,6:7) = zero !Vanishing tangential field
        uin(i,j,ku1+2,6:7) = zero !Vanishing tangential field
        !Tangential magnetic field
        !uin(i,j,ku1  ,6:7) = uin(i,j,ku1+5,6:7)
        !uin(i,j,ku1+1,6:7) = uin(i,j,ku1+4,6:7)
        !uin(i,j,ku1+2,6:7) = uin(i,j,ku1+3,6:7)
!        uin(i,j,ku1  ,6:7) = uin(i,j,ku1+3,6:7)*(uin(i,j,ku1  ,1)/uin(i,j,ku1+3,1))**half
!        uin(i,j,ku1+1,6:7) = uin(i,j,ku1+3,6:7)*(uin(i,j,ku1+1,1)/uin(i,j,ku1+3,1))**half
!        uin(i,j,ku1+2,6:7) = uin(i,j,ku1+3,6:7)*(uin(i,j,ku1+2,1)/uin(i,j,ku1+3,1))**half
     end do
  end do
  !$OMP END PARALLEL DO

  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dbxdx, dbydy)
  do j = ju1, ju2-1
     do i = iu1, iu2-1
        !Normal magnetic field
        dbxdx = (uin(i+1,j  ,ku1+2,6) - uin(i,j,ku1+2,6))/dx
        dbydy = (uin(i  ,j+1,ku1+2,7) - uin(i,j,ku1+2,7))/dy
        uin(i,j,ku1+2,8) = uin(i,j,ku1+3,8) + dz*(dbxdx + dbydy)
        dbxdx = (uin(i+1,j  ,ku1+1,6) - uin(i,j,ku1+1,6))/dx
        dbydy = (uin(i  ,j+1,ku1+1,7) - uin(i,j,ku1+1,7))/dy
        uin(i,j,ku1+1,8) = uin(i,j,ku1+2,8) + dz*(dbxdx + dbydy)
        dbxdx = (uin(i+1,j  ,ku1  ,6) - uin(i,j,ku1  ,6))/dx
        dbydy = (uin(i  ,j+1,ku1  ,7) - uin(i,j,ku1  ,7))/dy
        uin(i,j,ku1  ,8) = uin(i,j,ku1+1,8) + dz*(dbxdx + dbydy)
     end do
  end do

#if ISO == 0
  ! Zero gradient Temperature
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(rho, E, rhovx, rhovy, rhovz) &
  !$OMP PRIVATE(Bxc, Byc, Bzc)
  do j = ju1, ju2-1
     !$acc loop
     do i = iu1, iu2-1
        rho = uin(i,j,ku1+3,1) 
        E = uin(i,j,ku1+3,5)
        rhovx = uin(i,j,ku1+3,2) 
        rhovy = uin(i,j,ku1+3,3) 
        rhovz = uin(i,j,ku1+3,4)
        Bxc = half*(uin(i,j,ku1+3,6) + uin(i+1,j  ,ku1+3,6))
        Byc = half*(uin(i,j,ku1+3,7) + uin(i  ,j+1,ku1+3,7))
        Bzc = half*(uin(i,j,ku1+3,8) + uin(i  ,j  ,ku1+4,8))
        call get_csq(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,csq)
        do k = ku1+2, ku1, -1
           rho = uin(i,j,k,1)
           rhovx = uin(i,j,k,2) 
           rhovy = uin(i,j,k,3)
           rhovz = uin(i,j,k,4)
           Bxc = half*(uin(i,j,k,6)+uin(i+1,j  ,k  ,6))
           Byc = half*(uin(i,j,k,7)+uin(i  ,j+1,k  ,7))
           Bzc = half*(uin(i,j,k,8)+uin(i  ,j  ,k+1,8))
           call get_E(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,csq)
           uin(i,j,k,5) = E
        end do
     end do
  end do
  !$OMP END PARALLEL DO
#endif

#endif
  return
end subroutine zinner_ana
!===============================================================================
!> Compute boundary conditions in the z-direction, inner edge
!===============================================================================
subroutine zouter_ana
#if NDIM > 1
  use params
  use const
  use variables
  use stratified
  implicit none

  integer :: i, j, k
  real(dp):: H, factor, ratio_nyp1, ratio_nyp2, ratio_nyp3, dbxdx, dbydy
  real(dp):: rho, rhovx, rhovy, rhovz, E, Bxc, Byc, Bzc, csq

  H = ciso/Omega0
  factor = -dz/two/H/H
  if (floor) then
     ratio_nyp1 = one
     ratio_nyp2 = one
     ratio_nyp3 = one
  else
     ratio_nyp1 = exp(factor*(two*(zmax-half*dz) +      dz))
     ratio_nyp2 = exp(factor*(two*(zmax-half*dz) + 3.d0*dz))
     ratio_nyp3 = exp(factor*(two*(zmax-half*dz) + 5.d0*dz))
  endif

  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(rho, E, rhovx, rhovy, rhovz) &
  !$OMP PRIVATE(Bxc, Byc, Bzc) &
  !$OMP FIRSTPRIVATE(H, factor, ratio_nyp1, ratio_nyp2, ratio_nyp3)
  do j = ju1, ju2-1
     do i = iu1, iu2-1
#if ISO == 0
        rho = uin(i,j,ku2-3,1) 
        E = uin(i,j,ku2-3,5)
        rhovx = uin(i,j,ku2-3,2) 
        rhovy = uin(i,j,ku2-3,3) 
        rhovz = uin(i,j,ku2-3,4)
        Bxc = half*(uin(i,j,ku2-3,6) + uin(i+1,j  ,ku2-3,6))
        Byc = half*(uin(i,j,ku2-3,7) + uin(i  ,j+1,ku2-3,7))
        Bzc = half*(uin(i,j,ku2-3,8) + uin(i  ,j  ,ku2-2,8))
        call get_csq(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,csq)
        H = sqrt(csq)/Omega0/sqrt(gamma)
        factor = -dz/two/H/H
        ratio_nyp1 = one !exp(factor*(2*(zmax-half*dz)+     dz))
        ratio_nyp2 = one !exp(factor*(2*(zmax-half*dz)+3.d0*dz))
        ratio_nyp3 = one !exp(factor*(2*(zmax-half*dz)+5.d0*dz))
#endif
        !extrapolating density
        !factor=(log(uin(i,j,ku2-3,1))-log(uin(i,j,ku2-4,1)))/log(uin(i,j,ku2-3,1))
        !uin(i,j,ku2-2,1) = uin(i,j,ku2-3,1)**(1.d0+factor)
        !factor=(log(uin(i,j,ku2-2,1))-log(uin(i,j,ku2-3,1)))/log(uin(i,j,ku2-2,1))
        !uin(i,j,ku2-1,1) = uin(i,j,ku2-2,1)**(1.d0+factor)
        !factor=(log(uin(i,j,ku2-1,1))-log(uin(i,j,ku2-2,1)))/log(uin(i,j,ku2-1,1))
        !uin(i,j,ku2  ,1) = uin(i,j,ku2-1,1)**(1.d0+factor)
        uin(i,j,ku2-2,1) = max(uin(i,j,ku2-3,1)*ratio_nyp1,smallr)
        uin(i,j,ku2-1,1) = max(uin(i,j,ku2-2,1)*ratio_nyp2,smallr)
        uin(i,j,ku2  ,1) = max(uin(i,j,ku2-1,1)*ratio_nyp3,smallr)
     end do
  end do
  !$OMP END PARALLEL DO

  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do j = ju1, ju2
     do i = iu1, iu2
        !zero gradient BC for velocities
        uin(i,j,ku2-2,2:3) = uin(i,j,ku2-3,2:3)/uin(i,j,ku2-3,1)
        uin(i,j,ku2-1,2:3) = uin(i,j,ku2-3,2:3)/uin(i,j,ku2-3,1)
        uin(i,j,ku2  ,2:3) = uin(i,j,ku2-3,2:3)/uin(i,j,ku2-3,1)
        uin(i,j,ku2-2,2:3) = uin(i,j,ku2-2,2:3)*uin(i,j,ku2-2,1)
        uin(i,j,ku2-1,2:3) = uin(i,j,ku2-1,2:3)*uin(i,j,ku2-1,1)
        uin(i,j,ku2  ,2:3) = uin(i,j,ku2  ,2:3)*uin(i,j,ku2  ,1)
        !uin(i,j,ku2-2,4) = max(uin(i,j,ku2-2,4)*uin(i,j,ku2-2,1),0.d0)
        !uin(i,j,ku2-1,4) = max(uin(i,j,ku2-1,4)*uin(i,j,ku2-1,1),0.d0)
        !uin(i,j,ku2  ,4) = max(uin(i,j,ku2  ,4)*uin(i,j,ku2  ,1),0.d0)
        uin(i,j,ku2-2,4) = max(uin(i,j,ku2-3,4), zero)
        uin(i,j,ku2-1,4) = max(uin(i,j,ku2-3,4), zero)
        uin(i,j,ku2  ,4) = max(uin(i,j,ku2-3,4), zero)
        !uin(i,j,ku2-2,4) = uin(i,j,ku2-3,4)
        !uin(i,j,ku2-1,4) = uin(i,j,ku2-3,4)
        !uin(i,j,ku2  ,4) = uin(i,j,ku2-3,4)
        !Tangential magnetic field
        !uin(i,j,ku2-2,6:7) = uin(i,j,ku2-3,6:7)
        !uin(i,j,ku2-1,6:7) = uin(i,j,ku2-3,6:7)
        !uin(i,j,ku2  ,6:7) = uin(i,j,ku2-3,6:7)
        !Tangential magnetic field
        uin(i,j,ku2-2,6:7) = zero !- uin(i,j,ku2-3,6:7)
        uin(i,j,ku2-1,6:7) = zero !- uin(i,j,ku2-4,6:7)
        uin(i,j,ku2  ,6:7) = zero !- uin(i,j,ku2-5,6:7)
        !Tangential magnetic field
        !uin(i,j,ku2-2,6:7) = uin(i,j,ku2-3,6:7)*(uin(i,j,ku2-2,1)/uin(i,j,ku2-3,1))**0.5
        !uin(i,j,ku2-1,6:7) = uin(i,j,ku2-3,6:7)*(uin(i,j,ku2-1,1)/uin(i,j,ku2-3,1))**0.5
        !uin(i,j,ku2  ,6:7) = uin(i,j,ku2-3,6:7)*(uin(i,j,ku2  ,1)/uin(i,j,ku2-3,1))**0.5
     end do
  end do
  !$OMP END PARALLEL DO

  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dbxdx, dbydy)
  do j = ju1, ju2-1
     do i = iu1, iu2-1
        !Normal magnetic field
        dbxdx = (uin(i+1,j,ku2-2,6) - uin(i,j,ku2-2,6))/dx
        dbydy = (uin(i  ,j,ku2-2,7) - uin(i,j,ku2-2,7))/dy
        uin(i,j,ku2-1,8) = uin(i,j,ku2-2,8) - dz*(dbxdx + dbydy)
        dbxdx = (uin(i+1,j,ku2-1,6) - uin(i,j,ku2-1,6))/dx
        dbydy = (uin(i  ,j,ku2-1,7) - uin(i,j,ku2-1,7))/dy
        uin(i,j,ku2  ,8) = uin(i,j,ku2-1,8) - dz*(dbxdx + dbydy)
     end do
  end do
  !$OMP END PARALLEL DO

#if ISO == 0
  ! Zero gradient Temperature
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(rho, E, rhovx, rhovy, rhovz) &
  !$OMP PRIVATE(Bxc, Byc, Bzc)
  do j = ju1, ju2-1
     !$acc loop
     do i = iu1, iu2-1
        rho = uin(i,j,ku2-3,1) 
        E = uin(i,j,ku2-3,5)
        rhovx = uin(i,j,ku2-3,2) 
        rhovy = uin(i,j,ku2-3,3) 
        rhovz = uin(i,j,ku2-3,4)
        Bxc = half*(uin(i,j,ku2-3,6) + uin(i+1,j  ,ku2-3,6))
        Byc = half*(uin(i,j,ku2-3,7) + uin(i  ,j+1,ku2-3,7))
        Bzc = half*(uin(i,j,ku2-3,8) + uin(i  ,j  ,ku2-2,8))
        call get_csq(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,csq)
        do k = ku2-2, ku2-1
           rho = uin(i,j,k,1)
           rhovx = uin(i,j,k,2) 
           rhovy = uin(i,j,k,3) 
           rhovz = uin(i,j,k,4)
           Bxc = half*(uin(i,j,k,6) + uin(i+1,j  ,k  ,6))
           Byc = half*(uin(i,j,k,7) + uin(i  ,j+1,k  ,7))
           Bzc = half*(uin(i,j,k,8) + uin(i  ,j  ,k+1,8))
           call get_E(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,csq)
           uin(i,j,k,5) = E !uin(i,j,ku2-3,5)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
#endif

#endif
  return
end subroutine zouter_ana
