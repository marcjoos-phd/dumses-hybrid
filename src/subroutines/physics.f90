!===============================================================================
!> \file physics.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is physics subroutines.
!! \details
!! Contains compute_dt(), comp_speed(), comp_speed_fast(), comp_mhd_flux(),
!! dissipation(), resistivity(), viscosity()
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
!! \b last \b modified: 05-06-2015
!<
!===============================================================================
!> Compute timestep \c dt
!===============================================================================
subroutine compute_dt(dt)
#if MPI == 1
  use mpi
  use mpi_var
#endif
  use params
  use variables, only: uin, qin, x, y
  implicit none

  real(dp), intent(out) :: dt
  real(dp) :: dt_proc
  real(dp) :: rho, vx, vy, vz, p, bx, by, bz
  real(dp) :: c2, b2, d2, cf
  real(dp) :: velx, vely, velz
  integer  :: i, j, k
  integer  :: ierr

  !$py begin_statement

  call primitive(zero)

  dt = courant*dx/smallc

  !$acc kernels loop !private(rho, vx, vy, vz, p, bx, by, bz, c2, b2, d2, cf, velx, vely, velz)
  !$OMP PARALLEL DO REDUCTION(MIN: dt) SCHEDULE(RUNTIME) PRIVATE(c2, b2, d2) &
  !$OMP PRIVATE(cf, velx, vely, velz, rho, vx, vy, vz, p, bx, by, bz)
  do k = 1, nz
     do j = 1, ny
        do i = 1, nx
           rho = qin(i,j,k,ir)
           vx  = qin(i,j,k,iu)
           vy  = qin(i,j,k,iv)
           vz  = qin(i,j,k,iw)
           p   = qin(i,j,k,ip)
           bx  = qin(i,j,k,iA)
           by  = qin(i,j,k,iB)
           bz  = qin(i,j,k,iC)

#if OACC == 1
           c2   = gamma*p/rho
           b2   = bx*bx + by*by + bz*bz
           d2   = half*(c2 + b2/rho)
           cf   = sqrt(d2 + sqrt(d2**2 - c2*bx*bx/rho))
           velx = cf + abs(vx)
           if (nx == 1) velx = zero
#if NDIM > 1
           cf   = sqrt(d2 + sqrt(d2**2 - c2*by*by/rho))
           vely = cf + abs(vy)
#if GEOM == CARTESIAN
           if ((Omega0 > zero) .and. (.not. fargo)) then
              vely = vely + 1.5d0*Omega0*(xmax-xmin)*half
           endif
#endif
           if (ny == 1) vely = zero
#endif
#if NDIM == 3
           cf   = sqrt(d2 + sqrt(d2**2 - c2*bz*bz/rho))
           velz = cf + abs(vz)
           if (nz == 1) velz = zero
#endif
#else
           call comp_speed(rho, p, vx, bx, by, bz, velx)
           if(nx == 1) velx = zero
#if NDIM > 1
           call comp_speed(rho, p, vy, by, bz, bx, vely)
#if GEOM == CARTESIAN
           if ((Omega0 > zero) .and. (.not. fargo)) then
              vely = vely + 1.5d0*Omega0*(xmax-xmin)*half
           endif
           if(ny == 1) vely = zero
#endif
#endif
#if NDIM == 3
           call comp_speed(rho, p, vz, bz, bx, by, velz)
           if(nz == 1) velz = zero
#endif
#endif

#if NDIM == 1
           dt = min(dt, dx/velx)
#endif
#if NDIM == 2
#if GEOM == CARTESIAN
           dt = min(dt, (one/(velx/dx + vely/dy)))
#else
           dt = min(dt, (one/(velx/dx + vely/dy/x(i))))
#endif
#endif
#if NDIM == 3
#if GEOM == CARTESIAN
           dt = min(dt, (one/(velx/dx + vely/dy + velz/dz)))
#endif
#if GEOM == CYLINDRICAL
           dt = min(dt, (one/(velx/dx + vely/dy/x(i) + velz/dz)))
#endif
#if GEOM == SPHERICAL
           dt = min(dt, (one/(velx/dx + vely/dy/x(i) + velz/dz/x(i)/sin(y(j)))))
#endif
#endif
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

#if MPI == 1
  dt_proc = courant*dt
  call MPI_Allreduce(dt_proc, dt, 1, MPI_DOUBLE_PRECISION, MPI_MIN &
       , MPI_COMM_WORLD, ierr)
  !$acc update device(dt)
#else
  dt = courant*dt
#endif

  return
end subroutine compute_dt
!===============================================================================
!> Compute the fastest velocity at which information is exchanged at the 
!! interface
!! Input:
!!  - \c rho: density
!!  - \c p: pressure
!!  - \c v: velocity in the normal direction to the interface
!!  - \c bn: magnetic field in the normal direction to the interface
!!  - \c bt1: magnetic field in the first tangential direction
!!  - \c bt2: magnetic field in the second tangential direction
!!
!! Output:
!!  - \c vel: fastest velocity at the interface
!<
!===============================================================================
subroutine comp_speed(rho, p, v, bn, bt1, bt2, vel)
  use params
  implicit none

  real(dp), intent(in)  :: rho, p, v, bn, bt1, bt2
  real(dp), intent(out) :: vel
  real(dp) :: c2, b2, d2, cf

  c2  = gamma*p/rho
  b2  = bn*bn + bt1*bt1 + bt2*bt2
  d2  = half*(c2 + b2/rho)
  cf  = sqrt(d2 + sqrt(d2**2 - c2*bn*bn/rho))
  vel = cf + abs(v)

  return
end subroutine comp_speed
!===============================================================================
!> Compute the fast magnetosonic velocity
!! Input:
!!  - \c rho: density
!!  - \c p: pressure
!!  - \c bn: magnetic field in the normal direction to the interface
!!  - \c bt1: magnetic field in the first tangential direction
!!  - \c bt2: magnetic field in the second tangential direction
!!
!! Output:
!!  - \c vel: fast magnetosonic velocity
!<
!===============================================================================
subroutine comp_speed_fast(rho, p, bn, bt1, bt2, vel)
  use params
  implicit none

  real(dp), intent(in)  :: rho, p, bn, bt1, bt2
  real(dp), intent(out) :: vel
  real(dp) :: c2, b2, d2, cf

  c2  = gamma*p/rho
  b2  = bn*bn + bt1*bt1 + bt2*bt2
  d2  = half*(c2 + b2/rho)
  cf  = sqrt(d2 + sqrt(d2**2 - c2*bn*bn/rho))
  vel = cf

  return
end subroutine comp_speed_fast
!===============================================================================
!> Compute the fast magnetosonic velocity
!! Input:
!!  - \c rho: density
!!  - \c bn: magnetic field in the normal direction to the interface
!!
!! Output:
!!  - \c vel: Alfven velocity
!<
!===============================================================================
subroutine comp_alfven_speed(rho, bn, vel)
  use params
  implicit none

  real(dp), intent(in)  :: rho, bn
  real(dp), intent(out) :: vel

  vel = sqrt(bn*bn/rho)

  return
end subroutine comp_alfven_speed
!===============================================================================
!> Compute MHD fluxes from conservative variables
!! Input:
!!  - \c rho: density
!!  - \c p: pressure
!!  - \c vn: velocity in the normal direction to the interface
!!  - \c bn: magnetic field in the normal direction to the interface
!!  - \c bt1: magnetic field in the first tangential direction
!!  - \c bt2: magnetic field in the second tangential direction
!!
!! Output:
!!  - \c consvar: conservative variables
!!  - \c flux: MHD flux
!<
!===============================================================================
subroutine comp_mhd_flux(rho, p, vn, vt1, vt2, bn, bt1, bt2, consvar, flux)
  use params
  implicit none

  real(dp), intent(in) :: rho, p, vn, vt1, vt2, bn, bt1, bt2
  real(dp) :: entho, ecin, emag, etot, ptot, ploc
  real(dp), dimension(nvar), intent(out) :: consvar, flux

#if ISO == 1
  ploc = rho*ciso**2
#else
  ploc = p
#endif

  ! Local variables
  entho = one/(gamma - one)
  ecin  = half*(vn*vn + vt1*vt1 + vt2*vt2)*rho
  emag  = half*(bn*bn + bt1*bt1 + bt2*bt2)
  etot  = ploc*entho + ecin + emag
  ptot  = ploc + emag
  
  ! Define conservative variables
  consvar = (/ rho, etot, rho*vn, bn, rho*vt1, bt1, rho*vt2, bt2 /)

  ! Compute fluxes
  flux(1) = rho*vn
  flux(2) = (etot + ptot)*vn - bn*(vn*bn + vt1*bt1 + vt2*bt2)
  flux(3) = rho*vn*vn + ptot - bn*bn
  flux(4) = zero
  flux(5) = rho*vn*vt1 - bn*bt1
  flux(6) = bt1*vn - bn*vt1
  flux(7) = rho*vn*vt2 - bn*bt2
  flux(8) = bt2*vn - bn*vt2
  
  return
end subroutine comp_mhd_flux
!===============================================================================
!> Dissipation processes
!===============================================================================
subroutine dissipation
  use params
  use variables
  implicit none

  ! Ohmic resistivity
  if (eta > zero) call resistivity

  ! Navier-Stokes viscosity
  if (nu > zero) then
     flux = zero
     call viscosity
  endif

  call boundary

  return
end subroutine dissipation
!===============================================================================
!> Resistivity update
!===============================================================================
subroutine resistivity
  use params
  use variables
  implicit none

  real(dp) :: dBxdy, dBxdz, dBydx, dBydz, dBzdx, dBzdy
  real(dp) :: jx, jy, jz
  integer :: i, j, k
#if ISO == 0
  real(dp) :: jxp, jyp, jzp, bx, by, bz
#endif
  
  !$py begin_statement
  
  emfx = zero; emfy = zero; emfz = zero
  dBxdy = zero; dBxdz = zero
  dBydx = zero; dBydz = zero
  dBzdx = zero; dBzdy = zero

  !!$acc kernels loop !-> can't be parallelized because of the subroutine call
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dBxdy, dBxdz, dBydx, dBydz) &
  !$OMP PRIVATE(dBzdx, dBzdy, jx, jy)
  do k = kf1, kf2
     do j = jf1, jf2
        do i = if1, if2
           dBydx = (uin(i,j,k,iB) - uin(i-1,j,k,iB))/dx
           dBzdx = (uin(i,j,k,iC) - uin(i-1,j,k,iC))/dx
#if NDIM > 1
           dBxdy = (uin(i,j,k,iA) - uin(i,j-1,k,iA))/dy
           dBzdy = (uin(i,j,k,iC) - uin(i,j-1,k,iC))/dy
#endif
#if NDIM == 3
           dBxdz = (uin(i,j,k,iA) - uin(i,j,k-1,iA))/dz
           dBydz = (uin(i,j,k,iB) - uin(i,j,k-1,iB))/dz
#endif
#if NDIM > 1
           jx = dBzdy - dBydz
#endif
           jy = dBxdz - dBzdx; jz = dBydx - dBxdy

#if NDIM > 1
           call get_eta(eta, x(i))
           emfx(i,j,k) = -eta*jx*dt
#endif
           
           call get_eta(eta, half*(x(i) + x(i-1)))
           emfy(i,j,k) = -eta*jy*dt
           emfz(i,j,k) = -eta*jz*dt
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

#if ISO == 0
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(bx, by, bz, jx, jxp, jy, jyp) &
  !$OMP PRIVATE(jz, jzp)
  do k = kf1, kf2
     do j = jf1, jf2
        do i = if1, if2
           
           ! First direction enery flux
#if NDIM == 1
           by = half*(uin(i,j,k,iB) + uin(i-1,j,k,iB))
#endif
#if NDIM > 1
           by = forth*(uin(i,j,k,iB) + uin(i-1,j,k,iB) & 
                     + uin(i,j+1,k,iB) + uin(i-1,j+1,k,iB))
#endif
#if NDIM < 3
           bz = half*(uin(i,j,k,iC) + uin(i-1,j,k,iC))
#endif
#if NDIM == 3
           bz = forth*(uin(i,j,k,iC) + uin(i-1,j,k,iC) &
                     + uin(i,j+1,k,iC) + uin(i-1,j+1,k,iC))
#endif
#if NDIM < 3
           jy = -(uin(i,j,k,iC) - uin(i-1,j,k,iC))/dx
#endif
#if NDIM == 3
           jy = (uin(i,j,k,iA) - uin(i,j,k-1,iA))/dz &
              - (uin(i,j,k,iC) - uin(i-1,j,k,iC))/dx
           jyp = (uin(i,j,k+1,iA) - uin(i,j,k,iA))/dz &
               - (uin(i,j,k+1,iC) - uin(i-1,j,k+1,iC))/dx
           jy = half*(jy + jyp)
#endif

#if NDIM == 1
           jz = (uin(i,j,k,iB) - uin(i-1,j,k,iB))/dx
#endif
#if NDIM > 1
           jz = (uin(i,j,k,iB) - uin(i-1,j,k,iB))/dx &
              - (uin(i,j,k,iA) - uin(i,j-1,k,iA))/dy
           jzp = (uin(i,j+1,k,iB) - uin(i-1,j+1,k,iB))/dx &
               - (uin(i,j+1,k,iA) - uin(i,j,k,iA))/dy
           jz = half*(jz + jzp)
#endif

           flux(i,j,k,ip,1) = flux(i,j,k,ip,1) - eta*(jy*bz - jz*by)*dt/dx

           ! Second direction energy flux
#if NDIM > 1
           bx = forth*(uin(i,j,k,iA) + uin(i,j-1,j,iA) &
                     + uin(i+1,j,k,iA) + uin(i+1,j-1,k,iA))
#if NDIM == 2
           bz = half*(uin(i,j,k,iC) + uin(i,j-1,k,iC))
#endif
#if NDIM == 3
           bz = forth*(uin(i,j,k,iC) + uin(i,j-1,j,iC) &
                     + uin(i,j,k+1,iC) + uin(i,j-1,k+1,iC))
#endif

#if NDIM == 2
           jx = (uin(i,j,k,iC) - uin(i,j-1,k,iC))/dy
#endif
#if NDIM == 3
           jx = (uin(i,j,k,iC) - uin(i,j-1,k,iC))/dy &
              - (uin(i,j,k,iB) - uin(i,j,k-1,iB))/dz
           jxp = (uin(i,j,k+1,iC) - uin(i,j-1,k+1,iC))/dy &
               - (uin(i,j,k+1,iB) - uin(i,j,k,iB))/dz
           jx = half*(jx + jxp)
#endif

           jz = (uin(i,j,k,iB) - uin(i-1,j,k,iB))/dx &
              + (uin(i,j,k,iA) - uin(i,j-1,k,iA))/dy
           jzp = (uin(i+1,j,k,iB) - uin(i,j,k,iB))/dx &
               + (uin(i+1,j,k,iA) - uin(i+1,j-1,k,iA))/dy
           jz = half*(jz + jzp)

           flux(i,j,k,ip,2) = flux(i,j,k,ip,2) - eta*(jz*bx - jx*bz)*dt/dy
#endif

           ! Thrid direction energy flux
#if NDIM == 3
           bx = forth*(uin(i,j,k,iA) + uin(i,j,k-1,iA) &
                     + uin(i+1,j,k,iA) + uin(i+1,j,k-1,iA))
           by = forth*(uin(i,j,k,iB) + uin(i,j,k-1,iB) &
                     + uin(i,j+1,k,iB) + uin(i,j+1,k-1,iB))

           jx = (uin(i,j,k,iC) - uin(i,j-1,k,iC))/dy &
              - (uin(i,j,k,iB) - uin(i,j,k-1,iB))/dz
           jxp = (uin(i,j+1,k,iC) - uin(i,j,k,iC))/dy &
               - (uin(i,j+1,k,iB) - uin(i,j+1,k-1,iB))/dz
           jx = half*(jx + jxp)

           jy = (uin(i,j,k,iA) - uin(i,j,k-1,iA))/dz &
              + (uin(i,j,k,iC) - uin(i-1,j,k,iC))/dx
           jyp = (uin(i+1,j,k,iA) - uin(i+1,j,k-1,iA))/dz &
               + (uin(i+1,j,k,iC) - uin(i,j,k,iC))/dx
           jy = half*(jy + jyp)

           flux(i,j,k,ip,3) = flux(i,j,k,ip,3) + eta*(jx*by - jy*bx)*dt/dz
#endif
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO
#endif

  ! Magnetic_field update
  call constrained_transport

  return
end subroutine resistivity
!===============================================================================
!> Viscosity update
!===============================================================================
subroutine viscosity
  use params
  use variables
  implicit none

  real(dp), dimension(3) :: dudx, dudy, dudz
  real(dp) :: rho, u, v, w, uR, uL, uRR, uRL, uLR, uLL
  real(dp) :: txx, tyy, tzz, txy, txz, tyz
  integer :: i, j, k, i0, j0, k0, idim

  !$py begin_statement
  
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(rho, u, v, w, uR, uL, uRR, uLL) &
  !$OMP PRIVATE(uRL, uLR, txx, tyy, tzz, txy, txz, tyz)
  do k = kf1, kf2
     do j = jf1, jf2
        do i = if1, if2

           ! First direction viscous flux
           rho = half*(uin(i,j,k,ir) + uin(i-1,j,k,ir))
#if ISO == 0
           u   = half*(uin(i,j,k,iu)/uin(i,j,k,ir) &
                     + uin(i-1,j,k,iu)/uin(i-1,j,k,ir))
           v   = half*(uin(i,j,k,iv)/uin(i,j,k,ir) &
                     + uin(i-1,j,k,iv)/uin(i-1,j,k,ir))
           w   = half*(uin(i,j,k,iw)/uin(i,j,k,ir) &
                     + uin(i-1,j,k,iw)/uin(i-1,j,k,ir))
#endif

           do idim = 1, 3
              uR = uin(i,j,k,idim+1)/uin(i,j,k,ir)
              uL = uin(i-1,j,k,idim+1)/uin(i,j,k,ir)
              dudx(idim) = (uR - uL)/dx
           enddo

#if NDIM > 1
           do idim = 1, 2
              uRR = uin(i,j+1,k,idim+1)/uin(i,j+1,k,ir)
              uRL = uin(i-1,j+1,k,idim+1)/uin(i-1,j+1,k,ir)
              uLR = uin(i,j-1,k,idim+1)/uin(i,j-1,k,ir)
              uLL = uin(i-1,j-1,k,idim+1)/uin(i-1,j-1,k,ir)
              uR = uRR + uRL; uL = uLR + uLL
              dudy(idim) = forth*(uR - uL)/dy
           enddo
#else
           dudy = zero
#endif

#if NDIM == 3
           do idim = 1, 3, 2
              uRR = uin(i,j,k+1,idim+1)/uin(i,j,k+1,ir)
              uRL = uin(i-1,j,k+1,idim+1)/uin(i-1,j,k+1,ir)
              uLR = uin(i,j,k-1,idim+1)/uin(i,j,k-1,ir)
              uLL = uin(i-1,j,k-1,idim+1)/uin(i-1,j,k-1,ir)
              uR = uRR + uRL; uL = uLR + uLL
              dudz(idim) = forth*(uR - uL)/dz
           enddo
#else
           dudz = zero
#endif

           txx = -two3rd*nu*rho*(two*dudx(1) - dudy(2) - dudz(3))
           txy = -nu*rho*(dudy(1) + dudx(2))
           txz = -nu*rho*(dudz(1) + dudx(3))
           flux(i,j,k,iu,1) = txx*dt/dx
           flux(i,j,k,iv,1) = txy*dt/dx
           flux(i,j,k,iw,1) = txz*dt/dx
#if ISO == 0
           flux(i,j,k,ip,1) = u*txx + v*txy + w*txz
           flux(i,j,k,ip,1) = flux(i,j,k,ip,1)*dt/dx
#endif

           ! Second direction viscous flux
#if NDIM > 1
           rho = half*(uin(i,j,k,ir) + uin(i,j-1,k,ir))
#if ISO == 0
           u   = half*(uin(i,j,k,iu)/uin(i,j,k,ir) &
                     + uin(i,j-1,k,iu)/uin(i,j-1,k,ir))
           v   = half*(uin(i,j,k,iv)/uin(i,j,k,ir) &
                     + uin(i,j-1,k,iv)/uin(i,j-1,k,ir))
           w   = half*(uin(i,j,k,iw)/uin(i,j,k,ir) &
                     + uin(i,j-1,k,iw)/uin(i,j-1,k,ir))
#endif

           do idim = 1, 3
              uR = uin(i,j,k,idim+1)/uin(i,j,k,ir)
              uL = uin(i,j-1,k,idim+1)/uin(i,j-1,k,ir)
              dudy(idim) = (uR - uL)/dy
           enddo

           do idim = 1, 2
              uRR = uin(i+1,j,k,idim+1)/uin(i+1,j,k,ir)
              uRL = uin(i+1,j-1,k,idim+1)/uin(i+1,j-1,k,ir)
              uLR = uin(i-1,j,k,idim+1)/uin(i-1,j,k,ir)
              uLL = uin(i-1,j-1,k,idim+1)/uin(i-1,j-1,k,ir)
              uR = uRR + uRL; uL = uLR + uLL
              dudx(idim) = forth*(uR - uL)/dx
           enddo

#if NDIM == 3
           do idim = 2, 3
              uRR = uin(i,j,k+1,idim+1)/uin(i,j,k+1,ir)
              uRL = uin(i,j-1,k+1,idim+1)/uin(i,j-1,k+1,ir)
              uLR = uin(i,j,k-1,idim+1)/uin(i,j,k-1,ir)
              uLL = uin(i,j-1,k-1,idim+1)/uin(i,j-1,k-1,ir)
              uR = uRR + uRL; uL = uLR + uLL
              dudz(idim) = forth*(uR - uL)/dz
           enddo
#else
           dudz = zero
#endif

           tyy = -two3rd*nu*rho*(two*dudy(2) - dudx(1) - dudz(3))
           txy = -nu*rho*(dudy(1) + dudx(2))
           tyz = -nu*rho*(dudz(2) + dudy(3))
           flux(i,j,k,iu,2) = txy*dt/dy
           flux(i,j,k,iv,2) = tyy*dt/dy
           flux(i,j,k,iw,2) = tyz*dt/dy
#if ISO == 0
           flux(i,j,k,ip,2) = u*txy + v*tyy + w*tyz
           flux(i,j,k,ip,2) = flux(i,j,k,ip,2)*dt/dy
#endif
#endif

           ! Third direction viscous flux
#if NDIM == 3
           rho = half*(uin(i,j,k,ir) + uin(i,j,k-1,ir))
#if ISO == 0
           u   = half*(uin(i,j,k,iu)/uin(i,j,k,ir) &
                     + uin(i,j,k-1,iu)/uin(i,j,k-1,ir))
           v   = half*(uin(i,j,k,iv)/uin(i,j,k,ir) &
                     + uin(i,j,k-1,iv)/uin(i,j,k-1,ir))
           w   = half*(uin(i,j,k,iw)/uin(i,j,k,ir) &
                     + uin(i,j,k-1,iw)/uin(i,j,k-1,ir))
#endif

           do idim = 1, 3
              uR = uin(i,j,k,idim+1)/uin(i,j,k,ir)
              uL = uin(i,j,k-1,idim+1)/uin(i,j,k-1,ir)
              dudz(idim) = (UR - uL)/dz
           enddo

           do idim = 1, 3, 2
              uRR = uin(i+1,j,k,idim+1)/uin(i+1,j,k,ir)
              uRL = uin(i+1,j,k-1,idim+1)/uin(i+1,j,k-1,ir)
              uLR = uin(i-1,j,k,idim+1)/uin(i-1,j,k,ir)
              uLL = uin(i-1,j,k-1,idim+1)/uin(i-1,j,k-1,ir)
              uR = uRR + uRL; uL = uLR + uLL
              dudx(idim) = forth*(uR - uL)/dx
           enddo

           do idim = 2, 3
              uRR = uin(i,j+1,k,idim+1)/uin(i,j+1,k,ir)
              uRL = uin(i,j+1,k-1,idim+1)/uin(i,j+1,k-1,ir)
              uLR = uin(i,j-1,k,idim+1)/uin(i,j-1,k,ir)
              uLL = uin(i,j-1,k-1,idim+1)/uin(i,j-1,k-1,ir)
              uR = uRR + uRL; uL = uLR + uLL
              dudy(idim) = forth*(uR - uL)/dy
           enddo

           tzz = -two3rd*nu*rho*(two*dudz(3) - dudx(1) - dudy(2))
           txz = -nu*rho*(dudz(1) + dudx(3))
           tyz = -nu*rho*(dudz(2) + dudy(3))
           flux(i,j,k,iu,3) = txz*dt/dz
           flux(i,j,k,iv,3) = tyz*dt/dz
           flux(i,j,k,iw,3) = tzz*dt/dz
#if ISO == 0
           flux(i,j,k,ip,3) = u*txz + v*tyz + w*tzz
           flux(i,j,k,ip,3) = flux(i,j,k,ip,3)*dt/dz
#endif
#endif
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

  ! Conservative update for momentum and enery
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(i0, j0, k0)
  do idim = 1, ndim
     i0 = 0; j0 = 0; k0 = 0
     if (idim == 1) i0 = 1
     if (idim == 2) j0 = 1
     if (idim == 3) k0 = 1
     
     do k = 1, nz
        do j = 1, ny
           do i = 1, nx
              uin(i,j,k,iu:ip) = uin(i,j,k,iu:ip) &
               + (flux(i,j,k,iu:ip,idim) - flux(i+i0,j+j0,k+k0,iu:ip,idim))
           enddo
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

  return
end subroutine viscosity
