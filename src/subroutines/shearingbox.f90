!===============================================================================
!> \file shearingbox.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is shearing box subroutines.
!! \details
!! Contains shearing_boundary(), shearing_flux(), shearing_emf(), 
!! shearing_slope(), shearing_boundary_old(), shearing_flux_old(), 
!! shearing_emf_old(), shearing_slope_old()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Romain Teyssier, 
!! Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          02-20-2014
!! \b last \b modified: 05-11-2015
!<
!===============================================================================
!> Compute shearing boundary conditions
!===============================================================================
subroutine shearing_boundary(itime)
  use params
  use variables
  use mpi_var
  implicit none

  real(dp), intent(in) :: itime
  real(dp), dimension(:,:,:,:), allocatable :: inner_dq, outer_dq
  real(dp), dimension(:,:,:,:), allocatable :: inner, outer
  real(dp) :: deltay, epsi
  real(dp) :: eps, lambda
  integer  :: jp
  integer  :: jremap, jremapp
  integer  :: i, j, k, ighost, ivar

  !$py begin_statement

  deltay = 1.5d0*Omega0*(xmax - xmin)*itime
  deltay = mod(deltay, (ymax - ymin))
  epsi   = mod(deltay, dy)
  jp     = int(deltay/dy)

  allocate(inner_dq(1:nghost,ju1:ju2,ku1:ku2,1:nvar+3))
  allocate(outer_dq(1:nghost,ju1:ju2,ku1:ku2,1:nvar+3))
  allocate(inner(1:nghost,ju1:ju2,ku1:ku2,1:nvar+3))
  allocate(outer(1:nghost,ju1:ju2,ku1:ku2,1:nvar+3))

  !$acc data pcreate(inner_dq, outer_dq, inner, outer)
  if (nxslice*nyslice*nzslice > 1) then
     !$acc update host(uin)
     call get_shear_quant(inner, outer, itime)
     !$acc update device(inner, outer)
   
     call shearing_slope(inner, inner_dq)
     call shearing_slope(outer, outer_dq)
  else
     !$acc kernels
     inner = uin(iu1+nghost:iu1+2*nghost-1,:,:,:)
     outer = uin(iu2-2*nghost+1:iu2-nghost,:,:,:)
     !$acc end kernels
     !$py start_timing shearing slope
     call shearing_slope(inner, inner_dq)
     call shearing_slope(outer, outer_dq)
     !$py end_timing shearing slope
  endif

  !$acc kernels loop private(jremap, jremapp)
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(jremap, jremapp, eps, lambda)
  do k = ku1, ku2
     do j = ju1+1, ju2-1
        ! Inner boundary conditions
        if (xposition == 0) then
           if (nxslice*nyslice*nzslice > 1) then
              jremap = j - 1
              jremapp = jremap + 1
           else
              jremap = j - jp - 1
              jremapp = jremap + 1
              if (jremap <= 0) jremap = jremap + ny
              if (jremapp <= 0) jremapp = jremapp + ny
           endif
           eps    = one - epsi/dy
           lambda = half*eps*(eps - one)
  
           do ivar = 1, nvar-2
              uin(iu1:iu1+2,j,k,ivar) = &
                   (one - eps)*outer(1:nghost,jremap,k,ivar) &
                 + eps*outer(1:nghost,jremapp,k,ivar) &
                 + lambda*(outer_dq(1:nghost,jremap,k,ivar) &
                         - outer_dq(1:nghost,jremapp,k,ivar))
           enddo
  
           uin(iu1:iu1+2,j,k,nvar-1) = outer(1:nghost,jremap,k,nvar-1) &
                                     + eps*outer_dq(1:nghost,jremap,k,nvar-1)
  
           uin(iu1:iu1+2,j,k,nvar) = &
                (one - eps)*outer(1:nghost,jremap,k,nvar) &
              + eps*outer(1:nghost,jremapp,k,nvar) &
              + lambda*(outer_dq(1:nghost,jremap,k,nvar) &
                      - outer_dq(1:nghost,jremapp,k,nvar))
        endif
  
        ! Outer boundary conditions
        if (xposition == nxslice - 1) then
           if (nxslice*nyslice*nzslice > 1) then
              jremap = j
              jremapp = jremap + 1
           else
              jremap = j + jp
              jremapp = jremap + 1
              if (jremap > ny) jremap = jremap - ny
              if (jremapp > ny) jremapp = jremapp - ny
           endif
           eps     = epsi/dy
           lambda = half*eps*(eps - one)
  
           do ivar = 1, nvar-3
              uin(iu2-2:iu2,j,k,ivar) = &
                   (one - eps)*inner(1:nghost,jremap,k,ivar) &
                 + eps*inner(1:nghost,jremapp,k,ivar) &
                 + lambda*(inner_dq(1:nghost,jremapp,k,ivar) &
                         - inner_dq(1:nghost,jremap,k,ivar))
           enddo
  
           uin(iu2-1:iu2,j,k,iA) = &
                (one - eps)*inner(2:nghost,jremap,k,iA) &
              + eps*inner(2:nghost,jremapp,k,iA) &
              + lambda*(inner_dq(2:nghost,jremapp,k,iA) &
                      - inner_dq(2:nghost,jremap,k,iA))
  
           uin(iu2-2:iu2,j,k,iB) = inner(1:nghost,jremap,k,iB) &
                                 + eps*inner_dq(1:nghost,jremap,k,iB)
  
           uin(iu2-2:iu2,j,k,iC) = &
                (one - eps)*inner(1:nghost,jremap,k,iC) &
              + eps*inner(1:nghost,jremapp,k,iC) &
              + lambda*(inner_dq(1:nghost,jremapp,k,iC) &
                      - inner_dq(1:nghost,jremap,k,iC))
  
           uin(iu2,j,k,iA+3) = &
                (one - eps)*inner(nghost,jremap,k,iA+3) &
              + eps*inner(nghost,jremapp,k,iA+3) &
              + lambda*(inner_dq(nghost,jremapp,k,iA+3) &
                      - inner_dq(nghost,jremap,k,iA+3))
        endif
     enddo
  enddo
  !$OMP END PARALLEL DO

  ! Inner radial boundary conditions, right face B-field
#if MPI == 1
  if (xposition == 0) then
#endif
     !$acc kernels loop
     do i = iu1+2, iu1, -1
        uin(i,:,:,iA+3) = uin(i+1,:,:,iA)
     enddo
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do j = ju1, ju2-1
        uin(iu1:iu1+2,j,:,iB+3) = uin(iu1:iu1+2,j+1,:,iB)
     enddo
     !$OMP END PARALLEL DO
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = ku1, ku2-1
        uin(iu1:iu1+2,:,k,iC+3) = uin(iu1:iu1+2,:,k+1,iC)
     enddo
     !$OMP END PARALLEL DO
#if MPI == 1
  endif
#endif

  ! Outer radial boundary conditions, right face B-field
#if MPI == 1
  if (xposition == nxslice - 1) then
#endif
     !$acc kernels loop
     do i = iu2-2, iu2-1
        uin(i,:,:,iA+3) = uin(i+1,:,:,iA)
     enddo
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do j = ju1, ju2-1
        uin(iu2-2:iu2,j,:,iB+3) = uin(iu2-2:iu2,j+1,:,iB)
     enddo
     !$OMP END PARALLEL DO
     !$acc kernels loop
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = ku1, ku2-1
        uin(iu2-2:iu2,:,k,iC+3) = uin(iu2-2:iu2,:,k+1,iC)
     enddo
     !$OMP END PARALLEL DO
#if MPI == 1
  endif
#endif

  !$acc end data

  deallocate(inner_dq, outer_dq)
  deallocate(inner, outer)

  return
end subroutine shearing_boundary
!===============================================================================
!> Compute fluxes at shearing boundary
!===============================================================================
subroutine shearing_flux(itime)
  use params
  use variables
  use mpi_var
  implicit none

  real(dp), intent(in) :: itime
  real(dp), dimension(:,:), allocatable :: inner, outer
  real(dp) :: deltay, epsi
  real(dp) :: eps
  integer  :: jp
  integer  :: jremap, jremapp
  integer  :: j, k

  !$py begin_statement

  deltay = 1.5d0*Omega0*(xmax - xmin)*itime
  deltay = mod(deltay, (ymax - ymin))
  epsi   = mod(deltay, dy)
  jp = int(deltay/dy)

  allocate(inner(jf1:jf2+1,kf1:kf2))
  allocate(outer(jf1-1:jf2,kf1:kf2))

  !$acc data pcreate(inner, outer)
  if (nxslice*nyslice*nzslice > 1) then
     !$acc update host(flux)
     call get_shear_flux(inner, outer, itime)
     !$acc update device(inner, outer)
  else
     !$acc kernels
     inner(jf1:jf2,kf1:kf2) = flux(if1,jf1:jf2,kf1:kf2,1,1)
     outer(jf1:jf2,kf1:kf2) = flux(if2,jf1:jf2,kf1:kf2,1,1)
     !$acc end kernels
  endif

  !$acc kernels loop private(jremap, jremapp)
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(jremap, jremapp, eps)
  do k = kf1, kf2
     do j = jf1, jf2
        ! Inner boundary conditions
        if (xposition == 0) then
           if (nxslice*nyslice*nzslice > 1) then
              jremap = j - 1
              jremapp = jremap + 1
           else
              jremap = j - jp - 1
              jremapp = jremap + 1
              if (jremap <= 0) jremap = jremap + ny
              if (jremapp <= 0) jremapp = jremapp + ny
           endif
           eps     = one - epsi/dy

           flux(iu1+nghost,j,k,1,1) = flux(iu1+nghost,j,k,1,1) &
                + (one - eps)*outer(jremap,k) + eps*outer(jremapp,k)
           flux(iu1+nghost,j,k,1,1) = half*flux(iu1+nghost,j,k,1,1)
        endif

        ! Outer boundary conditions
        if (xposition == nxslice - 1) then
           if (nxslice*nyslice*nzslice > 1) then
              jremap = j
              jremapp = jremap + 1
           else
              jremap = j + jp
              jremapp = jremap + 1
              if (jremap > ny) jremap = jremap - ny
              if (jremapp > ny) jremapp = jremapp - ny
           endif
           eps     = epsi/dy

           flux(iu2-nghost+1,j,k,1,1) = flux(iu2-nghost+1,j,k,1,1) &
                + (one - eps)*inner(jremap,k) + eps*inner(jremapp,k)
           flux(iu2-nghost+1,j,k,1,1) = half*flux(iu2-nghost+1,j,k,1,1)
        endif
     enddo
  enddo
  !$OMP END PARALLEL DO

  !$acc end data
  
  deallocate(inner, outer)

  call shearing_emf(itime)

  return
end subroutine shearing_flux
!===============================================================================
!> Compute fluxes at shearing boundary
!===============================================================================
subroutine shearing_emf(itime)
  use params
  use variables
  use mpi_var
  implicit none

  real(dp), intent(in) :: itime
  real(dp), dimension(:,:), allocatable :: inner, outer
  real(dp) :: deltay, epsi
  real(dp) :: eps
  integer  :: jp
  integer  :: jremap, jremapp
  integer  :: j, k

  !$py begin_statement

  deltay = 1.5d0*Omega0*(xmax - xmin)*itime
  deltay = mod(deltay, (ymax - ymin))
  epsi   = mod(deltay, dy)
  jp     = int(deltay/dy)

  allocate(inner(ju1:ju2,ku1:ku2))
  allocate(outer(ju1:ju2,ku1:ku2))

  !$acc data pcreate(inner, outer)
  if (nxslice*nyslice*nzslice > 1) then
     !$acc update host(emfy)
     call get_shear_emf(inner, outer, itime)
     !$acc update device(inner, outer)
  else
     !$acc kernels
     inner = emfy(iu1+nghost,ju1:ju2,ku1:ku2)
     outer = emfy(iu2-nghost+1,ju1:ju2,ku1:ku2)
     !$acc end kernels
  endif

  !$acc kernels loop private(jremap, jremapp)
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(jremap, jremapp, eps)
  do k = ku1, ku2
     do j = ju1+1, ju2
        ! Inner boundary conditions
        if (xposition == 0) then
           if (nxslice*nyslice*nzslice > 1) then
              jremap = j - 1
              jremapp = jremap + 1
           else
              jremap = j - jp - 1
              jremapp = jremap + 1
              if (jremap <= 0) jremap = jremap + ny
              if (jremapp <= 0) jremapp = jremapp + ny
           endif
           eps     = one - epsi/dy
           
           emfy(iu1+nghost,j,k) = emfy(iu1+nghost,j,k) &
                + (one - eps)*outer(jremap,k) + eps*outer(jremapp,k)
           emfy(iu1+nghost,j,k) = half*emfy(iu1+nghost,j,k)
        endif
     enddo
  enddo
  !$OMP END PARALLEL DO

  !$acc kernels loop private(jremap, jremapp)
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(jremap, jremapp, eps)
  do k = ku1, ku2
     do j = ju1, ju2-1
        ! Outer boundary conditions
        if (xposition == nxslice - 1) then
           if (nxslice*nyslice*nzslice > 1) then
              jremap = j
              jremapp = jremap + 1
           else
              jremap = j + jp
              jremapp = jremap + 1
              if (jremap > ny) jremap = jremap - ny
              if (jremapp > ny) jremapp = jremapp - ny
           endif
           eps     = epsi/dy

           emfy(iu2-nghost+1,j,k) = emfy(iu2-nghost+1,j,k) &
                + (one - eps)*inner(jremap,k) + eps*inner(jremapp,k)
           emfy(iu2-nghost+1,j,k) = half*emfy(iu2-nghost+1,j,k)
        endif
     enddo
  enddo
  !$OMP END PARALLEL DO

  !$acc end data

  deallocate(inner, outer)

  return
end subroutine shearing_emf
!===============================================================================
!> Compute slopes at shearing boundaries
!===============================================================================
subroutine shearing_slope(quant, dq)
  use params
  implicit none

  real(dp), dimension(nghost,ju1:ju2,ku1:ku2,nvar+3), intent(in)  :: quant
  real(dp), dimension(nghost,ju1:ju2,ku1:ku2,nvar+3), intent(out) :: dq
  integer :: jslope_type
  integer :: ighost, ivar, j, k
  real(dp) :: dsgn, dlim, dcen, dlft, drgt, slop

  jslope_type = min(slope_type, 2)
  if (jslope_type == 0) then
     dq = zero
     return
  endif

#if NDIM == 3
  if (jslope_type == 1 .or. jslope_type == 2) then
     !$acc kernels
     do ighost = 1, nghost
        !$acc loop vector(128)
        !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dsgn, dlim, dcen, dlft) &
        !$OMP PRIVATE(drgt, slop)
        do k = ku1, ku2
           do j = ju1+1, ju2-1
              do ivar = 1, nvar-2
                 dlft = jslope_type*(quant(ighost,j,k,ivar) &
                                   - quant(ighost,j-1,k,ivar))
                 drgt = jslope_type*(quant(ighost,j+1,k,ivar) &
                                   - quant(ighost,j,k,ivar))
                 dcen = half*(dlft + drgt)/jslope_type
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft), abs(drgt))
                 dlim = slop
                 if ((dlft*drgt) <= zero) dlim = zero
                 dq(ighost,j,k,ivar) = dsgn*min(dlim, abs(dcen))
              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO
     enddo
     !$acc end kernels

     !$acc kernels
     do ighost = 1, nghost
        !$acc loop vector(128)
        !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dsgn, dlim, dcen, dlft) &
        !$OMP PRIVATE(drgt, slop)
        do k = ku1, ku2
           do j = ju1+1, ju2-1
              dlft = jslope_type*(quant(ighost,j,k,nvar) &
                                - quant(ighost,j-1,k,nvar))
              drgt = jslope_type*(quant(ighost,j+1,k,nvar) &
                                - quant(ighost,j,k,nvar))
              dcen = half*(dlft + drgt)/jslope_type
              dsgn = sign(one, dcen)
              slop = min(abs(dlft), abs(drgt))
              dlim = slop
              if ((dlft*drgt) <= zero) dlim = zero
              dq(ighost,j,k,nvar) = dsgn*min(dlim, abs(dcen))
           enddo
        enddo
        !$OMP END PARALLEL DO
     enddo
     !$acc end kernels

     !$acc kernels
     do ighost = 1, nghost
        !$acc loop vector(128)
        !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dsgn, dlim, dcen, dlft) &
        !$OMP PRIVATE(drgt, slop)
        do k = ku1, ku2
           do j = ju1+1, ju2-1
              do ivar = nvar+1, nvar+3, 2
                 dlft = jslope_type*(quant(ighost,j,k,ivar) &
                                   - quant(ighost,j-1,k,ivar))
                 drgt = jslope_type*(quant(ighost,j+1,k,ivar) &
                                   - quant(ighost,j,k,ivar))
                 dcen = half*(dlft + drgt)/jslope_type
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft), abs(drgt))
                 dlim = slop
                 if ((dlft*drgt) <= zero) dlim = zero
                 dq(ighost,j,k,ivar) = dsgn*min(dlim, abs(dcen))
              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO
     enddo
     !$acc end kernels
  else
     write(*,*) "Unknown slope type!"
     stop
  endif

  !$acc kernels
  do ighost = 1, nghost
     !$acc loop vector(128)
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = ku1, ku2
        do j = ju1+1, ju2-1
          dq(ighost,j,k,nvar-1) = quant(ighost,j+1,k,nvar-1) &
                                   - quant(ighost,j,k,nvar-1)
           dq(ighost,j,k,nvar+2) = quant(ighost,j+1,k,nvar+2) &
                                   - quant(ighost,j,k,nvar+2)
        enddo
     enddo
     !$OMP END PARALLEL DO
  enddo
  !$acc end kernels
#endif

  return
end subroutine shearing_slope
!===============================================================================
!> Compute shearing boundary conditions - old version
!===============================================================================
subroutine shearing_boundary_old(itime)
  use params
  use variables
  use mpi_var
  implicit none

  real(dp) :: itime
  real(dp), dimension(:,:,:,:), allocatable :: inner_dq, outer_dq
  real(dp), dimension(:,:,:,:), allocatable :: inner_quant, outer_quant
  real(dp), save :: deltay, epsi
  real(dp) :: eps, lambda
  integer  :: jp, jremap, jremapp
  integer  :: i, j, k, ighost, ivar

  if (verbose) print*, '> Entering shearing_boundary'

  deltay = 1.5d0*Omega0*(xmax - xmin)*itime
  deltay = mod(deltay, (ymax - ymin))
  epsi   = mod(deltay, dy)
  jp = int(deltay/dy)

#if MPI == 1
  allocate(inner_dq(ju1:ju1+ny*nyslice+2*nghost-1,ku1:ku2,nvar+3,nghost))
  allocate(outer_dq(ju1:ju1+ny*nyslice+2*nghost-1,ku1:ku2,nvar+3,nghost))
  allocate(inner_quant(ju1:ju1+ny*nyslice+2*nghost-1,ku1:ku2,nvar+3,nghost))
  allocate(outer_quant(ju1:ju1+ny*nyslice+2*nghost-1,ku1:ku2,nvar+3,nghost))
  call get_shearing_quant_old(inner_quant, outer_quant, deltay)
  call shearing_slope_old(inner_quant, inner_dq)
  call shearing_slope_old(outer_quant, outer_dq)
#else
  allocate(inner_dq(ju1:ju2,ku1:ku2,nvar+3,nghost))
  allocate(outer_dq(ju1:ju2,ku1:ku2,nvar+3,nghost))
  allocate(inner_quant(ju1:ju2,ku1:ku2,nvar+3,nghost))
  allocate(outer_quant(ju1:ju2,ku1:ku2,nvar+3,nghost))

  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do ighost = 1, nghost
     inner_quant(ju1:ju2,ku1:ku2,1:nvar+3,ighost) = &
          uin(iu1+nghost+ighost-1,ju1:ju2,ku1:ku2,1:nvar+3)
     outer_quant(ju1:ju2,ku1:ku2,1:nvar+3,ighost) = &
          uin(iu2-2*nghost+ighost,ju1:ju2,ku1:ku2,1:nvar+3)
  enddo
  !$OMP END PARALLEL DO
  call shearing_slope_old(inner_quant, inner_dq)
  call shearing_slope_old(outer_quant, outer_dq)
#endif

  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(jremap, jremapp, eps, lambda)
  do k = ku1, ku2
     do j = 1, ny
        ! Inner boundary conditions
#if MPI == 1
        if (xposition == 0) then
           jremap = j + ny*yposition - jp - 1
#else
           jremap = j - jp - 1
#endif
           jremapp = jremap + 1
           eps     = one - epsi/dy

#if MPI == 1
           if (jremap <= 0) jremap = jremap + ny*nyslice
           if (jremapp <= 0) jremapp = jremapp + ny*nyslice
#else
           if (jremap <= 0) jremap = jremap + ny
           if (jremapp <= 0) jremapp = jremapp + ny
#endif
           lambda = half*eps*(eps - one)

           do ivar = 1, nvar-2
              uin(iu1:iu1+2,j,k,ivar) = &
                   (one - eps)*outer_quant(jremap,k,ivar,1:nghost) &
                 + eps*outer_quant(jremapp,k,ivar,1:nghost) &
                 + lambda*(outer_dq(jremap,k,ivar,1:nghost) &
                         - outer_dq(jremapp,k,ivar,1:nghost))
           enddo

           uin(iu1:iu1+2,j,k,nvar-1) = outer_quant(jremap,k,nvar-1,1:nghost) &
                                     + eps*outer_dq(jremap,k,nvar-1,1:nghost)

           uin(iu1:iu1+2,j,k,nvar) = &
                (one - eps)*outer_quant(jremap,k,nvar,1:nghost) &
              + eps*outer_quant(jremapp,k,nvar,1:nghost) &
              + lambda*(outer_dq(jremap,k,nvar,1:nghost) &
                      - outer_dq(jremapp,k,nvar,1:nghost))
#if MPI == 1
        endif
#endif

        ! Outer boundary conditions
#if MPI == 1
        if (xposition == nxslice - 1) then
           jremap = j + ny*yposition + jp
#else
           jremap = j + jp
#endif
           jremapp = jremap + 1
           eps     = epsi/dy

#if MPI == 1
           if (jremap >= ny*nyslice + 1) jremap = jremap - ny*nyslice
           if (jremapp >= ny*nyslice + 1) jremapp = jremapp - ny*nyslice
#else
           if (jremap >= ny + 1) jremap = jremap - ny
           if (jremapp >= ny + 1) jremapp = jremapp - ny
#endif
           lambda = half*eps*(eps - one)

           do ivar = 1, nvar-3
              uin(iu2-2:iu2,j,k,ivar) = &
                   (one - eps)*inner_quant(jremap,k,ivar,1:nghost) &
                 + eps*inner_quant(jremapp,k,ivar,1:nghost) &
                 + lambda*(inner_dq(jremapp,k,ivar,1:nghost) &
                         - inner_dq(jremap,k,ivar,1:nghost))
           enddo

           uin(iu2-1:iu2,j,k,iA) = &
                (one - eps)*inner_quant(jremap,k,iA,2:nghost) &
              + eps*inner_quant(jremapp,k,iA,2:nghost) &
              + lambda*(inner_dq(jremapp,k,iA,2:nghost) &
                      - inner_dq(jremap,k,iA,2:nghost))

           uin(iu2-2:iu2,j,k,iB) = inner_quant(jremap,k,iB,1:nghost) &
                                 + eps*inner_dq(jremap,k,iB,1:nghost)

           uin(iu2-2:iu2,j,k,iC) = &
                (one - eps)*inner_quant(jremap,k,iC,1:nghost) &
              + eps*inner_quant(jremapp,k,iC,1:nghost) &
              + lambda*(inner_dq(jremapp,k,iC,1:nghost) &
                      - inner_dq(jremap,k,iC,1:nghost))

           uin(iu2,j,k,iA+3) = &
                (one - eps)*inner_quant(jremap,k,iA+3,nghost) &
              + eps*inner_quant(jremapp,k,iA+3,nghost) &
              + lambda*(inner_dq(jremapp,k,iA+3,nghost) &
                      - inner_dq(jremap,k,iA+3,nghost))
#if MPI == 1
        endif
#endif
     enddo
  enddo
  !$OMP END PARALLEL DO

  ! Inner radial boundary conditions, right face B-field
#if MPI == 1
  if (xposition == 0) then
#endif
     do i = iu1+2, iu1, -1
        uin(i,:,:,iA+3) = uin(i+1,:,:,iA)
     enddo
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do j = ju1, ju2-1
        uin(iu1:iu1+2,j,:,iB+3) = uin(iu1:iu1+2,j+1,:,iB)
     enddo
     !$OMP END PARALLEL DO
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = ku1, ku2-1
        uin(iu1:iu1+2,:,k,iC+3) = uin(iu1:iu1+2,:,k+1,iC)
     enddo
     !$OMP END PARALLEL DO
#if MPI == 1
  endif
#endif

  ! Outer radial boundary conditions, right face B-field
#if MPI == 1
  if (xposition == nxslice - 1) then
#endif
     do i = iu2-2, iu2-1
        uin(i,:,:,iA+3) = uin(i+1,:,:,iA)
     enddo
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do j = ju1, ju2-1
        uin(iu2-2:iu2,j,:,iB+3) = uin(iu2-2:iu2,j+1,:,iB)
     enddo
     !$OMP END PARALLEL DO
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = ku1, ku2-1
        uin(iu2-2:iu2,:,k,iC+3) = uin(iu2-2:iu2,:,k+1,iC)
     enddo
     !$OMP END PARALLEL DO
#if MPI == 1
  endif
#endif

  deallocate(inner_dq, outer_dq, inner_quant, outer_quant)

  return
end subroutine shearing_boundary_old
!===============================================================================
!> Compute fluxes at shearing boundary - old version
!===============================================================================
subroutine shearing_flux_old(itime)
  use mpi_var
  use params
  use variables
  implicit none

  real(dp) :: itime
  real(dp), dimension(:,:), allocatable :: inner_flux, outer_flux
  real(dp) :: deltay, epsi
  real(dp) :: eps
  integer  :: jp, jremap, jremapp
  integer  :: j, k

  if (verbose) print*, '> Entering shearing_flux'

  deltay = 1.5d0*Omega0*(xmax - xmin)*itime
  deltay = mod(deltay, (ymax - ymin))
  epsi   = mod(deltay, dy)
  jp = int(deltay/dy)

#if MPI == 1
  allocate(inner_flux(jf1:jf1+ny*nyslice,kf1:kf2))
  allocate(outer_flux(jf1:jf1+ny*nyslice,kf1:kf2))
  call get_shearing_flux_old(inner_flux, outer_flux, deltay)
#else
  allocate(inner_flux(jf1:jf2,kf1:kf2))
  allocate(outer_flux(jf1:jf2,kf1:kf2))
  inner_flux(jf1:jf2,kf1:kf2) = flux(iu1+nghost,jf1:jf2,kf1:kf2,1,1)
  outer_flux(jf1:jf2,kf1:kf2) = flux(iu2-nghost+1,jf1:jf2,kf1:kf2,1,1)
#endif
  
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(jremap, jremapp, eps)
  do k = kf1, kf2
     do j = jf1, jf2
        ! Inner boundary conditions
#if MPI == 1
        if (xposition == 0) then
           jremap = j + ny*yposition - jp - 1
#else
           jremap = j - jp - 1
#endif
           jremapp = jremap + 1
           eps     = one - epsi/dy

#if MPI == 1
           if (jremap <= 0) jremap = jremap + ny*nyslice
           if (jremapp <= 0) jremapp = jremapp + ny*nyslice
#else
           if (jremap <= 0) jremap = jremap + ny
           if (jremapp <= 0) jremapp = jremapp + ny
#endif

           flux(iu1+nghost,j,k,1,1) = flux(iu1+nghost,j,k,1,1) &
                                    + (one - eps)*outer_flux(jremap,k) &
                                    + eps*outer_flux(jremapp,k)
           flux(iu1+nghost,j,k,1,1) = half*flux(iu1+nghost,j,k,1,1)
#if MPI == 1
        endif
#endif

        ! Outer boundary conditions
#if MPI == 1
        if (xposition == nxslice - 1) then
           jremap = j + ny*yposition + jp 
#else
           jremap = j + jp
#endif
           jremapp = jremap + 1
           eps     = epsi/dy

#if MPI == 1
           if (jremap >= ny*nyslice + 1) jremap = jremap - ny*nyslice
           if (jremapp >= ny*nyslice + 1) jremapp = jremapp - ny*nyslice
#else
           if (jremap >= ny + 1) jremap = jremap - ny
           if (jremapp >= ny + 1) jremapp = jremapp - ny
#endif

           flux(iu2-nghost+1,j,k,1,1) = flux(iu2-nghost+1,j,k,1,1) &
                                      + (one - eps)*inner_flux(jremap,k) &
                                      + eps*inner_flux(jremapp,k)
           flux(iu2-nghost+1,j,k,1,1) = half*flux(iu2-nghost+1,j,k,1,1)
#if MPI == 1
        endif
#endif
     enddo
  enddo
  !$OMP END PARALLEL DO

  deallocate(inner_flux, outer_flux)

  call shearing_emf_old(itime)

  return
end subroutine shearing_flux_old
!===============================================================================
!> Compute electromotive forces at shearing boundaries - old version
!===============================================================================
subroutine shearing_emf_old(itime)
  use mpi_var
  use params
  use variables
  implicit none

  real(dp) :: itime
  real(dp), dimension(:,:), allocatable :: inner_emf, outer_emf
  real(dp) :: deltay, epsi
  real(dp) :: eps
  integer  :: jp, jremap, jremapp
  integer  :: j, k

  if (verbose) print*, '> Entering shearing_emf'

  deltay = 1.5d0*Omega0*(xmax - xmin)*itime
  deltay = mod(deltay, (ymax - ymin))
  epsi   = mod(deltay, dy)
  jp = int(deltay/dy)

#if MPI == 1
  allocate(inner_emf(ju1:ju1+ny*nyslice+2*nghost-1,ku1:ku2))
  allocate(outer_emf(ju1:ju1+ny*nyslice+2*nghost-1,ku1:ku2))
  call get_shearing_emf_old(inner_emf, outer_emf, deltay)
#else
  allocate(inner_emf(ju1:ju2,ku1:ku2))
  allocate(outer_emf(ju1:ju2,ku1:ku2))
  inner_emf(ju1:ju2,ku1:ku2) = emfy(iu1+nghost,ju1:ju2,ku1:ku2)
  outer_emf(ju1:ju2,ku1:ku2) = emfy(iu2-nghost+1,ju1:ju2,ku1:ku2)
#endif
  
  !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(jremap, jremapp, eps)
  do k = ku1, ku2
     do j = ju1, ju2
        ! Inner boundary conditions
#if MPI == 1
        if (xposition == 0) then
           jremap = j + ny*yposition - jp - 1
#else
           jremap = j - jp - 1
#endif
           jremapp = jremap + 1
           eps     = one - epsi/dy

#if MPI == 1
           if (jremap <= 0) jremap = jremap + ny*nyslice
           if (jremapp <= 0) jremapp = jremapp + ny*nyslice
#else
           if (jremap <= 0) jremap = jremap + ny
           if (jremapp <= 0) jremapp = jremapp + ny
#endif

           emfy(iu1+nghost,j,k) = emfy(iu1+nghost,j,k) &
                               + (one - eps)*outer_emf(jremap,k) &
                               + eps*outer_emf(jremapp,k)
           emfy(iu1+nghost,j,k) = half*emfy(iu1+nghost,j,k)
#if MPI == 1
        endif
#endif

        ! Outer boundary conditions
#if MPI == 1
        if (xposition == nxslice - 1) then
           jremap = j + ny*yposition + jp 
#else
           jremap = j + jp
#endif
           jremapp = jremap + 1
           eps     = epsi/dy

#if MPI == 1
           if (jremap >= ny*nyslice + 1) jremap = jremap - ny*nyslice
           if (jremapp >= ny*nyslice + 1) jremapp = jremapp - ny*nyslice
#else
           if (jremap >= ny + 1) jremap = jremap - ny
           if (jremapp >= ny + 1) jremapp = jremapp - ny
#endif

           emfy(iu2-nghost+1,j,k) = emfy(iu2-nghost+1,j,k) &
                                  + (one - eps)*inner_emf(jremap,k) &
                                  + eps*inner_emf(jremapp,k)
           emfy(iu2-nghost+1,j,k) = half*emfy(iu2-nghost+1,j,k)
#if MPI == 1
        endif
#endif
     enddo
  enddo
  !$OMP END PARALLEL DO

  deallocate(inner_emf, outer_emf)

  return
end subroutine shearing_emf_old
!===============================================================================
!> Compute slopes at shearing boundaries - old version
!===============================================================================
subroutine shearing_slope_old(quant, dq)
  use params
  implicit none

#if MPI == 1
  real(dp), dimension(ju1:ju1+ny*nyslice+2*nghost-1,ku1:ku2,nvar+3,nghost) :: quant, dq
#else
  real(dp), dimension(ju1:ju2,ku1:ku2,nvar+3,nghost) :: quant, dq
#endif
  integer :: jslope_type
  integer :: ighost, ivar, j, k
  real(dp) :: dsgn, dlim, dcen, dlft, drgt, slop

  jslope_type = min(slope_type, 2)
  if (jslope_type == 0) then
     dq = zero
     return
  endif

#if NDIM == 3
  if (jslope_type == 1 .or. jslope_type == 2) then
     do ighost = 1, nghost
        do ivar = 1, nvar-2
           !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dsgn, dlim, dcen, dlft) &
           !$OMP PRIVATE(drgt, slop)
           do k = ku1, ku2
              do j = ju1 + 1, ju1 + ny*nyslice + 2*nghost - 2
                 dlft = jslope_type*(quant(j,k,ivar,ighost) &
                                   - quant(j-1,k,ivar,ighost))
                 drgt = jslope_type*(quant(j+1,k,ivar,ighost) &
                                   - quant(j,k,ivar,ighost))
                 dcen = half*(dlft + drgt)/jslope_type
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft), abs(drgt))
                 dlim = slop
                 if ((dlft*drgt) <= zero) dlim = zero
                 dq(j,k,ivar,ighost) = dsgn*min(dlim, abs(dcen))
              enddo
           enddo
           !$OMP END PARALLEL DO
        enddo
        !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dsgn, dlim, dcen, dlft) &
        !$OMP PRIVATE(drgt, slop)
        do k = ku1, ku2
           do j = ju1 + 1, ju1 + ny*nyslice + 2*nghost - 2
              dlft = jslope_type*(quant(j,k,nvar,ighost) &
                                - quant(j-1,k,nvar,ighost))
              drgt = jslope_type*(quant(j+1,k,nvar,ighost) &
                                - quant(j,k,nvar,ighost))
              dcen = half*(dlft + drgt)/jslope_type
              dsgn = sign(one, dcen)
              slop = min(abs(dlft), abs(drgt))
              dlim = slop
              if ((dlft*drgt) <= zero) dlim = zero
              dq(j,k,nvar,ighost) = dsgn*min(dlim, abs(dcen))
           enddo
        enddo
        !$OMP END PARALLEL DO
        do ivar = nvar+1, nvar+3, 2
           !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(dsgn, dlim, dcen, dlft) &
           !$OMP PRIVATE(drgt, slop)
           do k = ku1, ku2
              do j = ju1 + 1, ju1 + ny*nyslice + 2*nghost - 2
                 dlft = jslope_type*(quant(j,k,ivar,ighost) &
                                   - quant(j-1,k,ivar,ighost))
                 drgt = jslope_type*(quant(j+1,k,ivar,ighost) &
                                   - quant(j,k,ivar,ighost))
                 dcen = half*(dlft + drgt)/jslope_type
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft), abs(drgt))
                 dlim = slop
                 if ((dlft*drgt) <= zero) dlim = zero
                 dq(j,k,ivar,ighost) = dsgn*min(dlim, abs(dcen))
              enddo
           enddo
           !$OMP END PARALLEL DO
        enddo
     enddo
  else
     write(*,*) "Unknown slope type!"
     stop
  endif

  do ighost = 1, nghost
     !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     do k = ku1, ku2
        do j = ju1+1, ju1 + ny*nyslice + 2*nghost - 2
           dq(j,k,nvar-1,ighost) = quant(j+1,k,nvar-1,ighost) &
                                   - quant(j,k,nvar-1,ighost)
           dq(j,k,nvar+2,ighost) = quant(j+1,k,nvar+2,ighost) &
                                   - quant(j,k,nvar+2,ighost)
        enddo
     enddo
     !$OMP END PARALLEL DO
  enddo
#endif

  return
end subroutine shearing_slope_old
