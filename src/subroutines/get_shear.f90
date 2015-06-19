!===============================================================================
!> \file get_shear.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is part of the shearing box subroutines.
!! \details
!! Contains get_shearing_quant(), get_shearing_flux(), get_shearing_emf(),
!! get_shearing_quant_old(), get_shearing_flux_old(), get_shearing_emf_old()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Romain Teyssier, 
!! Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          07-24-2014
!! \b last \b modified: 05-11-2015
!<
!===============================================================================
!> Get shearing quantities
!===============================================================================
subroutine get_shear_quant(inner, outer, itime)
#if MPI == 1
  use params
  use variables
  use mpi_var
  use mpi
  implicit none

  real(dp), dimension(nghost,ju1:ju2,ku1:ku2,nvar+3), intent(out) :: inner
  real(dp), dimension(nghost,ju1:ju2,ku1:ku2,nvar+3), intent(out) :: outer
  real(dp), intent(in) :: itime
  real(dp) :: Lx, Ly, yshear, deltay, epsi
  integer  :: jshft
  integer  :: npes_shift, ncell_shift
  integer  :: rank_send1, rank_send2, rank_recv1, rank_recv2
  integer :: size
  real(dp), allocatable, dimension(:,:,:,:) :: slbound, srbound
  real(dp), allocatable, dimension(:,:,:,:) :: rlbound, rrbound
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer :: ierr

  Lx = xmax - xmin
  Ly = ymax - ymin
  yshear = 1.5d0*Omega0*Lx*itime
  deltay = mod(yshear, Ly)
  epsi   = mod(deltay, dy)
  jshft  = int(deltay/dy)

  allocate(slbound(nghost,ju1:ju2,ku1:ku2,1:nvar+3))
  allocate(srbound(nghost,ju1:ju2,ku1:ku2,1:nvar+3))
  allocate(rlbound(nghost,ju1:ju2,ku1:ku2,1:nvar+3))
  allocate(rrbound(nghost,ju1:ju2,ku1:ku2,1:nvar+3))

  ! First, perform a periodic boundary conditions-like communication
  size = (ju2 - ju1 + 1)*(ku2 - ku1 + 1)*(nvar + 3)*nghost
  !$OMP PARALLEL WORKSHARE
  slbound(:,:,:,:) = uin(iu1+nghost:iu1+2*nghost-1,:,:,:)
  srbound(:,:,:,:) = uin(iu2-2*nghost+1:iu2-nghost,:,:,:)
  !$OMP END PARALLEL WORKSHARE
  
  call MPI_Sendrecv(slbound, size, MPI_DOUBLE_PRECISION, xleft, 10 &
       , rlbound, size, MPI_DOUBLE_PRECISION, xright, 10 &
       , MPI_COMM_WORLD, status, ierr)
  
  call MPI_Sendrecv(srbound, size, MPI_DOUBLE_PRECISION, xright, 11 &
          , rrbound, size, MPI_DOUBLE_PRECISION, xleft, 11 &
          , MPI_COMM_WORLD, status, ierr)
  
  ! Second, shift boundaries
  ! Outer boundary
  if (xposition == 0) then
     if (jshft <= ny) then
        outer(:,ju1+jshft:ju2,:,:) = rrbound(:,ju1:ju2-jshft,:,:)
        rank_send1 = yright
        rank_recv1 = yleft

        size = nghost*jshft*(ku2-ku1+1)*(nvar+3)
        call MPI_Sendrecv(rrbound(:,ju2-jshft+1-2*nghost:ju2-2*nghost,:,:) &
             , size, MPI_DOUBLE_PRECISION, rank_send1, 10 &
             , outer(:,ju1:ju1+jshft-1,:,:), size, MPI_DOUBLE_PRECISION &
             , rank_recv1, 10, MPI_COMM_WORLD, status, ierr)
     else
        npes_shift  = jshft/ny
        ncell_shift = mod(jshft, ny)
        rank_send1 = mype + (npes_shift+1)*nxslice*nzslice
        rank_recv1 = mype - (npes_shift+1)*nxslice*nzslice
        rank_send2 = mype + (npes_shift)*nxslice*nzslice
        rank_recv2 = mype - (npes_shift)*nxslice*nzslice
        if (rank_send1 >= npes) rank_send1 = rank_send1 - npes
        if (rank_recv1 < 0) rank_recv1 = rank_recv1 + npes
        if (rank_send2 >= npes) rank_send2 = rank_send2 - npes
        if (rank_recv2 < 0) rank_recv2 = rank_recv2 + npes

        size = nghost*ncell_shift*(ku2-ku1+1)*(nvar+3)
        call MPI_Sendrecv(rrbound(:,ju2-ncell_shift+1-2*nghost:ju2-2*nghost,:,:) &
             , size, MPI_DOUBLE_PRECISION, rank_send1, 11 &
             , outer(:,ju1:ju1+ncell_shift-1,:,:), size &
             , MPI_DOUBLE_PRECISION, rank_recv1, 11, MPI_COMM_WORLD, status &
             , ierr)
     
        size = nghost*(ju2-ju1-ncell_shift+1)*(ku2-ku1+1)*(nvar+3)
        call MPI_Sendrecv(rrbound(:,ju1:ju2-ncell_shift,:,:) &
             , size, MPI_DOUBLE_PRECISION, rank_send2, 12 &
             , outer(:,ju1+ncell_shift:ju2,:,:), size &
             , MPI_DOUBLE_PRECISION, rank_recv2, 12, MPI_COMM_WORLD, status &
             , ierr)
     endif
  endif

  ! Inner boundary
  if (xposition == nxslice - 1) then
     if (jshft <= ny) then
        inner(:,ju1:ju2-jshft,:,:) = rlbound(:,ju1+jshft:ju2,:,:)
        rank_send1 = yleft
        rank_recv1 = yright

        size = nghost*jshft*(ku2-ku1+1)*(nvar+3)
        call MPI_Sendrecv(rlbound(:,ju1+2*nghost:ju1+jshft-1+2*nghost,:,:) &
             , size, MPI_DOUBLE_PRECISION, rank_send1, 10 &
             , inner(:,ju2-jshft+1:ju2,:,:), size, MPI_DOUBLE_PRECISION &
             , rank_recv1, 10, MPI_COMM_WORLD, status, ierr)
     else
        npes_shift  = jshft/ny
        ncell_shift = mod(jshft, ny)
        rank_send1 = mype - (npes_shift+1)*nxslice*nzslice
        rank_recv1 = mype + (npes_shift+1)*nxslice*nzslice
        rank_send2 = mype - (npes_shift)*nxslice*nzslice
        rank_recv2 = mype + (npes_shift)*nxslice*nzslice
        if (rank_send1 < 0) rank_send1 = rank_send1 + npes
        if (rank_recv1 >= npes) rank_recv1 = rank_recv1 - npes
        if (rank_send2 < 0) rank_send2 = rank_send2 + npes
        if (rank_recv2 >= npes) rank_recv2 = rank_recv2 - npes

        size = nghost*ncell_shift*(ku2-ku1+1)*(nvar+3)
        call MPI_Sendrecv(rlbound(:,ju1+2*nghost:ju1+ncell_shift-1+2*nghost,:,:) &
             , size, MPI_DOUBLE_PRECISION, rank_send1, 11 &
             , inner(:,ju2-ncell_shift+1:ju2,:,:), size &
             , MPI_DOUBLE_PRECISION, rank_recv1, 11, MPI_COMM_WORLD, status &
             , ierr)
        
        size = nghost*(ju2-ju1-ncell_shift+1)*(ku2-ku1+1)*(nvar+3)
        call MPI_Sendrecv(rlbound(:,ju1+ncell_shift:ju2,:,:) &
             , size, MPI_DOUBLE_PRECISION, rank_send2, 12 &
             , inner(:,ju1:ju2-ncell_shift,:,:), size &
             , MPI_DOUBLE_PRECISION, rank_recv2, 12, MPI_COMM_WORLD, status &
             , ierr)
     endif
  endif

  deallocate(slbound, srbound, rlbound, rrbound)
#endif

  return
end subroutine get_shear_quant
!===============================================================================
!> Get shearing flux
!===============================================================================
subroutine get_shear_flux(inner, outer, itime)
#if MPI == 1
  use params
  use variables
  use mpi_var
  use mpi
  implicit none

  real(dp), dimension(jf1-1:jf2,kf1:kf2), intent(out) :: outer
  real(dp), dimension(jf1:jf2+1,kf1:kf2), intent(out) :: inner
  real(dp), intent(in) :: itime
  real(dp) :: Lx, Ly, yshear, deltay, epsi
  integer  :: jshft
  integer  :: npes_shift, ncell_shift
  integer  :: rank_send1, rank_send2, rank_send3
  integer  :: rank_recv1, rank_recv2, rank_recv3
  integer :: size
  real(dp), allocatable, dimension(:,:) :: slbound, srbound
  real(dp), allocatable, dimension(:,:) :: rlbound, rrbound
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer :: ierr

  Lx = xmax - xmin
  Ly = ymax - ymin
  yshear = 1.5d0*Omega0*Lx*itime
  deltay = mod(yshear, Ly)
  epsi   = mod(deltay, dy)
  jshft  = int(deltay/dy)

  allocate(slbound(jf1:jf2,kf1:kf2))
  allocate(srbound(jf1:jf2,kf1:kf2))
  allocate(rlbound(jf1:jf2,kf1:kf2))
  allocate(rrbound(jf1:jf2,kf1:kf2))

  ! First, perform a periodic boundary conditions-like communication
  size = (jf2 - jf1 + 1)*(kf2 - kf1 + 1)
  !$OMP PARALLEL WORKSHARE
  slbound(:,:) = flux(if1,:,:,1,1)
  srbound(:,:) = flux(if2,:,:,1,1)
  !$OMP END PARALLEL WORKSHARE
  
  call MPI_Sendrecv(slbound, size, MPI_DOUBLE_PRECISION, xleft, 10 &
       , rlbound, size, MPI_DOUBLE_PRECISION, xright, 10 &
       , MPI_COMM_WORLD, status, ierr)
  
  call MPI_Sendrecv(srbound, size, MPI_DOUBLE_PRECISION, xright, 11 &
       , rrbound, size, MPI_DOUBLE_PRECISION, xleft, 11 &
       , MPI_COMM_WORLD, status, ierr)

  ! Second, shift boundaries
  ! Outer boundary
  if (xposition == 0) then
     if (jshft <= ny-1) then
        rank_send1 = yright
        rank_recv1 = yleft
        rank_send2 = yleft
        rank_recv2 = yright

        if (jshft == 0) then
           outer(jf1:jf2-1,:) = rrbound(jf1:jf2-1,:)
        
           size = kf2 - kf1 + 1
           call MPI_Sendrecv(rrbound(jf2-1,:), size, MPI_DOUBLE_PRECISION &
                & , rank_send1, 10, outer(jf1-1,:), size, MPI_DOUBLE_PRECISION &
                & , rank_recv1, 10, MPI_COMM_WORLD, status, ierr)

           size = kf2 - kf1 + 1
           call MPI_Sendrecv(rrbound(jf1,:), size, MPI_DOUBLE_PRECISION &
                & , rank_send2, 11, outer(jf2,:), size, MPI_DOUBLE_PRECISION &
                & , rank_recv2, 11, MPI_COMM_WORLD, status, ierr)
        else
           outer(jf1+jshft:jf2,:) = rrbound(jf1:jf2-jshft,:)
           
           size = (jshft + 1)*(kf2 - kf1 + 1)
           call MPI_Sendrecv(rrbound(jf2-jshft-1:jf2-1,:), size &
                & , MPI_DOUBLE_PRECISION, rank_send1, 10 &
                & , outer(jf1-1:jf1+jshft-1,:), size, MPI_DOUBLE_PRECISION &
                & , rank_recv1, 10, MPI_COMM_WORLD, status, ierr)
        endif
     else
        npes_shift  = jshft/ny
        ncell_shift = mod(jshft, ny)
        rank_send1 = mype + (npes_shift+1)*nxslice*nzslice
        rank_recv1 = mype - (npes_shift+1)*nxslice*nzslice
        rank_send2 = mype + (npes_shift)*nxslice*nzslice
        rank_recv2 = mype - (npes_shift)*nxslice*nzslice
        rank_send3 = mype + (npes_shift-1)*nxslice*nzslice
        rank_recv3 = mype - (npes_shift-1)*nxslice*nzslice
        if (rank_send1 >= npes) rank_send1 = rank_send1 - npes
        if (rank_recv1 < 0) rank_recv1 = rank_recv1 + npes
        if (rank_send2 >= npes) rank_send2 = rank_send2 - npes
        if (rank_recv2 < 0) rank_recv2 = rank_recv2 + npes
        if (rank_send3 >= npes) rank_send3 = rank_send3 - npes
        if (rank_recv3 < 0) rank_recv3 = rank_recv3 + npes

        if (ncell_shift == 0) then
           size = kf2 - kf1 + 1
           call MPI_Sendrecv(rrbound(jf2-1,:), size, MPI_DOUBLE_PRECISION &
                & , rank_send1, 10, outer(jf1-1,:), size, MPI_DOUBLE_PRECISION &
                & , rank_recv1, 10, MPI_COMM_WORLD, status, ierr)

           size = (jf2 - jf1)*(kf2 - kf1 + 1)
           call MPI_Sendrecv(rrbound(jf1:jf2-1,:), size, MPI_DOUBLE_PRECISION &
                & , rank_send2, 11, outer(jf1:jf2-1,:), size &
                & , MPI_DOUBLE_PRECISION, rank_recv2, 11, MPI_COMM_WORLD &
                & , status, ierr)

           size = kf2 - kf1 + 1
           call MPI_Sendrecv(rrbound(jf1,:), size, MPI_DOUBLE_PRECISION &
                & , rank_send3, 12, outer(jf2,:), size, MPI_DOUBLE_PRECISION &
                & , rank_recv3, 12, MPI_COMM_WORLD, status, ierr)
        else
           size = (ncell_shift + 1)*(kf2 - kf1 + 1)
           call MPI_Sendrecv(rrbound(jf2-ncell_shift-1:jf2-1,:), size &
                & , MPI_DOUBLE_PRECISION, rank_send1, 10 &
                & , outer(jf1-1:jf1+ncell_shift-1,:), size &
                & , MPI_DOUBLE_PRECISION, rank_recv1, 10, MPI_COMM_WORLD &
                & , status, ierr)

           size = (jf2 - jf1 + 1 - ncell_shift)*(kf2 - kf1 + 1)
           call MPI_Sendrecv(rrbound(jf1:jf2-ncell_shift,:), size &
                & , MPI_DOUBLE_PRECISION, rank_send2, 11 &
                & , outer(jf1+ncell_shift:jf2,:), size, MPI_DOUBLE_PRECISION &
                & , rank_recv2, 11, MPI_COMM_WORLD, status, ierr)
        endif
     endif
  endif

  ! Inner boundary
  if (xposition == nxslice - 1) then
     if (jshft <= ny-1) then
        rank_send1 = yleft
        rank_recv1 = yright

        inner(jf1:jf2-1-jshft,:) = rlbound(jf1+jshft:jf2-1,:)

        size = (jshft + 2)*(kf2 - kf1 + 1)
        call MPI_Sendrecv(rlbound(jf1:jf1+1+jshft,:), size &
             & , MPI_DOUBLE_PRECISION, rank_send1, 10 &
             & , inner(jf2-jshft:jf2+1,:), size, MPI_DOUBLE_PRECISION &
             & , rank_recv1, 10, MPI_COMM_WORLD, status, ierr)
     else
        npes_shift  = jshft/ny
        ncell_shift = mod(jshft, ny)
        rank_send1 = mype - (npes_shift+1)*nxslice*nzslice
        rank_recv1 = mype + (npes_shift+1)*nxslice*nzslice
        rank_send2 = mype - (npes_shift)*nxslice*nzslice
        rank_recv2 = mype + (npes_shift)*nxslice*nzslice
        if (rank_send1 < 0) rank_send1 = rank_send1 + npes
        if (rank_recv1 >= npes) rank_recv1 = rank_recv1 - npes
        if (rank_send2 < 0) rank_send2 = rank_send2 + npes
        if (rank_recv2 >= npes) rank_recv2 = rank_recv2 - npes

        size = (ncell_shift + 2)*(kf2 - kf1 + 1)
        call MPI_Sendrecv(rlbound(jf1:jf1+1+ncell_shift,:) &
             , size, MPI_DOUBLE_PRECISION, rank_send1, 11 &
             , inner(jf2-ncell_shift:jf2+1,:), size &
             , MPI_DOUBLE_PRECISION, rank_recv1, 11, MPI_COMM_WORLD, status &
             , ierr)
        
        size = (jf2 - jf1 - ncell_shift)*(kf2 - kf1 + 1)
        call MPI_Sendrecv(rlbound(jf1+ncell_shift:jf2-1,:) &
             , size, MPI_DOUBLE_PRECISION, rank_send2, 12 &
             , inner(jf1:jf2-1-ncell_shift,:), size &
             , MPI_DOUBLE_PRECISION, rank_recv2, 12, MPI_COMM_WORLD, status &
             , ierr)
     endif
  endif

  deallocate(slbound, srbound, rlbound, rrbound)
#endif

  return
end subroutine get_shear_flux
!===============================================================================
!> Get shearing EMF
!===============================================================================
subroutine get_shear_emf(inner, outer, itime)
#if MPI == 1
  use params
  use variables
  use mpi_var
  use mpi
  implicit none

  real(dp), dimension(ju1:ju2,ku1:ku2), intent(out) :: outer
  real(dp), dimension(ju1:ju2,ku1:ku2), intent(out) :: inner
  real(dp), intent(in) :: itime
  real(dp) :: Lx, Ly, yshear, deltay, epsi
  integer  :: jshft
  integer  :: npes_shift, ncell_shift
  integer  :: rank_send1, rank_send2, rank_send3
  integer  :: rank_recv1, rank_recv2, rank_recv3
  integer  :: size
  real(dp), allocatable, dimension(:,:) :: slbound, srbound
  real(dp), allocatable, dimension(:,:) :: rlbound, rrbound
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer :: ierr

  Lx = xmax - xmin
  Ly = ymax - ymin
  yshear = 1.5d0*Omega0*Lx*itime
  deltay = mod(yshear, Ly)
  epsi   = mod(deltay, dy)
  jshft  = int(deltay/dy)

  allocate(slbound(ju1:ju2,ku1:ku2))
  allocate(srbound(ju1:ju2,ku1:ku2))
  allocate(rlbound(ju1:ju2,ku1:ku2))
  allocate(rrbound(ju1:ju2,ku1:ku2))

  ! First, perform a periodic boundary conditions-like communication
  size = (ju2 - ju1 + 1)*(ku2 - ku1 + 1)
  !$OMP PARALLEL WORKSHARE
  slbound(:,:) = emfy(iu1+nghost,:,:)
  srbound(:,:) = emfy(iu2-nghost+1,:,:)
  !$OMP END PARALLEL WORKSHARE
  
  call MPI_Sendrecv(slbound, size, MPI_DOUBLE_PRECISION, xleft, 10 &
       , rlbound, size, MPI_DOUBLE_PRECISION, xright, 10 &
       , MPI_COMM_WORLD, status, ierr)
  
  call MPI_Sendrecv(srbound, size, MPI_DOUBLE_PRECISION, xright, 11 &
       , rrbound, size, MPI_DOUBLE_PRECISION, xleft, 11 &
       , MPI_COMM_WORLD, status, ierr)

  ! Second, shift boundaries
  ! Outer boundary
  if (xposition == 0) then
     if (jshft <= ny-1) then
        outer(ju1+nghost+jshft:ju2-nghost,:) &
             = rrbound(ju1+nghost:ju2-nghost-jshft,:)
        rank_send1 = yright
        rank_recv1 = yleft
        rank_send2 = yleft
        rank_recv2 = yright

        size = (jshft+nghost)*(ku2-ku1+1)
        call MPI_Sendrecv(rrbound(ju2-jshft-2*nghost+1:ju2-nghost,:) &
             , size, MPI_DOUBLE_PRECISION, rank_send1, 10 &
             , outer(ju1:ju1+nghost+jshft-1,:), size, MPI_DOUBLE_PRECISION &
             , rank_recv1, 10, MPI_COMM_WORLD, status, ierr)

        size = nghost*(ku2-ku1+1)
        call MPI_Sendrecv(rrbound(ju1+nghost+jshft:ju1+2*nghost-1+jshft,:) &
             , size, MPI_DOUBLE_PRECISION, rank_send2, 11 &
             , outer(ju2-nghost+1:ju2,:), size, MPI_DOUBLE_PRECISION &
             , rank_recv2, 11, MPI_COMM_WORLD, status, ierr)
     else
        npes_shift  = jshft/ny
        ncell_shift = mod(jshft, ny)
        rank_send1 = mype + (npes_shift+1)*nxslice*nzslice
        rank_recv1 = mype - (npes_shift+1)*nxslice*nzslice
        rank_send2 = mype + (npes_shift)*nxslice*nzslice
        rank_recv2 = mype - (npes_shift)*nxslice*nzslice
        rank_send3 = mype + (npes_shift-1)*nxslice*nzslice
        rank_recv3 = mype - (npes_shift-1)*nxslice*nzslice
        if (rank_send1 >= npes) rank_send1 = rank_send1 - npes
        if (rank_recv1 < 0) rank_recv1 = rank_recv1 + npes
        if (rank_send2 >= npes) rank_send2 = rank_send2 - npes
        if (rank_recv2 < 0) rank_recv2 = rank_recv2 + npes
        if (rank_send3 >= npes) rank_send3 = rank_send3 - npes
        if (rank_recv3 < 0) rank_recv3 = rank_recv3 + npes

        if (ncell_shift < nghost) then
           size = (nghost+ncell_shift)*(ku2-ku1+1)
           call MPI_Sendrecv(rrbound(ju1+nghost:ju1+2*nghost-1+ncell_shift,:) &
                , size, MPI_DOUBLE_PRECISION, rank_send1, 10 &
                , outer(ju2-nghost-ncell_shift+1:ju2,:), size &
                , MPI_DOUBLE_PRECISION, rank_recv1, 10, MPI_COMM_WORLD, status &
                , ierr)

           size = (ju2-ju1-2*nghost+1)*(ku2-ku1+1)
           call MPI_Sendrecv(rrbound(ju1+nghost:ju2-nghost,:), size &
                , MPI_DOUBLE_PRECISION, rank_send2, 11 &
                , outer(ju1+nghost-ncell_shift:ju2-nghost-ncell_shift,:), size &
                , MPI_DOUBLE_PRECISION, rank_recv2, 11, MPI_COMM_WORLD, status &
                , ierr)

           size = (nghost-ncell_shift)*(ku2-ku1+1)
           call MPI_Sendrecv(rrbound(ju2-2*nghost+ncell_shift+1:ju2-nghost,:) &
                , size, MPI_DOUBLE_PRECISION, rank_send3, 12 &
                , outer(ju1:ju1+nghost-ncell_shift-1,:), size &
                , MPI_DOUBLE_PRECISION, rank_recv3, 12, MPI_COMM_WORLD, status &
                , ierr)
        else
           size = (nghost+ncell_shift)*(ku2-ku1+1)
           call MPI_Sendrecv(rrbound(ju1+nghost:ju1+2*nghost-1+ncell_shift,:) &
                , size, MPI_DOUBLE_PRECISION, rank_send1, 10 &
                , outer(ju2-nghost-ncell_shift+1:ju2,:), size &
                , MPI_DOUBLE_PRECISION, rank_recv1, 10, MPI_COMM_WORLD, status &
                , ierr)

           size = (ju2-ju1-nghost-ncell_shift+1)*(ku2-ku1+1)
           call MPI_Sendrecv(rrbound(ju1+nghost:ju2-ncell_shift,:), size &
                , MPI_DOUBLE_PRECISION, rank_send2, 11 &
                , outer(ju1:ju2-nghost-ncell_shift,:), size &
                , MPI_DOUBLE_PRECISION, rank_recv2, 11, MPI_COMM_WORLD, status &
                , ierr)
        endif
     endif
  endif

  ! Inner boundary
  if (xposition == nxslice - 1) then
     if (jshft <= ny-1) then
        inner(ju1+nghost:ju2-nghost-jshft,:) &
             = rlbound(ju1+nghost+jshft:ju2-nghost,:)
        rank_send1 = yleft
        rank_recv1 = yright
        rank_send2 = yright
        rank_recv2 = yleft

        size = (jshft+nghost)*(ku2-ku1+1)
        call MPI_Sendrecv(rlbound(ju1+nghost:ju1+jshft+2*nghost-1,:) &
             , size, MPI_DOUBLE_PRECISION, rank_send1, 10 &
             , inner(ju2-jshft-nghost+1:ju2,:), size, MPI_DOUBLE_PRECISION &
             , rank_recv1, 10, MPI_COMM_WORLD, status, ierr)

        size = nghost*(ku2-ku1+1)
        call MPI_Sendrecv(rlbound(ju2-2*nghost+1-jshft:ju2-nghost-jshft,:) &
             , size, MPI_DOUBLE_PRECISION, rank_send2, 11 &
             , inner(ju1:ju1+nghost-1,:), size, MPI_DOUBLE_PRECISION &
             , rank_recv2, 11, MPI_COMM_WORLD, status, ierr)
     else
        npes_shift  = jshft/ny
        ncell_shift = mod(jshft, ny)
        rank_send1 = mype + (npes_shift-1)*nxslice*nzslice
        rank_recv1 = mype - (npes_shift-1)*nxslice*nzslice
        rank_send2 = mype + (npes_shift)*nxslice*nzslice
        rank_recv2 = mype - (npes_shift)*nxslice*nzslice
        rank_send3 = mype + (npes_shift+1)*nxslice*nzslice
        rank_recv3 = mype - (npes_shift+1)*nxslice*nzslice
        if (rank_send1 >= npes) rank_send1 = rank_send1 - npes
        if (rank_recv1 < 0) rank_recv1 = rank_recv1 + npes
        if (rank_send2 >= npes) rank_send2 = rank_send2 - npes
        if (rank_recv2 < 0) rank_recv2 = rank_recv2 + npes
        if (rank_send3 >= npes) rank_send3 = rank_send3 - npes
        if (rank_recv3 < 0) rank_recv3 = rank_recv3 + npes

        if (ncell_shift < nghost) then
           size = (nghost+ncell_shift)*(ku2-ku1+1)
           call MPI_Sendrecv(rlbound(ju2-2*nghost-ncell_shift+1:ju2-nghost,:) &
                , size, MPI_DOUBLE_PRECISION, rank_send1, 10 &
                , inner(ju1:ju1+nghost+ncell_shift-1,:), size &
                , MPI_DOUBLE_PRECISION, rank_recv1, 10, MPI_COMM_WORLD, status &
                , ierr)

           size = (ju2-ju1-2*nghost+1)*(ku2-ku1+1)
           call MPI_Sendrecv(rlbound(ju1+nghost:ju2-nghost,:), size &
                , MPI_DOUBLE_PRECISION, rank_send2, 11 &
                , inner(ju1+nghost+ncell_shift:ju2-nghost+ncell_shift,:), size &
                , MPI_DOUBLE_PRECISION, rank_recv2, 11, MPI_COMM_WORLD, status &
                , ierr)

           size = (nghost-ncell_shift)*(ku2-ku1+1)
           call MPI_Sendrecv(rlbound(ju1+nghost:ju1+2*nghost-ncell_shift-1,:) &
                , size, MPI_DOUBLE_PRECISION, rank_send3, 12 &
                , inner(ju2-nghost+ncell_shift+1:ju2,:), size &
                , MPI_DOUBLE_PRECISION, rank_recv3, 12, MPI_COMM_WORLD, status &
                , ierr)
        else
           size = (ju2-ju1-ncell_shift-nghost+1)*(ku2-ku1+1)
           call MPI_Sendrecv(rlbound(ju1+nghost:ju2-ncell_shift,:), size &
                , MPI_DOUBLE_PRECISION, rank_send2, 10 &
                , inner(ju1+nghost+ncell_shift:ju2,:), size &
                , MPI_DOUBLE_PRECISION, rank_recv2, 10, MPI_COMM_WORLD, status &
                , ierr)

           size = (nghost+ncell_shift)*(ku2-ku1+1)
           call MPI_Sendrecv(rlbound(ju2-2*nghost-ncell_shift+1:ju2-nghost,:) &
                , size, MPI_DOUBLE_PRECISION, rank_send1, 11 &
                , inner(ju1:ju1+nghost+ncell_shift-1,:), size &
                , MPI_DOUBLE_PRECISION, rank_recv1, 11, MPI_COMM_WORLD, status &
                , ierr)
        endif
     endif
  endif

  deallocate(slbound, srbound, rlbound, rrbound)
#endif

  return
end subroutine get_shear_emf
!===============================================================================
!> Get shearing quantities - old version
!===============================================================================
# if MPI == 1
subroutine get_shearing_quant_old(inner_quant, outer_quant, deltay)
  use mpi
  use params
  use variables
  use mpi_var
  implicit none

  real(dp), dimension(ju1:ju1+ny*nyslice+2*nghost-1,ku1:ku2,nvar+3,nghost) :: inner_quant, outer_quant
  real(dp) :: deltay, y1

  real(dp), dimension(:,:,:,:,:), allocatable :: sent_quant
  real(dp), dimension(:,:,:,:,:), allocatable :: recv_quant1, recv_quant2
  integer :: jksize
  integer :: ypos_comm, ypos_commp, jstart, jstartp
  integer :: send_cpu1, send_cpu2, recv_cpu1, recv_cpu2
  integer :: send_cpu1p, send_cpu2p, recv_cpu1p, recv_cpu2p
  integer :: send_cpu1m, send_cpu2m, recv_cpu1m, recv_cpu2m
  integer :: ighost
  integer :: nbound=2
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer :: ierr

  jksize = (ju2 - ju1 + 1)*(ku2 - ku1 + 1)*(nvar + 3)*nghost*nbound
  allocate(sent_quant(ju1:ju2,ku1:ku2,nvar+3,nghost,nbound))
  allocate(recv_quant1(ju1:ju2,ku1:ku2,nvar+3,nghost,nbound))
  allocate(recv_quant2(ju1:ju2,ku1:ku2,nvar+3,nghost,nbound))

  !$OMP PARALLEL WORKSHARE
  inner_quant(:,:,:,:) = 1000.d0
  outer_quant(:,:,:,:) = 1000.d0
  !$OMP END PARALLEL WORKSHARE

  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do ighost = 1, nghost
     sent_quant(ju1:ju2,ku1:ku2,:,ighost,1) &
          = uin(iu1-1+nghost+ighost,ju1:ju2,ku1:ku2,:)
     sent_quant(ju1:ju2,ku1:ku2,:,ighost,2) &
          = uin(iu2-2*nghost+ighost,ju1:ju2,ku1:ku2,:)
  enddo
  !$OMP END PARALLEL DO

  if ((xposition == 0) .and. (nxslice > 1)) then
     if (deltay == 0) then
        ypos_comm = yposition - 1
        if (ypos_comm < 0) ypos_comm = ypos_comm + nyslice
     else
        y1 = y(ju1+nghost) - half*dy - deltay
        if (y1 < ymin) y1 = y1 + (ymax - ymin)
        ypos_comm = floor((y1 - ymin)/(y(ju2-nghost) - y(ju1+nghost) + dy))
     endif
     ypos_commp = ypos_comm + 1
     if (ypos_commp > nyslice - 1) ypos_commp = ypos_commp - nyslice

     send_cpu1 = zposition*nxslice + (nxslice - 1) + nxslice*nzslice*ypos_comm
     send_cpu2 = zposition*nxslice + (nxslice - 1) + nxslice*nzslice*ypos_commp
     recv_cpu1 = send_cpu1
     recv_cpu2 = zposition*nxslice + (nxslice - 1) + nxslice*nzslice*ypos_commp

     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu2, 11 &
          , recv_quant1, jksize, MPI_DOUBLE_PRECISION, recv_cpu1, 10 &
          , MPI_COMM_WORLD, status, ierr)
     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu1, 13 &
          , recv_quant2, jksize, MPI_DOUBLE_PRECISION, recv_cpu2, 12 &
          , MPI_COMM_WORLD, status, ierr)

     jstart  = ju1 + nghost + ny*ypos_comm
     jstartp = ju1 + nghost + ny*ypos_commp

     if (abs(ypos_commp - ypos_comm) > 1) then
        outer_quant(jstart-nghost:jstart+ny+nghost-1,ku1:ku2,:,:) &
             = recv_quant1(ju1:ju2,ku1:ku2,:,:,2)
        outer_quant(jstartp-nghost:jstartp+ny+nghost-1,ku1:ku2,:,:) &
             = recv_quant2(ju1:ju2,ku1:ku2,:,:,2)
     endif

     if (abs(ypos_commp - ypos_comm) == 1) then
        if (ypos_commp > ypos_comm) then
           outer_quant(jstart-nghost:jstart+ny-1,ku1:ku2,:,:) &
                = recv_quant1(ju1:ju2-nghost,ku1:ku2,:,:,2)
           outer_quant(jstartp:jstartp+ny+nghost-1,ku1:ku2,:,:) &
                = recv_quant2(ju1+nghost:ju2,ku1:ku2,:,:,2)
        else
           outer_quant(jstart:jstart+ny+nghost-1,ku1:ku2,:,:) &
                = recv_quant1(ju1+nghost:ju2,ku1:ku2,:,:,2)
           outer_quant(jstartp-nghost:jstartp+ny-1,ku1:ku2,:,:) &
                = recv_quant2(ju1:ju2-nghost,ku1:ku2,:,:,2)
        endif
     endif

     if (abs(ypos_commp - ypos_comm) == 0) then
        outer_quant(ju1:ju2,ku1:ku2,:,:) = recv_quant1(ju1:ju2,ku1:ku2,:,:,2)
     endif
  endif

  if ((xposition == nxslice - 1) .and. (nxslice > 1)) then
     if (deltay == 0) then
        ypos_comm = yposition
     else
        y1 = y(ju1+nghost) - half*dy + deltay
        if (y1 > ymax) y1 = y1 - (ymax - ymin)
        ypos_comm = floor((y1 - ymin)/(y(ju2-nghost) - y(ju1+nghost) + dy))
     endif
     ypos_commp = ypos_comm + 1
     if (ypos_commp > nyslice - 1) ypos_commp = ypos_commp - nyslice

     send_cpu1 = zposition*nxslice + nxslice*nzslice*ypos_comm
     send_cpu2 = zposition*nxslice + nxslice*nzslice*ypos_commp
     recv_cpu1 = send_cpu1
     recv_cpu2 = zposition*nxslice + nxslice*nzslice*ypos_commp
  
     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu2, 10 &
          , recv_quant1, jksize, MPI_DOUBLE_PRECISION, recv_cpu1, 11 &
          , MPI_COMM_WORLD, status, ierr)
     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu1, 12 &
          , recv_quant2, jksize, MPI_DOUBLE_PRECISION, recv_cpu2, 13 &
          , MPI_COMM_WORLD, status, ierr)

     jstart  = ju1 + nghost + ny*ypos_comm
     jstartp = ju1 + nghost + ny*ypos_commp

     if (abs(ypos_commp - ypos_comm) > 1) then
        inner_quant(jstart-nghost:jstart+ny+nghost-1,ku1:ku2,:,:) &
             = recv_quant1(ju1:ju2,ku1:ku2,:,:,1)
        inner_quant(jstartp-nghost:jstartp+ny+nghost-1,ku1:ku2,:,:) &
             = recv_quant2(ju1:ju2,ku1:ku2,:,:,1)
     endif

     if (abs(ypos_commp - ypos_comm) == 1) then
        if (ypos_commp > ypos_comm) then
           inner_quant(jstart-nghost:jstart+ny-1,ku1:ku2,:,:) &
                = recv_quant1(ju1:ju2-nghost,ku1:ku2,:,:,1)
           inner_quant(jstartp:jstartp+ny+nghost-1,ku1:ku2,:,:) &
                = recv_quant2(ju1+nghost:ju2,ku1:ku2,:,:,1)
        else
           inner_quant(jstart:jstart+ny+nghost-1,ku1:ku2,:,:) &
                = recv_quant1(ju1+nghost:ju2,ku1:ku2,:,:,1)
           inner_quant(jstartp-nghost:jstartp+ny-1,ku1:ku2,:,:) &
                = recv_quant2(ju1:ju2-nghost,ku1:ku2,:,:,1)
        endif
     endif
        
     if (abs(ypos_commp - ypos_comm) == 0) then
        inner_quant(ju1:ju2,ku1:ku2,:,:) = recv_quant1(ju1:ju2,ku1:ku2,:,:,1)
     endif
  endif

  if ((nxslice == 1) .and. (nyslice > 1)) then
     if (deltay == 0) then
        ypos_comm = yposition - 1
        if (ypos_comm < 0) ypos_comm = ypos_comm + nyslice
     else
        y1 = y(ju1+nghost) - half*dy - deltay
        if (y1 < ymin) y1 = y1 + (ymax - ymin)
        ypos_comm = floor((y1 - ymin)/(y(ju2-nghost) - y(ju1+nghost) + dy))
     endif
     ypos_commp = ypos_comm + 1
     if (ypos_commp > nyslice - 1) ypos_commp = ypos_commp - nyslice

     send_cpu1p = zposition*nxslice + (nxslice - 1) + nxslice*nzslice*ypos_comm
     send_cpu2p = zposition*nxslice + (nxslice - 1) + nxslice*nzslice*ypos_commp
     recv_cpu1p = send_cpu1p
     recv_cpu2p = zposition*nxslice + (nxslice - 1) + nxslice*nzslice*ypos_commp

     if (deltay == 0) then
        ypos_comm = yposition
     else
        y1 = y(ju1+nghost) - half*dy + deltay
        if (y1 > ymax) y1 = y1 - (ymax - ymin)
        ypos_comm = floor((y1 - ymin)/(y(ju2-nghost) - y(ju1+nghost) + dy))
     endif
     ypos_commp = ypos_comm + 1
     if (ypos_commp > nyslice - 1) ypos_commp = ypos_commp - nyslice

     send_cpu1m = zposition*nxslice + nxslice*nzslice*ypos_comm
     send_cpu2m = zposition*nxslice + nxslice*nzslice*ypos_commp
     recv_cpu1m = send_cpu1m
     recv_cpu2m = zposition*nxslice + nxslice*nzslice*ypos_commp
  
     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu2p, 10&
          , recv_quant1, jksize, MPI_DOUBLE_PRECISION, recv_cpu1m, 10 &
          , MPI_COMM_WORLD, status, ierr)
     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu1p, 11&
          , recv_quant2, jksize, MPI_DOUBLE_PRECISION, recv_cpu2m, 11 &
          , MPI_COMM_WORLD, status, ierr)

     jstart  = ju1 + nghost + ny*ypos_comm
     jstartp = ju1 + nghost + ny*ypos_commp

     if (abs(ypos_commp - ypos_comm) > 1) then
        inner_quant(jstart-nghost:jstart+ny+nghost-1,ku1:ku2,:,:) &
             = recv_quant1(ju1:ju2,ku1:ku2,:,:,1)
         inner_quant(jstartp-nghost:jstartp+ny+nghost-1,ku1:ku2,:,:) &
              = recv_quant2(ju1:ju2,ku1:ku2,:,:,1)
     else
        if (ypos_commp > ypos_comm) then
           inner_quant(jstart-nghost:jstart+ny-1,ku1:ku2,:,:) &
                = recv_quant1(ju1:ju2-nghost,ku1:ku2,:,:,1)
           inner_quant(jstartp:jstartp+ny+nghost-1,ku1:ku2,:,:) &
                = recv_quant2(ju1+nghost:ju2,ku1:ku2,:,:,1)
        else
           inner_quant(jstart:jstart+ny+nghost-1,ku1:ku2,:,:) &
                = recv_quant1(ju1+nghost:ju2,ku1:ku2,:,:,1)
           inner_quant(jstartp-nghost:jstartp+ny-1,ku1:ku2,:,:) &
                = recv_quant2(ju1:ju2-nghost,ku1:ku2,:,:,1)
        endif
     endif

     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu1m, 10&
          , recv_quant2, jksize, MPI_DOUBLE_PRECISION, recv_cpu2p, 10 &
          , MPI_COMM_WORLD, status, ierr)
     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu2m, 11&
          , recv_quant1, jksize, MPI_DOUBLE_PRECISION, recv_cpu1p, 11 &
          , MPI_COMM_WORLD, status, ierr)

     if (deltay == 0) then
        ypos_comm = yposition - 1
        if (ypos_comm < 0) ypos_comm = ypos_comm + nyslice
     else
        y1 = y(ju1+nghost) - half*dy - deltay
        if (y1 < ymin) y1 = y1 + (ymax - ymin)
        ypos_comm = floor((y1 - ymin)/(y(ju2-nghost) - y(ju1+nghost) + dy))
     endif
     ypos_commp = ypos_comm + 1
     if (ypos_commp > nyslice - 1) ypos_commp = ypos_commp - nyslice

     jstart  = ju1 + nghost + ny*ypos_comm
     jstartp = ju1 + nghost + ny*ypos_commp

     if (abs(ypos_commp - ypos_comm) > 1) then
        outer_quant(jstart-nghost:jstart+ny+nghost-1,ku1:ku2,:,:) &
             = recv_quant1(ju1:ju2,ku1:ku2,:,:,2)
        outer_quant(jstartp-nghost:jstartp+ny+nghost-1,ku1:ku2,:,:) &
             = recv_quant2(ju1:ju2,ku1:ku2,:,:,2)
     endif

     if (abs(ypos_commp - ypos_comm) == 1) then
        if (ypos_commp > ypos_comm) then
           outer_quant(jstart-nghost:jstart+ny-1,ku1:ku2,:,:) &
                = recv_quant1(ju1:ju2-nghost,ku1:ku2,:,:,2)
           outer_quant(jstartp:jstartp+ny+nghost-1,ku1:ku2,:,:) &
                = recv_quant2(ju1+nghost:ju2,ku1:ku2,:,:,2)
        else
           outer_quant(jstart:jstart+ny+nghost-1,ku1:ku2,:,:) &
                = recv_quant1(ju1+nghost:ju2,ku1:ku2,:,:,2)
           outer_quant(jstartp-nghost:jstartp+ny-1,ku1:ku2,:,:) &
                = recv_quant2(ju1:ju2-nghost,ku1:ku2,:,:,2)
        endif
     endif
  endif

  if ((nxslice == 1) .and. (nyslice == 1)) then
     outer_quant(ju1:ju2,ku1:ku2,:,:) &
          = sent_quant(ju1:ju2,ku1:ku2,:,:,2)
     inner_quant(ju1:ju2,ku1:ku2,:,:) &
          = sent_quant(ju1:ju2,ku1:ku2,:,:,1)
  endif

  deallocate(sent_quant, recv_quant1, recv_quant2)

  return
end subroutine get_shearing_quant_old
!===============================================================================
!> Get shearing fluxes - old version
!===============================================================================
subroutine get_shearing_flux_old(inner_quant, outer_quant, deltay)
  use mpi
  use params
  use variables
  use mpi_var
  implicit none

  real(dp), dimension(jf1:jf1+ny*nyslice,kf1:kf2) :: inner_quant, outer_quant
  real(dp) :: deltay, y1

  real(dp), dimension(:,:,:), allocatable :: sent_quant
  real(dp), dimension(:,:,:), allocatable :: recv_quant1, recv_quant2
  integer :: jksize
  integer :: ypos_comm, ypos_commp, jstart, jstartp
  integer :: send_cpu1, send_cpu2, recv_cpu1, recv_cpu2
  integer :: send_cpu1p, send_cpu2p, recv_cpu1p, recv_cpu2p
  integer :: send_cpu1m, send_cpu2m, recv_cpu1m, recv_cpu2m
  integer :: ighost
  integer :: nbound=2
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer :: ierr

  jksize = (jf2 - jf1 + 1)*(kf2 - kf1 + 1)*nbound
  allocate(sent_quant(jf1:jf2,kf1:kf2,nbound))
  allocate(recv_quant1(jf1:jf2,kf1:kf2,nbound))
  allocate(recv_quant2(jf1:jf2,kf1:kf2,nbound))

  !$OMP PARALLEL WORKSHARE
  inner_quant(:,:) = 1000.d0
  outer_quant(:,:) = 1000.d0
  sent_quant(jf1:jf2,kf1:kf2,1) = flux(if1,jf1:jf2,kf1:kf2,1,1)
  sent_quant(jf1:jf2,kf1:kf2,2) = flux(if2,jf1:jf2,kf1:kf2,1,1)
  !$OMP END PARALLEL WORKSHARE

  if ((xposition == 0) .and. (nxslice > 1)) then
     if (deltay == 0) then
        ypos_comm = yposition - 1
        if (ypos_comm < 0) ypos_comm = ypos_comm + nyslice
     else
        y1 = y(ju1+nghost) - half*dy - deltay
        if (y1 < ymin) y1 = y1 + (ymax - ymin)
        ypos_comm = floor((y1 - ymin)/(y(ju2-nghost) - y(ju1+nghost) + dy))
     endif
     ypos_commp = ypos_comm + 1
     if (ypos_commp > nyslice - 1) ypos_commp = ypos_commp - nyslice

     send_cpu1 = zposition*nxslice + (nxslice - 1) + nxslice*nzslice*ypos_comm
     send_cpu2 = zposition*nxslice + (nxslice - 1) + nxslice*nzslice*ypos_commp
     recv_cpu1 = send_cpu1
     recv_cpu2 = zposition*nxslice + (nxslice - 1) + nxslice*nzslice*ypos_commp

     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu2, 11 &
          , recv_quant1, jksize, MPI_DOUBLE_PRECISION, recv_cpu1, 10 &
          , MPI_COMM_WORLD, status, ierr)
     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu1, 13 &
          , recv_quant2, jksize, MPI_DOUBLE_PRECISION, recv_cpu2, 12 &
          , MPI_COMM_WORLD, status, ierr)

     jstart  = jf1 + ny*ypos_comm
     jstartp = jf1 + ny*ypos_commp

     if (abs(ypos_commp - ypos_comm) > 1) then
        outer_quant(jstart:jstart+ny,kf1:kf2)   = recv_quant1(jf1:jf2,kf1:kf2,2)
        outer_quant(jstartp:jstartp+ny,kf1:kf2) = recv_quant2(jf1:jf2,kf1:kf2,2)
     endif

     if (abs(ypos_commp - ypos_comm) == 1) then
        if (ypos_commp > ypos_comm) then
           outer_quant(jstart:jstart+ny-1,kf1:kf2) &
                = recv_quant1(jf1:jf2-1,kf1:kf2,2)
           outer_quant(jstartp:jstartp+ny,kf1:kf2) &
                = recv_quant2(jf1:jf2,kf1:kf2,2)
        else
           outer_quant(jstart:jstart+ny,kf1:kf2) &
                = recv_quant1(jf1:jf2,kf1:kf2,2)
           outer_quant(jstartp:jstartp+ny-1,kf1:kf2) &
                = recv_quant2(jf1:jf2-1,kf1:kf2,2)
        endif
     endif

     if (abs(ypos_commp - ypos_comm) == 0) then
        outer_quant(jf1:jf2,kf1:kf2) = recv_quant1(jf1:jf2,kf1:kf2,2)
     endif
  endif

  if ((xposition == nxslice - 1) .and. (nxslice > 1)) then
     if (deltay == 0) then
        ypos_comm = yposition
     else
        y1 = y(ju1+nghost) - half*dy + deltay
        if (y1 > ymax) y1 = y1 - (ymax - ymin)
        ypos_comm = floor((y1 - ymin)/(y(ju2-nghost) - y(ju1+nghost) + dy))
     endif
     ypos_commp = ypos_comm + 1
     if (ypos_commp > nyslice - 1) ypos_commp = ypos_commp - nyslice

     send_cpu1 = zposition*nxslice + nxslice*nzslice*ypos_comm
     send_cpu2 = zposition*nxslice + nxslice*nzslice*ypos_commp
     recv_cpu1 = send_cpu1
     recv_cpu2 = zposition*nxslice + nxslice*nzslice*ypos_commp
  
     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu2, 10 &
          , recv_quant1, jksize, MPI_DOUBLE_PRECISION, recv_cpu1, 11 &
          , MPI_COMM_WORLD, status, ierr)
     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu1, 12 &
          , recv_quant2, jksize, MPI_DOUBLE_PRECISION, recv_cpu2, 13 &
          , MPI_COMM_WORLD, status, ierr)

     jstart  = jf1 + ny*ypos_comm
     jstartp = jf1 + ny*ypos_commp

     if (abs(ypos_commp - ypos_comm) > 1) then
        inner_quant(jstart:jstart+ny,kf1:kf2)   = recv_quant1(jf1:jf2,kf1:kf2,1)
        inner_quant(jstartp:jstartp+ny,kf1:kf2) = recv_quant2(jf1:jf2,kf1:kf2,1)
     endif

     if (abs(ypos_commp - ypos_comm) == 1) then
        if (ypos_commp > ypos_comm) then
           inner_quant(jstart:jstart+ny-1,kf1:kf2) &
                = recv_quant1(jf1:jf2-1,kf1:kf2,1)
           inner_quant(jstartp:jstartp+ny,kf1:kf2) &
                = recv_quant2(jf1:jf2,kf1:kf2,1)
        else
           inner_quant(jstart:jstart+ny,kf1:kf2) &
                = recv_quant1(jf1:jf2,kf1:kf2,1)
           inner_quant(jstartp:jstartp+ny-1,kf1:kf2) &
                = recv_quant2(jf1:jf2-1,kf1:kf2,1)
         endif
     endif
        
     if (abs(ypos_commp - ypos_comm) == 0) then
        inner_quant(jf1:jf2,kf1:kf2) = recv_quant1(jf1:jf2,kf1:kf2,1)
     endif
  endif

  if ((nxslice == 1) .and. (nyslice > 1)) then
     if (deltay == 0) then
        ypos_comm = yposition - 1
        if (ypos_comm < 0) ypos_comm = ypos_comm + nyslice
     else
        y1 = y(ju1+nghost) - half*dy - deltay
        if (y1 < ymin) y1 = y1 + (ymax - ymin)
        ypos_comm = floor((y1 - ymin)/(y(ju2-nghost) - y(ju1+nghost) + dy))
     endif
     ypos_commp = ypos_comm + 1
     if (ypos_commp > nyslice - 1) ypos_commp = ypos_commp - nyslice

     send_cpu1p = zposition*nxslice + (nxslice - 1) + nxslice*nzslice*ypos_comm
     send_cpu2p = zposition*nxslice + (nxslice - 1) + nxslice*nzslice*ypos_commp
     recv_cpu1p = send_cpu1p
     recv_cpu2p = zposition*nxslice + (nxslice - 1) + nxslice*nzslice*ypos_commp

     if (deltay == 0) then
        ypos_comm = yposition
     else
        y1 = y(ju1+nghost) - half*dy + deltay
        if (y1 > ymax) y1 = y1 - (ymax - ymin)
        ypos_comm = floor((y1 - ymin)/(y(ju2-nghost) - y(ju1+nghost) + dy))
     endif
     ypos_commp = ypos_comm + 1
     if (ypos_commp > nyslice - 1) ypos_commp = ypos_commp - nyslice

     send_cpu1m = zposition*nxslice + nxslice*nzslice*ypos_comm
     send_cpu2m = zposition*nxslice + nxslice*nzslice*ypos_commp
     recv_cpu1m = send_cpu1m
     recv_cpu2m = zposition*nxslice + nxslice*nzslice*ypos_commp
  
     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu2p, 10&
          , recv_quant1, jksize, MPI_DOUBLE_PRECISION, recv_cpu1m, 10 &
          , MPI_COMM_WORLD, status, ierr)
     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu1p, 11&
          , recv_quant2, jksize, MPI_DOUBLE_PRECISION, recv_cpu2m, 11 &
          , MPI_COMM_WORLD, status, ierr)

     jstart  = jf1 + ny*ypos_comm
     jstartp = jf1 + ny*ypos_commp

     if (abs(ypos_commp - ypos_comm) > 1) then
        inner_quant(jstart:jstart+ny,kf1:kf2)   = recv_quant1(jf1:jf2,kf1:kf2,1)
        inner_quant(jstartp:jstartp+ny,kf1:kf2) = recv_quant2(jf1:jf2,kf1:kf2,1)
     else
        if (ypos_commp > ypos_comm) then
           inner_quant(jstart:jstart+ny-1,kf1:kf2) &
                = recv_quant1(jf1:jf2-1,kf1:kf2,1)
           inner_quant(jstartp:jstartp+ny,kf1:kf2) &
                = recv_quant2(jf1:jf2,kf1:kf2,1)
        else
           inner_quant(jstart:jstart+ny,kf1:kf2) &
                = recv_quant1(jf1:jf2,kf1:kf2,1)
           inner_quant(jstartp:jstartp+ny-1,kf1:kf2) &
                = recv_quant2(jf1:jf2-1,kf1:kf2,1)
        endif
     endif

     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu1m, 10&
          , recv_quant2, jksize, MPI_DOUBLE_PRECISION, recv_cpu2p, 10 &
          , MPI_COMM_WORLD, status, ierr)
     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu2m, 11&
          , recv_quant1, jksize, MPI_DOUBLE_PRECISION, recv_cpu1p, 11 &
          , MPI_COMM_WORLD, status, ierr)

     if (deltay == 0) then
        ypos_comm = yposition - 1
        if (ypos_comm < 0) ypos_comm = ypos_comm + nyslice
     else
        y1 = y(ju1+nghost) - half*dy - deltay
        if (y1 < ymin) y1 = y1 + (ymax - ymin)
        ypos_comm = floor((y1 - ymin)/(y(ju2-nghost) - y(ju1+nghost) + dy))
     endif
     ypos_commp = ypos_comm + 1
     if (ypos_commp > nyslice - 1) ypos_commp = ypos_commp - nyslice

     jstart  = jf1 + ny*ypos_comm
     jstartp = jf1 + ny*ypos_commp

     if (abs(ypos_commp - ypos_comm) > 1) then
        outer_quant(jstart:jstart+ny,kf1:kf2)   = recv_quant1(jf1:jf2,kf1:kf2,2)
        outer_quant(jstartp:jstartp+ny,kf1:kf2) = recv_quant2(jf1:jf2,kf1:kf2,2)
     endif

     if (abs(ypos_commp - ypos_comm) == 1) then
        if (ypos_commp > ypos_comm) then
           outer_quant(jstart:jstart+ny-1,kf1:kf2) &
                = recv_quant1(jf1:jf2-1,kf1:kf2,2)
           outer_quant(jstartp:jstartp+ny,kf1:kf2) &
                = recv_quant2(jf1:jf2,kf1:kf2,2)
        else
           outer_quant(jstart:jstart+ny,kf1:kf2) &
                = recv_quant1(jf1:jf2,kf1:kf2,2)
           outer_quant(jstartp:jstartp+ny-1,kf1:kf2) &
                = recv_quant2(jf1:jf2-1,kf1:kf2,2)
        endif
     endif
  endif

  if ((nxslice == 1) .and. (nyslice == 1)) then
     outer_quant(jf1:jf2,kf1:kf2) = sent_quant(jf1:jf2,kf1:kf2,2)
     inner_quant(jf1:jf2,kf1:kf2) = sent_quant(jf1:jf2,kf1:kf2,1)
  endif

  deallocate(sent_quant, recv_quant1, recv_quant2)

  return
end subroutine get_shearing_flux_old
!===============================================================================
!> Get shearing EMF - old version
!===============================================================================
subroutine get_shearing_emf_old(inner_quant, outer_quant, deltay)
  use mpi
  use params
  use variables
  use mpi_var
  implicit none

  real(dp), dimension(ju1:ju1+ny*nyslice+2*nghost-1,ku1:ku2) :: inner_quant, outer_quant
  real(dp) :: deltay, y1

  real(dp), dimension(:,:,:), allocatable :: sent_quant
  real(dp), dimension(:,:,:), allocatable :: recv_quant1, recv_quant2
  integer :: jksize
  integer :: ypos_comm, ypos_commp, jstart, jstartp
  integer :: send_cpu1, send_cpu2, recv_cpu1, recv_cpu2
  integer :: send_cpu1p, send_cpu2p, recv_cpu1p, recv_cpu2p
  integer :: send_cpu1m, send_cpu2m, recv_cpu1m, recv_cpu2m
  integer :: ighost
  integer :: nbound=2
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer :: ierr

  jksize = (ju2 - ju1 + 1)*(ku2 - ku1 + 1)*nbound
  allocate(sent_quant(ju1:ju2,ku1:ku2,nbound))
  allocate(recv_quant1(ju1:ju2,ku1:ku2,nbound))
  allocate(recv_quant2(ju1:ju2,ku1:ku2,nbound))

  inner_quant(:,:) = 1000.d0
  outer_quant(:,:) = 1000.d0

  sent_quant(ju1:ju2,ku1:ku2,1) = emfy(iu1+nghost,ju1:ju2,ku1:ku2)
  sent_quant(ju1:ju2,ku1:ku2,2) = emfy(iu2-nghost+1,ju1:ju2,ku1:ku2)

  if ((xposition == 0) .and. (nxslice > 1)) then
     if (deltay == 0) then
        ypos_comm = yposition - 1
        if (ypos_comm < 0) ypos_comm = ypos_comm + nyslice
     else
        y1 = y(ju1+nghost) - half*dy - deltay
        if (y1 < ymin) y1 = y1 + (ymax - ymin)
        ypos_comm = floor((y1 - ymin)/(y(ju2-nghost) - y(ju1+nghost) + dy))
     endif
     ypos_commp = ypos_comm + 1
     if (ypos_commp > nyslice - 1) ypos_commp = ypos_commp - nyslice

     send_cpu1 = zposition*nxslice + (nxslice - 1) + nxslice*nzslice*ypos_comm
     send_cpu2 = zposition*nxslice + (nxslice - 1) + nxslice*nzslice*ypos_commp
     recv_cpu1 = send_cpu1
     recv_cpu2 = zposition*nxslice + (nxslice - 1) + nxslice*nzslice*ypos_commp

     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu2, 11 &
          , recv_quant1, jksize, MPI_DOUBLE_PRECISION, recv_cpu1, 10 &
          , MPI_COMM_WORLD, status, ierr)
     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu1, 13 &
          , recv_quant2, jksize, MPI_DOUBLE_PRECISION, recv_cpu2, 12 &
          , MPI_COMM_WORLD, status, ierr)

     jstart  = ju1 + nghost + ny*ypos_comm
     jstartp = ju1 + nghost + ny*ypos_commp

     if (abs(ypos_commp - ypos_comm) > 1) then
        outer_quant(jstart-nghost:jstart+ny+nghost-1,ku1:ku2) &
             = recv_quant1(ju1:ju2,ku1:ku2,2)
        outer_quant(jstartp-nghost:jstartp+ny+nghost-1,ku1:ku2) &
             = recv_quant2(ju1:ju2,ku1:ku2,2)
     endif

     if (abs(ypos_commp - ypos_comm) == 1) then
        if (ypos_commp > ypos_comm) then
           outer_quant(jstart-nghost:jstart+ny-1,ku1:ku2) &
                = recv_quant1(ju1:ju2-nghost,ku1:ku2,2)
           outer_quant(jstartp:jstartp+ny+nghost-1,ku1:ku2) &
                = recv_quant2(ju1+nghost:ju2,ku1:ku2,2)
        else
           outer_quant(jstart:jstart+ny+nghost-1,ku1:ku2) &
                = recv_quant1(ju1+nghost:ju2,ku1:ku2,2)
           outer_quant(jstartp-nghost:jstartp+ny-1,ku1:ku2) &
                = recv_quant2(ju1:ju2-nghost,ku1:ku2,2)
        endif
     endif

     if (abs(ypos_commp - ypos_comm) == 0) then
        outer_quant(ju1:ju2,ku1:ku2) = recv_quant1(ju1:ju2,ku1:ku2,2)
     endif
  endif

  if ((xposition == nxslice - 1) .and. (nxslice > 1)) then
     if (deltay == 0) then
        ypos_comm = yposition
     else
        y1 = y(ju1+nghost) - half*dy + deltay
        if (y1 > ymax) y1 = y1 - (ymax - ymin)
        ypos_comm = floor((y1 - ymin)/(y(ju2-nghost) - y(ju1+nghost) + dy))
     endif
     ypos_commp = ypos_comm + 1
     if (ypos_commp > nyslice - 1) ypos_commp = ypos_commp - nyslice

     send_cpu1 = zposition*nxslice + nxslice*nzslice*ypos_comm
     send_cpu2 = zposition*nxslice + nxslice*nzslice*ypos_commp
     recv_cpu1 = send_cpu1
     recv_cpu2 = zposition*nxslice + nxslice*nzslice*ypos_commp
  
     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu2, 10 &
          , recv_quant1, jksize, MPI_DOUBLE_PRECISION, recv_cpu1, 11 &
          , MPI_COMM_WORLD, status, ierr)
     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu1, 12 &
          , recv_quant2, jksize, MPI_DOUBLE_PRECISION, recv_cpu2, 13 &
          , MPI_COMM_WORLD, status, ierr)

     jstart  = ju1 + nghost + ny*ypos_comm
     jstartp = ju1 + nghost + ny*ypos_commp

     if (abs(ypos_commp - ypos_comm) > 1) then
        inner_quant(jstart-nghost:jstart+ny+nghost-1,ku1:ku2) &
             = recv_quant1(ju1:ju2,ku1:ku2,1)
        inner_quant(jstartp-nghost:jstartp+ny+nghost-1,ku1:ku2) &
             = recv_quant2(ju1:ju2,ku1:ku2,1)
     endif

     if (abs(ypos_commp - ypos_comm) == 1) then
        if (ypos_commp > ypos_comm) then
           inner_quant(jstart-nghost:jstart+ny-1,ku1:ku2) &
                = recv_quant1(ju1:ju2-nghost,ku1:ku2,1)
           inner_quant(jstartp:jstartp+ny+nghost-1,ku1:ku2) &
                = recv_quant2(ju1+nghost:ju2,ku1:ku2,1)
        else
           inner_quant(jstart:jstart+ny+nghost-1,ku1:ku2) &
                = recv_quant1(ju1+nghost:ju2,ku1:ku2,1)
           inner_quant(jstartp-nghost:jstartp+ny-1,ku1:ku2) &
                = recv_quant2(ju1:ju2-nghost,ku1:ku2,1)
        endif
     endif
        
     if (abs(ypos_commp - ypos_comm) == 0) then
        inner_quant(ju1:ju2,ku1:ku2) = recv_quant1(ju1:ju2,ku1:ku2,1)
     endif
  endif

  if ((nxslice == 1) .and. (nyslice > 1)) then
     if (deltay == 0) then
        ypos_comm = yposition - 1
        if (ypos_comm < 0) ypos_comm = ypos_comm + nyslice
     else
        y1 = y(ju1+nghost) - half*dy - deltay
        if (y1 < ymin) y1 = y1 + (ymax - ymin)
        ypos_comm = floor((y1 - ymin)/(y(ju2-nghost) - y(ju1+nghost) + dy))
     endif
     ypos_commp = ypos_comm + 1
     if (ypos_commp > nyslice - 1) ypos_commp = ypos_commp - nyslice

     send_cpu1p = zposition*nxslice + (nxslice - 1) + nxslice*nzslice*ypos_comm
     send_cpu2p = zposition*nxslice + (nxslice - 1) + nxslice*nzslice*ypos_commp
     recv_cpu1p = send_cpu1p
     recv_cpu2p = zposition*nxslice + (nxslice - 1) + nxslice*nzslice*ypos_commp

     if (deltay == 0) then
        ypos_comm = yposition
     else
        y1 = y(ju1+nghost) - half*dy + deltay
        if (y1 > ymax) y1 = y1 - (ymax - ymin)
        ypos_comm = floor((y1 - ymin)/(y(ju2-nghost) - y(ju1+nghost) + dy))
     endif
     ypos_commp = ypos_comm + 1
     if (ypos_commp > nyslice - 1) ypos_commp = ypos_commp - nyslice

     send_cpu1m = zposition*nxslice + nxslice*nzslice*ypos_comm
     send_cpu2m = zposition*nxslice + nxslice*nzslice*ypos_commp
     recv_cpu1m = send_cpu1m
     recv_cpu2m = zposition*nxslice + nxslice*nzslice*ypos_commp
  
     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu2p, 10&
          , recv_quant1, jksize, MPI_DOUBLE_PRECISION, recv_cpu1m, 10 &
          , MPI_COMM_WORLD, status, ierr)
     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu1p, 11&
          , recv_quant2, jksize, MPI_DOUBLE_PRECISION, recv_cpu2m, 11 &
          , MPI_COMM_WORLD, status, ierr)

     jstart  = ju1 + nghost + ny*ypos_comm
     jstartp = ju1 + nghost + ny*ypos_commp

     if (abs(ypos_commp - ypos_comm) > 1) then
        inner_quant(jstart-nghost:jstart+ny+nghost-1,ku1:ku2) &
             = recv_quant1(ju1:ju2,ku1:ku2,1)
         inner_quant(jstartp-nghost:jstartp+ny+nghost-1,ku1:ku2) &
              = recv_quant2(ju1:ju2,ku1:ku2,1)
     else
        if (ypos_commp > ypos_comm) then
           inner_quant(jstart-nghost:jstart+ny-1,ku1:ku2) &
                = recv_quant1(ju1:ju2-nghost,ku1:ku2,1)
           inner_quant(jstartp:jstartp+ny+nghost-1,ku1:ku2) &
                = recv_quant2(ju1+nghost:ju2,ku1:ku2,1)
        else
           inner_quant(jstart:jstart+ny+nghost-1,ku1:ku2) &
                = recv_quant1(ju1+nghost:ju2,ku1:ku2,1)
           inner_quant(jstartp-nghost:jstartp+ny-1,ku1:ku2) &
                = recv_quant2(ju1:ju2-nghost,ku1:ku2,1)
        endif
     endif

     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu1m, 10&
          , recv_quant2, jksize, MPI_DOUBLE_PRECISION, recv_cpu2p, 10 &
          , MPI_COMM_WORLD, status, ierr)
     call MPI_Sendrecv(sent_quant, jksize, MPI_DOUBLE_PRECISION, send_cpu2m, 11&
          , recv_quant1, jksize, MPI_DOUBLE_PRECISION, recv_cpu1p, 11 &
          , MPI_COMM_WORLD, status, ierr)

     if (deltay == 0) then
        ypos_comm = yposition - 1
        if (ypos_comm < 0) ypos_comm = ypos_comm + nyslice
     else
        y1 = y(ju1+nghost) - half*dy - deltay
        if (y1 < ymin) y1 = y1 + (ymax - ymin)
        ypos_comm = floor((y1 - ymin)/(y(ju2-nghost) - y(ju1+nghost) + dy))
     endif
     ypos_commp = ypos_comm + 1
     if (ypos_commp > nyslice - 1) ypos_commp = ypos_commp - nyslice

     jstart  = ju1 + nghost + ny*ypos_comm
     jstartp = ju1 + nghost + ny*ypos_commp

     if (abs(ypos_commp - ypos_comm) > 1) then
        outer_quant(jstart-nghost:jstart+ny+nghost-1,ku1:ku2) &
             = recv_quant1(ju1:ju2,ku1:ku2,2)
        outer_quant(jstartp-nghost:jstartp+ny+nghost-1,ku1:ku2) &
             = recv_quant2(ju1:ju2,ku1:ku2,2)
     endif

     if (abs(ypos_commp - ypos_comm) == 1) then
        if (ypos_commp > ypos_comm) then
           outer_quant(jstart-nghost:jstart+ny-1,ku1:ku2) &
                = recv_quant1(ju1:ju2-nghost,ku1:ku2,2)
           outer_quant(jstartp:jstartp+ny+nghost-1,ku1:ku2) &
                = recv_quant2(ju1+nghost:ju2,ku1:ku2,2)
        else
           outer_quant(jstart:jstart+ny+nghost-1,ku1:ku2) &
                = recv_quant1(ju1+nghost:ju2,ku1:ku2,2)
           outer_quant(jstartp-nghost:jstartp+ny-1,ku1:ku2) &
                = recv_quant2(ju1:ju2-nghost,ku1:ku2,2)
        endif
     endif
  endif

  if ((nxslice == 1) .and. (nyslice == 1)) then
     outer_quant(ju1:ju2,ku1:ku2) &
          = sent_quant(ju1:ju2,ku1:ku2,2)
     inner_quant(ju1:ju2,ku1:ku2) &
          = sent_quant(ju1:ju2,ku1:ku2,1)
  endif

  deallocate(sent_quant, recv_quant1, recv_quant2)

  return
end subroutine get_shearing_emf_old
#endif
