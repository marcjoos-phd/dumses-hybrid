!===============================================================================
!> \file history.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is history subroutines.
!! \details
!! Contains history(), compute_yz_mean(), compute_xz_mean(), compute_x_mean(),
!! compute_y_mean(), compute_z_mean()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Pierre Kestener, 
!! Romain Teyssier, Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          06-13-2014 
!! \b last \b modified: 04-29-2015
!<
!===============================================================================
!> \brief
!! This routine produces history file of the simulation.
!===============================================================================
subroutine history
  use params
  use mpi_var
  use variables
  implicit none

  character(LEN=80) :: filename="history.txt", hist_format
  real(dp) :: mass, maxwell, reynolds, dtau, magp, divB
  real(dp), dimension(ndim) :: meanB
  real(dp), dimension(:), allocatable       :: tmp
  real(dp), dimension(:,:,:,:), allocatable :: loc_var
  real(dp), dimension(:,:), allocatable     :: loc_var_yz_mean
  integer :: i, j, k, ipos, nvar_mean
  logical :: fexist

  if (verbose) print*, "> Writing history..."

#if NDIM == 1
  do ipos = 1, nx
     if (x(ipos) >= half) exit
  enddo

  if (mype == 0) then
     if (time == zero) then
        open(unit=2, file=filename, status="unknown")
        write(2, "('# ', A11, 1X, 3(A14,1X))") "time", "rho", "rhoux", "p"
     else
        inquire(file=filename, exist=fexist)
        if (fexist) then
           open(unit=2, file=filename, status="old", position="append")
        else
           open(unit=2, file=filename, status="unknown")
           write(2, "('# ', A11, 1X, 3(A14,1X))") "time", "rho", "rhoux", "p"
        endif
     endif
     hist_format = '(1X, 1PE12.5, 1X, 3(E14.5,1X))'
     write(2, hist_format) time, uin(ipos,1,1,ir), uin(ipos,1,1,iu) &
          , uin(ipos,1,1,ip)
     close(2)
  endif
#endif

#if NDIM == 3
  dtau = dx*dy*dz/(xmax - xmin)/(ymax - ymin)/(zmax - zmin)

  nvar_mean = 3
  allocate(loc_var(iu1:iu2,ju1:ju2,ku1:ku2,nvar_mean))
  allocate(loc_var_yz_mean(iu1:iu2,nvar_mean))

  !$OMP PARALLEL WORKSHARE
  loc_var(:,:,:,1) = uin(:,:,:,ir)
  loc_var(:,:,:,2) = uin(:,:,:,iu)/uin(:,:,:,ir)
  loc_var(:,:,:,3) = uin(:,:,:,iw)/uin(:,:,:,ir)
  !$OMP END PARALLEL WORKSHARE
  
  call compute_yz_mean(loc_var, loc_var_yz_mean, nvar_mean)
  deallocate(loc_var)

  mass    = zero; magp     = zero
  meanB   = zero; divB     = zero
  maxwell = zero; reynolds = zero
  !$OMP PARALLEL DO REDUCTION(+: mass,magp,maxwell,meanB,reynolds,divB) &
  !$OMP SCHEDULE(RUNTIME)
  do k = 1, nz
     do j = 1, ny
        do i = 1, nx
           mass = mass + uin(i,j,k,ir)
           magp = magp &
                + forth*(uin(i,j,k,iA) + uin(i+1,j,k,iA))**2 &
                + forth*(uin(i,j,k,iB) + uin(i,j+1,k,iB))**2 &
                + forth*(uin(i,j,k,iC) + uin(i,j,k+1,iC))**2
           maxwell = maxwell &
                - forth*(uin(i,j,k,iA) + uin(i+1,j,k,iA)) &
                *(uin(i,j,k,iC) + uin(i,j,k+1,iC))
           meanB(1) = meanB(1) + uin(i,j,k,iA)
           meanB(2) = meanB(2) + uin(i,j,k,iB)
           meanB(3) = meanB(3) + uin(i,j,k,iC)

           reynolds = reynolds + uin(i,j,k,ir)*dtau &
                *(uin(i,j,k,iu)/uin(i,j,k,ir) - loc_var_yz_mean(i,2)) &
                *(uin(i,j,k,iw)/uin(i,j,k,ir) - loc_var_yz_mean(i,3))

           divB = divB &
                + (uin(i+1,j,k,iA) - uin(i,j,k,iA))/dx &
                + (uin(i,j+1,k,iB) - uin(i,j,k,iB))/dy &
                + (uin(i,j,k+1,iC) - uin(i,j,k,iC))/dz
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

  deallocate(loc_var_yz_mean)

  magp    = half*magp*dtau
  mass    = mass*dtau
  maxwell = maxwell*dtau
  meanB   = meanB*dtau

#if MPI == 1
  allocate(tmp(1))

  tmp = (/ mass /)
  call vol_average(tmp, 1)
  mass = tmp(1)

  tmp = (/ magp /)
  call vol_average(tmp, 1)
  magp = tmp(1)

  tmp = (/ maxwell /)
  call vol_average(tmp, 1)
  maxwell = tmp(1)

  tmp = (/ reynolds /)
  call vol_average(tmp, 1)
  reynolds = tmp(1)

  tmp = (/ divB /)
  call vol_average(tmp, 1)
  divB = tmp(1)

  deallocate(tmp)

  call vol_average(meanB, 3)
#endif

  if (mype == 0) then
     if (time == zero) then
        open(unit=2, file=filename, status="unknown")
           write(2, "('# ', A11, 1X, A12, 1X, 9(A14,1X))") "time", "dt", "mass" &
                , "maxwell", "reynolds", "max+rey", "magp", "meanBx", "meanBy" &
                , "meanBz", "divB"
     else
        inquire(file=filename, exist=fexist)
        if (fexist) then
           open(unit=2, file=filename, status="old", position="append")
        else
           open(unit=2, file=filename, status="unknown")
           write(2, "('# ', A11, 1X, A12, 1X, 9(A14,1X))") "time", "dt", "mass" &
                , "maxwell", "reynolds", "max+rey", "magp", "meanBx", "meanBy" &
                , "meanBz", "divB"
        endif
     endif
     hist_format = '(1X, 1PE12.5, 1X, 1PE12.5, 1X, 9(E14.5,1X))'
     write(2, hist_format) time, dt, mass, maxwell, reynolds, maxwell+reynolds &
          , magp, meanB(1), meanB(2), meanB(3), divB
     close(2)
  endif
#endif

  return
end subroutine history
!===============================================================================
!> This routine computes volume average
!===============================================================================
subroutine vol_average(quant, nquant)
#if MPI == 1
  use params
  use mpi
  implicit none

  integer, intent(in) :: nquant
  real(dp), dimension(nquant), intent(inout) :: quant
  real(dp), dimension(nquant) :: vol_quant
  integer :: ierr

  call MPI_Reduce(quant, vol_quant, nquant, MPI_DOUBLE_PRECISION, MPI_SUM, 0 &
       , MPI_COMM_WORLD, ierr)
  quant = vol_quant

#endif
  return
end subroutine vol_average
!===============================================================================
!> This routine computes mean in the yz-plane
!===============================================================================
subroutine compute_yz_mean(quant, mean_quant, nquant)
#if MPI == 1
  use mpi
  use mpi_var
#endif
  use params

  integer, intent(in) :: nquant
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nquant), intent(in) :: quant
  real(dp), dimension(iu1:iu2,nquant), intent(out) :: mean_quant
  real(dp), dimension(:,:), allocatable :: mean_quant_local
  real(dp) :: mean_quant_loc
  integer :: ivar, i, j, k
#if MPI == 1
  integer :: ierr, commx
#endif

  allocate(mean_quant_local(iu1:iu2, nquant))

  do ivar = 1, nquant
     do i = iu1, iu2
        mean_quant_loc = zero
        !$OMP PARALLEL DO REDUCTION(+: mean_quant_loc) SCHEDULE(RUNTIME)
        do k = 1, nz
           do j = 1, ny
              mean_quant_loc = mean_quant_loc + quant(i,j,k,ivar)
           enddo
        enddo
        !$OMP END PARALLEL DO
        mean_quant_local(i,ivar) = mean_quant_loc
     enddo
  enddo

#if MPI == 1
  call create_comm_plane(commx, 1)
  mean_quant = zero
  call MPI_Reduce(mean_quant_local, mean_quant, (iu2-iu1+1)*nquant &
       , MPI_DOUBLE_PRECISION, MPI_SUM, 0, commx, ierr)
  mean_quant = mean_quant/ny/nz/nyslice/nzslice
  call MPI_Bcast(mean_quant, (iu2-iu1+1)*nquant, MPI_DOUBLE_PRECISION &
       , 0, commx, ierr)
  call MPI_Comm_Free(commx, ierr)
#else
  mean_quant = mean_quant_local/ny/nz
#endif
  
  deallocate(mean_quant_local)

  return
end subroutine compute_yz_mean
!===============================================================================
!> This routine computes mean in the xz-plane
!===============================================================================
subroutine compute_xz_mean(quant, mean_quant, nquant)
#if MPI == 1
  use mpi
  use mpi_var
#endif
  use params

  integer, intent(in) :: nquant
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nquant), intent(in) :: quant
  real(dp), dimension(ju1:ju2,nquant), intent(out) :: mean_quant
  real(dp), dimension(:,:), allocatable :: mean_quant_local
  real(dp) :: mean_quant_loc
  integer :: ivar, i, j, k
#if MPI == 1
  integer :: ierr, commy
#endif

  allocate(mean_quant_local(ju1:ju2,nquant))

  do ivar = 1, nquant
     do j = ju1, ju2
        mean_quant_loc = zero
        !$OMP PARALLEL DO REDUCTION(+: mean_quant_loc) SCHEDULE(RUNTIME)
        do k = 1, nz
           do i = 1, nx
              mean_quant_loc = mean_quant_loc + quant(i,j,k,ivar)
           enddo
        enddo
        !$OMP END PARALLEL DO
        mean_quant_local(j,ivar) = mean_quant_loc
     enddo
  enddo

#if MPI == 1
  call create_comm_plane(commy, 2)
  mean_quant = zero
  call MPI_Reduce(mean_quant_local, mean_quant, (ju2-ju1+1)*nquant &
       , MPI_DOUBLE_PRECISION, MPI_SUM, xposition, commy, ierr)
  mean_quant = mean_quant/nx/nz/nxslice/nzslice
  call MPI_Bcast(mean_quant, (ju2-ju1+1)*nquant, MPI_DOUBLE_PRECISION &
       , xposition, commy, ierr)
  call MPI_Comm_Free(commy, ierr)
#else
  mean_quant = mean_quant_local/nx/nz
#endif

  deallocate(mean_quant_local)

  return
end subroutine compute_xz_mean
!===============================================================================
!> This routine computes mean in the x direction
!===============================================================================
subroutine compute_x_mean(quant, mean_quant, nquant)
#if MPI == 1
  use mpi
  use mpi_var
#endif
  use params
  implicit none

  integer, intent(in) :: nquant
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nquant), intent(in) :: quant
  real(dp), dimension(ju1:ju2,ku1:ku2,nquant), intent(out) :: mean_quant
  real(dp), dimension(:,:,:), allocatable :: mean_quant_local
  real(dp) :: mean_quant_loc
  integer :: i, j, k, ivar
#if MPI == 1
  integer :: ierr, commx
#endif

  allocate(mean_quant_local(ju1:ju2,ku1:ku2,nquant))

  mean_quant_local = zero
  do ivar = 1, nquant
     do k = ku1, ku2
        do j = ju1, ju2
           !$OMP PARALLEL DO REDUCTION(+: mean_quant_loc) SCHEDULE(RUNTIME)
           do i = 1, nx
              mean_quant_loc = mean_quant_loc + quant(i,j,k,ivar)
           enddo
           !$OMP END PARALLEL DO
           mean_quant_local(j,k,ivar) = mean_quant_loc
        enddo
     enddo
  enddo
  
#if MPI == 1
  call create_comm_line(commx, 1)
  mean_quant = zero
  call MPI_Reduce(mean_quant_local, mean_quant, (ju2-ju1+1)*(ku2-ku1+1)*nquant &
       , MPI_DOUBLE_PRECISION, MPI_SUM, 0, commx, ierr)
  mean_quant = mean_quant/nx/nxslice
  call MPI_Bcast(mean_quant, (ju2-ju1+1)*(ku2-ku1+1)*nquant &
       , MPI_DOUBLE_PRECISION, 0, commx, ierr)
  call MPI_Comm_Free(commx, ierr)
#else
  mean_quant = mean_quant_local/nx
#endif

  deallocate(mean_quant_local)

  return
end subroutine compute_x_mean
!===============================================================================
!> This routine computes mean in the y direction
!===============================================================================
subroutine compute_y_mean(quant, mean_quant, nquant)
#if MPI == 1
  use mpi
  use mpi_var
#endif
  use params
  implicit none

  integer, intent(in) :: nquant
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nquant), intent(in) :: quant
  real(dp), dimension(ju1:ju2,ku1:ku2,nquant), intent(out) :: mean_quant
  real(dp), dimension(:,:,:), allocatable :: mean_quant_local
  real(dp) :: mean_quant_loc
  integer :: i, j, k, ivar
#if MPI == 1
  integer :: ierr, commy
#endif

  allocate(mean_quant_local(iu1:iu2,ku1:ku2,nquant))
  
  mean_quant_local = zero
  do ivar = 1, nquant
     do k = ku1, ku2
        do i = iu1, iu2
           !$OMP PARALLEL DO REDUCTION(+: mean_quant_loc) SCHEDULE(RUNTIME)
           do j = 1, ny
              mean_quant_loc = mean_quant_loc + quant(i,j,k,ivar)
           enddo
           !$OMP END PARALLEL DO
           mean_quant_local(i,k,ivar) = mean_quant_loc
        enddo
     enddo
  enddo

#if MPI == 1
  call create_comm_line(commy, 2)
  mean_quant = zero
  call MPI_Reduce(mean_quant_local, mean_quant, (iu2-iu1+1)*(ku2-ku1+1)*nquant &
       , MPI_DOUBLE_PRECISION, MPI_SUM, 0, commy, ierr)
  mean_quant = mean_quant/ny/nyslice
  call MPI_Bcast(mean_quant, (iu2-iu1+1)*(ku2-ku1+1)*nquant &
       , MPI_DOUBLE_PRECISION, 0, commy, ierr)
  call MPI_Comm_Free(commy, ierr)
#else
  mean_quant = mean_quant_local/ny
#endif

  deallocate(mean_quant_local)

  return
end subroutine compute_y_mean
!===============================================================================
!> This routine computes mean in the z direction
!===============================================================================
subroutine compute_z_mean(quant, mean_quant, nquant)
#if MPI == 1
  use mpi
  use mpi_var
#endif
  use params
  implicit none

  integer, intent(in) :: nquant
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nquant), intent(in) :: quant
  real(dp), dimension(ju1:ju2,ku1:ku2,nquant), intent(out) :: mean_quant
  real(dp), dimension(:,:,:), allocatable :: mean_quant_local
  real(dp) :: mean_quant_loc
  integer :: i, j, k, ivar
#if MPI == 1
  integer :: ierr, commz
#endif

  allocate(mean_quant_local(iu1:iu2,ju1:ju2,nquant))

  mean_quant_local = zero
  do ivar = 1, nquant
     do j = ju1, ju2
        do i = iu1, iu2
           !$OMP PARALLEL DO REDUCTION(+: mean_quant_loc) SCHEDULE(RUNTIME)
           do k = 1, nz
              mean_quant_loc = mean_quant_loc + quant(i,j,k,ivar)
           enddo
           !$OMP END PARALLEL DO
           mean_quant_local(i,j,ivar) = mean_quant_loc
        enddo
     enddo
  enddo
  
#if MPI == 1
  call create_comm_line(commz, 3)
  mean_quant = zero
  call MPI_Reduce(mean_quant_local, mean_quant, (iu2-iu1+1)*(ju2-ju1+1)*nquant &
       , MPI_DOUBLE_PRECISION, MPI_SUM, 0, commz, ierr)
  mean_quant = mean_quant/nz/nzslice
  call MPI_Bcast(mean_quant, (iu2-iu1+1)*(ju2-ju1+1)*nquant &
       , MPI_DOUBLE_PRECISION, 0, commz, ierr)
  call MPI_Comm_Free(commz, ierr)
#else
  mean_quant = mean_quant_local/nz
#endif

  deallocate(mean_quant_local)

  return
end subroutine compute_z_mean
