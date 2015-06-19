!===============================================================================
!> \file history.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is history subroutines.
!! \details
!! Contains history()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Pierre Kestener, 
!! Romain Teyssier, Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          07-08-2014 
!! \b last \b modified: 07-08-2014
!<
!===============================================================================
!> \brief
!! This routine produces history file of the simulation.
!===============================================================================
subroutine history
  use params
  use mpi_var
  use variables
#if MPI == 1
  use mpi
#endif
  implicit none

  integer :: ipos, kpos
  integer :: i, j, k
  character(LEN=80) :: filename="history.txt", hist_format
  real(dp) :: rho, dvx, dvy
  real(dp) :: mass
  logical  :: fexist

  if (verbose) print*, "> Writing history..."

#if MPI == 1
  if (mod(nxslice,2) == 0) then
     ipos = 1
  else
     ipos = nx/2
  endif
  if (mod(nzslice,2) == 0) then
     kpos = 1
  else
     kpos = nz/2
  endif
#else
  ipos = nx/2
  kpos = nz/2
#endif

  rho = uin(ipos,1,kpos,ir)
  dvx = uin(ipos,1,kpos,iu)/rho
  dvy = uin(ipos,1,kpos,iv)/rho

  mass = zero
  !$OMP PARALLEL DO REDUCTION(+: mass) SCHEDULE(RUNTIME)
  do k = 1, nz
     do j = 1, ny
        do i = 1, nx
           mass = mass + uin(i,j,k,ir)
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

#if MPI == 1
  if ((xposition == nxslice/2) .and. (zposition == nzslice/2)) then
#else
  if (mype == 0) then
#endif
     inquire(file=filename, exist=fexist)
     if (fexist) then
        open(unit=2, file=filename, status="old", position="append")
     else
        open(unit=2, file=filename, status="unknown")
        write(2, "('# ', A11, 1X, 5(A14,1X))") "time", "dt", "rho", "dvx/ciso" &
             , "dvy/ciso", "mass"
     endif
     hist_format = '(1X, 6(E14.5,1X))'
     write(2, hist_format) time, dt, rho, dvx/ciso, dvy/ciso, mass
     close(2)
  endif

  return
end subroutine history
