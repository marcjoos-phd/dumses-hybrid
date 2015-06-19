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
!! \b created:          07-07-2014 
!! \b last \b modified: 07-07-2014
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
  real(dp) :: rho, dvx, dvy
  logical :: fexist

  if (verbose) print*, "> Writing history..."

  rho = maxval(uin(1:nx,1:ny,1:nz,ir))
  dvx = maxval(uin(1:nx,1:ny,1:nz,iu)/uin(1:nx,1:ny,1:nz,ir))
  dvy = maxval(uin(1:nx,1:ny,1:nz,iv)/uin(1:nx,1:ny,1:nz,ir))

  if (mype == 0) then
     inquire(file=filename, exist=fexist)
     if (fexist) then
        open(unit=2, file=filename, status="old", position="append")
     else
        open(unit=2, file=filename, status="unknown")
        write(2, "('# ', A11, 1X, 4(A14,1X))") "time", "dt", "rho", "dvx/ciso" &
             , "dvy/ciso"
     endif
     hist_format = '(1X, 5(E14.5,1X))'
     write(2, hist_format) time, dt, rho, dvx/ciso, dvy/ciso
     close(2)
  endif

  return
end subroutine history
