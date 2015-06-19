!===============================================================================
!> \file merge.f90
!! \brief
!! \b This is a merging and converting tool for DUMSES outputs
!! \details
!! Contains ...
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang
!! \copyright
!! Copyrights 2014-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          08-26-2014
!! \b last \b modified: 05-28-2015
!<
!===============================================================================
!> Merge and convert DUMSES output files
!===============================================================================
program merco
#if MPI == 1
  use mpi
#endif
  use params
  use variables
  use mpi_var
  use parallel
  use file
  implicit none

  ! Print variable
  character(len=20) :: version=""

  ! Timing variables
  integer :: t0, t1, irate
  real(dp) :: tcpu0, tcpu1, cputime, elptime

  ! Initialization
  call init_parallel
  call init_IO_libraries

  ! Get parameters
  call get_parameters

  if (old) version = " (old output format)"

  if (verbose) then
     print '(" > Output to process: ", I6)', ndump
     print '(" > Convert DUMSES output", A20)', version 
     print '(" > From:   ", A9)', type_in
     print '(" > To:     ", A9)', type_out
     if (xscale*yscale*zscale > 1) then
        print '(" > Resize the grid with:")'
        print '("   - xscale: ", I3)', xscale
        print '("   - yscale: ", I3)', yscale
        print '("   - zscale: ", I3)', zscale
     endif
  endif

  ! Timing
  call cpu_time(tcpu0)
  call system_clock(count=t0, count_rate=irate)

  ! Get metadata
  call get_metadata
  if (verbose) then
     print '(" ")'
     print '("Metadata:")'
     print '(" > Size of a subdomain: ", I6, ", ", I6, ", ", I6)', nx, ny, nz
     print '(" > Computed on ", I6, " MPI process(es)")', npes
  endif
#if MPI == 1
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

  ! Write metadata
  call write_metadata

  ! Read & write data
  call readwrite_data

  ! Timing
#if MPI == 1
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
  if(myrank == 0) then
     call system_clock(count=t1, count_rate=irate)
     elptime = real(t1 - t0, kind=8)/real(irate, kind=8)
     call cpu_time(tcpu1)
     cputime = tcpu1 - tcpu0
     print '(/ 3X,"Elapsed time : ",1PE10.3," s",/ &
             &,3X,"CPU time     : ",1PE10.3," s", /)', elptime, cputime
  endif

  ! Finalization
  call finalize_IO_libraries
  call finalize_parallel

end program merco
