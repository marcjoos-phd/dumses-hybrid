!===============================================================================
!> \file parallel.f90
!! \brief
!! \b This is part of a merging program for DUMSES outputs:
!! Manage parallelism
!! \details
!! Contains init_parallel(), finalize_parallel()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang
!! \copyright
!! Copyrights 2014-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          27-08-2014
!! \b last \b modified: 29-08-2014
!<
!===============================================================================
module parallel

contains
!===============================================================================
!> Initialize parallel libraries
!===============================================================================
  subroutine init_parallel
    use params
#if MPI == 1
    use mpi
#endif
    use mpi_var
    implicit none
    
#if MPI == 1
    call MPI_Init(ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, nproc, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, myrank, ierr)
#endif

    if (myrank == 0) verbose = .true.
  
    if (verbose) then
       print '("This is MERCO, a MERging and COnverting tool for DUMSES outputs")'
       print '("===================================================================")'
       print '("Execution with ", I3, " MPI process(es)")', nproc
       print '("===================================================================")'
    endif
  
    return
  end subroutine init_parallel
!===============================================================================
!> Finalize parallel libraries
!===============================================================================
  subroutine finalize_parallel
#if MPI == 1
    use mpi
    use mpi_var
    implicit none
  
    call MPI_Finalize(ierr)
#endif
    
    return
  end subroutine finalize_parallel
end module parallel
