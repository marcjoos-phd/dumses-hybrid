!===============================================================================
!> \file file.f90
!! \brief
!! \b This is part of a merging program for DUMSES outputs:
!! Module file, contains subroutines for file handling
!! \details
!! Contains init_IO_libraries(), finalize_IO_libraries(), get_filename(),
!! convtoasc()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang
!! \copyright
!! Copyrights 2014, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          28-08-2014
!! \b last \b modified: 29-08-2014
!<
!===============================================================================
module file

contains
!===============================================================================
!> Initialize libraries
!===============================================================================
  subroutine init_IO_libraries
#if HDF5 == 1 || PHDF5 == 1
    use hdf5
#endif
    use mpi_var
    implicit none

#if HDF5 == 1 || PHDF5 == 1
    call H5open_f(ierr)
#endif

    return
  end subroutine init_IO_libraries
!===============================================================================
!> Finalize libraries
!===============================================================================
  subroutine finalize_IO_libraries
#if HDF5 == 1 || PHDF5 == 1
    use hdf5
#endif
    use mpi_var
    implicit none

#if HDF5 == 1 || PHDF5 == 1
    call H5close_f(ierr)
#endif

    return
  end subroutine finalize_IO_libraries
!===============================================================================
!> Get \c filename from process rank and number of the output
!===============================================================================
  subroutine get_filename(ndump, mype, prefix, filename)
    implicit none

    integer           :: ndump, mype
    character(LEN=6)  :: snumdir, snumfile
    character(LEN=*)  :: prefix
    character(LEN=80) :: filename, filedir, filecmd
    
    call convtoasc(ndump, snumdir)
    call convtoasc(mype, snumfile)
    filedir  = 'output_'//trim(snumdir)//'/'
    filename = trim(filedir)//trim(prefix)//'.'//trim(snumfile)
    
    return
  end subroutine get_filename
!===============================================================================
!> Convert a integer in a 6 characters string
!===============================================================================
  subroutine convtoasc(number, sstring)
    implicit none
    integer :: number, istring, num, nums10, i
    character(LEN=6) :: sstring
    character(LEN=10), parameter :: nstring="0123456789"
    
    num    = 1000000
    nums10 = num/10
    do i = 1, 6
       istring      = 1 + mod(number,num)/nums10
       sstring(i:i) = nstring(istring:istring)
       num    = num/10
       nums10 = nums10/10
    enddo
    
    return
  end subroutine convtoasc
end module file
