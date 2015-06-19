!===============================================================================
!> \file input.f90
!! \brief
!! \b This is part of a merging program for DUMSES outputs:
!! Get input parameters
!! \details
!! Contains get_parameters(), default_parameters()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang
!! \copyright
!! Copyrights 2014-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          08-27-2014
!! \b last \b modified: 05-27-2015
!<
!===============================================================================
!> Get parameters
!===============================================================================
subroutine get_parameters
  use params
  implicit none
  
#if INTER == 1
  integer :: i, n
  integer :: iargc
  character(len=9)   :: opt
  character(len=128) :: arg
  
  call default_parameters
  
  n = iargc()
  do i = 1, n, 2
     call getarg(i, opt)
     if (i == n) then
        if (verbose) print '("Option ", A9," has no argument")', opt
        stop 2
     endif
     call getarg(i+1,arg)
     select case (opt)
     case ('-ndump')
        read(arg,*) ndump
     case ('-old')
        read(arg,*) old
     case ('-type_in')
        read(arg,*) type_in
     case ('-type_out')
        read(arg,*) type_out
     case ('-xscale')
        read(arg,*) xscale
     case ('-yscale')
        read(arg,*) yscale
     case ('-zscale')
        read(arg,*) zscale
     case default
        if (verbose) print '("Unknown option ", A9," ignored")', opt
     end select
  end do
#else
  namelist /input_params/ ndump, old, type_in, type_out, xscale, yscale, zscale
  
  call default_parameters
  
  open(unit=1, file='merco.in', status='old')
  read(1, input_params)
  close(1)
#endif
    
  type_in  = trim(type_in)
  type_out = trim(type_out)
  
  if ((type_in /= 'binary') .and. (type_in /= 'pnetcdf') &
       & .and. (type_in /= 'hdf5') .and. (type_in /= 'phdf5')) then
     print '("Input format not recognized!")'
     stop
  endif
  
  if ((type_out /= 'binary') .and. (type_out /= 'pnetcdf') &
       & .and. (type_out /= 'hdf5') .and. (type_out /= 'phdf5')) then
     print '("Output format not recognized!")'
     stop
  endif

  return
end subroutine get_parameters
!===============================================================================
!> Defaults parameters
!===============================================================================
subroutine default_parameters
  use params
  implicit none
  
  ndump    = 1
  old      = .false.
  type_in  = 'binary'
  type_out = 'pnetcdf'
  xscale   = 1
  yscale   = 1
  zscale   = 1
  
  return
end subroutine default_parameters
