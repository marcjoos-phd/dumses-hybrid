!===============================================================================
!> \file commons.f90
!! \brief
!! \b This is part of a merging program for DUMSES outputs:
!! Commons module with variables declaration
!! \details
!! Contains precision(), const(), params(), variables(), mpi_var()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang
!! \copyright
!! Copyrights 2014-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          08-27-2014
!! \b last \b modified: 05-28-2015
!<
!===============================================================================
!> Precision module; define double precision
!===============================================================================
module precision

  integer, parameter :: dp=kind(1.0D0)

end module precision
!===============================================================================
!> Constant module; define constants of the code
!===============================================================================
module const
  use precision

#ifndef NDIM
  integer, parameter :: ndim=3
#else
  integer, parameter :: ndim=NDIM
#endif

end module const
!===============================================================================
!> Parameters modules; parameters of the problem
!===============================================================================
module params
  use precision
  use const

  ! Verbosity
  logical :: verbose=.false.

  ! Variable parameters
  integer, parameter :: nvar=8
  integer, parameter :: ir=1
  integer, parameter :: iu=2
  integer, parameter :: iv=3
  integer, parameter :: iw=4
  integer, parameter :: ip=5
  integer, parameter :: iA=6
  integer, parameter :: iB=7
  integer, parameter :: iC=8
  integer, parameter :: nghost=3

  ! Model parameters
  integer :: nx, ny, nz, n1, n2, n3
  integer :: n1s, n2s, n3s

  ! MPI parameters
  integer :: nxslice, nyslice, nzslice
  integer :: mype, npes
  integer :: xposition, yposition, zposition
  integer :: xleft, xright, yleft, yright, zleft, zright
  
  ! Output parameters
  integer, parameter :: nb_int=6, nb_real=13
  integer, parameter :: nb_int_old=9, nb_mpi_old=9, nb_real_old=5

  ! Input parameters
  integer :: ndump
  character(len=80) :: type_in, type_out
  logical :: old
  integer :: xscale, yscale, zscale

end module params
!===============================================================================
!> Variables module; variables of the problem
!===============================================================================
module variables
  use params

  real(dp), dimension(:,:,:,:), allocatable :: uin
  real(dp), dimension(:,:,:,:), allocatable :: uout
  ! Metadata
  integer, dimension(nb_int)   :: meta_int
  real(dp), dimension(nb_real) :: meta_real

  ! Old metadata format
  real(dp), dimension(:), allocatable :: x, y, z
  integer, dimension(nb_int_old)      :: para_int
  integer, dimension(nb_mpi_old)      :: para_mpi
  real(dp), dimension(nb_real_old)    :: para_real
  real(dp), dimension(6)              :: boxSize

end module variables
!===============================================================================
!> MPI variables module
!===============================================================================
module mpi_var
  use precision

  integer :: nproc=1, myrank=0
  integer :: ierr

end module mpi_var
