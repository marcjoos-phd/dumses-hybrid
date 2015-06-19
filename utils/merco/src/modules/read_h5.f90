!===============================================================================
!> \file read_h5.f90
!! \brief
!! \b This is part of a merging program for DUMSES outputs:
!! write_h5 module with subroutines read data in HDF5 format
!! \details
!! Contains get_1d_int(), get_1d_real(), get_2d_real(), get_3d_real(),
!! get_4d_real()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Pierre Kestener
!! \copyright
!! Copyrights 2014-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          29-08-2014
!! \b last \b modified: 26-03-2015
!<
!===============================================================================
module read_h5

contains
!===============================================================================
#if HDF5 == 1 || PHDF5 == 1
!> Read a 1D array of integers in HDF5
!===============================================================================
  subroutine get_1d_int(loc_id, dsetname, array, nx)
    use precision
    use hdf5
    implicit none
  
    integer(hid_t)         :: loc_id
    character(LEN=*)       :: dsetname
    integer                :: nx
    integer, dimension(nx) :: array
    
    ! HDF5 vars
    integer(hid_t) :: h5_dset
    integer(hsize_t), dimension(1) :: dims
    integer        :: ierr
    
    ! Read the array in dataset dsetname
    dims = (/ nx /)
    call H5Dopen_f(loc_id, trim(dsetname), h5_dset, ierr)
    call H5Dread_f(h5_dset, H5T_NATIVE_INTEGER, array, dims, ierr)
    call H5Dclose_f(h5_dset, ierr)
    
    return
  end subroutine get_1d_int
!===============================================================================
!> Read a 1D array of reals in HDF5
!===============================================================================
  subroutine get_1d_real(loc_id, dsetname, array, nx)
    use precision
    use hdf5
    implicit none
    
    integer(hid_t)          :: loc_id
    character(LEN=*)        :: dsetname
    integer                 :: nx
    real(dp), dimension(nx) :: array
    
    ! HDF5 vars
    integer(hid_t) :: h5_dset
    integer(hsize_t), dimension(1) :: dims
    integer :: ierr
    
    ! Read the array in dataset dsetname
    dims = (/ nx /)
    call H5Dopen_f(loc_id, trim(dsetname), h5_dset, ierr)
    call H5Dread_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims, ierr)
    call H5Dclose_f(h5_dset, ierr)
    
    return
  end subroutine get_1d_real
!===============================================================================
!> Read a 3D array of reals in HDF5
!===============================================================================
  subroutine get_3d_real(loc_id, dsetname, array)
    use params
    use mpi_var
    use hdf5
    implicit none
  
    integer(hid_t)                :: loc_id
    character(LEN=*)              :: dsetname
    real(dp), dimension(n1,n2,n3) :: array
    
    ! HDF5 vars
    integer(hid_t) :: h5_dspace, h5_dfile
    integer(hid_t) :: h5_dset, dxpl_id
    integer        :: rank=3
    integer(hsize_t), dimension(3) :: dims, dims_file
    integer(hsize_t), dimension(3) :: start, stride, count, blockSize
  
    dims         = (/ n1, n2, n3 /)
#if PHDF5 == 1
    if (type_in == 'phdf5') then
       ! Get the array in dataset dsetname
       dims_file    = dims
       dims_file(3) = dims_file(3)*nxslice*nyslice*nzslice
       
       call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
       call H5Screate_simple_f(rank, dims_file, h5_dfile, ierr)
       
       start     = (/ 0,0,0 /)
       stride    = (/ 1,1,1 /)
       count     = dims
       blockSize = (/ 1,1,1 /)
       if (mype <= npes-1) then
          call H5Sselect_hyperslab_f(h5_dspace, H5S_SELECT_SET_F, start &
               , count, ierr, stride, blockSize)
       else
          call H5Sselect_none_f(h5_dspace, ierr)
       endif

       start(1)  = 0
       start(2)  = 0
       start(3)  = dims(3)*mype
       stride    = (/ 1,1,1 /)
       count     = dims
       blockSize = (/ 1,1,1 /)
       if (mype <= npes-1) then
          call H5Sselect_hyperslab_f(h5_dfile, H5S_SELECT_SET_F, start &
               , count, ierr, stride, blockSize)
       else
          call H5Sselect_none_f(h5_dfile, ierr)
       endif
       
#if MPI == 1
       call H5Pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
       call H5Pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)
       
#endif
       call H5Dopen_f(loc_id, trim(dsetname), h5_dset, ierr)
       
#if MPI == 1
       call H5Dread_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims &
            , ierr, mem_space_id=h5_dspace, file_space_id=h5_dfile &
            , xfer_prp=dxpl_id)
#else
       call H5Dread_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims &
            , ierr, mem_space_id=h5_dspace, file_space_id=h5_dfile)
#endif

#if MPI == 1
       call H5Pclose_f(dxpl_id, ierr)
#endif
       call H5Dclose_f(h5_dset, ierr)
       call H5Sclose_f(h5_dspace, ierr)
       call H5Sclose_f(h5_dfile, ierr)
    else
       call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
       call H5Dopen_f(loc_id, trim(dsetname), h5_dset, ierr)
       call H5Dread_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims, ierr)
       call H5Dclose_f(h5_dset, ierr)
       call H5Sclose_f(h5_dspace, ierr)
    endif
#else
    call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
    call H5Dopen_f(loc_id, trim(dsetname), H5T_NATIVE_DOUBLE, h5_dset, ierr)
    call H5Dread_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims, ierr)
    call H5Dclose_f(h5_dset, ierr)
    call H5Sclose_f(h5_dspace, ierr)
#endif
    
    return
  end subroutine get_3d_real
!===============================================================================
!> Read a 4D array of reals in HDF5
!===============================================================================
  subroutine get_4d_real(loc_id, dsetname, array, n1, n2, n3, n4)
    use hdf5
    use precision
    implicit none
    
    integer(hid_t)                  :: loc_id
    character(LEN=*)                :: dsetname
    integer                         :: n1, n2, n3, n4
    real(dp), dimension(n1,n2,n3,n4) :: array
  
    ! HDF5 vars
    integer        :: ierr
    integer(hid_t) :: h5_dset
    
    integer(hsize_t), dimension(4) :: dims
    
    dims = (/ n1, n2, n3, n4 /)
    call H5Dopen_f(loc_id, trim(dsetname), h5_dset, ierr)
    call H5Dread_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims, ierr)
    call H5Dclose_f(h5_dset, ierr)
    
    return
  end subroutine get_4d_real
#endif
end module read_h5
