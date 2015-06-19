!===============================================================================
!> \file write_h5.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is write_h5 modules; subroutines to write HDF5 data.
!! \details
!! Contains dump_1d_int(), dump_1d_real(), dump_2d_real(), dump_3d_real(),
!! dump_4d_real()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Romain Teyssier, 
!! Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          12-15-2014 
!! \b last \b modified: 03-03-2015
!<
!===============================================================================
module write_h5

contains
#if PHDF5 == 1 || HDF5 == 1
!===============================================================================
!> Write a 1D array of reals in HDF5
!===============================================================================
  subroutine dump_1d_real(loc_id, dsetname, array, nx)
    use precision
    use hdf5
    implicit none
  
    integer(hid_t)          :: loc_id
    character(LEN=*)        :: dsetname
    integer                 :: nx
    real(dp), dimension(nx) :: array
    
    ! HDF5 vars
    integer        :: ierr
    integer(hid_t) :: h5_dspace
    integer(hid_t) :: h5_dset
    integer        :: rank=1
    integer(hsize_t), dimension(1) :: dims
    
    ! Dump the array in dataset dsetname
    dims = (/ nx /)
    call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
    call H5Dcreate_f(loc_id, trim(dsetname), H5T_NATIVE_DOUBLE, &
         & h5_dspace, h5_dset, ierr)
    call H5Dwrite_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims, ierr)
    call H5Dclose_f(h5_dset, ierr)
    call H5Sclose_f(h5_dspace, ierr)
    
    return
  end subroutine dump_1d_real
!===============================================================================
!> Write a 1D array of integers in HDF5
!===============================================================================
  subroutine dump_1d_int(loc_id, dsetname, array, nx)
    use precision
    use hdf5
    implicit none
  
    integer(hid_t)         :: loc_id
    character(LEN=*)       :: dsetname
    integer                :: nx
    integer, dimension(nx) :: array
    
    ! HDF5 vars
    integer        :: ierr
    integer(hid_t) :: h5_dspace
    integer(hid_t) :: h5_dset
    integer        :: rank=1
    integer(hsize_t), dimension(1) :: dims
    
    ! Dump the array in dataset dsetname
    dims = (/ nx /)
    call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
    call H5Dcreate_f(loc_id, trim(dsetname), H5T_NATIVE_INTEGER, &
         & h5_dspace, h5_dset, ierr)
    call H5Dwrite_f(h5_dset, H5T_NATIVE_INTEGER, array, dims, ierr)
    call H5Dclose_f(h5_dset, ierr)
    call H5Sclose_f(h5_dspace, ierr)
    
    return
  end subroutine dump_1d_int
!===============================================================================
!> Write a 2D array of reals in HDF5
!===============================================================================
  subroutine dump_2d_real(loc_id, dsetname, array, nx, ny)
    use precision
    use hdf5
    implicit none
  
    integer(hid_t)             :: loc_id
    character(LEN=*)           :: dsetname
    integer                    :: nx, ny
    real(dp), dimension(nx,ny) :: array
    
    ! HDF5 vars
    integer        :: ierr
    integer(hid_t) :: h5_dspace
    integer(hid_t) :: h5_dset
    integer        :: rank=2
    integer(hsize_t), dimension(2) :: dims
    
    ! Dump the array in dataset dsetname
    dims = (/ nx, ny /)
    call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
    call H5Dcreate_f(loc_id, trim(dsetname), H5T_NATIVE_DOUBLE, &
         & h5_dspace, h5_dset, ierr)
    call H5Dwrite_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims, ierr)
    call H5Dclose_f(h5_dset, ierr)
    call H5Sclose_f(h5_dspace, ierr)
    
    return
  end subroutine dump_2d_real
!===============================================================================
!> Write a 3D array of reals in HDF5
!===============================================================================
  subroutine dump_3d_real(loc_id, dsetname, array)
    use params
    use mpi_var
    use hdf5
    implicit none
  
    integer(hid_t)                :: loc_id
    character(LEN=*)              :: dsetname
    real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2) :: array
    
    ! HDF5 vars
    integer        :: ierr
    integer(hid_t) :: h5_dspace, h5_dfile
    integer(hid_t) :: h5_dset, dxpl_id
    integer        :: rank=3
    integer(hsize_t), dimension(3) :: dims, dims_file
    integer(hsize_t), dimension(3) :: start, stride, count, blockSize
  
    dims         = (/ nx+2*nghost, ny, nz /)
#if NDIM > 1
    dims(2)      = ny+2*nghost
#endif
#if NDIM == 3
    dims(3)      = nz+2*nghost
#endif
#if PHDF5 == 1
    if (io_type == 'phdf5') then
       ! Dump the array in dataset dsetname
       dims_file    = dims
       dims_file(3) = dims_file(3)*nxslice*nyslice*nzslice
       
       call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
       call H5Screate_simple_f(rank, dims_file, h5_dfile, ierr)
       
       start     = (/ 0,0,0 /)
       stride    = (/ 1,1,1 /)
       count     = dims
       blockSize = (/ 1,1,1 /)
       call H5Sselect_hyperslab_f(h5_dspace, H5S_SELECT_SET_F, start &
            , count, ierr, stride, blockSize)
       
       start(1)  = 0
       start(2)  = 0
       start(3)  = dims(3)*mype
       stride    = (/ 1,1,1 /)
       count     = dims
       blockSize = (/ 1,1,1 /)
       call H5Sselect_hyperslab_f(h5_dfile, H5S_SELECT_SET_F, start &
            , count, ierr, stride, blockSize)
       
       ! Enable parallel collective IO
       call H5Pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
       call H5Pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)
       
       call H5Dcreate_f(loc_id, trim(dsetname), H5T_NATIVE_DOUBLE &
            , h5_dfile, h5_dset, ierr, H5P_DEFAULT_F, H5P_DEFAULT_F &
            , H5P_DEFAULT_F)
       
       call H5Dwrite_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims &
            , ierr, mem_space_id=h5_dspace, file_space_id=h5_dfile &
            , xfer_prp=dxpl_id)
       
       call H5Pclose_f(dxpl_id, ierr)
       call H5Dclose_f(h5_dset, ierr)
       call H5Sclose_f(h5_dspace, ierr)
       call H5Sclose_f(h5_dfile, ierr)
    else
       call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
       call H5Dcreate_f(loc_id, trim(dsetname), H5T_NATIVE_DOUBLE, h5_dspace &
            , h5_dset, ierr)
       call H5Dwrite_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims, ierr)
       call H5Dclose_f(h5_dset, ierr)
       call H5Sclose_f(h5_dspace, ierr)
    endif
#else
    call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
    call H5Dcreate_f(loc_id, trim(dsetname), H5T_NATIVE_DOUBLE, h5_dspace &
         , h5_dset, ierr)
    call H5Dwrite_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims, ierr)
    call H5Dclose_f(h5_dset, ierr)
    call H5Sclose_f(h5_dspace, ierr)
#endif
  
    return
  end subroutine dump_3d_real
#endif
end module write_h5
