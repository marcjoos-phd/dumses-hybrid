!===============================================================================
!> \file toolbox.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is toolbox modules.
!! \details
!! Contains debugtools
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Romain Teyssier, 
!! Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          12-17-2013 
!! \b last \b modified: 06-25-2015
!<
!===============================================================================
!> Debugging module; tools to debug the application
!===============================================================================
module debugtools
  use file
#if HDF5 == 1 || PHDF5 == 1
  use hdf5
  use write_h5
#endif
#if PNCDF == 1
  use pnetcdf
#endif
#if MPI == 1
  use mpi
#endif

  contains
!===============================================================================
!> Write a debug file in HDF5 or Parallel NetCDF
!===============================================================================
    subroutine write_debug_file(ndump, filename, varname, array &
         , nxarr, nyarr, nzarr)
      use params
      use mpi_var
      implicit none

      integer           :: ndump, nxarr, nyarr, nzarr
      character(LEN=*)  :: varname, filename
      character(LEN=80) :: openedfile
      real(dp), dimension(nxarr,nyarr,nzarr) :: array
      logical           :: fexists
      integer           :: ierr
      
      integer, parameter :: nb_int=6
      integer, dimension(nb_int) :: meta_int

#if HDF5 == 1 || PHDF5 == 1
      integer(hid_t) :: file_id, fapl_id
#endif

#if PNCDF == 1
      integer :: nout, ncid, niid, piid, xdimid, ydimid, zdimid, varid
      integer, dimension(3) :: sdimid
      integer :: dimz
      integer(dp) :: nbint, onedp=1
      integer(kind=MPI_OFFSET_KIND) :: nxtot, nytot, nztot
      integer(kind=MPI_OFFSET_KIND), dimension(3) :: dims, start, count
#endif

      !$py begin_statement

      if (io_type .eq. 'phdf5' .or. io_type .eq. 'hdf5') then
#if HDF5 == 1 || PHDF5 == 1
#if PHDF5 == 1 && MPI == 1
         if (io_type .eq. 'phdf5') then
            call get_filename('output', ndump, 0, trim(filename), openedfile)
            openedfile = trim(openedfile)//'.h5'
            call H5open_f(ierr)

            call H5Pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
            call H5Pset_fapl_mpio_f(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL &
                 , ierr)

            inquire(file=openedfile, exist=fexists)
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            if (fexists) then
               call H5Fopen_f(openedfile, H5F_ACC_RDWR_F, file_id, ierr &
                    , access_prp=fapl_id)
            else
               if (mype .eq. 0) then
                  call H5Fcreate_f(openedfile, H5F_ACC_TRUNC_F, file_id, ierr)
                  meta_int = (/ nxarr, nyarr, nzarr, nxslice, nyslice, nzslice /)
                  call dump_1d_int(file_id, "meta_int", meta_int, nb_int)
                  call H5Fclose_f(file_id, ierr)
               endif
               call MPI_Barrier(MPI_COMM_WORLD, ierr)
               call H5Fopen_f(openedfile, H5F_ACC_RDWR_F, file_id, ierr &
                    , access_prp=fapl_id)
            endif
         else
#endif
            call get_filename('output', ndump, mype, trim(filename), openedfile)
            call H5open_f(ierr)

            inquire(file=openedfile, exist=fexists)
            if (fexists) then
               call H5Fopen_f(openedfile, H5F_ACC_RDWR_F, file_id, ierr)
            else
               call H5Fcreate_f(openedfile, H5F_ACC_TRUNC_F, file_id, ierr)
               meta_int = (/ nxarr, nyarr, nzarr, nxslice, nyslice, nzslice /)
               call dump_1d_int(file_id, "meta_int", meta_int, nb_int)
            endif
#if PHDF5 == 1 && MPI == 1
         endif
#endif

         call dump_3d_debug(file_id, varname, array, nxarr, nyarr, nzarr)

#if MPI == 1
         if (io_type .eq. 'phdf5') then
            call H5Pclose_f(fapl_id, ierr)
         endif
#endif
         call H5Fclose_f(file_id, ierr)
         call H5close_f(ierr)
#endif
      else if (io_type == 'pnetcdf') then
#if PNCDF == 1 && MPI == 1
         call get_filename('output', ndump, 0, trim(filename), openedfile)
         openedfile = trim(openedfile)//'.nc'

         inquire(file=openedfile, exist=fexists)
         call MPI_Barrier(MPI_COMM_WORLD, ierr)
         if (fexists) then
            nout = nfmpi_open(MPI_COMM_SELF, openedfile &
                 & , ior(NF_NOWRITE, NF_64BIT_OFFSET), MPI_INFO_NULL, ncid)
            nout = nfmpi_redef(ncid)
         else
            nout = nfmpi_create(MPI_COMM_WORLD, openedfile &
                 & , ior(NF_WRITE, NF_64BIT_OFFSET), MPI_INFO_NULL, ncid)
         endif

         nbint = nb_int
         nout  = nfmpi_def_dim(ncid, "nb_int", nbint, niid)
         nout  = nfmpi_def_var(ncid, "meta_int", NF_INT, 1, (/niid/), piid)

         nxtot = nxarr
         nytot = nyarr
         nztot = nzarr*nxslice*nyslice*nzslice
         nout  = nfmpi_def_dim(ncid, "nx", nxtot, xdimid)
         nout  = nfmpi_def_dim(ncid, "ny", nytot, ydimid)
         nout  = nfmpi_def_dim(ncid, "nz", nztot, zdimid)
         sdimid = (/ xdimid, ydimid, zdimid /)
         nout  = nfmpi_def_var(ncid, varname, NF_DOUBLE, 3, sdimid, varid)

         nout = nfmpi_enddef(ncid)

         nout = nfmpi_begin_indep_data(ncid)
         if (mype == 0) then
            meta_int = (/ nxarr, nyarr, nzarr, nxslice, nyslice, nzslice /)
            nout = nfmpi_put_vara_int(ncid, piid, (/onedp/), (/nbint/) &
                 & , meta_int)
         endif
         nout = nfmpi_end_indep_data(ncid)

         dims = (/ nxarr, nyarr, nzarr /)
         dimz  = mype*dims(3)+1
         start = (/ 1, 1, dimz /)
         count = dims
     
         nout = nfmpi_put_vara_double_all(ncid, varid, start, count, array)
         
         nout = nfmpi_close(ncid)
#endif
      endif

      return
    end subroutine write_debug_file
!===============================================================================
!> Write a 3D array of reals in HDF5
!===============================================================================
#if PHDF5 == 1 || HDF5 == 1
    subroutine dump_3d_debug(loc_id, dsetname, array, nxarr, nyarr, nzarr)
      use params
      use mpi_var
      use hdf5
      implicit none
    
      integer          :: nxarr, nyarr, nzarr
      integer(hid_t)   :: loc_id
      character(LEN=*) :: dsetname
      real(dp), dimension(nxarr,nyarr,nzarr) :: array
      
      ! HDF5 vars
      integer        :: ierr
      integer(hid_t) :: h5_dspace, h5_dfile
      integer(hid_t) :: h5_dset, dxpl_id
      integer        :: rank=3
      integer(hsize_t), dimension(3) :: dims, dims_file
      integer(hsize_t), dimension(3) :: start, stride, count, blockSize
    
      dims         = (/ nxarr, nyarr, nzarr /)
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
    end subroutine dump_3d_debug
#endif
  
  end module debugtools
