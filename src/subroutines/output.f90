!===============================================================================
!> \file output.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is output subroutines.
!! \details
!! Contains output(), dump_1d_real(), dump_1d_int(), dump_2d_real(), 
!! dump_3d_real()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Pierre Kestener, 
!! Romain Teyssier, Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          04-15-2013 
!! \b last \b modified: 03-27-2015
!<
!===============================================================================
!> \brief
!! This routine produces output of the simulation.
!! \details
!! This routine produces output; several format are available, depending on the
!! configuration chosen to generate Makefile, and on the \c io_type variable.
!! Available formats are:
!!  - Fortran binary
!!  - sequential HDF5
!!  - parallel HDF5
!!  - parallel NetCDF 
!<
!===============================================================================
subroutine output
#if MPI == 1
  use mpi
#endif
  use params
  use variables
  use mpi_var
  use file
#if PHDF5 == 1 || HDF5 == 1
  use hdf5
  use write_h5
#endif
#if PNCDF == 1
  use pnetcdf
#endif
  implicit none

  character(LEN=80) :: filename

#if PHDF5 == 1 || HDF5 == 1
  ! HDF5 variables
  integer(hid_t) :: file_id, fapl_id
#endif

#if PNCDF == 1 && MPI == 1
  ! Parallel NetCDF variables
  integer :: nout, ncid, niid, nrid, piid, prid, xdimid, ydimid, zdimid
  integer :: rhoid, Eid, vxid, vyid, vzid, Bxid, Byid, Bzid, Bxrid, Byrid, Bzrid
  integer, dimension(3) :: sdimid
  integer :: dimz
  integer(dp) :: nbint, nbreal, onedp=1
  integer(kind=MPI_OFFSET_KIND) :: nxtot, nytot, nztot
  integer(kind=MPI_OFFSET_KIND), dimension(3) :: dims, start, count
#endif

  ! Data variables
  integer, parameter :: nb_int=6, nb_real=13
  integer, dimension(nb_int)   :: meta_int
  real(dp), dimension(nb_real) :: meta_real
  integer :: ierr

  !$py begin_statement

#if PHDF5 == 1 || HDF5 == 1
  if (io_type == 'phdf5' .or. io_type == 'hdf5') then

     call H5open_f(ierr)
     if (io_type == 'phdf5') then
        call get_filename(ndump, 0, 'data', filename)
        filename = trim(filename)//'.h5'
        if (mype == 0) then
           call H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr)
           
           meta_int  = (/ nx, ny, nz, nxslice, nyslice, nzslice /)
           meta_real = (/ time, dt, dx, dy, dz, xmin, xmax, ymin, ymax &
                , zmin, zmax, gamma, ciso /)
           
           call dump_1d_int(file_id, "meta_int", meta_int, nb_int)
           call dump_1d_real(file_id, "meta_real", meta_real, nb_real)
        
           call H5Fclose_f(file_id, ierr)
        endif
#if MPI == 1
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
     else
        call get_filename(ndump, mype, 'data', filename)
        filename = trim(filename)//'.h5'
        call H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr)
        
        meta_int  = (/ nx, ny, nz, nxslice, nyslice, nzslice /)
        meta_real = (/ time, dt, dx, dy, dz, xmin, xmax, ymin, ymax &
             , zmin, zmax, gamma, ciso /)
        
        call dump_1d_int(file_id, "meta_int", meta_int, nb_int)
        call dump_1d_real(file_id, "meta_real", meta_real, nb_real)
        
        call H5Fclose_f(file_id, ierr)
     endif

#if PHDF5 == 1 && MPI == 1
  if (io_type == 'phdf5') then
     call H5Pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
     call H5Pset_fapl_mpio_f(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)
     call H5Fopen_f(filename, H5F_ACC_RDWR_F, file_id, ierr, access_prp=fapl_id)
  else
     call H5Fopen_f(filename, H5F_ACC_RDWR_F, file_id, ierr)
  endif
#else
  call H5Fopen_f(filename, H5F_ACC_RDWR_F, file_id, ierr)
#endif
  
  call dump_3d_real(file_id, "rho", uin(:,:,:,ir))
  call dump_3d_real(file_id, "vx", uin(:,:,:,iu))
  call dump_3d_real(file_id, "vy", uin(:,:,:,iv))
  call dump_3d_real(file_id, "vz", uin(:,:,:,iw))
  call dump_3d_real(file_id, "E", uin(:,:,:,ip))
  call dump_3d_real(file_id, "Bx", uin(:,:,:,iA))
  call dump_3d_real(file_id, "By", uin(:,:,:,iB))
  call dump_3d_real(file_id, "Bz", uin(:,:,:,iC))
  call dump_3d_real(file_id, "Bxr", uin(:,:,:,iA+3))
  call dump_3d_real(file_id, "Byr", uin(:,:,:,iB+3))
  call dump_3d_real(file_id, "Bzr", uin(:,:,:,iC+3))

  call H5Fclose_f(file_id, ierr)
#if PHDF5 == 1 && MPI == 1
  if (io_type == 'phdf5') then
     call H5Pclose_f(fapl_id, ierr)
  endif
#endif
  call H5close_f(ierr)
  endif
#endif
#if PNCDF == 1 && MPI == 1
  if (io_type == 'pnetcdf') then
     call get_filename(ndump, 0, 'data', filename)
     filename = trim(filename)//'.nc'

     nout = nfmpi_create(MPI_COMM_WORLD, trim(filename) &
          , ior(NF_WRITE, NF_64BIT_OFFSET), MPI_INFO_NULL, ncid)

     nbint = nb_int; nbreal = nb_real
     nout = nfmpi_def_dim(ncid, "nb_int", nbint, niid)
     nout = nfmpi_def_dim(ncid, "nb_real", nbreal, nrid)
     
     nout = nfmpi_def_var(ncid, "meta_int", NF_INT, 1, (/niid/), piid)
     nout = nfmpi_def_var(ncid, "meta_real", NF_DOUBLE, 1, (/nrid/), prid)

     nxtot = nx+2*nghost
     nytot = ny
     nztot = nz*nxslice*nyslice*nzslice
#if NDIM > 1
     nytot = ny+2*nghost
#endif
#if NDIM == 3
     nztot = (nz+2*nghost)*nxslice*nyslice*nzslice
#endif     

     nout = nfmpi_def_dim(ncid, "nx", nxtot, xdimid)
     nout = nfmpi_def_dim(ncid, "ny", nytot, ydimid)
     nout = nfmpi_def_dim(ncid, "nz", nztot, zdimid)
     sdimid = (/ xdimid, ydimid, zdimid /)
     
     nout = nfmpi_def_var(ncid, "rho", NF_DOUBLE, 3, sdimid, rhoid)
     nout = nfmpi_def_var(ncid, "E", NF_DOUBLE, 3, sdimid, Eid)
     nout = nfmpi_def_var(ncid, "vx", NF_DOUBLE, 3, sdimid, vxid)
     nout = nfmpi_def_var(ncid, "vy", NF_DOUBLE, 3, sdimid, vyid)
     nout = nfmpi_def_var(ncid, "vz", NF_DOUBLE, 3, sdimid, vzid)
     nout = nfmpi_def_var(ncid, "Bx", NF_DOUBLE, 3, sdimid, Bxid)
     nout = nfmpi_def_var(ncid, "By", NF_DOUBLE, 3, sdimid, Byid)
     nout = nfmpi_def_var(ncid, "Bz", NF_DOUBLE, 3, sdimid, Bzid)
     nout = nfmpi_def_var(ncid, "Bxr", NF_DOUBLE, 3, sdimid, Bxrid)
     nout = nfmpi_def_var(ncid, "Byr", NF_DOUBLE, 3, sdimid, Byrid)
     nout = nfmpi_def_var(ncid, "Bzr", NF_DOUBLE, 3, sdimid, Bzrid)
     
     nout = nfmpi_enddef(ncid)

     nout = nfmpi_begin_indep_data(ncid)
     if (mype == 0) then
        meta_int  = (/ nx, ny, nz, nxslice, nyslice, nzslice /)
        meta_real = (/ time, dt, dx, dy, dz, xmin, xmax, ymin, ymax, zmin, zmax &
             , gamma, ciso /)
        
        nout = nfmpi_put_vara_int(ncid, piid, (/onedp/), (/nbint/), meta_int)
        nout = nfmpi_put_vara_double(ncid, prid, (/onedp/), (/nbreal/) &
             , meta_real)
     endif
     nout = nfmpi_end_indep_data(ncid)

     dims = (/ nx+2*nghost, ny, nz /)
#if NDIM > 1
     dims(2) = ny+2*nghost
#endif
#if NDIM == 3
     dims(3) = nz+2*nghost
#endif
     dimz  = mype*dims(3)+1
     start = (/ 1, 1, dimz /)
     count = dims
     
     nout = nfmpi_put_vara_double_all(ncid, rhoid, start, count, uin(:,:,:,ir))
     nout = nfmpi_put_vara_double_all(ncid, Eid, start, count, uin(:,:,:,ip))
     nout = nfmpi_put_vara_double_all(ncid, vxid, start, count, uin(:,:,:,iu))
     nout = nfmpi_put_vara_double_all(ncid, vyid, start, count, uin(:,:,:,iv))
     nout = nfmpi_put_vara_double_all(ncid, vzid, start, count, uin(:,:,:,iw))
     nout = nfmpi_put_vara_double_all(ncid, Bxid, start, count, uin(:,:,:,iA))
     nout = nfmpi_put_vara_double_all(ncid, Byid, start, count, uin(:,:,:,iB))
     nout = nfmpi_put_vara_double_all(ncid, Bzid, start, count, uin(:,:,:,iC))
     nout = nfmpi_put_vara_double_all(ncid, Bxrid, start, count, uin(:,:,:,iA+3))
     nout = nfmpi_put_vara_double_all(ncid, Byrid, start, count, uin(:,:,:,iB+3))
     nout = nfmpi_put_vara_double_all(ncid, Bzrid, start, count, uin(:,:,:,iC+3))
     
     nout = nfmpi_close(ncid)
  endif
#endif

  if (io_type == 'binary') then
     call get_filename(ndump, mype, 'data', filename)
     filename = trim(filename)//'.bin'
     
     open(unit=10, file=trim(filename), status='unknown', form='unformatted')
     write(10) nx, ny, nz, nxslice, nyslice, nzslice
     write(10) time, dt, dx, dy, dz
     write(10) xmin, xmax, ymin, ymax, zmin, zmax
     write(10) gamma, ciso
     write(10) uin
     close(10)
  endif

  return
end subroutine output
