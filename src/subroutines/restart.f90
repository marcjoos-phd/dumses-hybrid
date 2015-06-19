!===============================================================================
!> \file restart.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is restart subroutines.
!! \details
!! Contains restart_run()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Pierre Kestener, 
!! Romain Teyssier, Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          04-30-2014 
!! \b last \b modified: 03-27-2015
!<
!===============================================================================
!> \brief
!! This routine restarts a simulation from a previous run.
!<
!===============================================================================
subroutine restart_run
#if MPI == 1
  use mpi
#endif
  use params
  use variables
  use mpi_var
  use file
#if PHDF5 == 1 || HDF5 == 1
  use hdf5
  use read_h5
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
  integer :: nout, ncid
  integer :: ndims, nvars, natts, nulms
  character(LEN=20) :: name
  integer(dp) :: nbint, nbreal, onedp=1
  integer :: dimz
  integer(kind=MPI_OFFSET_KIND), dimension(3) :: dims, start, count
  integer :: i, id, type, ndimvar, dimid(3)
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
        call get_filename_r(ndump, 0, 'data', filename)
        filename = trim(filename)//'.h5'

        if (mype == 0) then
           call H5Fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, ierr)
           call get_1d_int(file_id, "meta_int", meta_int, nb_int)
           call get_1d_real(file_id, "meta_real", meta_real, nb_int)

           call H5Fclose_f(file_id, ierr)
        endif
#if MPI == 1
        call MPI_Bcast(meta_int, nb_int, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(meta_real, nb_real, MPI_DOUBLE_PRECISION, 0 &
             , MPI_COMM_WORLD, ierr)
#endif

#if PHDF5 == 1 && MPI == 1
        call H5Pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
        call H5Pset_fapl_mpio_f(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)
        call H5Fopen_f(filename, H5F_ACC_RDWR_F, file_id, ierr &
             , access_prp=fapl_id)
#else
        call H5Fopen_f(filename, H5F_ACC_RDWR_F, file_id, ierr)
#endif
     else
        call get_filename_r(ndump, mype, 'data', filename)
        call H5Fopen_f(filename, H5F_ACC_RDONLY_F, file_id, ierr)
        
        call get_1d_int(file_id, "meta_int", meta_int, nb_int)
        call get_1d_real(file_id, "meta_real", meta_real, nb_int)
     endif

     nx      = meta_int(1); ny      = meta_int(2); nz      = meta_int(3)
     nxslice = meta_int(4); nyslice = meta_int(5); nzslice = meta_int(6)

     time  = meta_real(1);  dt   = meta_real(2)
     dx    = meta_real(3);  dy   = meta_real(4); dz = meta_real(5)
     xmin  = meta_real(6);  xmax = meta_real(7)
     ymin  = meta_real(8);  ymax = meta_real(9)
     zmin  = meta_real(10); zmax = meta_real(11)
     gamma = meta_real(12); ciso = meta_real(13)

     call get_3d_real(file_id, "rho", uin(:,:,:,ir))
     call get_3d_real(file_id, "vx", uin(:,:,:,iu))
     call get_3d_real(file_id, "vy", uin(:,:,:,iv))
     call get_3d_real(file_id, "vz", uin(:,:,:,iw))
     call get_3d_real(file_id, "E", uin(:,:,:,ip))
     call get_3d_real(file_id, "Bx", uin(:,:,:,iA))
     call get_3d_real(file_id, "By", uin(:,:,:,iB))
     call get_3d_real(file_id, "Bz", uin(:,:,:,iC))
     call get_3d_real(file_id, "Bxr", uin(:,:,:,iA+3))
     call get_3d_real(file_id, "Byr", uin(:,:,:,iB+3))
     call get_3d_real(file_id, "Bzr", uin(:,:,:,iC+3))

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
     call get_filename_r(ndump, 0, 'data', filename)
     filename = trim(filename)//'.nc'

     if (mype == 0) then
        nout = nfmpi_open(MPI_COMM_SELF, filename &
             , ior(NF_NOWRITE, NF_64BIT_OFFSET), MPI_INFO_NULL, ncid)
        nout = nfmpi_inq(ncid, ndims, nvars, natts, nulms)
        
        do i = 1, nvars
           nout = nfmpi_inq_var(ncid, i, name, type, ndimvar, dimid, natts)

           if (trim(name) == "meta_int") then
              nbint = nb_int
              nout  = nfmpi_get_vara_int_all(ncid, i, (/onedp/), (/nbint/) &
                   , meta_int)
           endif
           if (trim(name) == "meta_real") then
              nbreal = nb_real
              nout   = nfmpi_get_vara_double_all(ncid, i, (/onedp/), (/nbreal/)&
                   , meta_real)
           endif
        enddo
        nout = nfmpi_close(ncid)
     endif

     call MPI_Bcast(meta_int, nb_int, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(meta_real, nb_real, MPI_DOUBLE_PRECISION, 0 &
          , MPI_COMM_WORLD, ierr)

     nx      = meta_int(1); ny      = meta_int(2); nz      = meta_int(3)
     nxslice = meta_int(4); nyslice = meta_int(5); nzslice = meta_int(6)
     
     time  = meta_real(1);  dt   = meta_real(2)
     dx    = meta_real(3);  dy   = meta_real(4); dz = meta_real(5)
     xmin  = meta_real(6);  xmax = meta_real(7)
     ymin  = meta_real(8);  ymax = meta_real(9)
     zmin  = meta_real(10); zmax = meta_real(11)
     gamma = meta_real(12); ciso = meta_real(13)

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

     nout = nfmpi_open(MPI_COMM_WORLD, filename &
          , ior(NF_NOWRITE, NF_64BIT_OFFSET), MPI_INFO_NULL, ncid)
     nout = nfmpi_inq(ncid, ndims, nvars, natts, nulms)

     do i = 1, nvars
        nout = nfmpi_inq_var(ncid, i, name, type, ndimvar, dimid, natts)
     
        if (trim(name) == 'rho') id = ir
        if (trim(name) == 'vx')  id = iu
        if (trim(name) == 'vy')  id = iv
        if (trim(name) == 'vz')  id = iw
        if (trim(name) == 'E')   id = ip
        if (trim(name) == 'Bx')  id = iA
        if (trim(name) == 'By')  id = iB
        if (trim(name) == 'Bz')  id = iC
        if (trim(name) == 'Bxr') id = iA + 3 
        if (trim(name) == 'Byr') id = iB + 3 
        if (trim(name) == 'Bzr') id = iC + 3 

        if (ndimvar == 3) then
           nout = nfmpi_get_vara_double_all(ncid, i, start, count &
                , uin(:,:,:,id))
        endif
     enddo
     nout = nfmpi_close(ncid)
  endif
#endif
  
  if (io_type == 'binary') then
     call get_filename_r(ndump, mype, 'data', filename)
     filename = trim(filename)//'.bin'

     open(unit=10, file=trim(filename), status='old', form='unformatted')
     read(10) nx, ny, nz, nxslice, nyslice, nzslice
     read(10) time, dt, dx, dy, dz
     read(10) xmin, xmax, ymin, ymax, zmin, zmax
     read(10) gamma, ciso
     read(10) uin
     close(10)
  endif

  return
end subroutine restart_run
