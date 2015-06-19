!===============================================================================
!> \file metadata.f90
!! \brief
!! \b This is part of a merging program for DUMSES outputs:
!! subroutines to read & write metadata
!! \details
!! Contains get_metadata(), get_metadata_h5(), get_metadata_nc(),
!! get_metadata_bin(), get_1d_real(), get_1d_int(), write_metadata(), 
!! write_metadata_h5(), write_metadata_nc(), write_metadata_bin(), 
!! dump_1d_real(), dump_1d_int()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Pierre Kestener
!! \copyright
!! Copyrights 2014-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          08-28-2014
!! \b last \b modified: 05-28-2015
!<
!===============================================================================
!> Get metadata
!===============================================================================
subroutine get_metadata
  use params
  use variables
  implicit none

  if (type_in == 'hdf5' .or. type_in == 'phdf5') then
     call get_metadata_h5
  else if (type_in == 'pnetcdf') then
     call get_metadata_nc
  else if (type_in == 'binary') then
     call get_metadata_bin
  endif
  
  if (old) then
     meta_int(1) = para_int(4)
     meta_int(2) = para_int(5)
     meta_int(3) = para_int(6)
     meta_int(4) = para_int(7)
     meta_int(5) = para_int(8)
     meta_int(6) = para_int(9)

     meta_real(1)  = para_real(1)
     meta_real(2)  = para_real(2)
     meta_real(3)  = para_real(3)
     meta_real(4)  = para_real(4)
     meta_real(5)  = para_real(5)
     meta_real(6)  = boxSize(1)
     meta_real(7)  = boxSize(2)
     meta_real(8)  = boxSize(3)
     meta_real(9)  = boxSize(4)
     meta_real(10) = boxSize(5)
     meta_real(11) = boxSize(6)
     ! Since gamma and ciso are not written in the old DUMSES output format,
     ! they are set to -1.d0
     meta_real(12) = -1.d0
     meta_real(13) = -1.d0
  endif
  nx      = meta_int(1)
  ny      = meta_int(2)
  nz      = meta_int(3)
  nxslice = meta_int(4)
  nyslice = meta_int(5)
  nzslice = meta_int(6)
  npes    = nxslice*nyslice*nzslice
  n1      = nx + 2*nghost
#if NDIM > 1
  n2      = ny + 2*nghost
#else 
  n2      = 1
#endif
#if NDIM == 3
  n3      = nz + 2*nghost
#else
  n3      = 1
#endif

  if (xscale*yscale*zscale > 1) then
     meta_int(1)  = meta_int(1)*xscale
     meta_int(2)  = meta_int(2)*yscale
     meta_int(3)  = meta_int(3)*zscale
     meta_real(3) = meta_real(3)/xscale
     meta_real(4) = meta_real(4)/yscale
     meta_real(5) = meta_real(5)/zscale
  endif

  n1s     = meta_int(1) + 2*nghost
#if NDIM > 1
  n2s     = meta_int(2) + 2*nghost
#else 
  n2s     = 1
#endif
#if NDIM == 3
  n3s     = meta_int(3) + 2*nghost
#else
  n3s     = 1
#endif

  return
end subroutine get_metadata
!===============================================================================
!> Write metadata
!===============================================================================
subroutine write_metadata
  use params
  use variables
  implicit none

  if (type_out == 'hdf5' .or. type_out == 'phdf5') then
     call write_metadata_h5
  else if (type_out == 'pnetcdf') then
     call write_metadata_nc
  else if (type_out == 'binary') then
     call write_metadata_bin
  endif

  return
end subroutine write_metadata
!===============================================================================
!> Get metadata - HDF5 version
!===============================================================================
subroutine get_metadata_h5
#if HDF5 == 1 || PHDF5 == 1
  use params
  use variables
  use file
  use hdf5
  use read_h5
#if MPI == 1
  use mpi
#endif
  use mpi_var
  implicit none

  character(len=80) :: fileread
  integer(hid_t)    :: filer_id

  if (myrank == 0) then
     if (old) then
        call get_filename(ndump, 0, 'slices', fileread)
        call H5Fopen_f(trim(fileread), H5F_ACC_RDONLY_F, filer_id, ierr)
        call get_1d_int(filer_id, "para_int", para_int, nb_int_old)
        call get_1d_real(filer_id, "para_real", para_real, nb_real_old)
        call get_1d_real(filer_id, "boxSize", boxSize, 6)
        call get_1d_int(filer_id, "para_mpi", para_mpi, nb_mpi_old)
        call H5Fclose_f(filer_id, ierr)
     else
        call get_filename(ndump, 0, 'data', fileread)
        fileread = trim(fileread)//".h5"
        call H5Fopen_f(trim(fileread), H5F_ACC_RDONLY_F, filer_id, ierr)
        call get_1d_int(filer_id, "meta_int", meta_int, nb_int)
        call get_1d_real(filer_id, "meta_real", meta_real, nb_real)
        call H5Fclose_f(filer_id, ierr)
     endif
  endif

#if MPI == 1
  if (old) then
     call MPI_Bcast(para_int, nb_int_old, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(para_real, nb_real_old, MPI_DOUBLE_PRECISION, 0 &
          & , MPI_COMM_WORLD, ierr)
     call MPI_Bcast(boxSize, 6, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  else
     call MPI_Bcast(meta_int, nb_int, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(meta_real, nb_real, MPI_DOUBLE_PRECISION, 0 &
          & , MPI_COMM_WORLD, ierr)
  endif
#endif
#endif

  return
end subroutine get_metadata_h5
!===============================================================================
!> Get metadata - ParallelNetCDF version
!===============================================================================
subroutine get_metadata_nc
#if PNCDF == 1
  use params
  use variables
  use file
  use pnetcdf
#if MPI == 1
  use mpi
#endif
  use mpi_var
  implicit none

  character(len=80) :: fileread
  integer :: nout, ncid
  integer :: ndims, nvars, natts, nulms
  character(LEN=20) :: name
  integer(dp) :: nbint, nbreal, onedp=1, sixdp=6
  integer :: i, type, ndimvar, dimid(3)

  if (old) then
     call get_filename(ndump, 0, 'slices', fileread)
  else
     call get_filename(ndump, 0, 'data', fileread)
     fileread = trim(fileread)//".nc"
  endif
  if (myrank == 0) then
     nout = nfmpi_open(MPI_COMM_SELF, fileread &
          & , ior(NF_NOWRITE, NF_64BIT_OFFSET), MPI_INFO_NULL, ncid)
     nout = nfmpi_inq(ncid, ndims, nvars, natts, nulms)
     
     if (old) then
        do i = 1, nvars
           nout = nfmpi_inq_var(ncid, i, name, type, ndimvar, dimid, natts)
           
           if (trim(name) == "para_int") then
              nbint = nb_int_old
              nout  = nfmpi_get_vara_int_all(ncid, i, (/onedp/), (/nbint/) &
                   & , para_int)
           endif
           if (trim(name) == "para_real") then
              nbreal = nb_real_old
              nout   = nfmpi_get_vara_double_all(ncid, i, (/onedp/), (/nbreal/)&
                   & , para_real)
           endif
           if (trim(name) == "boxSize") then
              nout = nfmpi_get_vara_double_all(ncid, i, (/onedp/), (/sixdp/) &
                   & , boxSize)
           endif
        enddo
     else
        do i = 1, nvars
           nout = nfmpi_inq_var(ncid, i, name, type, ndimvar, dimid, natts)
           
           if (trim(name) == "meta_int") then
              nbint = nb_int
              nout  = nfmpi_get_vara_int_all(ncid, i, (/onedp/), (/nbint/) &
                   & , meta_int)
           endif
           if (trim(name) == "meta_real") then
              nbreal = nb_real_old
              nout   = nfmpi_get_vara_double_all(ncid, i, (/onedp/), (/nbreal/)&
                   & , meta_real)
           endif
        enddo
     endif
     nout = nfmpi_close(ncid)
  endif

#if MPI == 1
  if (old) then
     call MPI_Bcast(para_int, nb_int_old, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(para_real, nb_real_old, MPI_DOUBLE_PRECISION, 0 &
          & , MPI_COMM_WORLD, ierr)
     call MPI_Bcast(boxSize, 6, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  else
     call MPI_Bcast(meta_int, nb_int, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(meta_real, nb_real, MPI_DOUBLE_PRECISION, 0 &
          & , MPI_COMM_WORLD, ierr)
  endif
#endif
#endif

  return
end subroutine get_metadata_nc
!===============================================================================
!> Get metadata - Fortran binary version
!===============================================================================
subroutine get_metadata_bin
  use params
  use variables
  use file
#if MPI == 1
  use mpi
#endif
  use mpi_var
  implicit none
  
  character(len=80) :: fileread
  real(dp) :: time, dt, dx, dy, dz
  integer :: nhist, nspec
  real(dp) :: xmin, xmax, ymin, ymax, zmin, zmax
  real(dp) :: gamma, ciso
  
  if (old) then
     call get_filename(ndump, 0, 'slices', fileread)
     if (myrank == 0) then
        open(unit=10, file=trim(fileread), status='old', action='read' &
             , form='unformatted')
        read(10) time, dt, dx, dy, dz
        read(10) ndump, nhist, nspec
        read(10) nx, ny, nz
        read(10) nxslice, nyslice, nzslice
        
        npes = nxslice*nyslice*nzslice
        n1   = nx + 2*nghost
#if NDIM > 1
        n2   = ny + 2*nghost
#else
        n2   = 1
#endif
#if NDIM == 3
        n3   = nz + 2*nghost
#else
        n3   = 1
#endif
        allocate(x(n1), y(n2), z(n3))
        
        read(10) x, y, z
        
        xmin = x(1); xmax = x(n1)
        ymin = y(1); ymax = y(n2)
        zmin = z(1); zmax = z(n3)
        boxSize   = (/ xmin, xmax, ymin, ymax, zmin, zmax /)
        para_int  = (/ ndump, nhist, nspec, nx, ny, nz, nxslice, nyslice &
             &, nzslice /)
        para_real = (/ time, dt, dx, dy, dz /)
        
        deallocate(x, y, z)
     endif
   
#if MPI == 1
     call MPI_Bcast(para_int, nb_int_old, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(para_real, nb_real_old, MPI_DOUBLE_PRECISION, 0 &
          & , MPI_COMM_WORLD, ierr)
     call MPI_Bcast(boxSize, 6, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif
  else
     call get_filename(ndump, 0, 'data', fileread)
     fileread = trim(fileread)//'.bin'

     open(unit=10, file=trim(fileread), status='old', form='unformatted')
     read(10) nx, ny, nz, nxslice, nyslice, nzslice
     read(10) time, dt, dx, dy, dz
     read(10) xmin, xmax, ymin, ymax, zmin, zmax
     read(10) gamma, ciso
     close(10)

     meta_int  = (/ nx, ny, nz, nxslice, nyslice, nzslice /)
     meta_real = (/ time, dt, dx, dy, dz, xmin, xmax, ymin, ymax &
          & , zmin, zmax, gamma, ciso /)

#if MPI == 1
     call MPI_Bcast(meta_int, nb_int, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(meta_real, nb_real, MPI_DOUBLE_PRECISION, 0 &
          & , MPI_COMM_WORLD, ierr)
#endif
  endif

  return
end subroutine get_metadata_bin
!===============================================================================
!> Write metadata - HDF5 version
!===============================================================================
subroutine write_metadata_h5
#if HDF5 == 1 || PHDF5 == 1
  use params
  use variables
  use file
  use hdf5
  use write_h5
  use mpi_var
  implicit none

  character(len=80) :: filewrite
  integer(hid_t)    :: filew_id
  integer :: nfp, i

  if (type_out == 'phdf5') then
     call get_filename(ndump, 0, 'converted', filewrite)
        filewrite = trim(filewrite)//'.h5'
     if (myrank == 0) then
        call H5Fcreate_f(filewrite, H5F_ACC_TRUNC_F, filew_id, ierr)
   
        call dump_1d_int(filew_id, "meta_int", meta_int, nb_int)
        call dump_1d_real(filew_id, "meta_real", meta_real, nb_real)

        call H5Fclose_f(filew_id, ierr)
      endif
   else if (type_out == 'hdf5') then
#if MPI == 1
      nfp = npes/nproc
      if (mod(npes, nproc) /= 0) nfp = nfp + 1
#else
      nfp = nproc
#endif
        
      do i = 0, nfp-1
#if MPI == 1
         mype = i*nproc + myrank
#else
         mype = i
#endif
         if (mype <= npes-1) then
            call get_filename(ndump, mype, 'converted', filewrite)
            filewrite = trim(filewrite)//'.h5'
            call H5Fcreate_f(filewrite, H5F_ACC_TRUNC_F, filew_id, ierr)
            
            call dump_1d_int(filew_id, "meta_int", meta_int, nb_int)
            call dump_1d_real(filew_id, "meta_real", meta_real, nb_real)
            
            call H5Fclose_f(filew_id, ierr)
         endif
      enddo
   endif
#endif
   
   return
 end subroutine write_metadata_h5
!===============================================================================
!> Write metadata - ParallelNetCDF version
!===============================================================================
 subroutine write_metadata_nc
#if PNCDF == 1
   use params
   use variables
   use file
   use pnetcdf
#if MPI == 1
   use mpi
#endif
   use mpi_var
   implicit none

  character(len=80) :: filewrite
  integer :: nout, ndid
  integer :: ndims, nvars, natts, nulms
  integer :: niid, nrid, bsid, nmid
  integer :: prid, piid, nsid, pmid
  integer(dp) :: nbint, nbreal, nbmpi, onedp=1, sixdp=6
  integer :: ndimvar, dimid(3)

  if (myrank == 0) then
     ! Create file
     call get_filename(ndump, 0, 'converted', filewrite)
        filewrite = trim(filewrite)//'.nc'
     nout = nfmpi_create(MPI_COMM_SELF, filewrite &
          , ior(NF_CLOBBER,NF_64BIT_OFFSET), MPI_INFO_NULL, ndid)
     
     ! Define dimensions
     nbint = nb_int; nbreal = nb_real
     nout = nfmpi_def_dim(ndid, "nb_int", nbint, niid)
     nout = nfmpi_def_dim(ndid, "nb_real", nbreal, nrid)
     
     ! Create variables
     nout = nfmpi_def_var(ndid, "meta_int", NF_INT, 1, (/niid/), piid)
     nout = nfmpi_def_var(ndid, "meta_real", NF_DOUBLE, 1, (/nrid/), prid)
     
     ! End of definitions
     nout = nfmpi_enddef(ndid)
     
     ! Write metadata
     nout = nfmpi_put_vara_int_all(ndid, piid, (/onedp/), (/nbint/), meta_int)
     nout = nfmpi_put_vara_double_all(ndid, prid, (/onedp/), (/nbreal/) &
          &, meta_real)
     
     ! Close file
     nout = nfmpi_close(ndid)
  endif
#endif

  return
end subroutine write_metadata_nc
!===============================================================================
!> Write metadata - Fortran binary version
!===============================================================================
subroutine write_metadata_bin
  use params
  use variables
  use file
  use mpi_var
  implicit none

  character(len=80) :: filewrite
  real(dp) :: time, dt, dx, dy, dz
  real(dp) :: xmin, xmax, ymin, ymax, zmin, zmax
  real(dp) :: gamma, ciso
  integer :: ix, iy, iz
  integer :: i, nfp

#if MPI == 1
  nfp = npes/nproc
  if (mod(npes, nproc) /= 0) nfp = nfp + 1
#else
  nfp = nproc
#endif
        
  do i = 0, nfp-1
#if MPI == 1
     mype = i*nproc + myrank
#else
     mype = i
#endif
     if (mype <= npes-1) then
        time = meta_real(1); gamma = meta_real(12)
        dt   = meta_real(2); ciso  = meta_real(13)
        dx   = meta_real(3); xmin = meta_real(6) ; xmax = meta_real(7) 
        dy   = meta_real(4); ymin = meta_real(8) ; ymax = meta_real(9) 
        dz   = meta_real(5); zmin = meta_real(10); zmax = meta_real(11)

        nx = meta_int(1); nxslice = meta_int(4)
        ny = meta_int(2); nyslice = meta_int(5)
        nz = meta_int(3); nzslice = meta_int(6)
        
        call get_filename(ndump, mype, 'converted', filewrite)
        filewrite = trim(filewrite)//'.bin'
        open(unit=10, file=trim(filewrite), status='unknown', form='unformatted')
        write(10) nx, ny, nz, nxslice, nyslice, nzslice
        write(10) time, dt, dx, dy, dz
        write(10) xmin, xmax, ymin, ymax, zmin, zmax
        write(10) gamma, ciso
        close(10)
     endif
  enddo
     
  return
end subroutine write_metadata_bin
