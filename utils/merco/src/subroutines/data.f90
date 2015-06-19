!===============================================================================
!> \file data.f90
!! \brief
!! \b This is part of a merging program for DUMSES outputs:
!! subroutines to read & write data
!! \details
!! Contains readwrite_data()
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
!> Read and write data
!===============================================================================
subroutine readwrite_data
  use params
  use variables
  use file
#if HDF5 == 1 || PHDF5 == 1
  use hdf5
  use read_h5
  use write_h5
#endif
#if PNCDF == 1
  use pnetcdf
#endif
#if MPI == 1
  use mpi
#endif
  use mpi_var
  implicit none
  
  character(len=80) :: fileread, filewrite, filename
  integer :: iu1, iu2, ju1, ju2, ku1, ku2, ix, iy, iz
  integer :: is1, is2, js1, js2, ks1, ks2, xpos, ypos, zpos
  real(dp) :: dx, dy, dz
  integer :: nfp, i, is, js, ks

#if HDF5 == 1 || PHDF5 == 1
  integer(hid_t) :: filer_id, filew_id, file_id, fapl_id
#endif

#if PNCDF == 1
  integer :: nout, nciid, ncoid
  integer :: ndims, nvars, natts, nulms
  character(LEN=20) :: name
  integer(kind=MPI_OFFSET_KIND), dimension(3) :: dims, start, count
  integer :: type, ndimvar, dimid(3), id, diminline, ivar
  integer(kind=MPI_OFFSET_KIND) :: nxglob, nyglob, nzglob
  integer :: ndid, xdimid, ydimid, zdimid
  integer :: niid, nrid, bsid, nmid
  integer, dimension(3) :: sdimid
  integer :: prid, piid, nsid, pmid
  integer :: rhoid, Eid, vxid, vyid, vzid, Bxid, Byid, Bzid, Bxrid, Byrid, Bzrid
#endif

  ! Compute state vector sizes
  iu1 = 1-nghost; iu2 = nx+nghost
  ju1 = 1       ; ju2 = 1
  ku1 = 1       ; ku2 = 1
#if NDIM > 1
  ju1 = 1-nghost; ju2 = ny+nghost  
#endif
#if NDIM == 3
  ku1 = 1-nghost; ku2 = nz+nghost  
#endif

  is1 = 1-nghost; is2 = meta_int(1)+nghost
  js1 = 1       ; js2 = 1
  ks1 = 1       ; ks2 = 1
#if NDIM > 1
  js1 = 1-nghost; js2 = meta_int(2)+nghost  
#endif
#if NDIM == 3
  ks1 = 1-nghost; ks2 = meta_int(3)+nghost  
#endif

  if (type_in == 'phdf5' .or. type_out == 'phdf5') then
#if PHDF5 == 1
     if (old) then
        call get_filename(ndump, 0, 'slices', fileread)
        call get_filename(ndump, 0, 'converted', filewrite)
     else
        call get_filename(ndump, 0, 'data', fileread)
        fileread = trim(fileread)//'.h5'
        call get_filename(ndump, 0, 'converted', filewrite)
        filewrite = trim(filewrite)//'.h5'
     endif
     if (type_in == 'phdf5') then
# if MPI == 1
        call H5Pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
        call H5Pset_fapl_mpio_f(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL &
             , ierr)
        call H5Fopen_f(fileread, H5F_ACC_RDWR_F, filer_id, ierr &
             , access_prp=fapl_id)
#else
        call H5Fopen_f(trim(fileread), H5F_ACC_RDWR_F, filer_id, ierr)
#endif
     endif
     if (type_out == 'phdf5') then
# if MPI == 1
        call H5Pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
        call H5Pset_fapl_mpio_f(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL &
             , ierr)
        call H5Fopen_f(filewrite, H5F_ACC_RDWR_F, filew_id, ierr &
             , access_prp=fapl_id)
#else
        call H5Fopen_f(trim(filewrite), H5F_ACC_RDWR_F, filew_id, ierr)
#endif
     endif
#endif
  endif

  if (type_out == 'pnetcdf') then
#if PNCDF == 1
     call get_filename(ndump, 0, 'converted', filewrite)
     filewrite = trim(filewrite)//'.nc'

     nout = nfmpi_open(MPI_COMM_WORLD, filewrite &
          , ior(NF_WRITE, NF_64BIT_OFFSET), MPI_INFO_NULL, ncoid)
     nout = nfmpi_redef(ncoid)

     nxglob = n1s
     nyglob = n2s
     nzglob = n3s*nxslice*nyslice*nzslice
     nout = nfmpi_def_dim(ncoid, "nx", nxglob, xdimid)
     nout = nfmpi_def_dim(ncoid, "ny", nyglob, ydimid)
     nout = nfmpi_def_dim(ncoid, "nz", nzglob, zdimid)
     sdimid = (/ xdimid, ydimid, zdimid /)
     
     nout = nfmpi_def_var(ncoid, "rho", NF_DOUBLE, 3, sdimid, rhoid)
     nout = nfmpi_def_var(ncoid, "E", NF_DOUBLE, 3, sdimid, Eid)
     nout = nfmpi_def_var(ncoid, "vx", NF_DOUBLE, 3, sdimid, vxid)
     nout = nfmpi_def_var(ncoid, "vy", NF_DOUBLE, 3, sdimid, vyid)
     nout = nfmpi_def_var(ncoid, "vz", NF_DOUBLE, 3, sdimid, vzid)
     nout = nfmpi_def_var(ncoid, "Bx", NF_DOUBLE, 3, sdimid, Bxid)
     nout = nfmpi_def_var(ncoid, "By", NF_DOUBLE, 3, sdimid, Byid)
     nout = nfmpi_def_var(ncoid, "Bz", NF_DOUBLE, 3, sdimid, Bzid)
     nout = nfmpi_def_var(ncoid, "Bxr", NF_DOUBLE, 3, sdimid, Bxrid)
     nout = nfmpi_def_var(ncoid, "Byr", NF_DOUBLE, 3, sdimid, Byrid)
     nout = nfmpi_def_var(ncoid, "Bzr", NF_DOUBLE, 3, sdimid, Bzrid)
     
     nout = nfmpi_enddef(ncoid)
#endif  
  endif

#if MPI == 1
  nfp = npes/nproc
  if (mod(npes, nproc) /= 0) nfp = nfp + 1
#else
  nfp = nproc
#endif
    
  allocate(uin(n1,n2,n3,nvar+3))
  allocate(uout(n1s,n2s,n3s,nvar+3))
        
  do i = 0, nfp-1
#if MPI == 1
     mype = i*nproc + myrank
#else
     mype = i
#endif
     
     ! Read data
     if (type_in == 'hdf5' .or. type_in == 'binary') then
        if (old) then
           call get_filename(ndump, mype, 'slices', fileread)
        else
           call get_filename(ndump, mype, 'data', fileread)
           if (type_in == 'hdf5') then
              fileread = trim(fileread)//'.h5'
           else
              fileread = trim(fileread)//'.bin'
           endif
        endif
     else
        if (old) then
           call get_filename(ndump, 0, 'slices', fileread)
        else
           call get_filename(ndump, 0, 'data', fileread)
           if (type_in == 'phdf5') then
              fileread = trim(fileread)//'.h5'
           else
              fileread = trim(fileread)//'.nc'
           endif
        endif
     endif

     ! First, recreate MPI variables and (x,y,z) from metadata
     xposition = mod(mype, nxslice)
     yposition = mod(mype/(nxslice*nzslice), nyslice)
     zposition = mod(mype/nxslice, nzslice)
     
     ! Retrieve uin depending on the input IO type
     if (type_in == 'hdf5') then
#if HDF5 == 1 || PHDF5 == 1
        if (mype <= npes-1) then
           call H5Fopen_f(trim(fileread), H5F_ACC_RDONLY_F, filer_id, ierr)
           if (old) then
              call get_4d_real(filer_id, "uin", uin, n1, n2, n3, nvar+3)
           else
              call get_3d_real(filer_id, "rho", uin(:,:,:,ir))
              call get_3d_real(filer_id, "vx", uin(:,:,:,iu))
              call get_3d_real(filer_id, "vy", uin(:,:,:,iv))
              call get_3d_real(filer_id, "vz", uin(:,:,:,iw))
              call get_3d_real(filer_id, "E", uin(:,:,:,ip))
              call get_3d_real(filer_id, "Bx", uin(:,:,:,iA))
              call get_3d_real(filer_id, "By", uin(:,:,:,iB))
              call get_3d_real(filer_id, "Bz", uin(:,:,:,iC))
              call get_3d_real(filer_id, "Bxr", uin(:,:,:,iA+3))
              call get_3d_real(filer_id, "Byr", uin(:,:,:,iB+3))
              call get_3d_real(filer_id, "Bzr", uin(:,:,:,iC+3))
           endif
           call h5Fclose_f(filer_id, ierr)
        endif
#endif
     else if (type_in == 'phdf5') then
#if PHDF5 == 1
        call get_3d_real(filer_id, "rho", uin(:,:,:,ir))
        if (old) then
           call get_3d_real(filer_id, "rho_vx", uin(:,:,:,iu))
           call get_3d_real(filer_id, "rho_vy", uin(:,:,:,iv))
           call get_3d_real(filer_id, "rho_vz", uin(:,:,:,iw))
        else
           call get_3d_real(filer_id, "vx", uin(:,:,:,iu))
           call get_3d_real(filer_id, "vy", uin(:,:,:,iv))
           call get_3d_real(filer_id, "vz", uin(:,:,:,iw))
        endif
        call get_3d_real(filer_id, "E", uin(:,:,:,ip))
        call get_3d_real(filer_id, "Bx", uin(:,:,:,iA))
        call get_3d_real(filer_id, "By", uin(:,:,:,iB))
        call get_3d_real(filer_id, "Bz", uin(:,:,:,iC))
        call get_3d_real(filer_id, "Bxr", uin(:,:,:,iA+3))
        call get_3d_real(filer_id, "Byr", uin(:,:,:,iB+3))
        call get_3d_real(filer_id, "Bzr", uin(:,:,:,iC+3))
#endif
     else if (type_in == 'pnetcdf') then
#if PNCDF == 1
        nout = nfmpi_open(MPI_COMM_WORLD, fileread &
             , ior(NF_NOWRITE, NF_64BIT_OFFSET), MPI_INFO_NULL, nciid)
        nout = nfmpi_inq(nciid, ndims, nvars, natts, nulms)

        do ivar = 1, nvars
           nout = nfmpi_inq_var(nciid, ivar, name, type, ndimvar, dimid, natts)
           
           if (trim(name) == 'rho') id = ir
           if (trim(name) == 'E') id = ip
           if (old) then
              if (trim(name) == 'rho_vx') id = iu
              if (trim(name) == 'rho_vy') id = iv
              if (trim(name) == 'rho_vz') id = iw
           else
              if (trim(name) == 'vx') id = iu
              if (trim(name) == 'vy') id = iv
              if (trim(name) == 'vz') id = iw
           endif
           if (trim(name) == 'Bx') id = iA
           if (trim(name) == 'By') id = iB
           if (trim(name) == 'Bz') id = iC
           if (trim(name) == 'Bxr') id = iA+3
           if (trim(name) == 'Byr') id = iB+3
           if (trim(name) == 'Bzr') id = iC+3
           
           dims      = (/ n1, n2, n3 /)
           diminline = mype*dims(3) + 1
           start = (/ 1, 1, diminline /)
           count = dims

           if(ndimvar.eq.3) then
              nout = nfmpi_get_vara_double_all(nciid, ivar, start, count &
                   , uin(:,:,:,id))
           endif
        enddo

        nout = nfmpi_close(nciid)
#endif
     else if (type_in == 'binary') then
        if (mype <= npes-1) then
           rewind(unit=10)
           open(unit=10, file=trim(fileread), status='old', form='unformatted', action='read')
           read(10) 
           read(10)
           read(10)
           read(10)
           if (old) then
              read(10)
           endif
           read(10) uin
           close(10)
        endif
     endif

     ! Resize array if needed
     if (xscale*yscale*zscale > 1) then
        do ks = ks1, ks2
           if (ks <= 0) zpos = ks/zscale
           if (ks > 0) zpos = 1 + (ks - 1)/zscale
           do js = js1, js2
              if (js <= 0) ypos = js/yscale
              if (js > 0) ypos = 1 + (js - 1)/yscale
              do is = is1, is2
                 if (is <= 0) xpos = is/xscale
                 if (is > 0) xpos = 1 + (is - 1)/xscale

                 uout(is-is1+1,js-js1+1,ks-ks1+1,:) = uin(xpos-is1+1,ypos-js1+1,zpos-ks1+1,:)
              enddo
           enddo
        enddo
     else
        uout = uin
     endif

     ! Write data
     if (type_out == 'hdf5' .or. type_out == 'binary') then
        call get_filename(ndump, mype, 'converted', filewrite)
        if (type_out == 'hdf5') then
           filewrite = trim(filewrite)//'.h5'
        else
           filewrite = trim(filewrite)//'.bin'
        endif
     else
        call get_filename(ndump, 0, 'converted', filewrite)
        if (type_out == 'phdf5') then
           filewrite = trim(filewrite)//'.h5'
        else
           filewrite = trim(filewrite)//'.nc'
        endif
     endif

     if (type_out == 'hdf5') then
#if HDF5 == 1 || PHDF5 == 1
        if (mype <= npes-1) then
           call H5Fopen_f(filewrite, H5F_ACC_RDWR_F, filew_id, ierr)

           call dump_3d_real(filew_id, "rho", uout(:,:,:,ir), i)
           call dump_3d_real(filew_id, "vx", uout(:,:,:,iu), i)
           call dump_3d_real(filew_id, "vy", uout(:,:,:,iv), i)
           call dump_3d_real(filew_id, "vz", uout(:,:,:,iw), i)
           call dump_3d_real(filew_id, "E", uout(:,:,:,ip), i)
           call dump_3d_real(filew_id, "Bx", uout(:,:,:,iA), i)
           call dump_3d_real(filew_id, "By", uout(:,:,:,iB), i)
           call dump_3d_real(filew_id, "Bz", uout(:,:,:,iC), i)
           call dump_3d_real(filew_id, "Bxr", uout(:,:,:,iA+3), i)
           call dump_3d_real(filew_id, "Byr", uout(:,:,:,iB+3), i)
           call dump_3d_real(filew_id, "Bzr", uout(:,:,:,iC+3), i)

           call H5Fclose_f(filew_id, ierr)
        endif
#endif
     else if (type_out == 'phdf5') then
#if PHDF5 == 1
        call dump_3d_real(filew_id, "rho", uout(:,:,:,ir), i)
        call dump_3d_real(filew_id, "vx", uout(:,:,:,iu), i)
        call dump_3d_real(filew_id, "vy", uout(:,:,:,iv), i)
        call dump_3d_real(filew_id, "vz", uout(:,:,:,iw), i)
        call dump_3d_real(filew_id, "E", uout(:,:,:,ip), i)
        call dump_3d_real(filew_id, "Bx", uout(:,:,:,iA), i)
        call dump_3d_real(filew_id, "By", uout(:,:,:,iB), i)
        call dump_3d_real(filew_id, "Bz", uout(:,:,:,iC), i)
        call dump_3d_real(filew_id, "Bxr", uout(:,:,:,iA+3), i)
        call dump_3d_real(filew_id, "Byr", uout(:,:,:,iB+3), i)
        call dump_3d_real(filew_id, "Bzr", uout(:,:,:,iC+3), i)
#endif
     else if (type_out == 'pnetcdf') then
#if PNCDF == 1
        dims      = (/ n1s, n2s, n3s /)
        diminline = mype*dims(3)+1 
        start = (/ 1, 1, diminline /)
        
        count = dims
        
        nout = nfmpi_put_vara_double_all(ncoid, rhoid, start, count, &
             uout(:,:,:,ir))
        nout = nfmpi_put_vara_double_all(ncoid, Eid, start, count, &
             uout(:,:,:,ip))
        nout = nfmpi_put_vara_double_all(ncoid, vxid, start, count, &
             uout(:,:,:,iu))
        nout = nfmpi_put_vara_double_all(ncoid, vyid, start, count, &
             uout(:,:,:,iv))
        nout = nfmpi_put_vara_double_all(ncoid, vzid, start, count, &
             uout(:,:,:,iw))
        nout = nfmpi_put_vara_double_all(ncoid, Bxid, start, count, &
             uout(:,:,:,iA))
        nout = nfmpi_put_vara_double_all(ncoid, Byid, start, count, &
             uout(:,:,:,iB))
        nout = nfmpi_put_vara_double_all(ncoid, Bzid, start, count, &
             uout(:,:,:,iC))
        nout = nfmpi_put_vara_double_all(ncoid, Bxrid, start, count, &
             uout(:,:,:,iA+3))
        nout = nfmpi_put_vara_double_all(ncoid, Byrid, start, count, &
             uout(:,:,:,iB+3))
        nout = nfmpi_put_vara_double_all(ncoid, Bzrid, start, count, &
             uout(:,:,:,iC+3))
#endif
     else if (type_out == 'binary') then
        if (mype <= npes-1) then
           open(unit=10, file=trim(filewrite), status='old', access='append' &
                , form='unformatted')
           write(10) uout
           close(10)
        endif
     endif
  enddo

  deallocate(uin, uout)

#if PHDF5 == 1
  if (type_in == 'phdf5') call H5Fclose_f(filer_id, ierr)
  if (type_out == 'phdf5') call H5Fclose_f(filew_id, ierr)
#if MPI == 1
  if (type_in == 'phdf5' .or. type_out == 'phdf5') call H5Pclose_f(fapl_id, ierr)
#endif
#endif

#if PNCDF == 1
  if (type_in == 'pnetcdf') nout = nfmpi_close(nciid)
  if (type_out == 'pnetcdf') nout = nfmpi_close(ncoid)
#endif
  
  return
end subroutine readwrite_data
