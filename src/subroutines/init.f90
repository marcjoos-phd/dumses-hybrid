!===============================================================================
!> \file init.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is initialization subroutines.
!! \details
!! Contains init_param(), init(), read_input(), default_parameters()
!! init_solver(), allocate_workspace(), init_grid()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Romain Teyssier, 
!! Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          04-15-2013 
!! \b last \b modified: 06-24-2015
!<
!===============================================================================
!> Initialize parameters, calling read_input(); allocate arrays, calling allocate_workspace()
!===============================================================================
subroutine init_param
  use variables
  use params
  use mpi_var
  implicit none

  !$py begin_statement

  call read_input(restart, tlim, verbose, debug, bdtypex, bdtypey, bdtypez &
     & , boundary_type, iriemann, iriemann2d, slope_type, courant, fargo &
     & , nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, Omega0, ciso, gamma &
     & , nu, eta, rhs, dtdump, dthist, dtspec, io_type, nxslice, nyslice &
     & , nzslice, nxglob, nyglob, nzglob)

  !calculate state vector sizes
  iu1 = 1-nghost; iu2 = nx+nghost
  ju1 = 1       ; ju2 = 1
  ku1 = 1       ; ku2 = 1
#if NDIM > 1
  ju1 = 1-nghost; ju2 = ny+nghost  
#endif
#if NDIM == 3
  ku1 = 1-nghost; ku2 = nz+nghost  
#endif

  ! calculate flux vector sizes
  if1 = 1; if2 = nx+1
  jf1 = 1; jf2 = 1
  kf1 = 1; kf2 = 1
#if NDIM > 1
  jf1 = 1; jf2 = ny+1
#endif
#if NDIM == 3
  kf1 = 1; kf2 = nz+1
#endif

  if (mype == 0) print*, "Size with ghost cells: ", iu2-iu1+1, ju2-ju1+1 &
       & , ku2-ku1+1

  if (nxslice*nyslice*nzslice > npes) then
     if (verbose) print*, " > Domain decomposition does not fit number of MPI processes!"
     stop
  endif
  
  call allocate_workspace
  !$OMP PARALLEL WORKSHARE
  uin    = zero
  qin    = zero; emfx         = zero
  gravin = zero; emfy         = zero
  flux   = zero; emfz         = zero
  dv     = zero; ds           = zero
  !$OMP END PARALLEL WORKSHARE

  return
end subroutine init_param
!===============================================================================
!> Initialize grid, calling init_grid(); initialize problem, calling condinit()
!===============================================================================
subroutine init
  use const
  use variables
  use params, only: verbose, restart, rhs, iu1, iu2, ju1, ju2, ku1, ku2 &
       & , nvar, fargo
  implicit none

  integer :: c0, c1, rate

  time = zero
  dt   = zero

  call init_grid
  !$acc update host(x, y, z)
  call condinit
  !$acc update device(uin)
  call boundary

  ! User define quantities or initial values
  call user_init

  ! Allocate and initialize uin_old if source terms in the equations
  if (rhs) then
     allocate(uin_old(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3))
     uin_old = zero
  endif

  ! Check if fargo is used, if it is and uin_old not initialized, initialize it
  if (fargo) then
     if (.not. allocated(uin_old)) then
        allocate(uin_old(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3))
        uin_old = zero
     endif
  endif

  if (restart > 0) then
     if (verbose) print "('Restart from output #', I6)", restart
     ndump = restart
     call restart_run
     !$acc update device(uin)
     ndump = ndump + 1
  endif

  return
end subroutine init
!===============================================================================
!> Read input parameters from namelist; default parameters come from default_parameters(); initialize solver, calling init_solver()
!===============================================================================
subroutine read_input(restart, tlim, verbose, debug, bdtypex, bdtypey, bdtypez &
     & , boundary_type, iriemann, iriemann2d, slope_type, courant, fargo &
     & , nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, Omega0, ciso, gamma &
     & , nu, eta, rhs, dtdump, dthist, dtspec, io_type, nxslice, nyslice &
     & , nzslice, nxglob, nyglob, nzglob)
  use mpi_var
  use const
  implicit none

  integer, intent(out)  :: restart
  real(dp), intent(out) :: tlim, courant, xmin, xmax, ymin, ymax, zmin, zmax
  real(dp), intent(out) :: Omega0, ciso, gamma, nu, eta
  real(dp), intent(out) :: dtdump, dthist, dtspec
  logical, intent(out)  :: verbose, debug, fargo, rhs
  character(LEN=20)     :: bdtypex, bdtypey, bdtypez
  character(LEN=20), dimension(ndim), intent(out) :: boundary_type
  character(LEN=10)     :: riemann, riemann2d
  character(LEN=10), intent(out) :: io_type
  integer, intent(out)  :: iriemann, iriemann2d
  integer, intent(out)  :: slope_type
  integer, intent(out)  :: nx, ny, nz
  integer, intent(out)  :: nxslice, nyslice, nzslice, nxglob, nyglob, nzglob
  character(LEN=80)     :: filename
  logical               :: fexist

  
  namelist /start_params/  restart, tlim, verbose, debug
  namelist /scheme_params/ bdtypex, bdtypey, bdtypez, riemann, riemann2d &
                         & , slope_type, courant, fargo
  namelist /model_params/  nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax &
                         & , Omega0, ciso, gamma, nu, eta, rhs
  namelist /output_params/ dtdump, dthist, dtspec, io_type
  namelist /mpi_params/    nxslice, nyslice, nzslice

  call default_parameters(restart, tlim, verbose, debug, bdtypex, bdtypey &
     & , bdtypez, riemann, riemann2d, slope_type, courant, fargo &
     & , nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, Omega0, ciso, gamma &
     & , nu, eta, rhs, dtdump, dthist, dtspec, io_type, nxslice, nyslice &
     & , nzslice)

  open(unit=1, file='input', status='old')
  read(1, start_params)
  read(1, scheme_params)
  read(1, model_params)
  read(1, output_params)
  read(1, mpi_params)
  close(1)
  if (verbose) then
     if(mype /= 0) verbose = .false.
  endif
  if (verbose) print*, "> Reading input"

#if NDIM == 1
  boundary_type = (/ bdtypex /)
#endif
#if NDIM == 2
  boundary_type = (/ bdtypex, bdtypey /)
#endif
#if NDIM == 3
  boundary_type = (/ bdtypex, bdtypey, bdtypez /)
#endif

  nxglob = nx
  nx     = nxglob/nxslice
#if NDIM > 1
  nyglob = ny
  ny     = nyglob/nyslice
#endif
#if NDIM == 3
  nzglob = nz
  nz     = nzglob/nzslice
#endif

  if ((nx == 1 .and. ny == 1) .or. (nx == 1 .and. nz == 1) &
       & .or. (ny == 1 .and. nz == 1)) then
     ndim_act = 1
  else if ((nx == 1 .and. ny > 1 .and. nz > 1) &
       & .or. (nx > 1 .and. ny == 1 .and. nz > 1) &
       & .or. (nx > 1 .and. ny > 1 .and. nz == 1)) then
     ndim_act = 2
  else
     ndim_act = 3
  endif

  call init_solver(riemann, riemann2d, iriemann, iriemann2d)

  if (mype == 0) then
     filename = "run_param.log"
     inquire(file=filename, exist=fexist)
     if (fexist) then
        open(unit=3, file=filename, status="old", position="append")
     else
        open(unit=3, file=filename, status="unknown")
        write(3, "('# ', 34(A13,1X))") "restart", "tlim", "verbose", "debug" &
             & , "bdtypex", "bdtypey", "bdtypez", "riemann", "riemann2d" &
             & , "slope_type", "courant", "fargo", "nx", "ny", "nz" &
             & , "xmin", "xmax", "ymin", "ymax", "zmin", "zmax", "Omega0" &
             & , "ciso", "gamma", "nu", "eta", "rhs", "dtdump", "dthist" &
             & , "dtspec", "io_type", "nxslice", "nyslice", "nzslice"
     endif
     
     write(3, "(2X, I13, 1X, E13.5, 2(1X, L13), 5(1X, A13), 1X, I13, 1X, E13.5 &
          & , 1X, L13, 3(1X, I13), 11(1X, E13.5), 1X, L13, 3(1X, E13.5) &
          & , 1X, A13, 3(1X, I13))") restart, tlim, verbose, debug, bdtypex &
          & , bdtypey, bdtypez, riemann, riemann2d, slope_type, courant, fargo &
          & , nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, Omega0, ciso &
          & , gamma, nu, eta, rhs, dtdump, dthist, dtspec, io_type, nxslice &
          & , nyslice, nzslice
     
     close(3)
  endif

  return
end subroutine read_input
!===============================================================================
!> Define default parameters
!===============================================================================
subroutine default_parameters(restart, tlim, verbose, debug, bdtypex, bdtypey &
     & , bdtypez, riemann, riemann2d, slope_type, courant, fargo &
     & , nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, Omega0, ciso, gamma &
     & , nu, eta, rhs, dtdump, dthist, dtspec, io_type, nxslice, nyslice &
     & , nzslice)
  use const
  implicit none

  integer, intent(out)  :: restart
  real(dp), intent(out) :: tlim, courant, xmin, xmax, ymin, ymax, zmin, zmax
  real(dp), intent(out) :: Omega0, ciso, gamma, nu, eta
  real(dp), intent(out) :: dtdump, dthist, dtspec
  logical, intent(out)  :: verbose, debug, fargo, rhs
  character(LEN=20), intent(out) :: bdtypex, bdtypey, bdtypez
  character(LEN=10), intent(out) :: riemann, riemann2d
  character(LEN=10), intent(out) :: io_type
  integer, intent(out)  :: slope_type
  integer, intent(out)  :: nx, ny, nz
  integer, intent(out)  :: nxslice, nyslice, nzslice

  restart = 0
  tlim    = 0.d0
  verbose = .false.
  debug   = .false.
  bdtypex = 'periodic'
  bdtypey = 'periodic'
  bdtypez = 'periodic'

  riemann   = 'hlld'
  riemann2d = 'hlld'

  slope_type = 2
  courant    = 0.7d0
  fargo      = .false.

  nx = 256
  ny = 256
  nz = 256

  xmin = 0.d0  ; xmax = 1.d0
  ymin = -0.5d0; ymax = 0.5d0
  zmin = 0.d0  ; zmax = 1.d0
  ciso   = 0.d0
  Omega0 = 0.d0
  gamma  = 5.d0/3.d0
  nu     = 0.d0
  eta    = 0.d0
  rhs    = .false.

  tlim    = 1.d0
  dtdump  = -1.d0
  dthist  = -1.d0
  dtspec  = -1.d0
  io_type = "binary"

  nxslice = 1; nyslice = 1; nzslice = 1

  return
end subroutine default_parameters
!===============================================================================
!> Initialize solver; attribute an index to 1D and 2D Riemann solvers
!===============================================================================
subroutine init_solver(riemann, riemann2d, iriemann, iriemann2d)
  use params, only: iroe, illf, ihll, ihlld, iupwind, iacoustic, ihllf, ihlla
  implicit none

  character(LEN=10), intent(in) :: riemann, riemann2d
  integer, intent(out) :: iriemann, iriemann2d

  select case(riemann)
  case('roe')
     iriemann = iroe
  case('llf')
     iriemann = illf
  case('hll')
     iriemann = ihll
  case('hlld')
     iriemann = ihlld
  case('upwind')
     iriemann = iupwind
  case('hydro')
     iriemann = iacoustic
  end select

  select case(riemann2d)
  case('roe')
     iriemann2d = iroe
  case('llf')
     iriemann2d = illf
  case('hll')     
     iriemann2d = ihll
  case('hlld')
     iriemann2d = ihlld
  case('upwind')
     iriemann2d = iupwind
  case('hydro')
     iriemann2d = iacoustic
  case('hllf')
     iriemann2d = ihllf
  case('hlla')
     iriemann2d = ihlla
  end select

  return
end subroutine init_solver
!===============================================================================
!> Allocate arrays
!===============================================================================
subroutine allocate_workspace
  use variables
  use params
  implicit none

  !$py begin_statement

  allocate(uin(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3))
  allocate(qin(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar))
  allocate(gravin(iu1:iu2,ju1:ju2,ku1:ku2,1:ndim))
  allocate(flux(if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim))
  allocate(emfx(iu1:iu2,ju1:ju2,ku1:ku2))
  allocate(emfy(iu1:iu2,ju1:ju2,ku1:ku2))
  allocate(emfz(iu1:iu2,ju1:ju2,ku1:ku2))
  allocate(x(iu1:iu2))
  allocate(y(ju1:ju2))
  allocate(z(ku1:ku2))
  allocate(dv(iu1:iu2,ju1:ju2,ku1:ku2))       
  allocate(ds(iu1:iu2,ju1:ju2,ku1:ku2,1:ndim))
#if NDIM == 3
  allocate(Ex(iu1:iu2,ju1:ju2,ku1:ku2))
  allocate(Ey(iu1:iu2,ju1:ju2,ku1:ku2))
#endif
#if NDIM > 1
  allocate(Ez(iu1:iu2,ju1:ju2,ku1:ku2))
#endif
  allocate(bfc(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3))
  allocate(dq(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:nvar,1:ndim))
  allocate(dbfc(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3,1:ndim))
  allocate(qm(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim))
  allocate(qp(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim))
  allocate(qRT(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3))
  allocate(qRB(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3))
  allocate(qLT(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3))
  allocate(qLB(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3))
  allocate(fgodunov(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim))
  allocate(fgodunov_pre(iu1:iu2,ju1:ju2,ku1:ku2,1:ndim))

  return
end subroutine allocate_workspace
!===============================================================================
!> Initialize grid (resolution, coordinates areas and volumes)
!===============================================================================
subroutine init_grid
  use params
  use variables, only: x, y, z, ds, dv
  use mpi_var
#if OACC == 1
  use openacc
#endif
  implicit none

  integer :: i, j, k, idim
  integer :: ilo, ihi, jlo, jhi, klo, khi

#if MPI == 1
  call grid_structure
#endif

  dx = (xmax-xmin)/nxglob
  dy = (ymax-ymin)/nyglob
  dz = (zmax-zmin)/nzglob

  if (mype == 0) then
     if(verbose) print*, "> Initialize grid"
     print*, 'Mesh size:'
     print*, 'dx: ', dx
#if NDIM > 1
     print*, 'dy: ', dy
#endif
#if NDIM == 3
     print*, 'dz: ', dz
#endif
  endif

  !$acc kernels loop
#if MPI == 1
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do i = iu1, iu2
     x(i) = xmin + dx*nx*xposition + dx/2 + (i - 1)*dx
  enddo
  !$OMP END PARALLEL DO
#else
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do i = iu1, iu2
     x(i) = xmin + dx/2 + (i - 1)*dx
  enddo
  !$OMP END PARALLEL DO
#endif
#if NDIM > 1
  !$acc kernels loop
#if MPI == 1
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do j = ju1, ju2
     y(j) = ymin + dy*ny*yposition + dy/2 + (j - 1)*dy
  enddo
  !$OMP END PARALLEL DO
#else
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do j = ju1, ju2
     y(j) = ymin + dy/2 + (j - 1)*dy
  enddo
  !$OMP END PARALLEL DO
#endif
#endif
#if NDIM == 3
  !$acc kernels loop
#if MPI == 1
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do k = ku1, ku2
     z(k) = zmin + dz*nz*zposition + dz/2 + (k - 1)*dz
  enddo
  !$OMP END PARALLEL DO
#else
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do k = ku1, ku2
     z(k) = zmin + dz/2 + (k - 1)*dz
  enddo
  !$OMP END PARALLEL DO
#endif
#endif

  ilo = min(1,iu1+1); ihi = max(1,iu2-1)
  jlo = min(1,ju1+1); jhi = max(1,ju2-1)
  klo = min(1,ku1+1); khi = max(1,ku2-1)

#if GEOM == CARTESIAN
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do k = klo,khi
     do j = jlo,jhi
        do i = ilo,ihi
           dv(i,j,k) = dx*dy*dz
           do idim = 1,ndim
              if (idim == 1) ds(i,j,k,idim) = dy*dz
              if (idim == 2) ds(i,j,k,idim) = dx*dz
              if (idim == 3) ds(i,j,k,idim) = dx*dy
           enddo
        end do
     end do
  end do
  !$OMP END PARALLEL DO
#endif
#if GEOM == CYLINDRICAL
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do k = klo,khi
     do j = jlo,jhi
        do i = ilo,ihi
           dv(i,j,k) = dx*x(i)*dy*dz
           do idim = 1,ndim
              if (idim == 1) ds(i,j,k,idim) = half*(x(i) + x(i-1))*dy*dz
              if (idim == 2) ds(i,j,k,idim) = dx*dz
              if (idim == 3) ds(i,j,k,idim) = x(i)*dx*dy
           enddo
        end do
     end do
  end do
  !$OMP END PARALLEL DO
#endif
#if GEOM == SPHERICAL
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do k = klo,khi
     do j = jlo,jhi
        do i = ilo,ihi
           dv(i,j,k) = dx*x(i)*dy**x(i)*sin(y(j))*dz
           do idim = 1,ndim
              if (idim == 1) ds(i,j,k,idim) = half*(x(i) + x(i-1))*dy*dz
              if (idim == 2) ds(i,j,k,idim) = half*(x(i) + x(i-1))*dx*dz
              if (idim == 3) ds(i,j,k,idim) = x(i)*sin(z(k))*dx*dy
           enddo
        end do
     end do
  end do
  !$OMP END PARALLEL DO
#endif

  return
end subroutine init_grid
