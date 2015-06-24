!===============================================================================
!> \file commons.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is common modules.
!! \details
!! Contains precision, const, params and variables
!! WARNING: You should never modify this file for a given problem; if you want
!! to do so, create another module instead in $PROBLEM_PATH/modules/.
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
!> Precision module; define single and double precisions
!===============================================================================
module precision

  integer, parameter :: sp=kind(1.0E0) !< single precision
#ifndef NPRE                            
  integer, parameter :: dp=kind(1.0D0) !< default double precision
#else                                   
#if NPRE==4                            
  integer, parameter :: dp=kind(1.0E0) !< double precision if NPRE==4
#else                                   
  integer, parameter :: dp=kind(1.0D0) !< double precision is NPRE!=4
#endif
#endif

end module precision
!===============================================================================
!> Constant module; define constants of the code
!===============================================================================
module const
  use precision

  real(dp), parameter :: bigreal = 1.0d+30               !< large real number
  real(dp), parameter :: zero    = 0.0d0                 !< 0
  real(dp), parameter :: one     = 1.0d0                 !< 1
  real(dp), parameter :: two     = 2.0d0                 !< 2
  real(dp), parameter :: three   = 3.0d0                 !< 3
  real(dp), parameter :: four    = 4.0d0                 !< 4
  real(dp), parameter :: two3rd  = 0.66666666666666667d0 !< 2/3
  real(dp), parameter :: half    = 0.5d0                 !< 1/2
  real(dp), parameter :: third   = 0.33333333333333333d0 !< 1/3
  real(dp), parameter :: forth   = 0.25d0                !< 1/4
  real(dp), parameter :: sixth   = 0.16666666666666667d0 !< 1/6
  real(dp), parameter :: pi      = 2.d0*asin(1.d0)       !< \f$ \pi \f$
  real(dp), parameter :: twopi   = 4.d0*asin(1.d0)       !< \f$ 2 \pi \f$

  ! Number of dimensions
#ifndef NDIM
  integer, parameter :: ndim=3    !< number of dimensions of the problem
#else      
  integer, parameter :: ndim=NDIM !< number of dimensions of the problem
#endif     
  integer :: ndim_act

end module const
!===============================================================================
!> Parameters modules; parameters of the problem
!===============================================================================
module params
  use precision
  use const

  ! Number of independant variables
  integer, parameter :: nvar=8 !< number of variables
  integer, parameter :: ir=1   !< density index
  integer, parameter :: iu=2   !< velocity index (first component) 
  integer, parameter :: iv=3   !< velocity index (second component)
  integer, parameter :: iw=4   !< velocity index (third component) 
  integer, parameter :: ip=5   !< pressure index 
  integer, parameter :: iA=6   !< magnetic field index (first component) 
  integer, parameter :: iB=7   !< magnetic field index (second component)
  integer, parameter :: iC=8   !< magnetic field index (third component) 

  ! Size of hydro kernel
  integer, parameter :: nghost=3 !< number of ghost cells
  integer :: iu1                 !< left edge index, x-direction 
  integer :: iu2                 !< right edge index, x-direction
  integer :: ju1                 !< left edge index, y-direction 
  integer :: ju2                 !< right edge index, y-direction
  integer :: ku1                 !< left edge index, z-direction 
  integer :: ku2                 !< right edge index, z-direction
  integer :: if1                 !< flux left edge index, x-direction 
  integer :: if2                 !< flux right edge index, x-direction
  integer :: jf1                 !< flux left edge index, y-direction 
  integer :: jf2                 !< flux right edge index, y-direction
  integer :: kf1                 !< flux left edge index, z-direction 
  integer :: kf2                 !< flux right edge index, z-direction

  ! Start parameters
  integer  :: restart         !< output index for restart
  real(dp) :: tlim=1000.      !< limit time
  logical  :: verbose=.false. !< verbosity
  logical  :: debug=.false.   !< turn-on debugging features

  ! Scheme parameters
  character(LEN=20) :: bdtypex='periodic' !< boundary conditions, x-direction
  character(LEN=20) :: bdtypey='periodic' !< boundary conditions, y-direction
  character(LEN=20) :: bdtypez='periodic' !< boundary conditions, z-direction
  character(LEN=20), dimension(ndim) :: boundary_type !< boundary conditions
  character(LEN=10) :: riemann  ='hlld'   !< 1D Riemann solver
  character(LEN=10) :: riemann2d='hlld'   !< 2D Riemann solver
  integer :: iriemann                     !< 1D Riemann solver index
  integer :: iriemann2d                   !< 2D Riemann solver index
  integer, parameter :: iroe=1            !< Roe solver index
  integer, parameter :: illf=2            !< Lax-Friedrich solver index
  integer, parameter :: ihll=3            !< HLL solver index
  integer, parameter :: ihlld=4           !< HLLD solver index
  integer, parameter :: iupwind=5         !< Upwind solver index
  integer, parameter :: iacoustic=6       !< Hydro solver index
  integer, parameter :: ihllf=7           !< HLLF solver index
  integer, parameter :: ihlla=8           !< HLLA solver index
  integer  :: slope_type=2                !< slope limiter type for hydro & MHD
  real(dp) :: courant                     !< Courant factor
  logical  :: fargo                       !< Use the FARGO algorithm
  real(dp), parameter :: smallr=1.d-10    !< Dimensioned small constant
  real(dp), parameter ::  smallc=1.d-10   !< Dimensioned small constant

  ! Model parameters
  integer  :: nx=1              !< number of cells in a subdomain, x-direction
  integer  :: ny=1              !< number of cells in a subdomain, y-direction
  integer  :: nz=1              !< number of cells in a subdomain, z-direction
  integer  :: nxglob=1          !< total number of cells, x-direction
  integer  :: nyglob=1          !< total number of cells, y-direction
  integer  :: nzglob=1          !< total number of cells, z-direction
  real(dp) :: nu=0.d0           !< viscosity
  real(dp) :: eta=0.d0          !< resistivity
  real(dp) :: xmin              !< minimum x coordinate
  real(dp) :: xmax              !< maximum x coordinate
  real(dp) :: ymin              !< minimum y coordinate
  real(dp) :: ymax              !< maximum y coordinate
  real(dp) :: zmin              !< minimum z coordinate
  real(dp) :: zmax              !< maximum z coordinate
  real(dp) :: Omega0=0.d0       !< angular velocity
  real(dp) :: ciso=0.d0         !< isothermal sound speed
  real(dp) :: gamma=5.d0/3.d0   !< ratio of specific heats \f$ \gamma \f$
  logical  :: rhs=.false.       !< .true. if there is a non-zero source term
   
  ! Mesh parameters
  real(dp) :: dx !< resolution in x-direction
  real(dp) :: dy !< resolution in y-direction
  real(dp) :: dz !< resolution in z-direction

  ! Output parameters
  real(dp) :: dtdump           !< elapsed time between two outputs
  real(dp) :: dthist           !< elapsed time between two history
  real(dp) :: dtspec           !< elapsed time between two special
  character(LEN=10) :: io_type !< output format type

  ! MPI parameters
  integer :: nxslice !< MPI domain decomposition: # of process in x direction
  integer :: nyslice !< MPI domain decomposition: # of process in y direction
  integer :: nzslice !< MPI domain decomposition: # of process in z direction

end module params
!===============================================================================
!> Variables module; variables of the problem
!===============================================================================
module variables
  use precision

  real(dp), dimension(:,:,:,:), allocatable   :: uin          !< MHD variables array
  real(dp), dimension(:,:,:,:), allocatable   :: uin_old      !< MHD variables array at the previous timestep
  real(dp), dimension(:,:,:,:), allocatable   :: qin            !< primitive array
  real(dp), dimension(:,:,:,:), allocatable   :: gravin       !< grav. array
  real(dp), dimension(:,:,:,:,:), allocatable :: flux         !< flux array
  real(dp), dimension(:,:,:), allocatable     :: emfx         !< EMF array, x-dir.
  real(dp), dimension(:,:,:), allocatable     :: emfy         !< EMF array, y-dir.
  real(dp), dimension(:,:,:), allocatable     :: emfz         !< EMF array, z-dir.
  real(dp), dimension(:), allocatable         :: x            !< x-coordinates array
  real(dp), dimension(:), allocatable         :: y            !< y-coordinates array
  real(dp), dimension(:), allocatable         :: z            !< z-coordinates array
#if NDIM == 3
  real(dp), dimension(:,:,:), allocatable     :: Ex           !< Elec. current, x-dir.
  real(dp), dimension(:,:,:), allocatable     :: Ey           !< Elec. current, y-dir.
#endif
#if NDIM > 1
  real(dp), dimension(:,:,:), allocatable     :: Ez           !< Elec. current, z-dir.
#endif
  real(dp), allocatable, dimension(:,:,:,:)   :: bfc          !< face-centered B field
  real(dp), allocatable, dimension(:,:,:,:,:) :: dq           !< slopes
  real(dp), allocatable, dimension(:,:,:,:,:) :: dbfc         !< magnetic slopes 
  real(dp), allocatable, dimension(:,:,:,:,:) :: qm           !< face averaged state
  real(dp), allocatable, dimension(:,:,:,:,:) :: qp           !< face averaged state
  real(dp), allocatable, dimension(:,:,:,:,:) :: qRT          !< face averaged corner state
  real(dp), allocatable, dimension(:,:,:,:,:) :: qRB          !< face averaged corner state
  real(dp), allocatable, dimension(:,:,:,:,:) :: qLT          !< face averaged corner state
  real(dp), allocatable, dimension(:,:,:,:,:) :: qLB          !< face averaged corner state
  real(dp), allocatable, dimension(:,:,:,:,:) :: fgodunov     !< Godunov flux
  real(dp), allocatable, dimension(:,:,:,:)   :: fgodunov_pre !< 

  ! Curvilinear coordinates parameters
  real(dp), dimension(:,:,:), allocatable   :: dv       !< cells volume array
  real(dp), dimension(:,:,:,:), allocatable :: ds       !< cells area array

  ! Time variables
  real(dp) :: dt   !< timestep
  real(dp) :: time !< time

  ! Output variables
  integer :: ndump !< # of the current output
  integer :: nspec !< # of the current special

end module variables
!===============================================================================
!> MPI variables module
!===============================================================================
module mpi_var
  use precision

  integer :: mype=0    !< # of the current MPI process
  integer :: npes=1    !< # of MPI processes
  integer :: xposition !< x-position of the current process in the MPI grid
  integer :: yposition !< y-position of the current process in the MPI grid
  integer :: zposition !< z-position of the current process in the MPI grid
  integer :: xleft     !< left neighbor of the current process, x-direction 
  integer :: xright    !< right neighbor of the current process, x-direction
  integer :: yleft     !< left neighbor of the current process, y-direction 
  integer :: yright    !< right neighbor of the current process, y-direction
  integer :: zleft     !< left neighbor of the current process, z-direction 
  integer :: zright    !< right neighbor of the current process, z-direction

end module mpi_var
!===============================================================================
!> OpenACC kernels configuration parameters
!===============================================================================
module oacc_params

  integer, parameter :: blockx_solver=32
  integer, parameter :: blocky_solver=8
  integer, parameter :: blockx_solver_mag=16
  integer, parameter :: blocky_solver_mag=8
  integer, parameter :: blockx_trace=32
  integer, parameter :: blocky_trace=8
  integer, parameter :: blockx_update=32
  integer, parameter :: blocky_update=4
  integer, parameter :: blockx_slope=16
  integer, parameter :: blocky_slope=32
  integer, parameter :: blockx_primitive=32
  integer, parameter :: blocky_primitive=4

end module oacc_params
