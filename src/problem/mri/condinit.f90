!===============================================================================
!> \file mri/condinit.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is the initial conditions subroutine for MRI problem.
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
!! \b last \b modified: 06-22-2015
!<
!===============================================================================
subroutine condinit
  use params
  use variables
  use mpi_var
  implicit none

  real(dp) :: p0, B0, amp
  real(dp) :: rvalue
  integer  :: i, j, k
  real(dp) :: d0, beta
  character(len=10) :: type

  integer, dimension(8) :: vtime
  integer, dimension(:), allocatable :: seed
  integer :: size

  namelist /init_params/ d0, beta, type, amp

  if(verbose) print*, 'Dimension (ndim) of the problem: ', ndim

  ! MRI initial conditions
  if (verbose) print*, "> Define initial conditions"
  d0   = 1.d0
  beta = 400.d0
  p0   = d0*ciso*ciso
  amp  = 1.d-2

  open(unit=1, file='input', status='old')
  read(1, init_params)
  close(1)

  if (type /= 'pyl') then
     B0 = sqrt(2.d0*p0/beta)
  else
     B0 = 3.d0*half*sqrt(d0*Omega0**2*(zmax - zmin)**2/beta)
  endif
 
  ! Initialize seed for Random Number Generator
  call date_and_time(values=vtime)
#if PGIFC == 1
  allocate(seed(2))
  seed(1) = (mype + 1)*(vtime(4)*(360000*vtime(5) + 6000*vtime(6) + 100*vtime(7) + vtime(8)))
  call random_seed(put=seed)
#else
   call random_seed(size=size)
  allocate(seed(size))
  seed = (mype + 1)*(vtime(4)*(360000*vtime(5) + 6000*vtime(6) + 100*vtime(7) + vtime(8)))
  call random_seed(put=seed)
#endif

  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2
           ! density initialization
           uin(i,j,k,ir) = d0
  
           if (type /= 'Radial') then
              ! rho*vx
              call random_number(rvalue)
              uin(i,j,k,iu) = uin(i,j,k,ir)*amp*(rvalue - half)*sqrt(p0)
              if (type == "Test") then
                 uin(i,j,k,iu) = 1.d-2*1.d-3*cos(twopi*10.d0*z(k))*cos(twopi*5.d0*x(i))
              endif

              ! rho*vy
              call random_number(rvalue)
              uin(i,j,k,iv) = uin(i,j,k,ir)*amp*(rvalue - half)*sqrt(p0)
              
              ! rho*vz
              call random_number(rvalue)
              uin(i,j,k,iw) = uin(i,j,k,ir)*amp*(rvalue - half)*sqrt(p0)
           endif

           ! Bx
           if (type == 'Radial') then
              uin(i,j,k,iA) = B0*sin(y(j)/(ymax - ymin)*twopi)
           endif

           ! By
           if (type == 'Toroidal') then
              uin(i,j,k,iB) = B0
           endif
  
           ! Bz
           if ((type == "noflux") .or. (type == "Test")) then
              uin(i,j,k,iC)   = B0*sin(twopi*x(i))
              uin(i,j,k,iC+3) = B0*sin(twopi*x(i))
           else if ((type == 'pyl') .or. (type == 'fluxZ')) then
              uin(i,j,k,iC)   = B0
              uin(i,j,k,iC+3) = B0
           endif
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  deallocate(seed)

  return
end subroutine condinit
