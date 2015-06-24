!===============================================================================
!> \file dumses.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is the main program.
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
!>
!! \mainpage
!! <center>
!! <b>DUMSES-Hybrid</b>\n
!! A versatile MHD grid code for astrophysics\n
!! Copyrights 2013-2015, CEA, <a href="mailto:marc.joos@cea.fr">Marc Joos</a>, Sébastien Fromang, Patrick Hennebelle, Romain Teyssier\n
!! This software is distributed under the CeCILL-A & GNU/GPL licences (see 
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>)
!! </center>
!!
!! Main Contributors to the code:
!!  - <b>Code architecture</b>: <a href="mailto:marc.joos@cea.fr">Marc Joos</a>, Sébastien Fromang, Patrick Hennebelle, Romain Teyssier
!!  - <b>Parallelization  </b>: Sébastien Fromang, Patrick Hennebelle, Romain Teyssier
!!  - <b>Hybridation      </b>: Marc Joos
!!  - <b>MHD              </b>: Sébastien Fromang, Patrick Hennebelle, Romain Teyssier
!!  - <b>Parallel I/O     </b>: Marc Joos, Pierre Kestener
!!
!! \c DUMSES is a 3D MPI/OpenMP & MPI/OpenACC Eulerian second-order Godunov (magneto)hydrodynamic simulation code in cartesian, spherical and cylindrical coordinates.
!!
!! This documentation is under the Creative Commons attribution-NonCommercial-ShareAlike license 4.0 International <http://creativecommons.org/licenses/by-nc-sa/4.0/>.
!<
!===============================================================================
program dumses
  use variables
  use params
  use mpi_var
  implicit none

  real(dp) :: thist, tdump, tspec

  ! Timing variables
  real(dp) :: tcpu0, tcpu1, cputime, elptime
  integer  :: ti0, ti1, tg0, tg1, rate
  real(dp) :: tcpug0, tcpug1, globcputime, globelptime
  logical  :: inter

  ndump = 0
  nspec = 0
  inter = .true.

  call init_parallel

  call init_param
  !$acc data copy(uin, qin, x, y, z, dv, ds) create(dt)
  call init

  !$acc update host(uin)
  if (restart == 0) then
     call output
     ndump = ndump + 1
     call history
  endif
  thist = time; tdump = time; tspec = time
  
  ! CPU time
  call cpu_time(tcpug0)
  ! Elapsed time
  call system_clock(count=tg0, count_rate=rate)

  do
     if (inter) then
        call cpu_time(tcpu0)
     endif
  
     !$py start_timing Timestep
     call compute_dt(dt)
     !$py end_timing Timestep
     if (rhs) then
        !$acc update host(uin)
        uin_old = uin
     endif
     call godunov
     if (rhs) then
        !$py start_timing Source term
        !$acc update host(uin)
        call source_term
        !$acc update device(uin)
        !$py end_timing Source term
     endif
     !$py start_timing FARGO
     if (fargo) call fargo_update
     !$py end_timing FARGO

     !$py start_timing Boundary
     call boundary
     !$py end_timing Boundary
     !$py start_timing Dissipation processes
     if ((nu > zero) .or. (eta > zero)) call dissipation
     !$py end_timing Dissipation processes

     if (mype == 0) print '("Time/dt: ", E14.8, ", ", E14.8)', time, dt
  
     if (inter) then
        if (mype == 0) then
           call cpu_time(tcpu1)
           cputime = tcpu1 - tcpu0
           print '("CPU time required (per timestep): ", E12.6E2)', cputime
        endif
     endif
  
     time = time + dt
     if (debug) then
        !$acc update host(uin)
        call history
        thist = thist + dthist
        call special
        nspec = nspec + 1
        tspec = tspec + dtspec
        call output
        ndump = ndump + 1
        tdump = tdump + dtdump
     else
        if (((time - dt) <= (thist + dthist)) .and. (time > (thist + dthist))) then
           !$acc update host(uin)
           call history
           thist = thist + dthist
        endif
        if (((time - dt) <= (tdump + dtdump)) .and. (time > (tdump + dtdump))) then
           !$acc update host(uin)
           call output
           ndump = ndump + 1
           tdump = tdump + dtdump
        endif
        if (((time - dt) <= (tspec + dtspec)) .and. (time > (tspec + dtspec))) then
           !$acc update host(uin)
           call special
           nspec = nspec + 1
           tspec = tspec + dtspec
        endif
     endif

     if (debug) then
        if (ndump > 5) then
           call deallocate_workspace
           exit
        endif
     else
        if((time-dt) > tlim) then
           call deallocate_workspace
           exit
        endif
     endif
  end do
  !$acc end data
  
  if(mype == 0) then
     ! Elapsed time - end
     call system_clock(count=tg1, count_rate=rate)
     globelptime = real(tg1 - tg0, kind=8)/real(rate, kind=8)
     ! CPU time - end
     call cpu_time(tcpug1)
     globcputime = tcpug1 - tcpug0
     print '(/ 3X,"Elapsed time : ",1PE10.3," s",/ &
          &,3X,"CPU time     : ",1PE10.3," s", /)', globelptime, globcputime
  endif

#if MPI == 1
  call finalize_mpi
#endif

end program dumses
