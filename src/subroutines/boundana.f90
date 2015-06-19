!===============================================================================
!> \file boundana.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This subroutine contains the analytical boundary conditions. This version
!! contains dummy subroutines only.
!! \details
!! Contains xinner_ana(), xouter_ana(), yinner_ana(), youter_ana(),
!! zinner_ana(), zouter_ana()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Romain Teyssier, 
!! Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          07-05-2015 
!! \b last \b modified: 
!<
!===============================================================================
!> Compute boundary conditions in the x-direction, inner edge
!===============================================================================
subroutine xinner_ana
  implicit none

  return
end subroutine xinner_ana
!===============================================================================
!> Compute boundary conditions in the x-direction, outer edge
!===============================================================================
subroutine xouter_ana
  implicit none

  return
end subroutine xouter_ana
!===============================================================================
!> Compute boundary conditions in the y-direction, inner edge
!===============================================================================
subroutine yinner_ana
#if NDIM>1
  implicit none

#endif
  return
end subroutine yinner_ana
!===============================================================================
!> Compute boundary conditions in the y-direction, outer edge
!===============================================================================
subroutine youter_ana
#if NDIM>1
  implicit none

#endif
  return
end subroutine youter_ana
!===============================================================================
!> Compute boundary conditions in the z-direction, inner edge
!===============================================================================
subroutine zinner_ana
#if NDIM>1
  implicit none

#endif
  return
end subroutine zinner_ana
!===============================================================================
!> Compute boundary conditions in the z-direction, inner edge
!===============================================================================
subroutine zouter_ana
#if NDIM>1
  implicit none

#endif
  return
end subroutine zouter_ana
