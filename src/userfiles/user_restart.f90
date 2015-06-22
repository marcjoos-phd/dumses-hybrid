!===============================================================================
!> \file user_restart.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is user's restart subroutines.
!! \details
!! Contains user_restart()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Romain Teyssier, 
!! Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          06-22-2015
!! \b last \b modified: 06-22-2015
!<
!===============================================================================
!> User customization at restart
!===============================================================================
subroutine user_restart
  use precision
  use params, only: nvar
  use variables, only: uin
  implicit none

  return
end subroutine user_restart
