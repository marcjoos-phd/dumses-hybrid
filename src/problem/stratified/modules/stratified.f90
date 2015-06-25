!===============================================================================
!> \file stratified.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is stratified module.
!! \details
!! Contains stratified
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Romain Teyssier, 
!! Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          05-18-2015
!! \b last \b modified: 05-26-2015
!<
!===============================================================================
!> Stratified module; define variables for stratified problem
!===============================================================================
module stratified
  use precision
  implicit none

  logical  :: floor, smooth, addmass, massDiff
  real(dp) :: zfloor, dfloor, d0, mass0

end module stratified
