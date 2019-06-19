!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!--------------------------------------------------------------------------------------------------
!
! 1D model profile for Cascadia region
!
! by Carene Larmat
!--------------------------------------------------------------------------------------------------

subroutine model_1D_EH45TcoldCrust1rq(xmesh,ymesh,zmesh,rho,vp,vs,qmu_atten,qkappa_atten)

!given a GLL point, return velocity model values of EH45TcoldCrust1rq

 use constants, only : NX_TOPO_FILE,NY_TOPO_FILE,MODELING_ATMO
 use generate_databases_par, only: nspec => NSPEC_AB,ibool,HUGEVAL,itopo_bathy
 use create_regions_mesh_ext_par
 implicit none

!GLL point location
 double  precision, intent (in) :: xmesh,ymesh,zmesh

 ! density, Vp and Vs,Attenuation
 real(kind=CUSTOM_REAL),intent(inout) :: vp,vs,rho,qmu_atten,qkappa_atten

  ! local parameters
  real(kind=CUSTOM_REAL) :: x,y,z
  real(kind=CUSTOM_REAL) :: depth
  real(kind=CUSTOM_REAL) :: elevation,distmin

  integer :: layer

  ! converts GLL point location to real
  x = xmesh
  y = ymesh
  z = zmesh

  ! get approximate topography elevation at target coordinates
  distmin = HUGEVAL
  elevation = 0.0
  if (MODELING_ATMO) then 

     call get_topo_bathy_elevation(x,y,elevation, &
              itopo_bathy,NX_TOPO_FILE,NY_TOPO_FILE) 
     !! NEED TO THINK ABOUT THIS - KEEP SIGNATURE TOPOGRAPHY TO DEEP DEPTHS
     depth = elevation - z 

  else
     call get_topo_elevation_free_closest(x,y,elevation,distmin, &
                    nspec,nglob_dummy,ibool,xstore_dummy,ystore_dummy,zstore_dummy, &
                    num_free_surface_faces,free_surface_ispec,free_surface_ijk)

     ! depth in Z-direction in m 
     if (distmin < HUGEVAL) then
       depth = elevation - z
     else
       depth = - z
     endif
  endif


 if (depth>=0.0.and.depth<=80.0) then 
   !vp in m/s
   vp = 263.484848485 + 18.4765512266*depth -0.321645021645*depth**2+ 0.0179292929293*depth**3
   !vs in m/s
   vs = 156.666666667 + 11.731962482*depth -0.222943722944*depth**2 + 0.00133838383838*depth**3
   !rho in kg/m**3
   rho = 1666.11111111 + 4.97847522848*depth -0.0927128427128*depth**2 + 0.000547138047138*depth**3
   !qkappa 
   qkappa_atten = 57822.0
   !qmu
   qmu_atten = 100.0 + 2.5 *depth 
   layer = 1 
 elseif (depth<=1000.0) then 
   !vp in m/s
   vp = 2700.0 
   !vs in m/s 
   vs = 1500.0 
   !rho in kg/m**3
   rho = 2300.0
   !qkappa
   qkappa_atten = 57822.0
   !qmu
   qmu_atten = 600.0 
   layer = 2 
 elseif (depth<=47222.0) then 
   !vp in m/s
   vp = 5394.01992828 + 0.00786710452582*depth -3.43848271848e-08*depth**2+ 4.01677053376e-13*depth**3
   !vs in m/s
   vs = 3114.2381052 + 0.00454234644957*depth -1.98676628832e-08*depth**2+ 2.32121139816e-13*depth**3
   !rho in kg/m3
   rho = 2600.0 -1.60048845473e-16*depth
   !qkappa
   qkappa_atten = 57822.0 -3.57416041531e-15*depth+ 2.35079889496e-19*depth**2
   !qmu
   qmu_atten = 600.0
   layer = 3
 elseif (depth<=85000.0) then 
   !vp in m/s 
   vp = 6600.3194992 + 0.00498453559379*depth +  2.42826929045e-10*depth**2 -1.23902129576e-15*depth**3
   !vs in m/s 
   vs = 3810.34863767 + 0.00289449221658*depth -1.20645549884e-10*depth**2 + 6.16395530201e-16*depth**3
   !rho in kg/m3
   rho = 2900.0
   !qkappa 
   qkappa_atten = 57822.0
   !qmu
   qmu_atten = 600.0
   layer = 4
 else 
   write(*,*) 'Error ',zmesh,depth,elevation
   stop'Need to implement deeper depth'
 endif

!Debug 
  if (layer==1) then 
    write(*,*) zmesh,depth,elevation,layer,vp
  endif


end subroutine model_1D_EH45TcoldCrust1rq
