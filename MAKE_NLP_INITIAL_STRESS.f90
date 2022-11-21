!    Copyright (C) 2012 The SPEED FOUNDATION
!    Author: Ilario Mazzieri
!
!    This file is part of SPEED.
!
!    SPEED is free software; you can redistribute it and/or modify it
!    under the terms of the GNU Affero General Public License as
!    published by the Free Software Foundation, either version 3 of the
!    License, or (at your option) any later version.
!
!    SPEED is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Affero General Public License for more details.
!
!    You should have received a copy of the GNU Affero General Public License
!    along with SPEED.  If not, see <http://www.gnu.org/licenses/>.


!> @brief NonLinear - Plasticity Material based on MPII - Iwan springs Model (ref. Fabian's Noah)
!
!! @author Srihari Sangaraju
!> @date August, 2022 
!> @version 1.0
!
! stress = initial stress (shear stresses are 0, principal stresses from soil overburden pressure ); nodal vector (6*nnode_loc)
! P0stress = hydrostatic/volumetric stress;
! SS1 = stress - POstress;                  nodal vector (6*nnod_loc)

     subroutine   MAKE_NLP_INITIAL_STRESS(nnloc, nelem, nn_cs_loc, cs_loc, &
                                        nsurf, nmat, tag_mat, sdeg_mat, prop_mat, &
                                        xs_loc, ys_loc, zs_loc, pressdep_flag, effstress_flag,&
                                        stress)

        implicit none

        logical, intent(in) :: pressdep_flag, effstress_flag

        integer*4, intent(in) :: nnloc, nn_cs_loc, nelem, nsurf, nmat
        integer*4, dimension(nmat), intent(in) :: tag_mat, sdeg_mat
        integer*4, dimension(0:nn_cs_loc), intent(in) :: cs_loc

        real*8, dimension(nmat,4), intent(in) :: prop_mat
        real*8, dimension(nnloc), intent(in) :: xs_loc, ys_loc,zs_loc
        real*8, dimension(nnloc*6), intent(out) :: stress !, SS1

        character*70 :: file_vel_surf
        
        integer*4 :: inode, isurf, iaz
        integer*4 :: nnode_file, ntria_file
        integer*4, dimension(:), allocatable :: node1_tria, node2_tria, node3_tria
        
        real*8 ::  max_spacing, tol, Klateral, gfact, Peff0
        real*8, dimension(nnloc) :: zs_depth, zs_allu, prev_lay_zz, soil_press
        real*8, dimension(:), allocatable :: x_tria, y_tria, z_tria
        
        ! This Coeeficient is Hardcoded for now
        tol = 2.d0
        gfact = 9.8067;
        
        stress = 0;
        soil_press = 0.d0;
        !SS1 = 0.001;

        ! Should we change this? (in Oral's code, Klateral is 1, for pressure independednt formulation)
        Klateral = 0.3d0;

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Finding the depth of LGL nodes from the velocity Initerfaces. 
        ! This is to calculate the overburden soil pressure (initial stress conditions)
        
        ! Step -1 :Reading the depth of node below the Topography Surface
        file_vel_surf = 'XYZ.out'
        zs_depth = -1.0e+30

        call READ_DIME_FILEXYZ(file_vel_surf, nnode_file, ntria_file)
        allocate(x_tria(nnode_file),y_tria(nnode_file),z_tria(nnode_file))
        allocate(node1_tria(ntria_file), node2_tria(ntria_file), node3_tria(ntria_file))

        call READ_FILEXYZ(file_vel_surf, nnode_file, ntria_file, &
                    x_tria, y_tria, z_tria, node1_tria, node2_tria, node3_tria, &
                    max_spacing)

        call GET_NODE_DEPTH_FROM_TOPO_INISTRESS(nnode_file, ntria_file, &					
                    x_tria, y_tria, z_tria, node1_tria, node2_tria, node3_tria, &
                    nnloc, xs_loc, ys_loc, zs_loc, zs_depth, max_spacing, tol)
     
        deallocate(x_tria, y_tria, z_tria, node1_tria, node2_tria, node3_tria)
        max_spacing = 0.d0

        ! Step 2 : Reading the Elevation of node above the Velocity Interfaces and calculating Overburden Pressure
        prev_lay_zz = 0.d0;
        do isurf = 1,nsurf
            zs_allu = -1.0e+30

            write(file_vel_surf,'(A3, I0, A4)') 'ALL', isurf, '.out'
            call READ_DIME_FILEXYZ(file_vel_surf, nnode_file, ntria_file)
            allocate(x_tria(nnode_file), y_tria(nnode_file), z_tria(nnode_file))
            allocate(node1_tria(ntria_file), node2_tria(ntria_file), node3_tria(ntria_file))
            
            call READ_FILEXYZ(file_vel_surf, nnode_file, ntria_file, &
                    x_tria, y_tria, z_tria, node1_tria, node2_tria, node3_tria, &
                    max_spacing)

            call GET_NODE_ELEVATION_FROM_VEL_SURF_INISTRESS(nnode_file, ntria_file, &					
                    x_tria, y_tria, z_tria, node1_tria, node2_tria, node3_tria, &
                    nnloc, xs_loc, ys_loc, zs_loc, zs_allu, max_spacing, tol)

            deallocate(x_tria, y_tria, z_tria, node1_tria, node2_tria, node3_tria)
            max_spacing = 0.d0

            ! Overburden soil Pressure
            if (isurf.eq.1) then
                do inode = 1,nnloc
                     ! Checking if the GLL point is lying above the 'isurf'-th velocity surface (or below)
                     ! if (zs_allu(inode).le.0.d0) then
                     !    press_contrib = 0.d0
                     ! else
                     !    thick_
                     !    pressure_contrib = prop_mat(isurf,1)*gfact*zs_depth(inode)
                     ! endif

                    if (zs_allu(inode).eq.0.d0) then
                        soil_press(inode) = prop_mat(isurf,1)*gfact*zs_depth(inode)
                    elseif (zs_allu(inode).gt.0.d0) then
                        soil_press(inode) = prop_mat(isurf,1)*gfact*zs_depth(inode) + (prop_mat(isurf+1,1) - prop_mat(isurf,1))*gfact*zs_allu(inode)
                    endif
                enddo
            else
                do inode = 1,nnloc
                    soil_press(inode) = soil_press(inode) + (prop_mat(isurf+1,1) - prop_mat(isurf,1))*gfact*zs_allu(inode)
                enddo
            endif
        enddo


        ! Stress Tensor under Initial Conditions
        ! i.e. vertical soil pressure in zz direction;
        ! lateral soil pressure in xx, yy directions
        ! Shear stresses (xy, yz, zx) are zero.
        do inode = 1,nnloc
            iaz = 6*(inode-1)
            stress(iaz +1) = Klateral*soil_press(inode)
            stress(iaz +2) = Klateral*soil_press(inode)
            stress(iaz +3) = soil_press(inode)

            ! Deviatoric Stress for Von-Mises Criterian
            !Peff0 = (stress(iaz +1) + stress(iaz +2) + stress(iaz +3))/3.d0
            !SS1(iaz +1) = max( (stress(iaz+1)- Peff0), 0.001d0)
            !SS1(iaz +2) = max( (stress(iaz+2)- Peff0), 0.001d0)
            !SS1(iaz +3) = max( (stress(iaz+3)- Peff0), 0.001d0)
        enddo

    end subroutine MAKE_NLP_INITIAL_STRESS




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> @brief Computes Depth of GLL node from topography (XYZ.out).
! If GLL node is below Topography surface -> depth is +ve
!! @author Ilario Mazzieri, Srihari Sangaraju
!> @date November 2022

!> @param[in] nn_elev number nodes in the triangular grid
!> @param[in] nn_elem number of triangular elements
!> @param[in] x_elev elevation values of local nodes
!> @param[in] y_elev elevation values of local nodes
!> @param[in] z_elev elevation values of local nodes
!> @param[in] node1_elem index triangle vertex 
!> @param[in] node2_elem index triangle vertex
!> @param[in] node3_elem index triangle vertex
!> @param[in] nn_s number of local GLL nodes in Partition
!> @param[in] xx_s vertex x- coordinate  of local nodes
!> @param[in] yy_s vertex y- coordinate  of local nodes
!> @param[in] zz_s vertex z- coordinate  of local nodes
!> @param[in] max_es  max topography spacing
!> @param[in] tol  tolerance given in CASE option
!> @param[out] zz_depth Depth of the nodes from complex topography

      subroutine GET_NODE_DEPTH_FROM_TOPO_INISTRESS(nn_elev, nn_elem, &
                                            xx_elev, yy_elev, zz_elev, &
                                            node1_elem, node2_elem, node3_elem, &
                                            nn_s, xx_s, yy_s, zz_s, &
                                            zz_depth, max_es, tol)      

      implicit none
      
      integer*4 :: ic, h
      integer*4 :: nn_elev, nn_elem, nn_s
      integer*4, dimension(nn_elem) :: node1_elem,node2_elem,node3_elem

      real*8 :: X1,Y1,Z1                                
      real*8 :: X2,Y2,Z2                                
      real*8 :: X3,Y3,Z3                                
      real*8 :: ux,uy,uz,vx,vy,vz                
      real*8 :: a,b,c
      real*8 :: zz_interp                                
      real*8 :: v0x,v0y,v1x,v1y,v2x,v2y                        
      real*8 :: dot00,dot01,dot02,dot11,dot12        
      real*8 :: invDenom,u,v                 
      real*8 :: max_es, tol, d2min

      real*8, dimension(nn_elev) :: xx_elev,yy_elev,zz_elev
      real*8, dimension(nn_s) :: xx_s,yy_s,zz_s
      real*8, dimension(nn_s) :: zz_depth

      d2min = (5 * max_es)**2 
      
      do ic = 1,nn_s
         if (zz_depth(ic).eq.-1.0e+30) then
               do h = 1,nn_elem
                     X1 = xx_elev(node1_elem(h)) 
                     Y1 = yy_elev(node1_elem(h)) 
                     Z1 = zz_elev(node1_elem(h)) 
                     ! If the (x,y) coordinates of triangle are too far from node coordinates, skip that triangle                        
                     if (((X1 - xx_s(ic))**2 + (Y1 - yy_s(ic))**2).le.d2min) then 

                           X2 = xx_elev(node2_elem(h)) 
                           Y2 = yy_elev(node2_elem(h)) 
                           Z2 = zz_elev(node2_elem(h)) 

                           X3 = xx_elev(node3_elem(h)) 
                           Y3 = yy_elev(node3_elem(h)) 
                           Z3 = zz_elev(node3_elem(h)) 

                           ! Compute vectors 
                           ! v0 = C - A
                           v0x=(X3 - X1) 
                           v0y=(Y3 - Y1) 

                           ! v1 = B - A
                           v1x=(X2 - X1) 
                           v1y=(Y2 - Y1) 

                           ! v2 = P - A
                           v2x=(xx_s(ic) - X1) 
                           v2y=(yy_s(ic) - Y1) 
                                       
                           ! Compute dot products
                           ! [u].[v] = ux * vx + uy * vy
                           ! dot([u],[v])

                           !dot00 = dot(v0, v0)
                           dot00 = v0x * v0x + v0y * v0y

                           !dot01 = dot(v0, v1)
                           dot01 = v0x * v1x + v0y * v1y

                           !dot02 = dot(v0, v2)
                           dot02 = v0x * v2x + v0y * v2y

                           !dot11 = dot(v1, v1)
                           dot11 = v1x * v1x + v1y * v1y

                           !dot12 = dot(v1, v2)
                           dot12 = v1x * v2x + v1y * v2y

                           ! Compute barycentric coordinates
                           invDenom = 1 / (dot00 * dot11 - dot01 * dot01)
                           u = (dot11 * dot02 - dot01 * dot12) * invDenom
                           v = (dot00 * dot12 - dot01 * dot02) * invDenom

                           !Point in triangle test 
                           !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                           if ( (u.ge.0.0d0).and.(v.ge.0.0d0).and.((u + v).le.1.0d0) ) then
                                                                     
                                    ! Build up the plane passing through the points P1, P2 and P3
                                    ux=(X1-X2)/sqrt((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2) 
                                    uy=(Y1-Y2)/sqrt((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2) 
                                    uz=(Z1-Z2)/sqrt((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2) 
                                    vx=(X3-X2)/sqrt((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2) 
                                    vy=(Y3-Y2)/sqrt((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2) 
                                    vz=(Z3-Z2)/sqrt((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2) 

                                    a = uy * vz - uz * vy 
                                    b = uz * vx - ux * vz 
                                    c = ux * vy - uy * vx 
                                                            
                                    zz_interp = -a/c * (xx_s(ic)-X1) -b/c * (yy_s(ic)-Y1) + Z1
                                    zz_depth(ic) = ( zz_interp - zz_s(ic) )
                                    
                                    if (abs(zz_depth(ic)).lt.tol) then
                                             zz_depth(ic) = 0.0d0
                                    endif
                           endif
                                                            
                           if ( (u.ge.0.0d0).and.(v.ge.0.0d0).and.((u + v).le.1.0d0) ) exit

                     endif !if (((X1 - xx_s(ic))**2 + (Y1 - yy_s(ic))**2).le.d2min) then 
               enddo ! do h = 1,nn_elem
         endif !if (zz_depth(ic).eq.-1.0e+30)
      enddo ! do ic = 1,nn_s
      
      return
      
      end subroutine GET_NODE_DEPTH_FROM_TOPO_INISTRESS




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> @brief Computes Elevation of GLL node above the velocity/alluvial basin surface (ALL.out).
! If GLG node is above the velocity surface -> elevation is +ve
!! @author Ilario Mazzieri, Srihari Sangaraju
!> @date November 2022

!> @param[in] n_elev number nodes in the triangular grid
!> @param[in] nn_elem number of triangular elements
!> @param[in] x_elev elevation values of local nodes
!> @param[in] y_elev elevation values of local nodes
!> @param[in] z_elev elevation values of local nodes
!> @param[in] node1_elem index triangle vertex 
!> @param[in] node2_elem index triangle vertex
!> @param[in] node3_elem index triangle vertex
!> @param[in] nn_s number of local nodes
!> @param[in] xx_s vertex x- coordinate  of local nodes
!> @param[in] yy_s vertex y- coordinate  of local nodes
!> @param[in] zz_s vertex z- coordinate  of local nodes
!> @param[in] max_as  max alluvial spacing
!> @param[in] tol  tolerance given in CASE option
!> @param[out] zz_elevation elevation of the nodes from alluvial

      subroutine GET_NODE_ELEVATION_FROM_VEL_SURF_INISTRESS(n_elev, nn_elem, &
                                              x_elev, y_elev, z_elev, &
                                              node1_elem, node2_elem, node3_elem, &
                                              nn_s, xx_s, yy_s, zz_s, zz_elevation, max_as, tol)
      
      implicit none
      
      integer*4 :: n_elev, nn_elem,nn_s
      integer*4 :: i, ic, h
      integer*4, dimension(nn_elem) :: node1_elem,node2_elem,node3_elem

      real*8 :: X1,Y1,Z1
      real*8 :: X2,Y2,Z2
      real*8 :: X3,Y3,Z3
      real*8 :: ux,uy,uz,vx,vy,vz
      real*8 :: a,b,c
      real*8 :: zz_interp
      real*8 :: v0x,v0y,v1x,v1y,v2x,v2y
      real*8 :: dot00,dot01,dot02,dot11,dot12
      real*8 :: invDenom,u,v
      real*8 :: tol, max_as, d2min 
      real*8 :: z_surf_min, z_surf_max

      real*8, dimension(n_elev) :: x_elev,y_elev,z_elev
      real*8, dimension(nn_s) :: xx_s,yy_s,zz_s
      real*8, dimension(nn_s) :: zz_elevation


      d2min = (5 * max_as)**2 
         
      z_surf_min = z_elev(1)
      z_surf_max = z_elev(1)
      do i = 1,n_elev
         if (z_elev(i).lt.z_surf_min) z_surf_min = z_elev(i)
         if (z_elev(i).gt.z_surf_max) z_surf_max = z_elev(i)
      enddo
         
      do ic = 1,nn_s
         if (zz_elevation(ic) .eq. -1.0e+30) then
               do h = 1,nn_elem
                  if ((zz_s(ic).lt.z_surf_min) .or. (zz_s(ic).gt.z_surf_max)) exit 

                  X1 = x_elev(node1_elem(h)) 
                  Y1 = y_elev(node1_elem(h)) 
                  Z1 = z_elev(node1_elem(h)) 
                                                      
                  if (((X1 - xx_s(ic))**2 + (Y1 - yy_s(ic))**2).le.d2min) then 

                        X2 = x_elev(node2_elem(h)) 
                        Y2 = y_elev(node2_elem(h)) 
                        Z2 = z_elev(node2_elem(h)) 

                        X3 = x_elev(node3_elem(h)) 
                        Y3 = y_elev(node3_elem(h)) 
                        Z3 = z_elev(node3_elem(h)) 

                        if (zz_s(ic).ge.min(Z1,Z2,Z3)) then 
                              ! Compute vectors 
                              ! v0 = C - A
                              v0x=(X3 - X1) 
                              v0y=(Y3 - Y1) 

                              ! v1 = B - A
                              v1x=(X2 - X1) 
                              v1y=(Y2 - Y1) 

                              ! v2 = P - A
                              v2x=(xx_s(ic) - X1) 
                              v2y=(yy_s(ic) - Y1) 

                              ! Compute dot products
                              ! [u].[v] = ux * vx + uy * vy
                              ! dot([u],[v])

                              !dot00 = dot(v0, v0)
                              dot00 = v0x * v0x + v0y * v0y

                              !dot01 = dot(v0, v1)
                              dot01 = v0x * v1x + v0y * v1y

                              !dot02 = dot(v0, v2)
                              dot02 = v0x * v2x + v0y * v2y

                              !dot11 = dot(v1, v1)
                              dot11 = v1x * v1x + v1y * v1y

                              !dot12 = dot(v1, v2)
                              dot12 = v1x * v2x + v1y * v2y

                              ! Compute barycentric coordinates
                              invDenom = 1 / (dot00 * dot11 - dot01 * dot01)
                              u = (dot11 * dot02 - dot01 * dot12) * invDenom
                              v = (dot00 * dot12 - dot01 * dot02) * invDenom

                              !Point in triangle test
                              !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                              if ( (u.ge.0.0d0).and.(v.ge.0.0d0).and.((u + v).le.1.0d0) ) then                                    
                                    ! Build up the plane passing through the points P1, P2 and P3
                                    ux=(X1-X2)/sqrt((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2) 
                                    uy=(Y1-Y2)/sqrt((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2) 
                                    uz=(Z1-Z2)/sqrt((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2) 
                                    vx=(X3-X2)/sqrt((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2) 
                                    vy=(Y3-Y2)/sqrt((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2) 
                                    vz=(Z3-Z2)/sqrt((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2) 

                                    a = uy * vz - uz * vy 
                                    b = uz * vx - ux * vz 
                                    c = ux * vy - uy * vx 
                                             
                                    zz_interp = -a/c * (xx_s(ic)-X1) -b/c * (yy_s(ic)-Y1) + Z1
                                    zz_elevation(ic) = -1.d0 * ( zz_interp - zz_s(ic) )
                                    if (abs(zz_elevation(ic)).lt.tol) then
                                          zz_elevation(ic) = 0.0d0
                                    endif
                              endif
                                       
                              if ( (u.ge.0.0d0).and.(v.ge.0.0d0).and.((u + v).le.1.0d0) ) exit
                              
                        endif !if (zz_s(ic).ge.min(Z1,Z2,Z3)) then 
                  endif !if (((X1 - xx_s(ic))**2 + (Y1 - yy_s(ic))**2).le.d2min) then 
            enddo !do h = 1,nn_elem
         endif !if (zz_elevation(ic).le.0.0d0) then
      enddo !do ic = 1,nn_s 
      
      return
      
      end subroutine GET_NODE_ELEVATION_FROM_VEL_SURF_INISTRESS