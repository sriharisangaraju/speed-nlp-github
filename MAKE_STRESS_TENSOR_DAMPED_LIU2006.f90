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

!> @param[in] nn number of 1-D Legendre nodes
!> @param[in] duxdx nodal values for spatial derivatives of the displacement 
!> @param[in] duydx nodal values for spatial derivatives of the displacement
!> @param[in] duzdx nodal values for spatial derivatives of the displacement
!> @param[in] duxdy nodal values for spatial derivatives of the displacement
!> @param[in] duydy nodal values for spatial derivatives of the displacement
!> @param[in] duzdy nodal values for spatial derivatives of the displacement
!> @param[in] duxdz nodal values for spatial derivatives of the displacement
!> @param[in] duydz nodal values for spatial derivatives of the displacement      
!> @param[in] duzdz nodal values for spatial derivatives of the displacement
!> @param[in] nsls number of relaxation times for SLS elements
!> @param[in] Mp nodal values of (Unrelaxed lambda + 2*Mu)
!> @param[in] Ms nodal values of Unrelaxed Shear Modulus
!> @param[in] wt_p Weights for Memory Variables, corresponding to respective Relaxation times (p-wave modulus)
!> @param[in] wt_s Weights for Memory Variables, corresponding to respective Relaxation times (s-wave modulus)
!> @param[in] zeta_node_t0 Relaxation functions / Memory Variables at Nodes (calculated in previous time step)
!> @param[inout] zeta_node_t1 Relaxation functions / Memory Variables at Nodes (For currrent time step)
!> @param[out] stress_xx nodal values for the damped stress tensor
!> @param[out] stress_yy nodal values for the damped stress tensor
!> @param[out] stress_zz nodal values for the damped stress tensor
!> @param[out] stress_yz nodal values for the damped stress tensor
!> @param[out] stress_zx nodal values for the damped stress tensor
!> @param[out] stress_xy nodal values for the damped stress tensor

!> @brief Calculating Stress Tensor after damping
!         Frequency Independent damping extended to 3D (based on Liu et al., 2006; Day at al., 2001; Day et al., 1998)
!! @author Srihari Sangaraju
!> @date November, 2022 
!> @version 1.0
!

      subroutine MAKE_STRESS_TENSOR_DAMPED_LIU2006(nnod_loc, nn_cs_loc, cs_loc, node_counter, &
                               nn, ie, nsls, exp_Trelax, &
                               Mp, Ms, wt_p, wt_s, &
                               duxdx,duydx,duzdx,&
                               duxdy,duydy,duzdy,&
                               duxdz,duydz,duzdz,&
                               zeta_node_t0, zeta_node_t1, &
                               stress_xx,stress_yy,stress_zz,&
                               stress_yz,stress_zx,stress_xy)
      
      
      implicit none
      
      integer*4, intent(in) :: nn, nsls, ie, nnod_loc, nn_cs_loc
      integer*4, intent(in), dimension(nnod_loc) :: node_counter
      integer*4, intent(in), dimension(0:nn_cs_loc) :: cs_loc
      
      real*8, intent(in), dimension(nsls)     :: exp_Trelax, wt_s, wt_p
      real*8, intent(in), dimension(nn,nn,nn) :: Mp, Ms
      real*8, intent(in), dimension(nn,nn,nn) :: duxdx, duxdy, duxdz, duydx, duydy, duydz, duzdx, duzdy, duzdz
      real*8, intent(in),    dimension(nnod_loc*6,nsls) :: zeta_node_t0
      real*8, intent(inout), dimension(nnod_loc*6,nsls) :: zeta_node_t1
      real*8, intent(out),   dimension(nn,nn,nn) :: stress_xx,stress_xy,stress_zx
      real*8, intent(out),   dimension(nn,nn,nn) :: stress_yy,stress_yz,stress_zz

      integer*4 :: p, q, r, isls, is, in, iaz

      real*8 :: lambda_comp_relaxfcn, lambda_comp_stress, dum_relaxfcn, dum1, dum2, dum3, dum4, dum5
      real*8 :: sum_relaxfcn_xx, sum_relaxfcn_yy, sum_relaxfcn_zz
      real*8 :: sum_relaxfcn_xy, sum_relaxfcn_yz, sum_relaxfcn_xz
      real*8, dimension(nn,nn,nn) :: strain_xx,strain_xy,strain_xz
      real*8, dimension(nn,nn,nn) :: strain_yy,strain_yz,strain_zz, strain_vol
      real*8, dimension(nsls) :: xz_relaxfcn

      
      ! Total Strain - Elastic (at all GLL nodes in 1 element)
      strain_xx = 0.d0!duxdx;      !strain_xy = 0.d0!(duxdy + duydx)/2.d0
      strain_yy = 0.d0!duydy;      !strain_yz = 0.d0!(duydz + duzdy)/2.d0
      strain_zz = 0.d0!duzdz;      !strain_xz = (duzdx + duxdz)/2.d0
      strain_xy = 0.d0!(duxdy + duydx)/2.d0
      strain_yz = 0.d0!(duydz + duzdy)/2.d0
      strain_xz = (duzdx + duxdz)/2.d0

      strain_vol = (strain_xx + strain_yy + strain_zz)

      ! Calculating Stress or Strain Relaxations because of viscosity 
      ! Based on relaxation functions - based on Liu et al., 2006;;;; extended for 3D case
        do r = 1,nn
            do q = 1,nn
                do p = 1,nn
                  
                    is = nn*nn*(r -1) + nn*(q -1) + p
                    in = cs_loc(cs_loc(ie -1) + is)

                    sum_relaxfcn_xx = 0.d0;   sum_relaxfcn_yy=0.d0;   sum_relaxfcn_zz=0.d0;
                    sum_relaxfcn_xy = 0.d0;   sum_relaxfcn_yz=0.d0;   sum_relaxfcn_xz=0.d0;

                    do isls = 1,nsls

                        ! For Principal Axes - dilatational related component of relaxation function
                        lambda_comp_relaxfcn = (wt_p(isls) * Mp(p,q,r) - 2.d0 * wt_s(isls) * Ms(p,q,r)) * &
                                                (1.d0-exp_Trelax(isls)) * strain_vol(p,q,r)

                        !XX - Calculating the stress relaxation here
                        iaz = 6*(in-1) + 1
                        dum_relaxfcn = exp_Trelax(isls)*zeta_node_t0(iaz,isls) +  &
                                       lambda_comp_relaxfcn + 2*Ms(p,q,r)*(1.d0-exp_Trelax(isls))*strain_xx(p,q,r)
                        zeta_node_t1(iaz,isls) = zeta_node_t1(iaz,isls) + (dum_relaxfcn/node_counter(in))
                        sum_relaxfcn_xx = sum_relaxfcn_xx + dum_relaxfcn

                        !YY
                        iaz = 6*(in-1) + 2
                        dum_relaxfcn = exp_Trelax(isls)*zeta_node_t0(iaz,isls) +  &
                                       lambda_comp_relaxfcn + 2*Ms(p,q,r)*(1.d0-exp_Trelax(isls))*strain_yy(p,q,r)
                        zeta_node_t1(iaz,isls) = zeta_node_t1(iaz,isls) + (dum_relaxfcn/node_counter(in))
                        sum_relaxfcn_yy = sum_relaxfcn_yy + dum_relaxfcn

                        !ZZ
                        iaz = 6*(in-1) + 3
                        dum_relaxfcn = exp_Trelax(isls)*zeta_node_t0(iaz,isls) +  &
                                       lambda_comp_relaxfcn + 2*Ms(p,q,r)*(1.d0-exp_Trelax(isls))*strain_zz(p,q,r)
                        zeta_node_t1(iaz,isls) = zeta_node_t1(iaz,isls) + (dum_relaxfcn/node_counter(in))
                        sum_relaxfcn_zz = sum_relaxfcn_zz + dum_relaxfcn


                        ! Shear Components (this is Absolute strain/ not shear strain)
                        !XY - Calculating the strain relaxation here
                        iaz = 6*(in-1) + 4
                        dum_relaxfcn = exp_Trelax(isls)*zeta_node_t0(iaz,isls) + wt_s(isls)*(1.d0-exp_Trelax(isls))* strain_xy(p,q,r)
                        zeta_node_t1(iaz,isls) = zeta_node_t1(iaz,isls) + (dum_relaxfcn/node_counter(in))
                        sum_relaxfcn_xy = sum_relaxfcn_xy + dum_relaxfcn

                        !YZ
                        iaz = 6*(in-1) + 5
                        dum_relaxfcn = exp_Trelax(isls)*zeta_node_t0(iaz,isls) + wt_s(isls)*(1.d0-exp_Trelax(isls))* strain_yz(p,q,r)
                        zeta_node_t1(iaz,isls) = zeta_node_t1(iaz,isls) + (dum_relaxfcn/node_counter(in))
                        sum_relaxfcn_yz = sum_relaxfcn_yz + dum_relaxfcn

                        !XZ
                        iaz = 6*(in-1) + 6
                        dum_relaxfcn = exp_Trelax(isls)*zeta_node_t0(iaz,isls) + wt_s(isls)*(1.d0-exp_Trelax(isls))* strain_xz(p,q,r)
                        zeta_node_t1(iaz,isls) = zeta_node_t1(iaz,isls) + (dum_relaxfcn/node_counter(in))
                        sum_relaxfcn_xz = sum_relaxfcn_xz + dum_relaxfcn

                        !xz_relaxfcn(isls) = dum_relaxfcn
                    enddo
                    
                    ! Calculating Total Stressess after Damping (or Relaxation)
                    lambda_comp_stress = (Mp(p,q,r) - 2*Ms(p,q,r)) * strain_vol(p,q,r)
                    stress_xx(p,q,r) = lambda_comp_stress + 2*Ms(p,q,r)*strain_xx(p,q,r) - sum_relaxfcn_xx;
                    stress_yy(p,q,r) = lambda_comp_stress + 2*Ms(p,q,r)*strain_yy(p,q,r) - sum_relaxfcn_yy;
                    stress_zz(p,q,r) = lambda_comp_stress + 2*Ms(p,q,r)*strain_zz(p,q,r) - sum_relaxfcn_zz;
                    stress_xy(p,q,r) = 2*Ms(p,q,r)* (strain_xy(p,q,r) - sum_relaxfcn_xy)
                    stress_yz(p,q,r) = 2*Ms(p,q,r)* (strain_yz(p,q,r) - sum_relaxfcn_yz)
                    stress_zx(p,q,r) = 2*Ms(p,q,r)* (strain_xz(p,q,r) - sum_relaxfcn_xz)

                    if ( (ie.eq.1) .and. ((p+q+r).eq.3) ) then
                        !dum5 = lambda_comp_stress + 2*Ms(p,q,r)*strain_xx(p,q,r);
                        !dum1 = max( abs(100.d0*sum_relaxfcn_xx/dum5), 0.00000001);

                        !dum5 = lambda_comp_stress + 2*Ms(p,q,r)*strain_zz(p,q,r);
                        !dum2 = max( abs(100.d0*sum_relaxfcn_zz/dum5), 0.00000001);

                        !dum3 = max( abs(100.d0*sum_relaxfcn_xy/strain_xy(p,q,r)), 0.00000001);
                        !dum4 = max( abs(100.d0*sum_relaxfcn_xz/strain_xz(p,q,r)), 0.00000001);
                        !write(*,*) '% strain Relaxation xz = ', strain_xz(p,q,r), sum_relaxfcn_xz, dum4
                        !write(*,*) '% strain Relaxation xx = ', dum1, 'zz = ', dum2, 'xy = ', dum3, 'xz = ', dum4
                        !write(*,*) 'Strain xx = ', strain_xx(p,q,r), 'Strain zz = ', strain_zz(p,q,r), 'xy = ', strain_xy(p,q,r), 'zx = ', strain_xz(p,q,r)
                    endif
                enddo
            enddo
        enddo

      return
      
      end subroutine MAKE_STRESS_TENSOR_DAMPED_LIU2006


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                        !!!
!!!     Making Strain Tensor from Stress Tensor for visco-elastic          !!!
!!!                                                                        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     subroutine MAKE_STRAIN_TENSOR_FROM_STRESSTENSOR(nn, Mu, Kbulk, &
                                        sxx, syy, szz, syz, szx, sxy, &
                                        strain_xx, strain_yy, strain_zz, &
                                        strain_yz, strain_xz, strain_xy)

      implicit none

      integer*4, intent(in) :: nn
      
      real*8, intent(in), dimension(nn,nn,nn) :: Mu, Kbulk
      real*8, intent(in), dimension(nn,nn,nn) :: sxx, syy, szz, syz, szx, sxy
      real*8, intent(out),   dimension(nn,nn,nn) :: strain_xx, strain_yy, strain_zz
      real*8, intent(out),   dimension(nn,nn,nn) :: strain_yz, strain_xz, strain_xy

      real*8, dimension(nn,nn,nn) :: strain_vol, stress_vol

        ! Bulk Modulus
        !Kbulk = Mp - 4.d0*Mu/3.d0
        !Kbulk = lambda + 2.d0*Mu/3.d0
        
        ! Strain = Volumetric Strain + Deviatoric Strain
        ! Stress = Volumetric Stress + Deiatoric Stress
        ! strain_vol = ex + ey + ez;  hydrostatic strain strain_hyd = strain_vol/3
        ! Volumetric/Hydrostatic Stress = Kbulk*strain_vol
        ! Deviatoric Stress = 2*Mu*Strain_deviatoric
        stress_vol = (sxx + syy + szz)/3.d0
        strain_vol = stress_vol / Kbulk         !component-wise matrix division

        strain_xx = 0.5d0* (sxx - stress_vol)/Mu + strain_vol
        strain_yy = 0.5d0* (syy - stress_vol)/Mu + strain_vol
        strain_zz = 0.5d0* (szz - stress_vol)/Mu + strain_vol
        ! strain_yz = 0.5d0* (syz/Mu)
        ! strain_xz = 0.5d0* (szx/Mu)
        ! strain_xy = 0.5d0* (sxy/Mu)
        strain_yz = (syz/Mu)
        strain_xz = (szx/Mu)
        strain_xy = (sxy/Mu)
        
     end subroutine MAKE_STRAIN_TENSOR_FROM_STRESSTENSOR