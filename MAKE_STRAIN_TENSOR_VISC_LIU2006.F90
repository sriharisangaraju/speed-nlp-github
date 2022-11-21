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
!> @param[in] nsls number of relaxation times for SLS elements
!> @param[in] Mp nodal values of (Unrelaxed lambda + 2*Mu)
!> @param[in] Ms nodal values of Unrelaxed Shear Modulus
!> @param[in] duxdx nodal values for spatial derivatives of the displacement 
!> @param[in] duydx nodal values for spatial derivatives of the displacement
!> @param[in] duzdx nodal values for spatial derivatives of the displacement
!> @param[in] duxdy nodal values for spatial derivatives of the displacement
!> @param[in] duydy nodal values for spatial derivatives of the displacement
!> @param[in] duzdy nodal values for spatial derivatives of the displacement
!> @param[in] duxdz nodal values for spatial derivatives of the displacement
!> @param[in] duydz nodal values for spatial derivatives of the displacement      
!> @param[in] duzdz nodal values for spatial derivatives of the displacement
!> @param[in] zeta_xx nodal values for relaxation function (calculated from linear combination of relaxation spectra values at NSLS points (relaxation times) )
!> @param[out] strain_visc_xx nodal values for the viscous strain relaxation tensor
!> @param[out] strain_visc_xx nodal values for the viscous strain relaxation tensor
!> @param[out] strain_visc_xx nodal values for the viscous strain relaxation tensor
!> @param[out] strain_visc_xx nodal values for the viscous strain relaxation tensor
!> @param[out] strain_visc_xx nodal values for the viscous strain relaxation tensor
!> @param[out] strain_visc_xx nodal values for the viscous strain relaxation tensor

!> @brief Frequency Independent damping 
!         extended to 3D (based on Liu et al., 2006; Day at al., 2001; Day et al., 1998)
!! @author Srihari Sangaraju
!> @date November, 2022 
!> @version 1.0
!

      subroutine MAKE_STRAIN_TENSOR_VISC_LIU2006(nnod_loc, nn_cs_loc, cs_loc, node_counter, &
                               nn, ie, nsls, exp_Trelax, &
                               Mp, Ms, wt_s, wt_p, &
                               duxdx,duydx,duzdx,&
                               duxdy,duydy,duzdy,&
                               duxdz,duydz,duzdz,&
                               zeta_node_t0, zeta_node_t1, &
                               stress_xx,stress_yy,stress_zz,&
                               stress_yz,stress_zx,stress_xy)
      
      
      implicit none
      
      integer*4, intent(in) :: nn, nsls, ie, nnod_loc, nn_cs_loc
      integer*4, intent(in), dimension(nnod_loc) :: node_counter
      integer*4, intent(in), dimension(nn_cs_loc) :: cs_loc
      
      real*8, intent(in), dimension(nsls)     :: exp_Trelax, wt_s, wt_p
      real*8, intent(in), dimension(nn,nn,nn) :: Mp, Ms
      real*8, intent(in), dimension(nn,nn,nn) :: duxdx, duxdy, duxdz, duydx, duydy, duydz, duzdx, duzdy, duzdz
      real*8, intent(in),    dimension(nnod_loc*6,nsls) :: zeta_node_t0
      real*8, intent(inout), dimension(nnod_loc*6,nsls) :: zeta_node_t1
      real*8, intent(out),   dimension(nn,nn,nn) :: stress_xx,stress_xy,stress_zx
      real*8, intent(out),   dimension(nn,nn,nn) :: stress_yy,stress_yz,stress_zz

      integer*4 :: p, q, r, isls, is, in, iaz

      real*8 :: lambda_comp_relaxfcn, lambda_comp_stress, dum_relaxfcn
      real*8 :: sum_relaxfcn_xx, sum_relaxfcn_yy, sum_relaxfcn_zz
      real*8 :: sum_relaxfcn_xy, sum_relaxfcn_yz, sum_relaxfcn_xz
      real*8, dimension(nn,nn,nn) :: strain_xx,strain_xy,strain_xz
      real*8, dimension(nn,nn,nn) :: strain_yy,strain_yz,strain_zz, strain_vol

      
      ! Total Strain - Elastic (at all GLL nodes in 1 element)
      strain_xx = duxdx;      strain_xy = (duxdy + duydx)/2.d0
      strain_yy = duydy;      strain_yz = (duydz + duzdy)/2.d0
      strain_zz = duzdz;      strain_xz = (duzdx + duxdz)/2.d0

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
                        iaz = 7*(in-1) + 1
                        dum_relaxfcn = exp_Trelax(isls)*zeta_node_t0(iaz,isls) +  &
                                       lambda_comp_relaxfcn + 2*Ms(p,q,r)*(1.d0-exp_Trelax(isls))*strain_xx(p,q,r)
                        zeta_node_t1(iaz,isls) = zeta_node_t1(iaz,isls) + (dum_relaxfcn/node_counter(in))
                        sum_relaxfcn_xx = sum_relaxfcn_xx + dum_relaxfcn

                        !YY
                        iaz = 7*(in-1) + 2
                        dum_relaxfcn = exp_Trelax(isls)*zeta_node_t0(iaz,isls) +  &
                                       lambda_comp_relaxfcn + 2*Ms(p,q,r)*(1.d0-exp_Trelax(isls))*strain_yy(p,q,r)
                        zeta_node_t1(iaz,isls) = zeta_node_t1(iaz,isls) + (dum_relaxfcn/node_counter(in))
                        sum_relaxfcn_yy = sum_relaxfcn_yy + dum_relaxfcn

                        !ZZ
                        iaz = 7*(in-1) + 3
                        dum_relaxfcn = exp_Trelax(isls)*zeta_node_t0(iaz,isls) +  &
                                       lambda_comp_relaxfcn + 2*Ms(p,q,r)*(1.d0-exp_Trelax(isls))*strain_zz(p,q,r)
                        zeta_node_t1(iaz,isls) = zeta_node_t1(iaz,isls) + (dum_relaxfcn/node_counter(in))
                        sum_relaxfcn_zz = sum_relaxfcn_zz + dum_relaxfcn


                        ! Shear Components (this is Absolute strain/ not shear strain)
                        !XY - Calculating the strain relaxation here
                        iaz = 7*(in-1) + 4
                        dum_relaxfcn = exp_Trelax(isls)*zeta_node_t0(iaz,isls) + wt_s(isls)*(1.d0-exp_Trelax(isls))* strain_xy(p,q,r)
                        zeta_node_t1(iaz,isls) = zeta_node_t1(iaz,isls) + (dum_relaxfcn/node_counter(in))
                        sum_relaxfcn_xy = sum_relaxfcn_xy + dum_relaxfcn

                        !YZ
                        iaz = 7*(in-1) + 5
                        dum_relaxfcn = exp_Trelax(isls)*zeta_node_t0(iaz,isls) + wt_s(isls)*(1.d0-exp_Trelax(isls))* strain_yz(p,q,r)
                        zeta_node_t1(iaz,isls) = zeta_node_t1(iaz,isls) + (dum_relaxfcn/node_counter(in))
                        sum_relaxfcn_yz = sum_relaxfcn_yz + dum_relaxfcn

                        !XZ
                        iaz = 7*(in-1) + 6
                        dum_relaxfcn = exp_Trelax(isls)*zeta_node_t0(iaz,isls) + wt_s(isls)*(1.d0-exp_Trelax(isls))* strain_xz(p,q,r)
                        zeta_node_t1(iaz,isls) = zeta_node_t1(iaz,isls) + (dum_relaxfcn/node_counter(in))
                        sum_relaxfcn_xz = sum_relaxfcn_xz + dum_relaxfcn
                    enddo
                    
                    ! Calculating Total Stressess after Damping (or Relaxation)
                    lambda_comp_stress = (Mp(p,q,r) - 2*Ms(p,q,r)) * strain_vol(p,q,r)
                    stress_xx(p,q,r) = lambda_comp_stress + 2*Ms(p,q,r)*strain_xx(p,q,r) - sum_relaxfcn_xx;
                    stress_yy(p,q,r) = lambda_comp_stress + 2*Ms(p,q,r)*strain_yy(p,q,r) - sum_relaxfcn_yy;
                    stress_zz(p,q,r) = lambda_comp_stress + 2*Ms(p,q,r)*strain_zz(p,q,r) - sum_relaxfcn_zz;
                    stress_xy(p,q,r) = 2*Ms(p,q,r)* (strain_xy(p,q,r) - sum_relaxfcn_xy)
                    stress_yz(p,q,r) = 2*Ms(p,q,r)* (strain_yz(p,q,r) - sum_relaxfcn_yz)
                    stress_zx(p,q,r) = 2*Ms(p,q,r)* (strain_xz(p,q,r) - sum_relaxfcn_xz)
                enddo
            enddo
        enddo

      return
      
      end subroutine MAKE_STRAIN_TENSOR_VISC_LIU2006




      !       subroutine MAKE_STRESS_TENSOR_DAMPED_LIU2006(nn, k_urlx, mu_urlx, &
      !                          strain_vis_dev_xx, strain_vis_dev_xy, strain_vis_dev_xz, & 
      !                          strain_vis_dev_yy, strain_vis_dev_yz, strain_vis_dev_zz, strain_vis_vol,&
      !                          sxx, syy, szz, syz, szx, sxy)
      
      
      ! implicit none
      
      ! integer*4, intent(in) :: nn

      ! real*8, intent(in),  dimension(nn,nn,nn) :: k_urlx, mu_urlx
      ! real*8, intent(in),  dimension(nn,nn,nn) :: strain_vis_dev_xx, strain_vis_dev_xy, strain_vis_dev_xz
      ! real*8, intent(in),  dimension(nn,nn,nn) :: strain_vis_dev_yy, strain_vis_dev_yz, strain_vis_dev_zz, strain_vis_vol
      ! real*8, intent(out), dimension(nn,nn,nn) :: sxx, syy, szz, syz, szx, sxy

      ! integer*4 :: p, q, r

      !   do r = 1,nn
      !       do q = 1,nn
      !           do p = 1,nn
                    
      !               sxx(p,q,r) = 3.d0*k_urlx(p,q,r)*strain_vis_vol(p,q,r) + 2.d0*mu_urlx(p,q,r)*strain_vis_dev_xx(p,q,r)
      !               syy(p,q,r) = 3.d0*k_urlx(p,q,r)*strain_vis_vol(p,q,r) + 2.d0*mu_urlx(p,q,r)*strain_vis_dev_yy(p,q,r)
      !               szz(p,q,r) = 3.d0*k_urlx(p,q,r)*strain_vis_vol(p,q,r) + 2.d0*mu_urlx(p,q,r)*strain_vis_dev_zz(p,q,r)

      !               syz(p,q,r) = 2.d0*mu_urlx(p,q,r)*strain_vis_dev_yz(p,q,r)
      !               szx(p,q,r) = 2.d0*mu_urlx(p,q,r)*strain_vis_dev_xz(p,q,r)
      !               sxy(p,q,r) = 2.d0*mu_urlx(p,q,r)*strain_vis_dev_xy(p,q,r)

      !           enddo
      !       enddo
      !   enddo


      ! return
      
      ! end subroutine MAKE_STRESS_TENSOR_DAMPED_LIU2006

