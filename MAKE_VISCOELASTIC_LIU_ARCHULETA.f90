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


!> @brief Initialise Damping Parameters - Frequency Independent damping (Liu & Archuleta 2006)
!! @author Srihari Sangaraju
!> @date November, 2022 
!> @version 1.0
!
! viscoelastic_Mp = unrelaxed P-wave modulus (based on lambda + 2*shear_modulus)
! viscoelastic_Ms = unrelaxed S-wave modulus (shear_modulus)
! 
    subroutine   MAKE_VISCOELASTIC_LIU_ARCHULETA(nmat, nsls, prop_mat, QS, QP, f_ref, &
                                            viselastic_Trelax, viselastic_wgt_s, viselastic_wgt_p, &
                                            viselastic_Ms, viselastic_Mp)

        implicit none

        integer*4, intent(in) :: nmat, nsls
        
        real*8, intent(in) :: f_ref
        real*8, intent(in), dimension(nmat)   :: QS, QP
        real*8, intent(in), dimension(nmat,4) :: prop_mat
        real*8, intent(out), dimension(nsls)  :: viselastic_Trelax
        real*8, intent(out), dimension(nmat)  :: viselastic_Ms, viselastic_Mp
        real*8, intent(out), dimension(nmat,nsls) :: viselastic_wgt_s, viselastic_wgt_p

        integer*4 :: imat, isls
        real*8 :: pi, chi_s, chi_p
        real*8, dimension(nsls)      :: viselastic_alpha, viselastic_beta
        complex*8 :: A1_s, A1_p, dum_cmplx

        pi = 4.d0*datan(1.d0);

        ! Discretisation of relaxation spectra using 8 relaxation times.
        ! Quality factor = Real(ViscoElasticModulus)/Imaginary(ViscoElasticModulus)
        ! The Goal here is to keep Quality Factor constant in frequency domain.
        ! Here Relaxation Modulus is calculated as summation over 8 different relaxation times (Equation 5, Liu 2006)
        ! Relaxation times, and regression coefficients, From Table 1 (Liu et al.2006)
        viselastic_Trelax = (/1.72333e-3,1.80701e-3,5.38887e-3,1.99322e-2,8.49833e-2,4.09335e-1,2.05951,13.2629/)
        viselastic_alpha  = (/1.66958e-2,3.81644e-2,9.84666e-3,-1.36803e-2,-2.85125e-2,-5.37309e-2,-6.65035e-2,-1.33696e-1/)
        viselastic_beta  = (/8.98758e-2,6.84635e-2,9.67052e-2,1.20172e-1,1.30728e-1,1.38746e-1,1.40705e-1,2.14647e-1/)


        ! Weights corresponding to each relaxation time, Equation 6a and 6b
        do imat = 1,nmat
            chi_s = (3.071+ 1.433* (QS(imat)**(-1.158))* log(QS(imat)/5.d0))/ (1.d0+ 0.415*QS(imat))
            chi_p = (3.071+ 1.433* (QP(imat)**(-1.158))* log(QP(imat)/5.d0))/ (1.d0+ 0.415*QP(imat))
            
            do isls = 1,nsls
                viselastic_wgt_s(imat, isls) = chi_s * (chi_s* viselastic_alpha(isls) + viselastic_beta(isls))        ! S-wave
                viselastic_wgt_p(imat, isls) = chi_p * (chi_p* viselastic_alpha(isls) + viselastic_beta(isls))        ! P-Wave
            enddo
        enddo

        ! Visco Elastic Modulus Calculation (Equation 1: Mu = A1*Mu_viselastic)
        do imat = 1,nmat
            A1_s =  dcmplx(1.d0,0.d0);
            A1_p =  dcmplx(1.d0,0.d0);

            do isls = 1,nsls
                dum_cmplx = 1.d0 +  DCMPLX(0.d0,1.d0)* (2.d0*pi*f_ref)* viselastic_Trelax(isls)
                A1_s = A1_s - (viselastic_wgt_s(imat, isls)/ dum_cmplx)
                A1_p = A1_p - (viselastic_wgt_p(imat, isls)/ dum_cmplx)
            enddo

            viselastic_Ms(imat) = prop_mat(imat,3) /cabs(A1_s)
            viselastic_Mp(imat) = (prop_mat(imat,2) + 2*prop_mat(imat,3)) /cabs(A1_p)
        enddo


    end subroutine MAKE_VISCOELASTIC_LIU_ARCHULETA