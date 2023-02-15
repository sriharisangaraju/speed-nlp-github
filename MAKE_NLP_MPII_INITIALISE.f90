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
! These functions, prepare the input data related to non-linearity (which are mentioned in *.mate file)
! The elastic/viscoelastic other parameters are assigned based on type of material (pressure dependency or total/effective stress formulation)
! For all Iwan springs of a given material, corresponding Yield stress and Strains are interpolated based on input G/Gmax curve
!! @author Srihari Sangaraju
!> @date August, 2022 
!> @version 1.0

     subroutine   MAKE_NLP_IWANSPRING_YIELDVALUES()

        use NLP_MPII
        use speed_par, only: nfunc, tag_func, func_type, func_indx, func_data, nfunc_data, &
                            nmat, tag_mat, sdeg_mat, prop_mat, mpi_id, mpi_np, nelem_loc, con_spx_loc
        use speed_exit_codes

        implicit none

        !include 'SPEED.MPI'

        integer*4 :: im, im_nlp, ifunc, ie, nn, dum_log, ispr

        real*8 :: ref_strain
        real*8 :: dum1, strain_min, strain_max


        if (mpi_id.eq.0) write(*,'(A)')
        if (mpi_id.eq.0) write(*,'(A)')'------- NLP Material :: Making IWANs Springs (MPII Model) -------' 

        allocate(mpii_mat(nmat_nlp))

        do im_nlp = 1, nmat_nlp
            do im = 1, nmat
                if (tag_mat(im) .eq. tag_mat_nlp(im_nlp)) then
                    mpii_mat(im_nlp)%mesh_block_id = tag_mat_nlp(im_nlp)

                    ! Elastic Modulii
                    mpii_mat(im_nlp)%rho = prop_mat(im,1)
                    mpii_mat(im_nlp)%lambda_elastic = prop_mat(im,2)
                    mpii_mat(im_nlp)%mu_elastic = prop_mat(im,3)   
                    mpii_mat(im_nlp)%mp_elastic = prop_mat(im,2) + 2.d0*prop_mat(im,3)   
                    mpii_mat(im_nlp)%VS = DSQRT(prop_mat(im,3)/prop_mat(im,1))
                    mpii_mat(im_nlp)%VP = DSQRT( (prop_mat(im,2) + 2*prop_mat(im,3))/prop_mat(im,1) )

                    dum1 =  (mpii_mat(im_nlp)%VP/mpii_mat(im_nlp)%VS)**2
                    mpii_mat(im_nlp)%Ni_elastic = (dum1- 2.d0)/ (2.0*(dum1- 1.d0))        !Poisson's Ratio

                    mpii_mat(im_nlp)%E_elastic = 2* mpii_mat(im_nlp)%mu_elastic * (1+mpii_mat(im_nlp)%Ni_elastic)
                    mpii_mat(im_nlp)%K_elastic   = mpii_mat(im_nlp)%E_elastic/(3.d0*(1.d0 - 2.d0*mpii_mat(im_nlp)%Ni_elastic)) 

                    if (damping_type .eq. 4) then
                        mpii_mat(im_nlp)%Mu_unrelax = visc_Mu(im)
                        mpii_mat(im_nlp)%Mp_unrelax = visc_Mp(im)
                        !mpii_mat(im_nlp)%lambda_viscoel = visc_Mp(im) - 2*visc_Mu(im)

                        ! For MPII, Only this is Implemented now : Viscoelastic Damping + Pressure Independednt + Total Stress Formulation
                        mpii_mat(im_nlp)%Gmax   = visc_Mu(im)
                        mpii_mat(im_nlp)%Emax   = 2.d0*visc_Mu(im)*(1+mpii_mat(im_nlp)%Ni_elastic)
                        mpii_mat(im_nlp)%Kmax   = mpii_mat(im_nlp)%Emax/(3.d0*(1.d0 - 2.d0*mpii_mat(im_nlp)%Ni_elastic)) 
                        mpii_mat(im_nlp)%lambdamax = visc_Mp(im) - 2*visc_Mu(im)
                        mpii_mat(im_nlp)%G_corr = mpii_mat(im_nlp)%Mu_unrelax                                           ! This changes for Effective stress formulation

                        if (nlp_pressdep_flag) then
                            ! Need to be updated in every time step
                            !mpii_mat(im_nlp)%G_corr = mpii_mat(im_nlp)%Mu_unrelax*abs(effective_mean_stress/effective_mean_stress_at_middle_of_soil_layer)
                            CALL EXIT()
                        endif
                        if (nlp_effstress_flag) then
                            call EXIT()
                        endif

                    elseif (damping_type.eq.2) then
                        mpii_mat(im_nlp)%Mu_unrelax = mpii_mat(im_nlp)%mu_elastic
                        mpii_mat(im_nlp)%Mp_unrelax = mpii_mat(im_nlp)%mp_elastic

                        ! For MPII, Only this is Implemented now : Viscoelastic Damping + Pressure Independednt + Total Stress Formulation
                        mpii_mat(im_nlp)%Gmax   = mpii_mat(im_nlp)%mu_elastic
                        mpii_mat(im_nlp)%Emax   = mpii_mat(im_nlp)%E_elastic
                        mpii_mat(im_nlp)%Kmax   = mpii_mat(im_nlp)%K_elastic
                        mpii_mat(im_nlp)%lambdamax = mpii_mat(im_nlp)%lambda_elastic
                        mpii_mat(im_nlp)%G_corr = mpii_mat(im_nlp)%mu_elastic                                           ! This changes for Effective stress formulation

                        if (nlp_pressdep_flag) then
                            ! Need to be updated in every time step
                            !mpii_mat(im_nlp)%G_corr = mpii_mat(im_nlp)%Mu_unrelax*abs(effective_mean_stress/effective_mean_stress_at_middle_of_soil_layer)
                            CALL EXIT()
                        endif
                        if (nlp_effstress_flag) then
                            call EXIT()
                        endif

                    else
                        CALL EXIT(EXIT_ANELASTIC)
                    endif

                endif
            enddo
        enddo

        max_nspr = 0
        ! Reading the Input G/Gmax curve data
        do im_nlp = 1,nmat_nlp
            do ifunc = 1,nfunc
                if (tag_func(ifunc) .eq. functag_mat_nlp(im_nlp)) then
                    if ((func_type(ifunc) .eq. 70) .or. (func_type(ifunc) .eq. 71)) then
                        mpii_mat(im_nlp)%nspring = int(func_data(func_indx(ifunc)))
                        allocate( mpii_mat(im_nlp)%spr_strain(mpii_mat(im_nlp)%nspring) )
                        allocate( mpii_mat(im_nlp)%spr_GbyGmax(mpii_mat(im_nlp)%nspring) )
                        allocate( mpii_mat(im_nlp)%spr_yldstress(mpii_mat(im_nlp)%nspring) )
                        allocate( mpii_mat(im_nlp)%spr_CNinv(mpii_mat(im_nlp)%nspring - 1) )
                    else
                        call EXIT(EXIT_FUNCTION_ERROR)
                    endif

                    if (mpii_mat(im_nlp)%nspring .gt. max_nspr) max_nspr = mpii_mat(im_nlp)%nspring

                    if (func_type(ifunc) .eq. 70) then

                        if (.NOT.nlp_pressdep_flag) then
                            ref_strain = func_data(func_indx(ifunc) + 1)

                            ! = ceiling(dlog10(ref_strain));
                            !strain_min = 10.d0**(-3+dum_log)
                            !strain_max = 10.d0**( 2+dum_log)

                            strain_min = 1.0d-6
                            strain_max = 1.0d-1                        
                            call nonlinear_mu_degradation_hyperbola(mpii_mat(im_nlp)%nspring, strain_min, strain_max, &
                                                            ref_strain, mpii_mat(im_nlp)%Gmax, &
                                                            mpii_mat(im_nlp)%spr_strain, mpii_mat(im_nlp)%spr_GbyGmax,&
                                                            mpii_mat(im_nlp)%spr_yldstress, mpii_mat(im_nlp)%spr_CNinv)
                        endif
                    elseif (func_type(ifunc) .eq. 71) then
                        call EXIT(EXIT_FUNCTION_ERROR)      ! Yet to be developed
                    endif
                    
                endif
            enddo
        enddo

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Allocating memory occupying variables for each element
        allocate(nlp_elem(nelem_loc))
        do ie = 1,nelem_loc
            im = con_spx_loc(con_spx_loc(ie -1)); 
            nn = sdeg_mat(im) +1

            do im_nlp = 1, nmat_nlp
                if (tag_mat(im) .eq. tag_mat_nlp(im_nlp)) then
                    allocate(nlp_elem(ie)%activefsur(nn,nn,nn))
                    allocate(nlp_elem(ie)%F2(nn,nn,nn,mpii_mat(im_nlp)%nspring))
                    allocate(nlp_elem(ie)%Sa1(nn,nn,nn,mpii_mat(im_nlp)%nspring,6))
                    nlp_elem(ie)%activefsur = 0;
                    nlp_elem(ie)%F2 = 0.d0;
                    nlp_elem(ie)%Sa1 = 0.000;
                endif
            enddo
        enddo

        do im_nlp = 1, nmat_nlp
            write(*,*) 'imat = ',im_nlp, 'Gmax = ',mpii_mat(im_nlp)%Gmax, 'ref_strain = ',ref_strain
            
            do ispr=1,(mpii_mat(im_nlp)%nspring-1)
                write(*,*) 'Strains = ',mpii_mat(im_nlp)%spr_strain(ispr) ,'G/Gmax = ', &
                    mpii_mat(im_nlp)%spr_GbyGmax(ispr), 'Yield Stress = ',mpii_mat(im_nlp)%spr_yldstress(ispr) ,&
                    'CNinv = ',  mpii_mat(im_nlp)%spr_CNinv(ispr)
            enddo
        enddo


        if (mpi_id.eq.0) write(*,'(A)')'Done.' 

    end subroutine MAKE_NLP_IWANSPRING_YIELDVALUES



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                        !!!
!!!  Interpolating G/Gmax Values with Strain at IWAN Spring Yield levels   !!!
!!!    (G/Gmax follows Hyperbolic Curve - Hardin and Drnevich (1972))      !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine   nonlinear_mu_degradation_hyperbola(nspr, stin_min, stin_max, ref_stin, mu_max, &
                                                strain, mu_deg_val, yldstress, cninv)

        implicit none

        integer*4 , intent(in) :: nspr
        integer*4 :: i

        real*8, intent(in) :: stin_min, stin_max, ref_stin, mu_max
        real*8, dimension(nspr), intent(out) :: strain, mu_deg_val, yldstress
        real*8, dimension(nspr-1), intent(out) :: cninv
        real*8 :: xl, xu, dx

        xl    = dlog10(stin_min)
        xu    = dlog10(stin_max)
        dx    = (xu-xl)/(nspr- 1)  

        ! Yield Stress of 1st spring is 0. Because this spring is made active in initial conditions (linear behaviour)
        do i = 1,nspr             
            strain(i) = 10.0d0**(xl + dx*(i- 1))
            mu_deg_val(i)     = ( 1.0d0/( 1.0d0 + dabs(strain(i)/ref_stin) ) )       !G/Gmax Value
            yldstress(i)     = (mu_deg_val(i)*mu_max) * strain(i)                     
        enddo

        call CN_INVERSE(nspr, mu_max, strain, yldstress, cninv)

     end subroutine nonlinear_mu_degradation_hyperbola


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                        !!!
!!!  Sum of Inverse of difference between Average 'G' and 'Gmax'           !!!
!!!             in range of each iwan spring                               !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     subroutine CN_INVERSE(nspr, mu_max, strain, yldstress, cninv)

        implicit none

        integer*4 , intent(in) :: nspr
        integer*4 :: i

        real*8, intent(in) :: mu_max
        real*8, dimension(nspr), intent(in) :: strain, yldstress
        real*8, dimension(nspr-1), intent(OUT) :: cninv

        real*8 :: sumdummy

        !
        cninv     = 0.0d0
        sumdummy  = 0.0d0
        do i = 1,nspr-1
            cninv(i) = ( 0.5d0*(strain(i+1)- strain(i)) / (yldstress(i+1)-yldstress(i)) ) &
                        - 0.5d0/mu_max - sumdummy
            sumdummy = sumdummy+ cninv(i)
        enddo

     end subroutine CN_INVERSE