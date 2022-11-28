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
!! @author Srihari Sangaraju
!> @date August, 2022 
!> @version 1.0
!> Inputs:
! its = timestep number
    ! This Following SubRoutine, calculates the Stress tensor at a LGL node, 
    ! based on strain tensor at that node and IWAN's MPII speing model.
     subroutine   MAKE_STRESS_TENSOR_IWAN(its, nnloc, nn_cs_loc, cs_loc, & 
                                        ie, nn, mpii_lay, nspr, nlp_pressdep_flag, &
                                        oldstrain, oldstress, &
                                        strainxx_el, strainyy_el, strainzz_el, &
                                        strainxy_el, strainyz_el, strainzx_el, &
                                        Sa1, F2, active_fsur, &
                                        sxx_el, syy_el, szz_el, &                                   ! Stress Tensor - increment (output)
                                        sxy_el, syz_el, szx_el)
                                        ! detxx_el, detyy_el, detzz_el, &                             ! increment in Total Strain
                                        ! detxy_el, detyz_el, detzx_el, &
                                        ! SS1xx_el, SS1yy_el, SS1zz_el, &
                                        ! SS1xy_el, SS1yz_el, SS1zx_el, &
                                        ! SSa1xx_el, SSa1yy_el, SSa1zz_el, &
                                        ! SSa1xy_el, SSa1yz_el, SSa1zx_el, &
                                        ! FF2_el, asurf_el, &                                         ! Iwan Springs (shifted) yield value, and active iwan spring 
                                        ! dpxx_el, dpyy_el, dpzz_el, dpxy_el, dpyz_el, dpzx_el, &     ! Increment in Plastic Strain

        use mpii_material_def

        implicit none

        type(mpii_material), intent(in) :: mpii_lay

        integer*4, intent(in)                         :: nnloc, nn_cs_loc
        integer*4, intent(in)                         :: ie, its, nn, nspr
        integer*4, dimension(0:nn_cs_loc), intent(in) :: cs_loc
        integer*4, dimension(nn,nn,nn), intent(inout) :: active_fsur

        logical, intent(in) :: nlp_pressdep_flag

        real*8, dimension(6*nnloc), intent(in)    :: oldstrain, oldstress
        real*8, dimension(nn,nn,nn), intent(in) ::  strainxx_el, strainyy_el, strainzz_el, &
                                                    strainxy_el, strainyz_el, strainzx_el
        
        real*8, dimension(nn,nn,nn),        intent(out)   :: sxx_el, syy_el, szz_el, &
                                                             sxy_el, syz_el, szx_el  
        real*8, dimension(nn,nn,nn,nspr),   intent(inout) :: F2
        real*8, dimension(nn,nn,nn,nspr,6), intent(inout) :: Sa1
        

        ! real*8, dimension(nn,nn,nn) ::  detxx_el, detyy_el, detzz_el, &
        !                                 detxy_el, detyz_el, detzx_el
        ! real*8, dimension(nn,nn,nn)::  SS1xx_el, SS1yy_el, SS1zz_el, &
        !                                SS1xy_el, SS1yz_el, SS1zx_el
                                       !dpxx_el, dpyy_el, dpzz_el, dpxy_el, dpyz_el, dpzx_el
        ! real*8, dimension(nn,nn,nn,nspr_max) ::  SSa1xx_el, SSa1yy_el, SSa1zz_el, &
        !                                          SSa1xy_el, SSa1yz_el, SSa1zx_el, FF2_el

        ! Variables Used only in this subroutine
        integer*4 :: i, j, k, is, in, iaz, k11

        real*8 :: node_epsm, node_dsigm
        real*8, dimension(6)      :: node_dstress, node_deps, node_S1 !eps = Total Strain; 
        real*8, dimension(nspr)   :: node_F2
        real*8, dimension(nspr, 6):: node_Sa1

        do k = 1,nn
            do j= 1,nn
                do i= 1,nn
                    is = nn*nn*(k -1) + nn*(j -1) + i
                    in = cs_loc(cs_loc(ie -1) + is)
                    iaz = 6*(in -1)


                    ! nodal stress/strain tensor alignment = [xx yy zz xy yz zx]
                    node_dstress = 0;

                    ! Increment in Total Strain
                    node_deps(1) = strainxx_el(i,j,k) - oldstrain(iaz+1);
                    node_deps(2) = strainyy_el(i,j,k) - oldstrain(iaz+2);
                    node_deps(3) = strainzz_el(i,j,k) - oldstrain(iaz+3);
                    node_deps(4) = strainxy_el(i,j,k) - oldstrain(iaz+4);
                    node_deps(5) = strainyz_el(i,j,k) - oldstrain(iaz+5);
                    node_deps(6) = strainzx_el(i,j,k) - oldstrain(iaz+6);
                    node_epsm   = (node_deps(1) + node_deps(2)+node_deps(3))/ 3.d0   ! Hydrostatic/ mean strain

                    ! Variables for Vonmises Criteria - Deviatoric Stress Tensor (previous Time step)
                    node_dsigm   = mpii_lay%Kmax * node_epsm* 3.d0      !Volumetric/Mean stress
                    node_S1(1) = oldstress(iaz+1) - node_dsigm; 
                    node_S1(2) = oldstress(iaz+2) - node_dsigm; 
                    node_S1(3) = oldstress(iaz+3) - node_dsigm; 
                    node_S1(4) = oldstress(iaz+4); 
                    node_S1(5) = oldstress(iaz+5); 
                    node_S1(6) = oldstress(iaz+6); 

                    ! Variables for Vonmises Criteria - Hardening rule - origin/back stress
                    node_Sa1(1:nspr,1:6) = Sa1(i,j,k,1:nspr,1:6);

                    ! Von Mises Stress
                    node_F2(1:nspr) = F2(i,j,k,1:nspr);
                    
                    ! if ( (ie.eq.5276) .and. (i.eq.3) .and. (j.eq.1) .and. (k.eq.3) )  then
                    !     write(*,*) 'input dstrain ', (node_deps(k11), k11=1,6)
                    ! endif

                    call NEOIWAN(its, nspr, mpii_lay%lambdamax, mpii_lay%Gmax, mpii_lay%G_corr, mpii_lay%Kmax,&
                                node_deps, node_epsm, node_S1, node_Sa1, node_F2, &
                                mpii_lay%spr_yldstress, mpii_lay%spr_CNinv, active_fsur(i,j,k), node_dstress, ie, i,j,k)

                    !Arrainging Total Stresses (Stress_new = Stress_old + stressIncrement)
                    sxx_el(i,j,k) = oldstress(iaz+1) + node_dstress(1);  sxy_el(i,j,k) = oldstress(iaz+4) + node_dstress(4);
                    syy_el(i,j,k) = oldstress(iaz+2) + node_dstress(2);  syz_el(i,j,k) = oldstress(iaz+5) + node_dstress(5);
                    szz_el(i,j,k) = oldstress(iaz+3) + node_dstress(3);  szx_el(i,j,k) = oldstress(iaz+6) + node_dstress(6);

                    ! if ( (ie.eq.5276) .and. (i.eq.3) .and. (j.eq.1) .and. (k.eq.3) )  then
                    !     write(*,*) 'output dstress = ', (node_dstress(k11), k11=1,6)
                    ! endif

                    ! Plastic strain increment array (gamma for shear not epsilon)
                    ! if (nlp_pressdep_flag) then
                    !     dpxx_el(i,j,k) = node_deps(1)- (node_dstress(1) - mpii_lay%lambdamax*3.d0*node_epsm)/ (2.d0*mpii_lay%G_corr)
                    !     dpyy_el(i,j,k) = node_deps(2)- (node_dstress(2) - mpii_lay%lambdamax*3.d0*node_epsm)/ (2.d0*mpii_lay%G_corr)
                    !     dpzz_el(i,j,k) = node_deps(3)- (node_dstress(3) - mpii_lay%lambdamax*3.d0*node_epsm)/ (2.d0*mpii_lay%G_corr)
                    !     dpxy_el(i,j,k) = node_deps(4)-  node_dstress(4)/ (2.d0*mpii_lay%Gmax)
                    !     dpyz_el(i,j,k) = node_deps(5)-  node_dstress(5)/ (2.d0*mpii_lay%Gmax)
                    !     dpzx_el(i,j,k) = node_deps(6)-  node_dstress(6)/ (2.d0*mpii_lay%Gmax)
                    ! endif

                enddo
            enddo
        enddo


end subroutine MAKE_STRESS_TENSOR_IWAN



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                        !!!
!!!                             NeoIWAN                                    !!!
!!!                                                                        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine NEOIWAN(its, nspr, lambdamax, Gmax, Gini, Kmax, deps,depsm, &
                        SS1, SSa1, FF2, R, CNinv, aktifsur, dSigma, ie1, i1, j1, k1)

        implicit none

        integer*4, intent(in) :: its, nspr
        integer*4, intent(inout) :: aktifsur

        real*8, intent(in)                      :: lambdamax, Gmax, Gini, Kmax
        real*8, dimension(6), intent(in)        :: deps
        real*8, dimension(nspr), intent(in)     :: R
        real*8, dimension(nspr-1), intent(in)   :: CNinv

        real*8, dimension(nspr), intent(inout)  :: FF2
        real*8, dimension(6), intent(inout)     :: dSigma, SS1
        real*8, dimension(nspr,6), intent(inout):: SSa1

        ! Other Variables used only inside the subroutine
        integer*4   :: i, j, m, ie1, i1, j1, k1, k11
        integer*4   :: start1, aktif, D, errorflag
        integer*4, dimension(6) :: INDX

        real*8      :: depsm, dsigm
        real*8, dimension(6) :: de, dS
        real*8, dimension(nspr) :: dF
        real*8, dimension(6,6) :: Ed    !, Esd


        !write(*,*) 'Element No, Node No: ',ie1, i1, j1 ,k1
        ! TotalStress = deviotoricStress + meanStress
        ! Mean Stress and Strain calculation
        !depsm   = (deps(1) + deps(2)+deps(3))/ 3.d0    
        dsigm   = Kmax* depsm* 3.d0                 

        ! Emax and Kmax are used only for principal directions 
        ! and held independent of saturation 
        
        ! For First Iteration - Finding Elastic Stresses
        ! Then Finding Deviatoric Stress Tensor (dS)
        ! For Pressure-Independent Formulation, initial value of SS1 is alomst eqial to zero.
        if (its.eq.0) then
            call IWAN_elastic(Gini, Gmax, lambdamax, deps, dSigma) 

            dS    = dSigma
            dS(1) = dSigma(1)- dsigm
            dS(2) = dSigma(2)- dsigm
            dS(3) = dSigma(3)- dsigm   

            SS1 = SS1 + dS 
            SSa1 = 0.0
            !dplast  = 0.0
            aktifsur  = 0
            ! if ( (ie1.eq.5276) .and. (i1.eq.3) .and. (j1.eq.1) .and. (k1.eq.3) )  write(*,*) 'aktif = ', aktif
            return
        endif

        ! Deviatoric Strain
        de = 0.d0
        de(1) = de(1) - depsm
        de(2) = de(2) - depsm
        de(3) = de(3) - depsm
        de(4) = deps(4)/2.d0
        de(5) = deps(5)/2.d0
        de(6) = deps(6)/2.d0
        !de(1:3) = deps(1:3) - depsm
        !de(4:6) = deps(4:6)
        

        ! If First spring is not yielded - Find VonMises Yield stress Condition (FF2(1)),
        ! Corresponding to Total Deviatoric Stress (SS1)??, and origin/back stress of kinematic hardening rule (SSa1)?? 
        if (aktifsur .eq. 0) then
            call IWAN_surface (SS1, SSa1(1,:), FF2(1))
            call IWAN_dsurface(SS1, SSa1(1,:), de, dF(1))
            ! if ( (ie1.eq.5276) .and. (i1.eq.3) .and. (j1.eq.1) .and. (k1.eq.3) )  then
            !     write(*,*) 'SS1 = ', (SS1(k11), k11=1,6)
            !     write(*,*) 'SS1 = ', (SSa1(1,k11), k11=1,6)
            !     write(*,*) 'FF2(1) = ', FF2(1), 'dF(1) = ', dF(1)
            ! endif
        endif

        aktif = 0
        if (aktifsur .GT. 0) then
            do j = 1,aktifsur 
                do m = 1,6
                    SSa1(j,m) = SS1(m)- R(j)/ sqrt(FF2(j))*(SS1(m)-SSa1(j,m))
                enddo

                call IWAN_surface (SS1,SSa1(j,:),FF2(j))
                call IWAN_dsurface(SS1,SSa1(j,:),de,dF(j))
                if ( (dF(j) .GE. 0.0)  .AND. (FF2(j) .GE. R(j)**2) ) aktif = aktif+ 1
            enddo
        endif
        
        ! Condition where, first IWAN spring has not yet Yielded
        if ( FF2(1) .lt. R(1)**2.d0  .AND. dF(1) .ge. 0.d0)   then
            call IWAN_elastic(Gini, Gmax, lambdamax, deps, dSigma) 
          
            dS    = dSigma
            dS(1) = dSigma(1)- dsigm
            dS(2) = dSigma(2)- dsigm
            dS(3) = dSigma(3)- dsigm   
          
            SS1 = SS1+ dS 
            !dplast = 0.0
          
            return    
        endif

        ! Coefficient Matrix (Ed) computation 
        ! Ed * dS (deviatoric stress increment) = de (deviatoric strain increment)
        ! dS is unknown here
        Ed    = 0.d0
        do i = 1,6
            Ed(i,i) = 0.5/Gmax
        enddo

        ! Finiding if next spring is yielded or not, and uptating Coefficient Matrix (Ed)
        start1 = 1
        if (aktif .GT. 0) call Ematris(start1, nspr, Ed, CNinv, SS1, SSa1, FF2, aktif)
        do j = aktif+1, nspr-1
            call IWAN_surface (SS1, SSa1(j,:), FF2(j) )
            call IWAN_dsurface(SS1, SSa1(j,:), de, dF(j) )

            if ( (dF(j) .GE. 0.d0)  .AND. (FF2(j) .GE. R(j)**2.d0) ) then 
                aktif = aktif+ 1
                start1 = aktif

                call Ematris(start1, nspr, Ed, CNinv, SS1, SSa1, FF2, aktif)
            else 
                EXIT
            endif
        enddo

        ! Solving for stresses (increment of stresses)
        ! Page No. 6-9 of Joyner 1975.
        call LUDCMP(Ed, 6, INDX, D, errorflag)
        if (errorflag .ne. 0) stop 'not invertible matrix'
        call LUBKSB(Ed, 6, INDX, de) ! solve Ed·x = de (de is used as input/ouput.... during the output X values are stored in de variable)
        dS = de

        SS1 = SS1+ dS 
        aktifsur = max(aktif,1)

        ! Total Stresses
        dSigma(1:3) = dS(1:3) + dsigm
        dSigma(4:6) = dS(4:6)

    end subroutine NEOIWAN


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                        !!!
!!!                            IWAN ELASTIC                                !!!
!!!                                                                        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine IWAN_elastic (Ginitial,Gmax,lambda,deps,dsigma)  

        real*8, intent(in) :: Ginitial, Gmax, lambda
        real*8, intent(in) :: deps(6)
        real*8, intent(out):: dsigma(6)
      
        real*8 :: depsvol
      
        depsvol = deps(1)+ deps(2)+ deps(6)
        
        ! Stress - Strain Relationship for Linear Isotropic Material
        dsigma(1) = 2.d0*Ginitial*deps(1)+ lambda*depsvol
        dsigma(2) = 2.d0*Ginitial*deps(2)+ lambda*depsvol
        dsigma(3) = 2.d0*Ginitial*deps(3)+ lambda*depsvol
        ! dsigma(4) = 2.d0*Gmax * deps(4)
        ! dsigma(5) = 2.d0*Gmax * deps(5)
        ! dsigma(6) = 2.d0*Gmax * deps(6)
        dsigma(4) = Gmax * deps(4)
        dsigma(5) = Gmax * deps(5)
        dsigma(6) = Gmax * deps(6)
        
      
    end subroutine IWAN_elastic  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Von mises yield criteria + Preger kinematic hardening 
    ! Yield surface calculation corrwsponding to a Iwan Spring computation 
    ! Vonmises stress of deviatoric part of (cauchy stress tensor - hardening parameter stress)
    ! Formula for Tn - Page No. 8 of Joyner 1975
    subroutine IWAN_surface(S, Sa, F)  
  
        real*8, intent(in) :: S(6)
        real*8, intent(in) :: Sa(6)
        real*8, intent(inout) :: F
       
        F = 0.d0
        F = 0.5d0* ( (S(1)-Sa(1))**2.d0 + (S(2)-Sa(2))**2.d0+ (S(3)-Sa(3))**2.d0 + &
                        2.d0* (S(4)-Sa(4))**2.d0 + &
                        2.d0* (S(5)-Sa(5))**2.d0 + &
                        2.d0* (S(6)-Sa(6))**2.d0)
      
    end subroutine IWAN_surface  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Iwan surface(s) movement computation 
    ! Not So clear - Ask Fabian
    ! dF = 0.5*Strain Work?? (to know, if the current yield surface is still active or not?)
    subroutine IWAN_dsurface(S,Sa,de,dF) 
  
        real*8, intent(in) :: S(6)
        real*8, intent(in) :: Sa(6)
        real*8, intent(in) :: de(6)
        real*8, intent(inout) :: dF
       
        dF = 0.0d0
        dF = 0.5d0*( (S(1)-Sa(1))*de(1) + (S(2)-Sa(2))*de(2) + (S(3)-Sa(3))*de(3) + &
                        2.d0* (S(4)-Sa(4))*de(4) + &
                        2.d0* (S(5)-Sa(5))*de(5) + &
                        2.d0* (S(6)-Sa(6))*de(6) )
      
      end subroutine IWAN_dsurface  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Ematris(start1,Nspr,Ed,CNinv,S1,Sa1,F1,aktif)

        integer*4, intent(IN)     :: start1
        integer*4, intent(IN)     :: Nspr
        real*8,  INTENT(INOUT)  :: Ed(6,6)
        real*8,    INTENT(IN)     :: CNinv(Nspr-1)
        real*8,    INTENT(IN)     :: S1(6)
        real*8,    INTENT(IN)     :: Sa1(Nspr,6)
        real*8,    INTENT(IN)     :: F1 (Nspr)
        integer*4, intent(IN)     :: aktif
        
        integer*4 :: j,m,k
        real*8    :: ss(6)
        
        ss(1) = 1.d0;
        ss(2) = 1.d0;
        ss(3) = 1.d0;
        ss(4) = 2.d0;
        ss(5) = 2.d0;
        ss(6) = 2.d0;
        
        j = start1
        do while (j .lt. aktif+1)
            do m = 1,6
                do k = 1,6    
                    Ed(m,k) = Ed(m,k)+ CNinv(j)* ss(k)* (S1(m)-Sa1(j,m))&
                                *(S1(k)-Sa1(j,k))/ (2.d0* F1(j))
                enddo
            enddo
            j = j+1
        enddo
    
    end subroutine Ematris

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !  ***************************************************************
    !  * Given an N x N matrix A, this routine replaces it by the LU *
    !  * decomposition of a rowwise permutation of itself. A and N   *
    !  * are input. INDX is an output vector which records the row   *
    !  * permutation effected by the partial pivoting; D is output   *
    !  * as -1 or 1, depending on whether the number of row inter-   *
    !  * changes was even or odd, respectively. This routine is used *
    !  * in combination with LUBKSB to solve linear equations or to  *
    !  * invert a matrix. Return code is 1, if matrix is singular.   *
    !  ***************************************************************
     
     Subroutine LUDCMP(A,N,INDX,D,CODE)
        IMPLICIT NONE
        integer*4, parameter :: nmax = 100
        real*8, parameter :: tiny = 1.5D-16
   
        real*8, intent(inout), dimension(N,N) :: A
        integer*4, intent(in) :: N
        integer*4, intent(out) :: D, CODE
        integer*4, intent(out), dimension(N) :: INDX
        !f2py depend(N) A, indx
   
        REAL*8  :: AMAX, DUM, SUMM, VV(NMAX)
        INTEGER*4 :: I,J,K,IMAX
   
        D=1; CODE=0
   
        DO I=1,N
          AMAX=0.d0
          DO J=1,N
            IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
          END DO ! j loop
          IF(AMAX.LT.TINY) THEN
            CODE = 1
            RETURN
          END IF
          VV(I) = 1.d0 / AMAX
        END DO ! i loop
   
        DO J=1,N
          DO I=1,J-1
            SUMM = A(I,J)
            DO K=1,I-1
              SUMM = SUMM - A(I,K)*A(K,J) 
            END DO ! k loop
            A(I,J) = SUMM
          END DO ! i loop
          AMAX = 0.d0
          DO I=J,N
            SUMM = A(I,J)
            DO K=1,J-1
              SUMM = SUMM - A(I,K)*A(K,J) 
            END DO ! k loop
            A(I,J) = SUMM
            DUM = VV(I)*DABS(SUMM)
            IF(DUM.GE.AMAX) THEN
              IMAX = I
              AMAX = DUM
            END IF
          END DO ! i loop  
          
          IF(J.NE.IMAX) THEN
            DO K=1,N
              DUM = A(IMAX,K)
              A(IMAX,K) = A(J,K)
              A(J,K) = DUM
            END DO ! k loop
            D = -D
            VV(IMAX) = VV(J)
          END IF
   
          INDX(J) = IMAX
          IF(DABS(A(J,J)) < TINY) A(J,J) = TINY
   
          IF(J.NE.N) THEN
            DUM = 1.d0 / A(J,J)
            DO I=J+1,N
              A(I,J) = A(I,J)*DUM
            END DO ! i loop
          END IF 
        END DO ! j loop
   
        RETURN
    END subroutine LUDCMP
   
   
   !  *****************************************************************
   !  Solves the set of N linear equations A . X = B.  Here A is     *
   !  input, not as the matrix A but rather as its LU decomposition, *
   !  determined by the routine LUDCMP. INDX is input as the permuta-*
   !  tion vector returned by LUDCMP. B is input as the right-hand   *
   !  side vector B, and returns with the solution vector X. A, N and*
   !  INDX are not modified by this routine and can be used for suc- *
   !  cessive calls with different right-hand sides. This routine is *
   !  also efficient for plain matrix inversion.                     *
   !  *****************************************************************
   
    Subroutine LUBKSB(A, N, INDX, B)
        integer*4, intent(in) :: N 
        real*8, intent(in), dimension(N,N) :: A
        integer*4, intent(in), dimension(N) :: INDX
        real*8, intent(inout), dimension(N) :: B
        !f2py depend(N) A, INDX, B
    
        REAL*8  SUMM
        integer*4 :: II, LL,J,I
    
        II = 0
    
        DO I=1,N
        LL = INDX(I)
        SUMM = B(LL)
        B(LL) = B(I)
        IF(II.NE.0) THEN
            DO J=II,I-1
            SUMM = SUMM - A(I,J)*B(J)
            END DO ! j loop
        ELSE IF(SUMM.NE.0.d0) THEN
            II = I
        END IF
        B(I) = SUMM
        END DO ! i loop
    
        DO I=N,1,-1
        SUMM = B(I)
        IF(I < N) THEN
            DO J=I+1,N
            SUMM = SUMM - A(I,J)*B(J)
            END DO ! j loop
        END IF
        B(I) = SUMM / A(I,I)
        END DO ! i loop
    
        RETURN
    END subroutine LUBKSB