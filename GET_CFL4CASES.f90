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

!> @brief Computes CFL condition when mechanical properties are given
!! node by node.
!! @author Ilario Mazzieri
!> @date September, 2013 - Creation
!> @version 1.0
!> @param[in,out]  time_step  deltat used for the time integration
!> @param[in] nn_loc number of local nodes
!> @param[in] loc_n_num local node numeration
!> @param[in] nm number of materials
!> @param[in] sdeg   polynomial degree vector 
!> @param[in] pm  material property vector (1-rho, 2-lambda, 3-mu, 4-zeta)
!> @param[in] xx_loc x-coordinate of local nodes
!> @param[in] yy_loc y-coordinate of local nodes
!> @param[in] zz_loc z-coordinate of local nodes
!> @param[in] cs_nnz_loc length of cs_loc
!> @param[in] cs_loc local connectivity vector
!> @param[in] time_step_cfl  maximum deltat available in order to avoid instability
!> @param[in] fmax  maximum frequency of the wave
!> @param[in] deltat_fixed  'yes' deltat fixed even if is extremely small for the time integration
!                     'not' deltat could be changed if is extremely small for the time integration
!> @param[in] mpi_comm  mpi common world
!> @param[in] mpi_np  total mpi processor number
!> @param[in] mpi_id  number of mpi process 
!> @param[in] vcase  case value 
!> @param[in] tcase  tag case value
!> @param[in] zs_elev local nodes elevation (z) from the surface
!> @param[in] zs_all local nodes elevation (z) from the alluvial basin
!> @param[in] sub_tag_all  see MAKE_NOT_HONORING.f90
!> @param[in] b_failCFL  fail if CFL does not hold
!> @param[in] damping_type 1-Kosloff&Kosloff 2-SLS 3-Rayleigh 
!> @param[in] QS quality factors for S-waves  
!> @param[in] QP quality factors for P-waves  
!> @note Compute the CFL condition (deltat <= deltat_max = h_min/(vp_max*N^2)) when the not honoring 
!!  strategy is applied (material properties given node by node) and print on the screen what is 
!!  the percent of deltat_max you are using for the time integration

      subroutine GET_CFL4CASES(time_step, nn_loc, loc_n_num, nm, tm, pm, sdeg, &
                             xx_loc, yy_loc, zz_loc, cs_nnz_loc, cs_loc, & 
                             time_step_cfl, fmax, deltat_fixed, &
                             mpi_comm, mpi_np, mpi_id, &
                             vcase, tcase, zs_elev, zs_all, vs_tria, thick_tria, sub_tag_all, b_failCFL, &
                             damping_type, QS, QP)
          
        use speed_exit_codes

        implicit none

        include 'SPEED.MPI'

        character*3 :: deltat_fixed

        integer*4 :: nm, im, nn_loc, cs_nnz_loc, mpi_comm, mpi_np, mpi_err, mpi_id
        integer*4 :: ie, i, j, k, nn, mcode, smcode, sdeg_deltat_cfl, sdeg_npoints
        integer*4 :: tcase, vcase, damping_type
        integer*4 :: nel_loc, ic1, ic2, n1, n2, istart, ifin

       integer*4, dimension(1) :: pos
        integer*4, dimension(nm) :: tm, sdeg
        integer*4, dimension(nn_loc) :: loc_n_num
        integer*4, dimension(nn_loc) :: sub_tag_all        
        integer*4, dimension(0:cs_nnz_loc) :: cs_loc

        real*8 :: time_step, time_step_cfl, fmax, qs_loc, qp_loc
        real*8 :: length_min,length,vs_length_min,vs_length,length_vp_min,length_vp
        real*8 :: percent_deltat
        real*8 :: num_of_points,vs_npoints,vp_deltat_cfl
        real*8 :: rho,lambda,mu,vs,vp
        real*8 :: x1_length_vp_min,y1_length_vp_min,z1_length_vp_min
        real*8 :: x2_length_vp_min,y2_length_vp_min,z2_length_vp_min

        real*8, dimension(:), allocatable :: ct,ww
        real*8, dimension(:), allocatable :: time_step_glo
        real*8, dimension(nn_loc) :: xx_loc, yy_loc, zz_loc
        real*8, dimension(nn_loc) :: zs_elev
        real*8, dimension(nn_loc) :: zs_all
        real*8, dimension(nn_loc) :: vs_tria, thick_tria
        

        real*8, dimension(nm,4) :: pm
        real*8, dimension(:,:), allocatable :: dd
        real*8, dimension(nm) :: QS, QP

        real*8, dimension(:,:,:), allocatable :: rho_el, lambda_el, mu_el, gamma_el 

        logical :: b_failCFL, b_CFL_failure = .FALSE.

        nel_loc = cs_loc(0) - 1

        length_vp_min = 1.d30
        vs_length_min = 1.d30
       

    
          
         do ie = 1, nel_loc
                im = cs_loc(cs_loc(ie -1) +0)
            
                nn = sdeg(im) +1
                allocate(rho_el(nn,nn,nn), lambda_el(nn,nn,nn),mu_el(nn,nn,nn),gamma_el(nn,nn,nn)) 
        
                rho_el = pm(im,1) 
                lambda_el = pm(im,2)
                mu_el = pm(im,3)
                gamma_el = pm(im,4)        
                QS(im) = qs_loc
                QP(im) = qp_loc
            
                if (vcase.eq.tm(im)) then
                  call MAKE_ELTENSOR_FOR_CASES(tcase, vcase,&
                                              nn, rho_el, lambda_el, mu_el, gamma_el,&  
                                              nn_loc, zs_elev, zs_all, vs_tria, thick_tria,&
                                              cs_nnz_loc, cs_loc, ie,&
                                              sub_tag_all, zz_loc, mpi_id, loc_n_num, &
                                              damping_type, qs_loc, qp_loc, &
                                              xx_loc, yy_loc,0)
                                                                                                                     
                                              
                                              
                endif
      
      

                smcode=sdeg(im) 
                rho = 0.d0
                lambda = 0.d0
                mu = 0.d0
                      
                do i = 1, nn
                  do j = 1, nn
                     do k = 1, nn
                        rho = rho + rho_el(i,j,k)   
                        lambda = lambda + lambda_el(i,j,k) 
                        mu = mu + mu_el(i,j,k)
                     enddo
                  enddo
                enddo
               
                rho = rho/(nn**3)
                lambda = lambda/(nn**3)
                mu = mu/(nn**3)
                   
                vs = dsqrt(mu/rho)
                vp = dsqrt((lambda+2*mu)/rho)

                n1 = nn*nn*(1 -1) +nn*(1 -1) +1
                n2 = nn*nn*(nn -1) +nn*(nn -1) +nn
                istart = cs_loc(ie -1) +n1
                ifin = cs_loc(ie -1) +n2

               
               do ic1 = istart, ifin 
                  do ic2 = ic1 + 1, ifin               
               
                  if(cs_loc(ic1) .ne. 0 .and. cs_loc(ic2) .ne. 0) then
               
                     length=dsqrt((xx_loc(cs_loc(ic1))-xx_loc(cs_loc(ic2)))**2 + &
                                  (yy_loc(cs_loc(ic1))-yy_loc(cs_loc(ic2)))**2 + &
                                  (zz_loc(cs_loc(ic1))-zz_loc(cs_loc(ic2)))**2)
                        
                     length_vp = length/vp
                    if (length_vp_min .gt. length_vp) then
                        length_vp_min = length_vp
                        x1_length_vp_min = xx_loc(cs_loc(ic1))
                        y1_length_vp_min = yy_loc(cs_loc(ic1))
                        z1_length_vp_min = zz_loc(cs_loc(ic1))
                        x2_length_vp_min = xx_loc(cs_loc(ic2))
                        y2_length_vp_min = yy_loc(cs_loc(ic2))
                        z2_length_vp_min = zz_loc(cs_loc(ic2))
                        vp_deltat_cfl = vp
                        sdeg_deltat_cfl = smcode
                        length_min = length
                     endif


                     vs_length = vs/length

                     if (vs_length_min.gt.vs_length) then
                        vs_length_min = vs_length
                        vs_npoints = vs
                        sdeg_npoints = smcode
                     endif
                                                    
                endif
 
               enddo
             enddo   
      
         deallocate(lambda_el, gamma_el, mu_el, rho_el)
         
        enddo
        
        nn = sdeg_deltat_cfl + 1

        time_step_cfl=length_vp_min !/((nn-1)**2)

     allocate(time_step_glo(mpi_np))
        
        call MPI_BARRIER(mpi_comm, mpi_err)
        call MPI_ALLGATHER(time_step_cfl, 1, SPEED_DOUBLE, time_step_glo, 1, SPEED_DOUBLE, mpi_comm, mpi_err)

            
        pos = minloc(time_step_glo)
         
       deallocate(time_step_glo)

        if((mpi_id + 1) .eq. pos(1)) then

           write(*,'(A)')' '
           write(*,'(A)') '--------Stability for time advancing algorithm --------'
           write(*,'(A,E12.4)') 'Min. el. length :', length_vp_min*vp_deltat_cfl
           write(*,'(A,E12.4)') 'Min. Vp         :', vp_deltat_cfl
           write(*,'(A,E12.4,E12.4,E12.4)')'Min. node 1     : ',x1_length_vp_min,y1_length_vp_min,z1_length_vp_min
           write(*,'(A,E12.4,E12.4,E12.4)')'Min. node 2     : ',x2_length_vp_min,y2_length_vp_min,z2_length_vp_min
           write(*,'(A)') '-------------------------------------------------------'
           if (time_step.le.time_step_cfl) then
                percent_deltat=time_step/time_step_cfl*100
                if (percent_deltat.le.1) then
                        write(*,'(A,E12.4)')'WARNING!!! deltat CFL =',&
                        time_step_cfl
                        write(*,'(A,F6.2,A)')'You are using ',percent_deltat,&
                                '% of critical time step.'
                        write(*,'(A)')'This time step is excedingly lower the deltat CFL'
                        if (deltat_fixed.eq.'not') then  
                                write(*,'(A)')'deltat chosen will be substituted with 10% of deltat CFL'
                                time_step=time_step_cfl*0.1
                                write(*,'(A,E12.4)')'deltat chosen :',time_step
                        endif
                elseif (percent_deltat.le.25) then
                        write(*,'(A,E12.4)')'OK!!! deltat CFL =',&
                        time_step_cfl
                        write(*,'(A,F6.2,A)')'You are using ',percent_deltat,&
                                '% of critical time step.'
                else
                        write(*,'(A,E12.4)')'WARNING!!! deltat CFL =',&
                        time_step_cfl
                        write(*,'(A,F6.2,A)')'You are using ',percent_deltat,&
                                '% of critical time step.'
                        write(*,'(A)')'This could be not enough for the correct working of ABC'
                        if (deltat_fixed.eq.'not') then
                                write(*,'(A)')'deltat chosen will be substituted with 20% of deltat CFL'
                                time_step=time_step_cfl*0.2
                                write(*,'(A,E12.4)')'deltat chosen :',time_step
                        endif
                endif
           elseif (time_step.gt.time_step_cfl) then
                write(*,'(2X,A,E12.4)')'ERROR!!! deltat CFL = ',&
                time_step_cfl,' < deltat = ',time_step
                write(*,'(A)')'The time advancing scheme will be unstable!'
                if (deltat_fixed.eq.'not') then
                        write(*,'(A)')'deltat chosen will be substituted with 20% of deltat CFL'
                        time_step=time_step_cfl*0.2
                        write(*,'(A,E12.4)')'deltat chosen :',time_step
                else
                    b_CFL_failure = .TRUE.
                endif
           endif
           write(*,'(A)')'-------------------------------------------------------'
           write(*,'(A)')' '
        
           if (b_failCFL .and. b_CFL_failure) then
              write(*,*) 'CFL does NOT hold, asked to quit. Program finished.'
              call MPI_FINALIZE(mpi_err)
              call EXIT(EXIT_CFL)
           endif

        !Writing of the results for the maximum number of points per wave length
        
          allocate(ct(nn),ww(nn),dd(nn,nn))
          call MAKE_LGL_NW(nn,ct,ww,dd)

          num_of_points = vs_length_min*(ct(1)-ct(nn))/(ct(int(nn/2))-ct(int(nn/2)+1))/fmax

          write(*,'(A)')'-----------Number of points per wave length------------'
          write(*,'(A,E12.4)')'Max. el. length       :',1/vs_length_min*vs_npoints
          write(*,'(A,E12.4)')'Max. Vs               :',vs_npoints
          write(*,'(A,E12.4)')'Points per wavelength :',num_of_points
          write(*,'(A)') ' '

          deallocate(ct,ww,dd)
        endif   
        
        call MPI_BARRIER(mpi_comm, mpi_err)
        
        return
               
        end subroutine GET_CFL4CASES

