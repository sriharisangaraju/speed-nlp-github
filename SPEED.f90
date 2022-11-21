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

!> @brief SPEED (SPectral Elements in Elastodynamics with discontinuous Galerkin) 
!! is an open-source code for the simulation of seismic wave propagation in 
!! three-dimensional complex media. SPEED is jointly developed by MOX (The Laboratory for Modeling and Scientific 
!! Computing, Department of Mathematics) and DICA (Department of Civil and Environmental Engineering)
!! at Politecnico di Milano.
!> @see Website http://speed.mox.polimi.it
!> @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0

! Here starts the code SPEED

      program SPEED

      use max_var
      use str_mesh 
      use str_mesh_scratch                   
      use DGJUMP
      use speed_par
      use speed_par_dg

      implicit none
 
      include 'SPEED.MPI'
      

!*****************************************************************************************************************
      
      allocate(mpi_stat(SPEED_STATUS_SIZE))
      mpi_comm = SPEED_COMM
      
      call INITIALIZATION(mpi_comm,mpi_np,mpi_id,mpi_ierr)

      start = MPI_WTIME()
      
      if (mpi_id.eq.0) then
         write(*,'(A)')''
         write(*,'(A)')''
         write(*,'(A)')'*******************************************************'
         write(*,'(A)')'*                                                     *'
         write(*,'(A)')'*                        SPEED                        *'
         write(*,'(A)')'*          SPectral Elements in Elastodynamics        *'
         write(*,'(A)')'*              with Discontinuous Galerkin            *'
         write(*,'(A)')'*                                                     *'
         write(*,'(A)')'*          © PoliMi, 2012, All Rights Reserved        *'
         write(*,'(A)')'*                                                     *'
         write(*,'(A)')'*******************************************************'
         write(*,'(A)')''
         write(*,'(A)')''
      endif

!***********************************************************************************
!     READ INPUT FILES AND ALLOCATE VARIABLES
!***********************************************************************************


     call READ_INPUT_FILES()

     
!***********************************************************************************
!     PARTITION OF THE GRID AND GENERATION OF LOCAL CONNECTIVITY 
!***********************************************************************************


     allocate(elem_domain(nelem));  elem_domain = mpi_id;
     call MAKE_PARTITION_AND_MPI_FILES()


!***********************************************************************************
!    START SETUP 4 MPI COMUNICATION (ONLY BETWEEN CONFORMING REGIONS)
!***********************************************************************************


     call MAKE_SETUP_MPI_CONFORMING()

                   
!***********************************************************************************
!    NEW SPECTRAL LOCAL CONNECTIVITY
!***********************************************************************************


     call MAKE_SPX_GRID_WITH_LOC_NUMERATION()

                                      
!*****************************************************************************************
!    BUILD SEISMIC MOMENT
!*****************************************************************************************

      if (srcmodflag.eq.0) then
          call MAKE_SEISMIC_MOMENT_OR_EXPLOSIVE_SOURCE()
      elseif (srcmodflag.eq.1) then
          call MAKE_SEISMIC_MOMENT_POINTSOURCES_NHFAULT()
      endif


!*****************************************************************************************
!    DAMPING 
!*****************************************************************************************
      
     make_damping_yes_or_not = 0
      
     if (damping_type .eq. 1) then !Kosloff & Kosloff damping
        do im = 1,nmat        
           if (abs(prop_mat(im,4)) .gt. 10e-5) then 
              make_damping_yes_or_not = 1                
           endif                                        
        enddo                                
     
        if(mpi_id.eq.0 .and. make_damping_yes_or_not .eq. 1) then        
                write(*,'(A)')                         
                write(*,'(A)') '----------------Building the load matrix---------------'
        endif
        
        if(mpi_id.eq.0 .and. make_damping_yes_or_not .eq. 0) then        
                write(*,'(A)')                                                        
                write(*,'(A)')'----------------Building the load matrix---------------'                        
                write(*,'(A)')'ATT: There are no materials with damping defined on'        
        endif
        
     elseif(damping_type .eq. 2) then !Standard Linear Solid damping (Frequency Independent damping)
        if(mpi_id.eq.0) then   
                write(*,'(A)')                     
                write(*,'(A)') '------------Computing anelastic coefficients-----------'
        endif
        !Assume 3 SLS for a range ~ [fmax/10, fmax*10] Hz
        N_SLS = 3
        allocate(Y_lambda(nmat,N_SLS),Y_mu(nmat,N_SLS),frequency_range(N_SLS))
        
        Y_lambda=0.d0;  Y_mu=0.d0;         
                 
        if (fmax .ne. 0.d0) &                
            call MAKE_ANELASTIC_COEFFICIENTS(nmat, N_SLS, prop_mat, QS, QP, fmax, &
                                         Y_lambda, Y_mu, frequency_range, mpi_id)
           
            ! Check coefficient signs
            if (b_failoncoeffs) then
              do im = 1,  nmat
                 do i = 1, N_SLS
                   if (Y_lambda(im,i) .lt. 0.d0 .or. Y_mu(im,i) .lt. 0) then
                      call MPI_FINALIZE(mpi_ierr)
                      if (mpi_id .eq. 0) write(*,'(A)') 'NEGATIVE anelastic coefficient found, program ended.'
                      if (mpi_ierr.ne.0) write(*,'(A,I6)')'MPI Finalization error - proc : ',mpi_id
                      call EXIT(EXIT_ANELASTIC)
                   endif
                 enddo

              enddo

            endif
     elseif(damping_type .eq. 3)  then  !Rayleigh damping
        if(mpi_id.eq.0) then   
                write(*,'(A)')                     
                write(*,'(A)') '------------Computing Rayleigh coefficients-----------'
        endif
        allocate(A0_ray(nmat),A1_ray(nmat))
        A0_ray = 1.d0; A1_ray = 1.d0;
        
              
        if (fmax .ne. 0.d0) &                
            call MAKE_RAYLEIGH_COEFFICIENTS(nmat, A0_ray, A1_ray, prop_mat, QS, fmax,mpi_id)

      elseif(damping_type.eq.4) then !Standard Linear Solid damping (Frequency Independent damping based on Liu & Archuleta 2006)
         if(mpi_id.eq.0) then   
            write(*,'(A)')                     
            write(*,'(A)') '------Computing Visco-Elatic Parameters (Liu & Archuleta 2006)-------'
            if ((ncase .gt. 0) .or. (nmat_nhe .gt. 0)) then
               write(*,'(A)')                     
               write(*,'(A)') 'Not-Honoring Approaches are not implemented for Damping Type 4'
               call exit(EXIT_SETUP)
            endif
            if (fref.le.0.d0) then
               write(*,'(A)')                     
               write(*,'(A)') 'Reference Frequency is not defined in *.mate file'
               call exit(EXIT_SETUP)
            endif
         endif

         N_SLS = 8;
         allocate(Trelax(N_SLS), visc_Mu(nmat), visc_Mp(nmat), &
                  visc_wgt_s(nmat,N_SLS), visc_wgt_p(nmat,N_SLS))
         call MAKE_VISCOELASTIC_LIU_ARCHULETA(nmat, N_SLS, prop_mat, QS, QP, fref, &
                                             Trelax, visc_wgt_s, visc_wgt_p, &
                                             visc_Mu, visc_Mp)

         if (nmat_nlp .gt. 0) then
            call MAKE_NLP_IWANSPRING_YIELDVALUES()
         endif
     endif   
      

!*****************************************************************************************
!   Not Honoring Enhanced (Reading Material properties from Tomo File)          
!*****************************************************************************************
    if (nmat_nhe .gt. 0) then
        call MAKE_NH_Enhanced()
        call MPI_BARRIER(mpi_comm, mpi_ierr)
    endif


!*****************************************************************************************
!   BUILD LOAD MATRIX AND SEISMIC MOMENT TENSOR OR EXPLOSIVE SOURCE
!*****************************************************************************************


      call MAKE_LOAD_MATRIX()

      
!**************************************************************************************            
!   FIND MONITORED NODES  (PGM & MLST input files) + deallocate(elem_domain)
!**************************************************************************************

 
     call FIND_MONITOR_POSITION()
     deallocate(elem_domain);


!*****************************************************************************************
!   CASES 4 NOT HONORING          
!*****************************************************************************************

   
    if(n_case .gt. 0)  then
                                      
       if (mpi_id.eq.0) write(*,'(A)')'----------------Making not-honoring case---------------'                  

       allocate(zs_elev(nnod_loc), zs_all(nnod_loc), sub_tag_all(nnod_loc), vs_tria(nnod_loc), thick(nnod_loc))  
       
       if(mpi_id .eq.0) then
          start1=MPI_WTIME()  
       endif 
       
       call MAKE_NOTHONORING(local_node_num, nnod_loc,&
                             n_case, tag_case, val_case, tol_case, &
                             con_nnz_loc, con_spx_loc, nmat, tag_mat, sdeg_mat, &                        
                             xx_spx_loc, yy_spx_loc, zz_spx_loc, zs_elev, zs_all, vs_tria, thick, &
                             sub_tag_all, mpi_id)
 
  
 
                                                                            
      ncase = 1
      CALL MPI_BARRIER(mpi_comm, mpi_ierr) 
      
      if (mpi_id.eq.0) then
         start2=MPI_WTIME()
         write(*,'(A)')
         write(*,'(A,F20.3,A)')'Set-up time CASE = ',start2-start1,' s'
         write(*,'(A)')'Made.'                  
      endif                        
                                    
    else
    
      ncase = 0

    endif

!*****************************************************************************************
!   RANDOMIZED PARAMETERS          
!*****************************************************************************************
    if (nmat_rnd .gt. 0) then
        allocate(lambda_rnd(nnod_loc),mu_rnd(nnod_loc),rho_rnd(nnod_loc))
        
        call MAKE_RANDOM_PARAM(local_node_num, nnod_loc, &
                               nmat_rnd, rand_mat,  &
                               xx_spx_loc, yy_spx_loc, zz_spx_loc, &
                               lambda_rnd, mu_rnd, rho_rnd, mpi_id)
 
    endif


!*****************************************************************************************
!    SETUP FOR BOUNDARY CONDITIONS    
!*****************************************************************************************

     
     call MAKE_BOUNDARY_CONDITIONS()


!*****************************************************************************************
!    SETUP FOR DG - INTERFACE CONDITIONS          
!*****************************************************************************************

     
     call MAKE_DG_INTERFACE_CONDITIONS()


!***********************************************************************************************      
!    COMPUTE CFL CONDITION
!***********************************************************************************************

      if (mpi_id.eq.0) write(*,'(A)') 
      if (mpi_id.eq.0) write(*,'(A)')'------------Setting calculation parameters-------------'


      
      if(n_case .gt. 0) then 

        
          call GET_CFL4CASES(deltat, nnod_loc, local_node_num, &
                      nmat, tag_mat, prop_mat, sdeg_mat,&
                      xx_spx_loc, yy_spx_loc, zz_spx_loc, &
                      con_nnz_loc, con_spx_loc, &
                      deltat_cfl, fmax, deltat_fixed, mpi_comm, mpi_np, mpi_id, &
                      n_case,val_case, tag_case, zs_elev, zs_all, vs_tria, thick, sub_tag_all, b_failCFL, &
                      damping_type, QS, QP)

      else if (nmat_nhe.gt.0) then
          call GET_CFL4NHE(deltat, nnod_loc, local_node_num, &
                             nmat, tag_mat, prop_mat, sdeg_mat, &
                             xx_spx_loc, yy_spx_loc, zz_spx_loc, &
                             con_nnz_loc, con_spx_loc, & 
                             rho_nhe, lambda_nhe, mu_nhe, &
                             deltat_cfl, fmax, deltat_fixed, &
                             mpi_comm, mpi_np, mpi_id, b_failCFL)
      else
       
          call GET_CFL(deltat, nnod_loc, local_node_num, &
                      nmat, tag_mat, prop_mat, sdeg_mat,&
                      xx_spx_loc, yy_spx_loc, zz_spx_loc, &
                      con_nnz_loc, con_spx_loc, &
                      deltat_cfl, fmax, deltat_fixed, mpi_comm, mpi_np, mpi_id, b_failCFL) 
       endif


      nts = int((tstop -tstart) / deltat)
      tstop = tstart + dfloat(nts)*deltat
      if (mpi_id.eq.0) then
         write(*,'(A,I15)')  'Number of time-steps : ',nts
         write(*,'(A,E14.5)')'Start time           : ',tstart
         write(*,'(A,E14.5)')'Final time           : ',tstop
      endif
      
!**************************************************************************************      
!     SETUP FOR SNAPSHOTS
!**************************************************************************************      
      
!      if (nsnaps.gt.0) then
!        do i = 1, nsnaps
!          itersnap(i) = int((tsnap(i) - tstart) / deltat)
!          if (itersnap(i) .lt. 0) itersnap(i) = 0
!          if (itersnap(i) .gt. nts) itersnap(i) = nts
!        enddo
        
!        initial_snap = 1  
!        if(tstart .gt. 0.d0) then
!           allocate(set_initial_snap(nsnaps))
!           set_initial_snap = tstart - tsnap 
           
!           do i = 2, nsnaps
!                   if(set_initial_snap(i) .lt. set_initial_snap(i-1) .and. set_initial_snap(i) .ge. 0) then
!                      initial_snap = i
!                   endif   
!           enddo        
!           deallocate(set_initial_snap)
!         endif                 
!      endif

!**************************************************************************************      
!    END SETUP 
!**************************************************************************************      

      finish = MPI_WTIME()
      time_in_seconds = finish-start
      
      if (mpi_id.eq.0) then
         write(*,'(A)') 
         write(*,'(A)')'-------------------------------------------------------'
         write(*,'(A,F20.4,A)')'Set-up time = ',time_in_seconds,' s'
         write(*,'(A)')'-------------------------------------------------------'
         write(*,'(A)')
      endif

      if (b_setuponly) then
        call MPI_FINALIZE(mpi_ierr)
        if (mpi_ierr.ne.0) write(*,'(A,I6)')'MPI Finalization error - proc : ',mpi_id
        if (mpi_id .eq. 0) then
            write(*, '(A)') 'Set-up ended. Program finished.'
        endif
        call EXIT(EXIT_SETUP)
      endif

!**************************************************************************************      
!    TIME INTEGRATION  
!**************************************************************************************      

      if (mpi_id.eq.0) then
         write(*,'(A)')
         write(*,'(A)')'------------Beginning of the time-loop-----------------'
         write(*,'(A)')
      endif
           
     
     call TIME_LOOP(el_new)
 
!**************************************************************************************      
!     FINALIZE AND DEALLOCATE MEMORY
!**************************************************************************************      

      call MPI_FINALIZE(mpi_ierr)
      
      if (mpi_ierr.ne.0) write(*,'(A,I6)')'MPI Finalization error - proc : ',mpi_id
                     
      if ((mpi_id.eq.0).and.(opt_out_data.gt.0)) then
         write(*,'(A)');  write(*,'(A)') 'Print output'      
      endif
            
      call DEALLOCATE_VARIABLES()

      if (mpi_id.eq.0) write(*,'(A)')
      if (mpi_id.eq.0) write(*,'(A)')'Bye.'
      
      call EXIT(EXIT_NORMAL)
            
      end program SPEED
      
      
      
      
