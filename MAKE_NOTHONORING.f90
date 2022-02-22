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


!> @brief Makes not-honoring technique. Mechanical properties given node by node.
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] nnloc number of local nodes
!> @param[in] loc_n_num local node numeration
!> @param[in] tag_case label for CASE 
!> @param[in] val_case depth given in CASE keyword
!> @param[in] tolerance  tolerance for identifying nodes
!> @param[in] cs_nnz_loc length of cs_loc
!> @param[in] cs_loc  spectral connectivity vector
!> @param[in] nm number of materials
!> @param[in] tag_mat labels for materials
!> @param[in] sdeg_mat polynomial degree vector
!> @param[in] xs_loc x-coordinate of spectral nodes
!> @param[in] ys_loc y-coordinate of spectral nodes
!> @param[in] zs_loc z-coordinate of spectral nodes
!> @param[in] sub_tag_all  labels for multi not honoring 
!> @param[in] mpi_id  mpi process id
!> @param[out] zs_elev elevation from topography 
!> @param[out] zs_all elevation form alluvial basin
!> @param[out] vs vs30 values assigned to the spectral nodes
!> @param[out] thick thickness of soft sediments

     subroutine MAKE_NOTHONORING(loc_n_num, nn_loc,&
                             n_case, tag_case, val_case, tol_case, &
               		         cs_nnz_loc, cs_loc, nm, tag_mat, sdeg_mat, &			
               		         xs_loc, ys_loc, zs_loc, zs_elev, zs_all, vs, thick, sub_tag_all,mpi_id)

!                                 SUBROUTINE FOR CASEs
!                            COMPUTATION OF THE ELEVATION
!
!                                (FILE -> 'XYZ.out')
!                                       z_elev 
!                     +--------------------x----------------------+  +
!		              |    - vvvvvvvvvvvvvv|vvvvvvvvvvvvvv-       |  | zz_spx_elevation (depth) == zs_elev
!                     |     --vvvvvvvvvvvvv| z_spx vvvvv--        |  |
!                     |        ---vvvvvvvv x vvvvvvvv---          |  +
!                     |           -----vvvv|vvvv-----             |  | zz_spx_alluvial  (all_depth) == zs_all
!                     |                ----x----                  |  +
!                     |                      z_alluvial           |
!                     |                      (FILE -> 'ALL.out')  |
!                     |                                           |
!                     | tag_mat = val_case                        | 
!                     +-------------------------------------------+
!
!                         v = ALLUVIAL BASIN MATERIAL
!                         depth =     z_elev - z_spx     -> zz_spx_elevation == zs_elev
!                         all_depth = z_spx - z_alluvial -> zz_spx_alluvial  == zs_all
!

  
     character*70 :: file_case_xyz						
     character*70 :: file_case_all	
     character*70 :: file_case_vs

     !character*70 :: file_nhe_proc, file_nhe_new
     !integer*4 :: inode, unit_mpi
     !real*8 :: rho,lambda,mu,gamma,qs,qp
  
     integer*4 :: n_case, nn_loc, cs_nnz_loc, nm, mpi_id
     integer*4 :: ncase,vcase,tcase						
     integer*4 :: n_elev,n_tria_elev						
     integer*4 :: start,finish							
     integer*4 :: n_all,n_tria_all, ival, icase
     !integer*4 :: tag_case, val_case					

     integer*4, dimension (:), allocatable :: node1_all,node2_all,node3_all	
     integer*4, dimension (:), allocatable :: node1_elev,node2_elev,node3_elev	 
     integer*4, dimension(3)  :: clock					
     integer*4, dimension(0:cs_nnz_loc) :: cs_loc
     integer*4, dimension(nm) :: tag_mat, sdeg_mat
     integer*4, dimension(nn_loc) :: loc_n_num
     integer*4, dimension(n_case) :: tag_case					
     integer*4, dimension(n_case) :: val_case				
     integer*4, dimension(nn_loc), intent(inout) :: sub_tag_all			

     real*8 :: tolerance							
     real*8 :: max_elev_spacing,max_all_spacing

     real*8, dimension(n_case) :: tol_case				
     real*8, dimension (:), allocatable :: x_elev,y_elev,z_elev, vs_elev, sedim		
     real*8, dimension (:), allocatable :: x_all,y_all,z_all			
     real*8, dimension(nn_loc) :: xs_loc, ys_loc, zs_loc 							
     real*8, dimension(nn_loc), intent(inout) :: zs_elev, zs_all, vs, thick				
		

     !Initialization for vs for the new not-honoring strategy				
     vs = 0.d0
     thick = 0.d0
     							   
      							   
!*************************************************************************************************
!                                   GRENOBLE - honoring
!*************************************************************************************************
    ncase = n_case					
        
        						
    tcase = tag_case(1) !tag_case(ncase)										
    


    if (tcase.eq.1) then										

	 if (mpi_id.eq.0) then									
	     write(*,'(A)')									
   	     write(*,'(A)')'CASE 1: GRENOBLE honoring'					
	     write(*,'(A)')'Reading Topography...'						
	 endif											

	 file_case_xyz ='XYZ.out'								

  	 zs_elev = -1.0e+30								


	 call READ_DIME_FILEXYZ(file_case_xyz,n_elev,n_tria_elev)				

	 allocate(x_elev(n_elev),y_elev(n_elev),z_elev(n_elev))					
	 allocate(node1_elev(n_tria_elev))							
	 allocate(node2_elev(n_tria_elev))							
	 allocate(node3_elev(n_tria_elev))							

	 call READ_FILEXYZ(file_case_xyz,n_elev,n_tria_elev,&					
 			   x_elev,y_elev,z_elev,&	     
			   node1_elev,node2_elev,node3_elev,&
			   max_elev_spacing)		     
														
 	 do icase = 1, ncase
 	   		call GET_NODE_DEPTH_FROM_SIMPLE(loc_n_num, n_elev, n_tria_elev,&					
		                        			x_elev, y_elev, z_elev,&				
					                        node1_elev, node2_elev, node3_elev,&			
               				                cs_nnz_loc, cs_loc, nm, tag_mat, sdeg_mat,&			
               				                nn_loc, xs_loc, ys_loc, zs_loc,&					
					                        zs_elev, val_case(icase), max_elev_spacing, tol_case(icase))		
     enddo
 	 deallocate(x_elev,y_elev,z_elev,node1_elev,node2_elev,node3_elev)							

	 if (mpi_id.eq.0) then									
	  	 write(*,'(A)')'Done'								
		 write(*,'(A)')									
 	 endif														


!*************************************************************************************************
!                          Calssical not honoring (1 ALL and 1 XYZ file)
!*************************************************************************************************


	 elseif (tcase .eq. 2 .or. tcase .eq. 3 .or. tcase .eq. 4 .or. tcase .eq. 6 &
	        .or. tcase .eq. 7 .or. tcase .eq. 8 .or. tcase .eq. 11 .or. tcase .eq. 12 &
	        .or. tcase .eq. 13 .or. tcase .eq. 14 .or. tcase .eq. 15 .or. tcase .eq. 18 &
	        .or. tcase .eq. 22  .or. tcase .eq. 27 .or. tcase .eq. 28 .or. tcase .eq. 40 &
	        .or. tcase .eq. 33 .or. tcase .eq. 38 .or. tcase .eq. 46 .or. tcase .eq. 60) then									
	        
		if (mpi_id.eq. 0 .and. tcase .eq. 2) then									
			write(*,'(A)')									
			write(*,'(A)')'CASE 2: GRENOBLE'					
			
	    elseif(mpi_id .eq. 0 .and. tcase .eq. 3) then		
			write(*,'(A)')									
			write(*,'(A)')'CASE 3: GUBBIO'	      				

	    elseif(mpi_id .eq. 0 .and. tcase .eq. 4) then		
			write(*,'(A)')									
			write(*,'(A)')'CASE 4: SULMONA'	    				
								
	    elseif(mpi_id .eq. 0 .and. tcase .eq. 6) then							
			write(*,'(A)')									
			write(*,'(A)')'CASE 6: FRIULI'	        			
								
	    elseif(mpi_id .eq. 0 .and. tcase .eq. 7) then		
	        write(*,'(A)')									
			write(*,'(A)')'CASE 7: AQUILA'	     				

	    elseif(mpi_id .eq. 0 .and. tcase .eq. 8) then		
			write(*,'(A)')									
			write(*,'(A)')'CASE 8: SANTIAGO'	     				

	    elseif(mpi_id .eq. 0 .and. tcase .eq. 11) then		
			write(*,'(A)')									
			write(*,'(A)')'CASE 11: CHRISTCHURCH'               			

	    elseif(mpi_id .eq. 0 .and. tcase .eq. 12) then		
			write(*,'(A)')									
			write(*,'(A)')'CASE 12: PO PLAIN'               			

	    elseif(mpi_id .eq. 0 .and. tcase .eq. 13) then		
			write(*,'(A)')									
			write(*,'(A)')'CASE 13: PO PLAIN-BEDROCK'               			

	    elseif(mpi_id .eq. 0 .and. tcase .eq. 14) then		
			write(*,'(A)')									
			write(*,'(A)')'CASE 14: WELLINGTON'	     				

	    elseif(mpi_id .eq. 0 .and. tcase .eq. 15) then		
			write(*,'(A)')									
			write(*,'(A)')'CASE 15: MARSICA'	     				

	    elseif(mpi_id .eq. 0 .and. tcase .eq. 18) then		
			write(*,'(A)')									
			write(*,'(A)')'CASE 18: BEIJING-TUTORIAL'	     				

	    elseif(mpi_id .eq. 0 .and. tcase .eq. 22) then		
			write(*,'(A)')									
			write(*,'(A)')'CASE 22: NORCIA'	
	 
	    elseif(mpi_id .eq. 0 .and. tcase .eq. 23) then		
			write(*,'(A)')									
			write(*,'(A)')'CASE 33: GRONINGEN-ZE'	

	    elseif(mpi_id .eq. 0 .and. tcase .eq. 27) then		
			write(*,'(A)')									
			write(*,'(A)')'CASE 27: AQUILA-OB'	

	    elseif(mpi_id .eq. 0 .and. tcase .eq. 28) then		
			write(*,'(A)')									
			write(*,'(A)')'CASE 28: NORCIA-OB'	
			
        elseif(mpi_id .eq. 0 .and. tcase .eq. 38) then		
			write(*,'(A)')									
			write(*,'(A)')'CASE 38: MONTELIMAR'	
			
	    elseif(mpi_id .eq. 0 .and. tcase .eq. 40) then		
			write(*,'(A)')									
			write(*,'(A)')'CASE 40: KUTCH'	     				
		
		elseif(mpi_id .eq. 0 .and. tcase .eq. 46) then		
			write(*,'(A)')									
			write(*,'(A)')'CASE 46: KUMAMOTO'
	
		elseif(mpi_id .eq. 0 .and. tcase .eq. 50) then		
				write(*,'(A)')									
				write(*,'(A)')'CASE 60: JAKARTA'
			endif											


		if(mpi_id .eq. 0) write(*,'(A)')'Reading Topography&Alluvial...'
		file_case_xyz ='XYZ.out'								
		file_case_all ='ALL.out'								

		zs_elev = -1.0e+30								
		zs_all = -1.0e+30								

		call READ_DIME_FILEXYZ(file_case_xyz,n_elev,n_tria_elev)				
		call READ_DIME_FILEXYZ(file_case_all,n_all,n_tria_all)					

		allocate(x_elev(n_elev),y_elev(n_elev),z_elev(n_elev))					
		allocate(node1_elev(n_tria_elev), node2_elev(n_tria_elev), node3_elev(n_tria_elev))	

		allocate(x_all(n_all),y_all(n_all),z_all(n_all))					
		allocate(node1_all(n_tria_all),node2_all(n_tria_all),node3_all(n_tria_all))

		call READ_FILEXYZ(file_case_xyz,n_elev,n_tria_elev,&					
				  x_elev,y_elev,z_elev,&				
				  node1_elev,node2_elev,node3_elev,&			
				  max_elev_spacing)
				  					
		call READ_FILEXYZ(file_case_all,n_all,n_tria_all,&					
				  x_all,y_all,z_all,&					
				  node1_all,node2_all,node3_all,&			
				  max_all_spacing)					

		do icase = 1, ncase
		 	 
		 	 call GET_NODE_DEPTH_FROM_ALLUVIAL(loc_n_num, n_all, n_tria_all, &					
						                  x_all, y_all, z_all, &					
						                  node1_all, node2_all, node3_all,&			
			                              cs_nnz_loc, cs_loc, nm, tag_mat, sdeg_mat, &	
                        			      nn_loc, xs_loc, ys_loc, zs_loc, &					
						                  zs_all, val_case(icase), max_all_spacing, tol_case(icase))		

		    call GET_NODE_DEPTH_FROM_CMPLX(loc_n_num, n_elev, n_tria_elev, &					
	   				                   x_elev, y_elev, z_elev, &				
						               node1_elev, node2_elev, node3_elev,&			
                                       cs_nnz_loc, cs_loc, nm, tag_mat, sdeg_mat, &
                                       nn_loc, xs_loc, ys_loc, zs_loc, &
				                       zs_elev, zs_all, &					
						               val_case(icase), max_elev_spacing, tol_case(icase))			

        enddo

		deallocate(x_elev, y_elev, z_elev, node1_elev, node2_elev, node3_elev)	
		deallocate(x_all, y_all, z_all, node1_all, node2_all, node3_all)

		if (mpi_id.eq.0) then									
			write(*,'(A)')'Done'								
			write(*,'(A)')									
		endif                                                                                   

!*************************************************************************************************
!                            L'AQUILA- MULTI BASIN
!*************************************************************************************************

	 elseif (tcase.eq. 70) then									
		if (mpi_id.eq.0 ) then									
		 	write(*,'(A)')									
			write(*,'(A)')'CASE 70: Aquila-multibasin'	
				
		endif											
								
		write(*,'(A)')'Reading Topography&Alluvial...'		
		 
		file_case_xyz ='XYZ.out'								

		zs_elev = -1.0e+30
		zs_all = 0.d0;								

		call READ_DIME_FILEXYZ(file_case_xyz,n_elev,n_tria_elev)
		allocate(x_elev(n_elev),y_elev(n_elev),z_elev(n_elev))					
		allocate(node1_elev(n_tria_elev), node2_elev(n_tria_elev), node3_elev(n_tria_elev))	

		call READ_FILEXYZ(file_case_xyz,n_elev,n_tria_elev,&					
				  x_elev,y_elev,z_elev,&				
				  node1_elev,node2_elev,node3_elev,&			
				  max_elev_spacing)
				  					
		call GET_NODE_DEPTH_FROM_CMPLX(loc_n_num, n_elev, n_tria_elev, &					
			   x_elev, y_elev, z_elev, &				
			   node1_elev, node2_elev, node3_elev,&			
			   cs_nnz_loc, cs_loc, nm, tag_mat, sdeg_mat, &
			   nn_loc, xs_loc, ys_loc, zs_loc, &
			   zs_elev, zs_all, &					
			   val_case(1), max_elev_spacing, tol_case(1))		

   	    deallocate(x_elev, y_elev, z_elev, node1_elev, node2_elev, node3_elev)
				  
		sub_tag_all = 3	
		ival = 3							
											
		do j = 1,2										    
			if (j.eq.1) then								
				file_case_all ='ALL1.out'
			else										
				file_case_all ='ALL2.out'						
			endif										
	
			zs_all = -1.0e+30

			call READ_DIME_FILEXYZ(file_case_all,n_all,n_tria_all)


			allocate(x_all(n_all), y_all(n_all), z_all(n_all))
			allocate(node1_all(n_tria_all), node2_all(n_tria_all), node3_all(n_tria_all))
			
   		    call READ_FILEXYZ(file_case_all,n_all,n_tria_all,&					
			          		  x_all,y_all,z_all,&					
					          node1_all,node2_all,node3_all,&			
							  max_all_spacing)				
							  										  

		    do icase = 1, ncase
			 
			 		call GET_NODE_DEPTH_FROM_ALLUVIAL(loc_n_num, n_all, n_tria_all, &					
							   x_all, y_all, z_all, &					
							   node1_all, node2_all, node3_all,&			
				               cs_nnz_loc, cs_loc, nm, tag_mat, sdeg_mat, &	
	                           nn_loc, xs_loc, ys_loc, zs_loc, &	
							   zs_all, val_case(icase), max_all_spacing, tol_case(icase))		
		    enddo    					

		    call MAKE_SUBTAG_ALLUVIAL(nn_loc, zs_all, j, sub_tag_all, xs_loc, ival)
		        	
			deallocate(x_all, y_all, z_all, node1_all, node2_all, node3_all)
			
			if (mpi_id.eq.0) then	
				write(*,'(A)')	
				write(*,'(A,I8)') 'ALLUVIAL Layer # ',j	
			endif

		 enddo      
                                   
		if (mpi_id.eq.0) then
			write(*,'(A)') 'Done'
			write(*,'(A)')	
		endif

!*************************************************************************************************
!                             MULTI - NOT honoring
!*************************************************************************************************

	elseif (tcase.eq. 5 .or. tcase .eq. 50) then									
		if (mpi_id.eq.0 .and. tcase .eq. 5) then									
			write(*,'(A)')									
			write(*,'(A)')'CASE 5: VOLVI for CASHIMA benchmark'		

		elseif (mpi_id.eq.0 .and. tcase .eq. 50) then									
			write(*,'(A)')									
			write(*,'(A)')'CASE 50: PLANE-WAVE benchmark'	
				
		endif											
								
     	write(*,'(A)')'Reading Topography&Alluvial...'					
		sub_tag_all = 4	
		ival = 4							
											
		do j = 1,3										    
			if (j.eq.1) then								
				file_case_all ='ALL1.out'
			elseif (j.eq.2) then								
				file_case_all ='ALL2.out'						
			else										
				file_case_all ='ALL3.out'						
			endif										
	
			zs_all = -1.0e+30

			call READ_DIME_FILEXYZ(file_case_all,n_all,n_tria_all)

			allocate(x_all(n_all), y_all(n_all), z_all(n_all))
			allocate(node1_all(n_tria_all), node2_all(n_tria_all), node3_all(n_tria_all))
			
   		    call READ_FILEXYZ(file_case_all,n_all,n_tria_all,&					
			          		  x_all,y_all,z_all,&					
					          node1_all,node2_all,node3_all,&			
					          max_all_spacing)					

		    do icase = 1, ncase
			 
			 		call GET_NODE_DEPTH_FROM_ALLUVIAL(loc_n_num, n_all, n_tria_all, &					
							   x_all, y_all, z_all, &					
							   node1_all, node2_all, node3_all,&			
				               cs_nnz_loc, cs_loc, nm, tag_mat, sdeg_mat, &	
	                           nn_loc, xs_loc, ys_loc, zs_loc, &	
							   zs_all, val_case(icase), max_all_spacing, tol_case(icase))		
		    enddo    					

		    call MAKE_SUBTAG_ALLUVIAL(nn_loc, zs_all, j, sub_tag_all, xs_loc, ival)
		        	
			deallocate(x_all, y_all, z_all, node1_all, node2_all, node3_all)
			
			if (mpi_id.eq.0) then	
				write(*,'(A)')	
				write(*,'(A,I8)') 'ALLUVIAL Layer # ',j	
			endif

		enddo      
                                   
		if (mpi_id.eq.0) then
			write(*,'(A)') 'Done'
			write(*,'(A)')	
		endif

				
!*************************************************************************************************
!                             XYZ map - ALL map - VS30 map
!*************************************************************************************************

	elseif (tcase.eq. 16 .or. tcase.eq. 19 .or. tcase .eq. 20 .or. tcase .eq. 21 &
	        .or. tcase .eq. 29 .or. tcase .eq. 35) then	
									
		if (mpi_id.eq.0 .and. tcase .eq. 16) then	        
			write(*,'(A)')									
			write(*,'(A)')'CASE 16: ISTANBUL'	     				
		endif
		if (mpi_id.eq.0 .and. tcase .eq. 19) then	        
			write(*,'(A)')									
			write(*,'(A)')'CASE 19: THESSALONIKI'	     				
		endif
	        if (mpi_id.eq.0 .and. tcase .eq. 20) then        
			write(*,'(A)')									
			write(*,'(A)')'CASE 20: ATHENS'	     				
		endif													
	        if (mpi_id.eq.0 .and. tcase .eq. 21) then        
			write(*,'(A)')									
			write(*,'(A)')'CASE 21: BEIJING '	   
		endif													
	        if (mpi_id.eq.0 .and. tcase .eq. 29) then        
			write(*,'(A)')									
			write(*,'(A)')'CASE 29: THESS-BEDROCK'  				
		endif		
		 	if (mpi_id.eq.0 .and. tcase .eq. 35) then        
			write(*,'(A)')									
			write(*,'(A)')'CASE 35: THESS+MYGD-FINAL'  	
        endif 											
        	
        if (mpi_id.eq.0) write(*,'(A)')'Reading Topography&Alluvial...'					


		file_case_xyz ='XYZ.out'								
		if(tcase .eq. 19 .or. tcase .eq. 21 .or. tcase .eq. 29)  file_case_all ='ALL.out'
		file_case_vs = 'VS_RS.out'								
														
		zs_elev = 0.d0	
		zs_all = 1.d0							
		if(tcase .eq. 19 .or. tcase .eq. 21 .or. tcase .eq. 29) zs_all = -1.0e+30								

		call READ_DIME_FILEXYZ(file_case_xyz,n_elev,n_tria_elev)	
					
		if(tcase .eq. 19 .or. tcase .eq. 21 .or. tcase .eq. 29) call READ_DIME_FILEXYZ(file_case_all,n_all,n_tria_all)					

		allocate(x_elev(n_elev),y_elev(n_elev),z_elev(n_elev),&
		         vs_elev(n_tria_elev),sedim(n_tria_elev))					
		allocate(node1_elev(n_tria_elev), node2_elev(n_tria_elev), node3_elev(n_tria_elev))

		if(tcase .eq. 19 .or. tcase .eq. 21 .or. tcase .eq. 29) allocate(x_all(n_all),y_all(n_all),z_all(n_all))					
		if(tcase .eq. 19 .or. tcase .eq. 21 .or. tcase .eq. 29) allocate(node1_all(n_tria_all),node2_all(n_tria_all),node3_all(n_tria_all))

		call READ_FILEXYZ(file_case_xyz,n_elev,n_tria_elev,&					
				  x_elev,y_elev,z_elev,&				
				  node1_elev,node2_elev,node3_elev,&			
				  max_elev_spacing)
				  					
		if(tcase .eq. 19 .or. tcase .eq. 21 .or. tcase .eq. 29) &
		   call READ_FILEXYZ(file_case_all,n_all,n_tria_all,&					
			            	 x_all,y_all,z_all,&					
				             node1_all,node2_all,node3_all,&			
				             max_all_spacing)					

        
        call READ_FILEVS(file_case_vs, n_tria_elev, vs_elev, sedim)


        do icase = 1, ncase 
		     if(tcase .eq. 19 .or. tcase .eq. 21 .or. tcase .eq. 29) call GET_NODE_DEPTH_FROM_ALLUVIAL(loc_n_num, n_all, n_tria_all, &					
			     	                    		   x_all, y_all, z_all, &					
				    		                       node1_all, node2_all, node3_all,&			
			                                       cs_nnz_loc, cs_loc, nm, tag_mat, sdeg_mat, &	
                            		               nn_loc, xs_loc, ys_loc, zs_loc, &	
						                           zs_all, val_case(icase), max_all_spacing, tol_case(icase))		

		     call GET_NODE_DEPTH_AND_VS(loc_n_num, n_elev, n_tria_elev, &					
	   				        x_elev, y_elev, z_elev, vs_elev, sedim,&				
						    node1_elev, node2_elev, node3_elev,&			
                            cs_nnz_loc, cs_loc, nm, tag_mat, sdeg_mat, &
                            nn_loc, xs_loc, ys_loc, zs_loc, &	
				            zs_elev, zs_all, vs, thick, &					
						    val_case(icase), max_elev_spacing, tol_case(icase))			
         enddo


		deallocate(x_elev, y_elev, z_elev,vs_elev,sedim, node1_elev, node2_elev, node3_elev)
		if(tcase .eq. 19 .or. tcase .eq. 21 .or. tcase .eq. 29)  deallocate(x_all, y_all, z_all, node1_all, node2_all, node3_all)

		if (mpi_id.eq.0) then									
			write(*,'(A)')'Done'								
			write(*,'(A)')									
		endif 
				
 
				

!*************************************************************************************************
!                             Atene - Partenone
!*************************************************************************************************

	elseif (tcase.eq. 30) then									
		if (mpi_id.eq.0) then									
			write(*,'(A)')									
			write(*,'(A)')'CASE 30: ATHENS-PARTHENON'		
			write(*,'(A)')'Reading Topography&Alluvial...'					
		endif											
								

		sub_tag_all = 4		
		ival = 4						
											
		do j = 1,3										    
			if (j.eq.1) then								
				file_case_all ='ALL1.out'
			elseif (j.eq.2) then								
				file_case_all ='ALL2.out'
			else								
				file_case_all ='ALL3.out'						
			endif										
	
			zs_all = -1.0e+30

			call READ_DIME_FILEXYZ(file_case_all,n_all,n_tria_all)

			allocate(x_all(n_all), y_all(n_all), z_all(n_all))
			allocate(node1_all(n_tria_all), node2_all(n_tria_all), node3_all(n_tria_all))
			
   		    call READ_FILEXYZ(file_case_all,n_all,n_tria_all,&					
					  x_all,y_all,z_all,&					
					  node1_all,node2_all,node3_all,&			
					  max_all_spacing)					

			do icase = 1, ncase
						call GET_NODE_DEPTH_FROM_ALLUVIAL(loc_n_num, n_all, n_tria_all, &					
							                  x_all, y_all, z_all, &					
							                  node1_all, node2_all, node3_all,&			
				                              cs_nnz_loc, cs_loc, nm, tag_mat, sdeg_mat, &	
	                        			      nn_loc, xs_loc, ys_loc, zs_loc, &	
							                  zs_all, val_case(icase), max_all_spacing, tol_case(icase))	
            enddo

			call MAKE_SUBTAG_ALLUVIAL(nn_loc, zs_all, j, sub_tag_all, xs_loc, ival)				

			deallocate(x_all, y_all, z_all, node1_all, node2_all, node3_all)
			
			if (mpi_id.eq.0) then	
				write(*,'(A)')	
				write(*,'(A,I8)') 'ALLUVIAL Layer # ',j	
			endif

		enddo 	

   
                   
                                   
		if (mpi_id.eq.0) then
			write(*,'(A)') 'Done'
			write(*,'(A)')	
		endif
								
									
!*************************************************************************************************
!                             Groningen 
!*************************************************************************************************

	elseif (tcase.eq. 31 .or. tcase .eq. 32) then									
		if (mpi_id.eq.0) then									
			write(*,'(A)')
			if(tcase.eq. 31) write(*,'(A)')'CASE 31: GRONINGEN'		
			if(tcase.eq. 32) write(*,'(A)')'CASE 32: GRONINGEN'		 
			write(*,'(A)')'Reading Topography&Alluvial...'					
		endif											
								

		sub_tag_all = 7		
		ival = 7						
											
		do j = 1,6										    
			if (j.eq.1) then								
				file_case_all ='ALL1.out'
			elseif(j.eq.2) then								
				file_case_all ='ALL2.out'						
			elseif(j.eq.3) then								
				file_case_all ='ALL3.out'						
			elseif(j.eq.4) then								
				file_case_all ='ALL4.out'						
			elseif(j.eq.5) then								
				file_case_all ='ALL5.out'						
			else								
				file_case_all ='ALL6.out'						
			endif										
	
			zs_all = -1.0e+30

			call READ_DIME_FILEXYZ(file_case_all,n_all,n_tria_all)

			allocate(x_all(n_all), y_all(n_all), z_all(n_all))
			allocate(node1_all(n_tria_all), node2_all(n_tria_all), node3_all(n_tria_all))
			
   		    call READ_FILEXYZ(file_case_all,n_all,n_tria_all,&					
					  x_all,y_all,z_all,&					
					  node1_all,node2_all,node3_all,&			
					  max_all_spacing)					

			do icase = 1, ncase
					call GET_NODE_DEPTH_FROM_ALLUVIAL(loc_n_num, n_all, n_tria_all, &					
							                  x_all, y_all, z_all, &					
							                  node1_all, node2_all, node3_all,&			
				                              cs_nnz_loc, cs_loc, nm, tag_mat, sdeg_mat, &	
	                        			      nn_loc, xs_loc, ys_loc, zs_loc, &	
							                  zs_all, val_case(icase), max_all_spacing, tol_case(icase))		
            enddo

			call MAKE_SUBTAG_ALLUVIAL(nn_loc, zs_all, j, sub_tag_all, xs_loc, ival)				

			deallocate(x_all, y_all, z_all, node1_all, node2_all, node3_all)
			
			if (mpi_id.eq.0) then	
				write(*,'(A)')	
				write(*,'(A,I8)') 'ALLUVIAL Layer # ',j	
			endif

		enddo !do j = 1,3 	

                !do i = 1, nn_loc
                !   write(*,*) zs_loc(i), sub_tag_all(i)
                !enddo
                !read(*,*)   
                   
                                   
		if (mpi_id.eq.0) then
			write(*,'(A)') 'Done'
			write(*,'(A)')	
		endif

!*************************************************************************************************
!                             Groningen layered model
!*************************************************************************************************

	!elseif (tcase.eq. 32) then		
	
    ! 	zs_elev = -1.0e+30								
    !   zs_all = -1.0e+30								
							
!		if (mpi_id.eq.0) then									
!			write(*,'(A)')									
!			write(*,'(A)')'CASE 32: GRONINGEN'		
!			write(*,'(A)')'Layered Model...'					
!		endif												
!		if (mpi_id.eq.0) then
!			write(*,'(A)') 'Done'
!			write(*,'(A)')	
!		endif
									

!*************************************************************************************************
!                             TEST honoring -- only topgraphy surface
!*************************************************************************************************

	elseif (tcase.eq.98) then									
                if (mpi_id.eq.0) then									
			write(*,'(A)')									
			write(*,'(A)')'CASE 98: TEST honoring'	      					
			write(*,'(A)')'Reading Topography...'						
		endif											

		 file_case_xyz ='XYZ.out'								
	  	 zs_elev = -1.0e+30								

	
		 call READ_DIME_FILEXYZ(file_case_xyz,n_elev,n_tria_elev)				

		 allocate(x_elev(n_elev),y_elev(n_elev),z_elev(n_elev))					
		 allocate(node1_elev(n_tria_elev))							
		 allocate(node2_elev(n_tria_elev))							
		 allocate(node3_elev(n_tria_elev))							

		 call READ_FILEXYZ(file_case_xyz,n_elev,n_tria_elev,&					
	 			   x_elev,y_elev,z_elev,&	      			
				   node1_elev,node2_elev,node3_elev,& 		
				   max_elev_spacing)		      			
														
	 	 do icase = 1, ncase
	 	 		
	 	 			call GET_NODE_DEPTH_FROM_SIMPLE(loc_n_num, n_elev, n_tria_elev,&					
						   x_elev, y_elev, z_elev,&				
						   node1_elev, node2_elev, node3_elev,&			
	               	       cs_nnz_loc, cs_loc, nm, tag_mat, sdeg_mat,&			
	               		   nn_loc, xs_loc, ys_loc, zs_loc,&					
						   zs_elev, val_case(iacse), max_elev_spacing, tol_case(icase))		

         enddo
         
	 	 deallocate(x_elev,y_elev,z_elev,node1_elev,node2_elev,node3_elev)	

		if (mpi_id.eq.0) then									
			write(*,'(A)')'Done'								
			write(*,'(A)')									
		endif													

!*************************************************************************************************
!                             TEST - NOT honoring - topo and alluvial surfaces
!*************************************************************************************************

	elseif (tcase.eq.99) then
	       if (mpi_id.eq.0) then
	       		write(*,'(A)')											
			write(*,'(A)')'CASE 99: TEST not honoring'
			write(*,'(A)')'Reading Topography&Alluvial...'
		endif		

		file_case_xyz ='XYZ.out'								
		file_case_all ='ALL.out'								
									
		zs_elev = -1.0e+30								
		zs_all = -1.0e+30								

		call READ_DIME_FILEXYZ(file_case_xyz,n_elev,n_tria_elev)				
		call READ_DIME_FILEXYZ(file_case_all,n_all,n_tria_all)					

		allocate(x_elev(n_elev),y_elev(n_elev),z_elev(n_elev))					
		allocate(node1_elev(n_tria_elev), node2_elev(n_tria_elev), node3_elev(n_tria_elev))

		allocate(x_all(n_all),y_all(n_all),z_all(n_all))					
		allocate(node1_all(n_tria_all),node2_all(n_tria_all),node3_all(n_tria_all))

		call READ_FILEXYZ(file_case_xyz,n_elev,n_tria_elev,&					
				  x_elev,y_elev,z_elev,&				
				  node1_elev,node2_elev,node3_elev,&			
				  max_elev_spacing)
				  					
		call READ_FILEXYZ(file_case_all,n_all,n_tria_all,&					
				  x_all,y_all,z_all,&					
				  node1_all,node2_all,node3_all,&			
				  max_all_spacing)					


	 	do icase = 1, ncase
	 	 
		    call GET_NODE_DEPTH_FROM_ALLUVIAL(loc_n_num, n_all, n_tria_all, &					
						   x_all, y_all, z_all, &					
						   node1_all, node2_all, node3_all,&			
			               cs_nnz_loc, cs_loc, nm, tag_mat, sdeg_mat, &	
                           nn_loc, xs_loc, ys_loc, zs_loc, &	
						   zs_all, val_case(icase), max_all_spacing, tol_case(icase))		

		    call GET_NODE_DEPTH_FROM_CMPLX(loc_n_num, n_elev, n_tria_elev, &					
	   				       x_elev, y_elev, z_elev, &				
						   node1_elev, node2_elev, node3_elev,&			
                           cs_nnz_loc, cs_loc, nm, tag_mat, sdeg_mat, &
                           nn_loc, xs_loc, ys_loc, zs_loc, &
				           zs_elev, zs_all, &					
						   val_case(icase), max_elev_spacing, tol_case(icase))			
       enddo


		deallocate(x_elev, y_elev, z_elev, node1_elev, node2_elev, node3_elev)	
		deallocate(x_all, y_all, z_all, node1_all, node2_all, node3_all)

		if (mpi_id.eq.0) then									
			write(*,'(A)')'Done'								
			write(*,'(A)')									
		endif                                                                                   


	endif ! TCASE	

     
     
     end subroutine MAKE_NOTHONORING
     
