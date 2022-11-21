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

!> @brief Reads dimensions in filemate (*.mate)
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] filemate  file name (*.mate)
!> @param[out] nb_mat  number of blocks (materials)
!> @param[out] nb_mat_nle  number of non-linear elastic blocks
!> @param[out] nb_mat_nlp  number of non-linear - Plasticity blocks
!> @param[out] nb_load_dirX   number of Dirichlet b.c. (x-dir)
!> @param[out] nb_load_ndirY  number of Dirichlet b.c. (y-dir)
!> @param[out] nb_load_ndirZ  number of Dirichlet b.c. (z-dir)
!> @param[out] nb_load_nneuX  number of Neumann boundary loads (x-dir)
!> @param[out] nb_load_nneuY  number of Neumann boundary loads (y-dir)
!> @param[out] nb_load_nneuZ  number of Neumann boundary loads (z-dir)
!> @param[out] nb_load_nneuN  number of Neumann boundary loads (normal direction)
!> @param[out] nb_load_poiX   number of point loads (x-dir)
!> @param[out] nb_load_poiY   number of point loads (y-dir)
!> @param[out] nb_load_poiZ   number of point loads (z-dir)
!> @param[out] nb_load_traX   number of traveling loads (x-dir)
!> @param[out] nb_load_traY   number of traveling loads (y-dir)
!> @param[out] nb_load_traZ   number of traveling loads (z-dir)
!> @param[out] nb_load_plaX   number of plane wave loads (x-dir)
!> @param[out] nb_load_plaY   number of plane wave loads (y-dir)
!> @param[out] nb_load_plaZ   number of plane wave loads (z-dir)
!> @param[out] nb_load_forX   number of volume forces (x-dir)
!> @param[out] nb_load_forY   number of volume forces (y-dir)
!> @param[out] nb_load_forZ   number of volume forces (z-dir)  
!> @param[out] nb_load_forc   
!> @param[out] nb_load_pres   number of surface forces (pressure)
!> @param[out] nb_load_shea   number of surface forces (shear) 
!> @param[out] nb_load_abc   number of abc boundary conditions
!> @param[out] nb_load_dg   number of DG interfaces
!> @param[out] nb_func  number of time functions
!> @param[out] nb_func_data  number of data for time functions
!> @param[out] nb_load_sism  number of seismic loads
!> @param[out] nb_load_expl  number of explosive sources
!> @param[out] nb_case  number of not honoring cases
!> @param[out] n_test  1 for test case mode

      subroutine READ_DIME_FILEMATE(filemate,nb_mat,&
                                    nb_mat_nle, nb_mat_rnd, nb_mat_nlp, &
                                    nb_load_dirX,nb_load_dirY,nb_load_dirZ, &
                                    nb_load_neuX,nb_load_neuY,nb_load_neuZ, &
                                    nb_load_neuN, &  
                                    nb_load_poiX,nb_load_poiY,nb_load_poiZ, &
                                    nb_load_traX,nb_load_traY,nb_load_traZ, &
                                    nb_load_plaX,nb_load_plaY,nb_load_plaZ, &        
                                    nb_load_forX,nb_load_forY,nb_load_forZ, &
                                    nb_load_forc,nb_load_pres,nb_load_shea, &
                                    nb_load_abc,nb_load_dg,nb_func,nb_func_data,&
                                    nb_load_sism,nb_load_expl,nb_case,nb_nhee,&
                                    n_test,nb_frac,srcmodflag)  

                                                                                                                         
                                                                                            

  implicit none

  character*70   :: filemate, srcname, dum_filename   
  character*100000 :: inline
  character*4 :: keyword

  integer*4 :: nb_mat
  integer*4 :: nb_mat_nle, nb_mat_rnd, nb_mat_nlp
  integer*4 :: nb_load_dirX,nb_load_dirY,nb_load_dirZ
  integer*4 :: nb_load_neuX,nb_load_neuY,nb_load_neuZ
  integer*4 :: nb_load_neuN
  integer*4 :: nb_load_traX,nb_load_traY,nb_load_traZ
  integer*4 :: nb_load_poiX,nb_load_poiY,nb_load_poiZ
  integer*4 :: nb_load_plaX,nb_load_plaY,nb_load_plaZ
  integer*4 :: nb_load_forX,nb_load_forY,nb_load_forZ
  integer*4 :: nb_load_forc,nb_load_pres,nb_load_shea
  integer*4 :: nb_load_abc, nb_load_dg
  integer*4 :: nb_load_sism
  integer*4 :: nb_load_expl
  integer*4 :: nb_case, n_test, nb_frac,nb_nhee,srcmodflag 
  integer*4 :: nb_func,nb_func_data
  integer*4 :: status
  integer*4 :: lab_fnc, type_fnc, ndat_fnc

  real*8 :: dum1

  nb_mat = 0;            nb_mat_nle = 0;        nb_case = 0;           nb_mat_rnd = 0;     
  nb_mat_nlp = 0;   
  nb_load_dirX = 0;      nb_load_dirY = 0;      nb_load_dirZ = 0
  nb_load_neuX = 0;      nb_load_neuY = 0;      nb_load_neuZ = 0
  nb_load_neuN = 0;      nb_load_poiX = 0;      nb_load_poiY = 0 
  nb_load_poiZ = 0;      nb_load_plaX = 0;      nb_load_plaY = 0 
  nb_load_plaZ = 0;      nb_load_forX = 0;      nb_load_forY = 0   
  nb_load_forZ = 0;      nb_load_forc = 0;      nb_load_pres = 0 
  nb_load_shea = 0;      nb_load_abc = 0;       nb_load_dg = 0
  nb_load_sism = 0;      nb_load_expl = 0;       nb_frac = 0;
  nb_load_traX = 0;      nb_load_traY = 0;      nb_load_traZ = 0

  nb_func = 0;   nb_func_data = 0
  n_test = 0;

  open(40,file=filemate)

  do
    read(40,'(A)',IOSTAT = status) inline

    if (status.ne.0) exit

      keyword = inline(1:4)

      select case (keyword)

           case('MATE')
           nb_mat = nb_mat + 1
           case('MATN')    ! Nonlinear - Elastic
           nb_mat_nle = nb_mat_nle + 1     
           case('MATR')              
           nb_mat_rnd = nb_mat_rnd + 1                                    
          case('MATP')     ! Nonlinear - Plasticity Model using IWAN springs (Fabian's)    
            nb_mat_nlp = nb_mat_nlp + 1      
           case('DIRX')        
           nb_load_dirX = nb_load_dirX + 1
           case('DIRY')
           nb_load_dirY = nb_load_dirY + 1
           case('DIRZ')
           nb_load_dirZ = nb_load_dirZ + 1
           case('NEUX')
           nb_load_neuX = nb_load_neuX + 1
           case('NEUY')
           nb_load_neuY = nb_load_neuY + 1
           case('NEUZ')
           nb_load_neuZ = nb_load_neuZ + 1
           case('NEUN')         
           nb_load_neuN = nb_load_neuN + 1                 
           case('PLOX')
           nb_load_poiX = nb_load_poiX + 1
           case('PLOY')
           nb_load_poiY = nb_load_poiY + 1
           case('PLOZ')
           nb_load_poiZ = nb_load_poiZ + 1

!           case('TLOX')
!            nb_load_traX = nb_load_traX + 1    
!           case('TLOY')
!            nb_load_traY = nb_load_traY + 1    

           case('TLOZ')
            nb_load_traZ = nb_load_traZ + 1    
           case('PLAX')        
            nb_load_plaX = nb_load_plaX + 1                
           case('PLAY')        
            nb_load_plaY = nb_load_plaY + 1                
           case('PLAZ')        
            nb_load_plaZ = nb_load_plaZ + 1                
           case('FORX')
            nb_load_forX = nb_load_forX + 1
           case('FORY')
            nb_load_forY = nb_load_forY + 1
           case('FORZ')
            nb_load_forZ = nb_load_forZ + 1
           case('FORC')
            nb_load_forc = nb_load_forc + 1
           case('PRES')
            nb_load_pres = nb_load_pres + 1
           case('SHEA')
            nb_load_shea = nb_load_shea + 1
           case('ABSO')
            nb_load_abc = nb_load_abc + 1
           case('DGIC')
            nb_load_dg = nb_load_dg + 1
           case('DGIX')
            nb_load_dg = nb_load_dg + 1
           case('SISM')        
            nb_load_sism = nb_load_sism + 1                
           case('EXPL')        
            nb_load_expl = nb_load_expl + 1                
           case('CASE')        
            nb_case = nb_case + 1
           case('NHEE')        
            ! Not-honoring Enhanced - Elastic- damping type 2
            nb_nhee = nb_nhee + 1
           case('SLIP')        
            ! Not-honoring Fault Plane
            read(inline(5:),*) srcname
            if (srcname.eq.'LOAD-SRCMOD2') then
              srcmodflag = 1;
            endif
           case('TEST')
            n_test = n_test + 1         
           case('FRAC')
            nb_frac = 1                                     
           case('FUNC')
            nb_func = nb_func + 1         
            read(inline(5:),*) lab_fnc, type_fnc
            
            select case (type_fnc)
               case(0,32)
                ! Case 32 - Ramp Source Time Function
                 nb_func_data = nb_func_data + 0
               case(1) 
                 ! RICKER WAVELET
                 nb_func_data = nb_func_data + 2
               case(2) 
                 nb_func_data = nb_func_data + 2
               case(3,30,31,33) 
                   ! TIME SERIES
                   ! Case 31 - Text File with Source Time Function
                 read(inline(5:),*) lab_fnc, type_fnc, ndat_fnc
                 nb_func_data = nb_func_data + 2*ndat_fnc
                
               case(4) 
                 ! DERIVATIVE OF THE RICKER WAVELET
                 nb_func_data = nb_func_data + 2
               case(5) 
                 ! DERIVATIVE OF THE GAUSSIAN WAVELET
                 nb_func_data = nb_func_data + 2

               case(6) 
                 nb_func_data = nb_func_data + 2
               case(7) 
                 nb_func_data = nb_func_data + 2
               case(8,9) 
                 nb_func_data = nb_func_data + 1

               case(12) 
                 ! SIGMOIDAL FUNC
                 nb_func_data = nb_func_data + 3
               case(13) 
                 ! GRENOBLE BENCHMARK
                 nb_func_data = nb_func_data + 2
               case(14) 
                 ! SCEC BENCHMARK
                 nb_func_data = nb_func_data + 2
               case(15) 
                 ! EXPLOSION    
                 nb_func_data = nb_func_data + 4
               case(50,55) 
                 ! VARIABLE TAU 
                 nb_func_data = nb_func_data + 2  
               case(60,62) 
                 ! LINEAR EQUIVALENT
                 read(inline(5:),*) lab_fnc, type_fnc, ndat_fnc               
                 nb_func_data = nb_func_data + 2*ndat_fnc                            
               case(61,63)                                      
                 read(inline(5:),*) lab_fnc, type_fnc, ndat_fnc               
                 nb_func_data = nb_func_data + 2*ndat_fnc               
               case(70) 
                  ! Hyperbolic G/Gmax Curve - based on Reference Strain - Hardin and Drnevich (1972)
                  ! FunctionLabel;  FunctionType(i.e 70 or 71); No.of Iwan Springs; ref_strain
                  !read(inline(5:),*) lab_fnc, type_fnc, ndat_fnc, dum1     
                  nb_func_data = nb_func_data + 2            
               !case(71) 
                  ! User Defined G/Gmax Curve (External Ascii file with data in two columns)
                  ! First Column - Strain Values      Second Column - Corresponding G/Gmax value
                  !read(inline(5:),*) lab_fnc, type_fnc, ndat_fnc, dum_filename
                  !nb_func_data = nb_func_data + 2*ndat_fnc             
               case(99) 
                 ! CASHIMA   
                 nb_func_data = nb_func_data + 2    
               case(100)
                 nb_func_data = nb_func_data + 1   
               case(773) 
                   ! TIME SERIES
                 read(inline(5:),*) lab_fnc, type_fnc, ndat_fnc
                 nb_func_data = nb_func_data + ndat_fnc

            end select 
         
         end select

      enddo
      
      close(40)
      
      return
      end subroutine READ_DIME_FILEMATE



