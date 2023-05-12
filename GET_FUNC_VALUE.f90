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

!> @brief Computes time evolution function.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] nb_fnc number of functions
!> @param[in] type_fnc function type
!> @param[in] ind_fnc indices for the data 
!> @param[in] nb_data_fnc number of data for each function
!> @param[in] data_fnc data for the calculation (depending on type_fnc) 
!> @param[in] id_fnc  number of the function
!> @param[in] time  instant time
!> @param[in] dist  distance form source point (for travelling load)
!> @param[in] vel  (constant) velocity of the travelling load
!> @param[out] GET_FUNC_VALUE value of the time function

           

      real*8 function GET_FUNC_VALUE(nb_fnc, type_fnc, ind_fnc, &
                                     data_fnc, nb_data_fnc, id_fnc, time, dist,vel)



      use binarysearch
      
      implicit none
      
      integer*4 :: nb_fnc,id_fnc,i,nb_data_fnc, idx !nb_timeval 
      integer*4 :: ind_start, ind_end
      
      integer*4, dimension(nb_fnc) :: type_fnc
      integer*4, dimension(nb_fnc +1) :: ind_fnc
      
      real*8 :: PI,t_t0,t0,t1,v0,v1,omega, fp, fac
      real*8 :: TAU,scaling,HDUR 
      real*8 :: amp, ps0, tplus, alpha,time,beta2,dist,vel
      real*8 :: wp_w, ts_w
      
      real*8, dimension(nb_data_fnc) :: data_fnc
      real*8, dimension(1) :: valmax


      GET_FUNC_VALUE = 0.0d0

      PI = 4.0d0 * datan(1.0d0)


 !     write(*,*) type_fnc(id_fnc)
 !     write(*,*) nb_fnc
 !     write(*,*)  type_fnc
 !     write(*,*)  ind_fnc
 !     write(*,*)  data_fnc
 !     write(*,*)  nb_data_fnc
 !     write(*,*)  id_fnc
 !     read(*,*)

      select case (type_fnc(id_fnc))
         
         case(0)
           GET_FUNC_VALUE = 1.0d0

         case(1)
           t_t0 = time - data_fnc(ind_fnc(id_fnc) +1) 
           GET_FUNC_VALUE = (1.0d0 - 2.0d0*data_fnc(ind_fnc(id_fnc))*t_t0*t_t0) &
                 * dexp(-1.0d0*data_fnc(ind_fnc(id_fnc))*t_t0*t_t0)
         
         case(2)
           PI = 4.0d0 * datan(1.0d0)
           t_t0 = time - data_fnc(ind_fnc(id_fnc) +1)
           GET_FUNC_VALUE = dcos(PI*data_fnc(ind_fnc(id_fnc))*t_t0) &
                 * dexp(-0.5d0*data_fnc(ind_fnc(id_fnc)) &
                 * data_fnc(ind_fnc(id_fnc))*t_t0*t_t0)
      
         case(3)
           ind_start = ind_fnc(id_fnc); ind_end = ind_fnc(id_fnc+1)-3;

!           write(*,*) id_fnc, ind_fnc
!           read(*,*)
!           print*, ind_start, data_fnc(ind_start)
!           print*, ind_end, data_fnc(ind_end-1)
          ! read*
          ! print*, data_fnc(ind_start:ind_end:2)
          ! read*
           
        !   open(300,file='th.out',position='APPEND')
        !   do i = ind_start, ind_end, 2
        !       write(300,*) data_fnc(i), data_fnc(i+1)
        !   enddo
        !   close(300)  

           
           !valmax = maxval(data_fnc(ind_start:ind_end:2))
           if (time >= data_fnc(ind_end-1)) then
               v1 = data_fnc(ind_end + 2);
               GET_FUNC_VALUE = v1;
           else
              idx = binarySearch_real(data_fnc(ind_start:ind_end:2), time)
              
!              write(*,*) idx
              
              t0 = data_fnc(ind_start-1 + 2*idx-1);        t1 = data_fnc(ind_start-1 + 2*idx+1)
              v0 = data_fnc(ind_start-1 + 2*idx);          v1 = data_fnc(ind_start-1 + 2*idx+2)
              
              GET_FUNC_VALUE = (v1 - v0) / (t1 - t0) * (time - t0)  + v0

           endif

             
                      
!           do i = ind_fnc(id_fnc), ind_fnc(id_fnc+1) -3,2
!              t0 = data_fnc(i);    t1 = data_fnc(i +2)
!              v0 = data_fnc(i +1);  v1 = data_fnc(i +3)
!              if ((time.ge.t0) .and. (time .le. t1))  then
!               write(*,*) time, t0, t1, v0, v1, GET_FUNC_VALUE
             !  open(300,file='th.out',position='APPEND')
             !  write(300,*) time, GET_FUNC_VALUE
             !  close(300)  
             ! read*

!                 GET_FUNC_VALUE = (v1 - v0) / (t1 - t0) * (time - t0)  + v0
!                 val2 = (v1 - v0) / (t1 - t0) * (time - t0)  + v0
!                 if (abs(val2-val1) .ne. 0) then 
!                    print*, time, val1, val2 
!                    read*
!                 endif
!                 return
!              endif    
!           enddo

         case(4)
           t_t0 = time - data_fnc(ind_fnc(id_fnc) +1)
           beta2 = data_fnc(ind_fnc(id_fnc))
           GET_FUNC_VALUE = 2.0d0*beta2*t_t0 &
                 * (-3.0d0 + 2.0d0*beta2*t_t0*t_t0) &
                 * dexp(-beta2*t_t0*t_t0)
      
         case(5)
           t_t0 = time - data_fnc(ind_fnc(id_fnc) +1)
           fp = data_fnc(ind_fnc(id_fnc));
           alpha = 2.d0*(pi*fp)**2
           fac = 2*pi*fp*dsqrt(dexp(1.d0))
           GET_FUNC_VALUE = fac*t_t0*dexp(-alpha*t_t0**2);
           
!          write(*,*) time, data_fnc(ind_fnc(id_fnc) +1),  fp, alpha, fac, GET_FUNC_VALUE
!          read(*,*)
           

         case(6)
           PI = 4.0d0 * datan(1.0d0)
           t_t0 = time - data_fnc(ind_fnc(id_fnc) +1)
           GET_FUNC_VALUE = (2.d0*PI*data_fnc(ind_fnc(id_fnc))*t_t0) &
                 * dexp(-0.5d0*4.d0*PI*PI*data_fnc(ind_fnc(id_fnc)) &
                 * data_fnc(ind_fnc(id_fnc))*t_t0*t_t0)


         case(7)
           t_t0 = time - data_fnc(ind_fnc(id_fnc) +1) 
           GET_FUNC_VALUE = (6*data_fnc(ind_fnc(id_fnc)) &
                           -24*data_fnc(ind_fnc(id_fnc))**2 * (t_t0**2) &
                            +8.0d0*data_fnc(ind_fnc(id_fnc))**3 * (t_t0**4)) &
                 * dexp(-1.0d0*data_fnc(ind_fnc(id_fnc))*t_t0*t_t0)


         case(8)
           
           omega = data_fnc(ind_fnc(id_fnc))*2*pi
           GET_FUNC_VALUE = dcos(omega*time)

         case(9)
           
           omega = data_fnc(ind_fnc(id_fnc))*2*pi
           GET_FUNC_VALUE = dsin(omega*time)

         

         case(12)
           !------------------------------------------------
           ! 12 - sigmf(t,[a c]) = amp*(1/(1+exp(-a*(t-c))))
           ! 
           !
           !   
           !  |                                                                                                
           !  |............*************************......amp      
           !  |          ** : 
           !  |     a   *   :
           !  |        *    :
           !  |      **     :   
           !  0******----------------------------------> Time
           !       |        |
           !       |____c___|                          |
           !       |        |                          |
           !                                       n*(1/f)+t0
           ! 
           t_t0 = time - data_fnc(ind_fnc(id_fnc) +2)
           GET_FUNC_VALUE = data_fnc(ind_fnc(id_fnc)) &
                           * (1/(1+exp(-data_fnc(ind_fnc(id_fnc) +1)*(t_t0))))

         case(13)
           ! - GRENOBLE BENCHMARK
           TAU = data_fnc(ind_fnc(id_fnc));  scaling = data_fnc(ind_fnc(id_fnc) +1)
           t0 = 2.0 * TAU;  HDUR = TAU/2.0;  t_t0 = time - t0
           GET_FUNC_VALUE = 0.5d0 * ( 1.0d0 + erf(scaling*(t_t0)/HDUR) )
 
         case(14)
           ! - SCEC BENCHMARK
           TAU = data_fnc(ind_fnc(id_fnc));  amp = data_fnc(ind_fnc(id_fnc) +1)
           t_t0 = time
           if (t_t0 .lt. 0.0d0) then
             GET_FUNC_VALUE = 0.0d0
           elseif (t_t0 .eq. 0.0d0) then
             GET_FUNC_VALUE = 0.5d0
           else
             GET_FUNC_VALUE = amp * (1.0d0 - (1 + t_t0/TAU)*exp(-t_t0/TAU))
           endif

         case(15)
           ! - EXPLOSION
           ps0 = data_fnc(ind_fnc(id_fnc)) 
           tplus = data_fnc(ind_fnc(id_fnc) +1); alpha = data_fnc(ind_fnc(id_fnc) +2)
           t_t0 = time - data_fnc(ind_fnc(id_fnc) +3)
           if (t_t0 .lt. 0.0d0) then
             GET_FUNC_VALUE = 0.0d0
           else
             GET_FUNC_VALUE = ps0 * (1 - (1 + time/tplus)*exp(-alpha*time/tplus))
           endif

!         case(30)
!           t_t0 = time - dist/vel;
!           if (t_t0 .le. 0) then 
!               GET_FUNC_VALUE = 0.d0;
!           else 
!              GET_FUNC_VALUE = 2.0d0*beta2*t_t0 &
!                 * (-3.0d0 + 2.0d0*beta2*t_t0*t_t0) &
!                 * dexp(-beta2*t_t0*t_t0)
!           endif
           
         case(30)
         
 !          write(*,*) 'Sono qui'
 !          write(*,*) ind_fnc(id_fnc),ind_fnc(id_fnc+1) -3,2
 !          read(*,*)

         
           do i = ind_fnc(id_fnc),ind_fnc(id_fnc+1) -3,2
              t0 = data_fnc(i)   ;  t1 = data_fnc(i +2)
              v0 = data_fnc(i +1);  v1 = data_fnc(i +3)
              
 !             write(*,*) t0,t1,v0,v1
              if (((time -dist/vel) .ge. t0) .and. ((time -dist/vel) .le. t1))  then
                  GET_FUNC_VALUE = (v1 - v0) / (t1 - t0) * (time - t0)  + v0
                  return
              endif    
           enddo

 !          ind_start = ind_fnc(id_fnc); ind_end = ind_fnc(id_fnc+1)-3;

           
 !          if (time - dist/vel >= data_fnc(ind_end-1)) then
 !              v1 = data_fnc(ind_end + 2);
 !              GET_FUNC_VALUE = v1;
 !          else
 !             idx = binarySearch_real(data_fnc(ind_start:ind_end:2), time-dist/vel)
 !             t0 = data_fnc(2*idx-1);        t1 = data_fnc(2*idx+1)
 !             v0 = data_fnc(2*idx);      v1 = data_fnc(2*idx+2)
 !             
 !             GET_FUNC_VALUE = (v1 - v0) / (t1 - t0) * (time - t0)  + v0!
 !
 !          endif






         case(60,62)
           ! FUNCTION FOR G/G0
           do i = ind_fnc(id_fnc),ind_fnc(id_fnc+1) -3,2     
             
              t0 = data_fnc(i);     t1 = data_fnc(i +2)                    
              v0 = data_fnc(i +1);  v1 = data_fnc(i +3)                
    
              if (abs(time).le.data_fnc(ind_fnc(id_fnc))) then          
                GET_FUNC_VALUE = data_fnc(ind_fnc(id_fnc)+1)            
              elseif ((abs(time).ge.t0).and.(abs(time).le.t1)) then     
                GET_FUNC_VALUE = (v1 - v0) / (t1 - t0) * (abs(time) - t0)  + v0
             elseif (abs(time).ge.data_fnc(ind_fnc(id_fnc+1)-2)) then   
                GET_FUNC_VALUE = data_fnc(ind_fnc(id_fnc+1)-1)          
              endif                                                 
           enddo

         case(61,63)
           ! FUNCTION FOR DAMPING
           do i = ind_fnc(id_fnc),ind_fnc(id_fnc+1) -3,2      
              t0 = data_fnc(i);    t1 = data_fnc(i +2)                  
              v0 = data_fnc(i +1); v1 = data_fnc(i +3)                  

              if (abs(time).le.data_fnc(ind_fnc(id_fnc))) then          
                GET_FUNC_VALUE = data_fnc(ind_fnc(id_fnc)+1)            
              elseif ((abs(time).ge.t0).and.(abs(time).le.t1)) then    
                GET_FUNC_VALUE = (v1 - v0) / (t1 - t0) * (abs(time) - t0)  + v0 
              elseif (abs(time).ge.data_fnc(ind_fnc(id_fnc+1)-2)) then 
                GET_FUNC_VALUE = data_fnc(ind_fnc(id_fnc+1)-1)          
              endif                                                 
           enddo

         case(99)
           ! - CASHIMA1 BENCHMARK                                           
           TAU = data_fnc(ind_fnc(id_fnc));  amp = data_fnc(ind_fnc(id_fnc) +1)
           t_t0 = time                                                  
           GET_FUNC_VALUE = ( exp( - (((2.0d0 * 2.0d0*dasin(1.0d0))*1.5d0)&         
                            *(t_t0 - TAU)/amp)**2.0d0)&                
                            * cos(((2.0d0 * 2.0d0*dasin(1.0d0))*1.5d0)& 
                            *(t_t0 - TAU) + 2.0d0*dasin(1.0d0)/2.0d0))  
      
         case(100)
          ! - TESTMODE
          amp = data_fnc(ind_fnc(id_fnc))       
          GET_FUNC_VALUE = dsin(amp*PI*time)
          
          
         case(101) 
          GET_FUNC_VALUE = time

         case(47)   !Gabor Wavelet
           wp_w = (8.d0*datan(1.d0)) * data_fnc(ind_fnc(id_fnc));
           ts_w = 0.45 * data_fnc(ind_fnc(id_fnc)+1) / data_fnc(ind_fnc(id_fnc));
           GET_FUNC_VALUE = dexp(-(wp_w*(time-ts_w)/data_fnc(ind_fnc(id_fnc)+1))**2) * &
                            dcos(data_fnc(ind_fnc(id_fnc)+2) + wp_w*(time-ts_w))

         case default
           GET_FUNC_VALUE = 0.d0
      
      end select

     
      return
      
      end function GET_FUNC_VALUE

