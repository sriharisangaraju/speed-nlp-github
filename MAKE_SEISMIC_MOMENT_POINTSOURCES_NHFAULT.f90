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

!> @brief Generates extend seismic source as Points sources
!> @ Usually we need to do Triangular Mesh of fault in cubit, and then 3ptool, 
!>   creates a triangles using SISM command. But here SISM command can just give
!>   Hypocenter, Points where Sesimic moment has to be applied.
!>   KNNSearch is used to fing Nearest GLL point of each point source
!! @author Srihari Sangaraju
!> @date September, 2020 
!> @version 1.0

      subroutine MAKE_SEISMIC_MOMENT_POINTSOURCES_NHFAULT()

      use max_var
      use speed_par
      use kdtree2_module
      

      implicit none
 
      include 'SPEED.MPI'      

      integer*4 :: n_neighbours, ipt, indxL
      character*70 :: file_srcout

      real(kdkind) :: query_vec(3)
      real(kdkind), dimension(:,:), allocatable :: nodes_in_xyz
      type(kdtree2), pointer :: kd2_obj_locnode
      type(kdtree2_result) :: result_temp(1)

!*****************************************************************************************
!                                    SEISMIC MOMENT
!*****************************************************************************************
                                                                      
      if (mpi_id.eq.0 .and. nload_sism_el.gt.0) & 
           write(*,'(A)') '-----------Building the Seismic Moment vector----------'

      if (nload_sism_el.gt.0) then

         if (mpi_id.eq.0) then
               file_srcout = 'srcmod_points.out'
               open(444,file=file_srcout)
         endif

       
         allocate (num_node_sism(nload_sism_el))
         num_node_sism = 1
         max_num_node_sism = 1

         allocate (sour_node_sism(max_num_node_sism,nload_sism_el))
         allocate (dist_sour_node_sism(max_num_node_sism,nload_sism_el))
         allocate (pos_sour_node_x(max_num_node_sism,nload_sism_el),&
                   pos_sour_node_y(max_num_node_sism,nload_sism_el),&
                   pos_sour_node_z(max_num_node_sism,nload_sism_el))
         
         sour_node_sism = 0
         dist_sour_node_sism = 0.d0
        
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Defining Kd-tree Object and then search
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         allocate(nodes_in_xyz(3,nnod_loc))
         do ipt = 1, nnod_loc
            nodes_in_xyz(1,ipt) = xx_spx_loc(ipt)
            nodes_in_xyz(2,ipt) = yy_spx_loc(ipt)
            nodes_in_xyz(3,ipt) = zz_spx_loc(ipt)
         enddo

         kd2_obj_locnode => kdtree2_create(nodes_in_xyz,sort=.false.,rearrange=.true.)
         n_neighbours = 1



         !Searching the node 'id' in the global numeration for each seismic faults.
         !sour_node_sism = node id (global numeration) generating the 'i'th seismic fault
         do i = 1,nload_sism_el

            query_vec(1) = val_sism_el(i,4)
            query_vec(2) = val_sism_el(i,5)
            query_vec(3) = val_sism_el(i,6)
            call kdtree2_n_nearest(kd2_obj_locnode,query_vec,n_neighbours,result_temp)
            indxL = result_temp(1)%idx    ! Inder of Nearest Neighbor

            sour_node_sism(1,i) = local_node_num(indxL)
            pos_sour_node_x(1,i) = xx_spx_loc(indxL)
            pos_sour_node_y(1,i) = yy_spx_loc(indxL)
            pos_sour_node_z(1,i) = zz_spx_loc(indxL)

            ! Check dist_sour_node_sism - In Other parts of code its distance between Hypocenter to Node
            dist_sour_node_sism(1,i) = (val_sism_el(i,4)-pos_sour_node_x(1,i))*(val_sism_el(i,4)-pos_sour_node_x(1,i)) + &
                                       (val_sism_el(i,5)-pos_sour_node_y(1,i))*(val_sism_el(i,5)-pos_sour_node_y(1,i)) + &
                                       (val_sism_el(i,6)-pos_sour_node_z(1,i))*(val_sism_el(i,6)-pos_sour_node_z(1,i))


            ! Collecting Nearest Neighbours From all Processors and selecting the nearest one
            allocate(sism_el_glo(mpi_np), dist_el_glo(mpi_np))
            allocate(posx_el_glo(mpi_np), posy_el_glo(mpi_np),posz_el_glo(mpi_np))                
            
            call MPI_BARRIER(mpi_comm, mpi_ierr)                  
            call MPI_ALLGATHER(sour_node_sism(1,i), 1, SPEED_INTEGER, sism_el_glo, 1, SPEED_INTEGER, mpi_comm, mpi_ierr)
            call MPI_ALLGATHER(dist_sour_node_sism(1,i), 1, SPEED_DOUBLE, dist_el_glo, 1, SPEED_DOUBLE, mpi_comm, mpi_ierr)
            call MPI_ALLGATHER(pos_sour_node_x(1,i), 1, SPEED_DOUBLE, posx_el_glo, 1, SPEED_DOUBLE, mpi_comm, mpi_ierr)
            call MPI_ALLGATHER(pos_sour_node_y(1,i), 1, SPEED_DOUBLE, posy_el_glo, 1, SPEED_DOUBLE, mpi_comm, mpi_ierr)
            call MPI_ALLGATHER(pos_sour_node_z(1,i), 1, SPEED_DOUBLE, posz_el_glo, 1, SPEED_DOUBLE, mpi_comm, mpi_ierr)

            call GET_MINVALUES(sism_el_glo, dist_el_glo, mpi_np, k, 1, mpi_np) 


            sour_node_sism(1,i) =  sism_el_glo(k);         
            dist_sour_node_sism(1,i) =  dsqrt(dist_el_glo(k))
            pos_sour_node_x(1,i) = posx_el_glo(k);
            pos_sour_node_y(1,i) = posy_el_glo(k);
            pos_sour_node_z(1,i) = posz_el_glo(k);

            deallocate(sism_el_glo, dist_el_glo, posx_el_glo, posy_el_glo, posz_el_glo) 


            if (mpi_id.eq.0 .and. i .eq. 1) write(*,'(A)')
            if (mpi_id.eq.0 .and. i .eq. 1) write(*,'(A)')'Seismic faults'
            if (mpi_id.eq.0 .and. i .eq. 1) write(*,'(A)')
            j = 1;
            if (mpi_id.eq.0) write(*,'(A,I6,A,I6,A)')'Seismic Faults - ',i, ' is generated by ',j,' nodes'
            do k = 1, j
               if (mpi_id.eq.0) write(*,*) 'Node :', sour_node_sism(k,i)
               if (mpi_id.eq.0) write(*,*) 'Position :', pos_sour_node_x(k,i), pos_sour_node_y(k,i), pos_sour_node_z(k,i)  
               if (mpi_id.eq.0) write(*,*)  'Distance from Actual Node [m] = ', dist_sour_node_sism(k,i)
               if (mpi_id.eq.0) write(*,*)  'Seismic Moment [Nm] = ', val_sism_el(k,13), val_sism_el(k,14), val_sism_el(k,15)
            enddo     
            if (mpi_id.eq.0) write(*,'(A)')


            if (mpi_id.eq.0) then
               write(444,*)val_sism_el(i,4),val_sism_el(i,5),val_sism_el(i,6),dist_sour_node_sism(1,i),&
                           pos_sour_node_x(1,i), pos_sour_node_y(1,i), pos_sour_node_z(1,i)
            endif

 
         enddo !i = 1,nload_sism_el

         deallocate(nodes_in_xyz)
         call kdtree2_destroy(kd2_obj_locnode)


      endif !if (nload_sism_el.gt.0)
    
      if (nfunc.le.0) nfunc = 1
      
      if (nload_sism_el.gt.0) then                         
         allocate (factor_seismic_moment(nload_sism_el,6)) 
         allocate (tau_seismic_moment(nload_sism_el,1)) 
      endif                                                 

       if ((mpi_id.eq.0).and.(nload_sism_el.gt.0)) write(*,'(A)')'Seismic Moment vector built'
       call MPI_BARRIER(mpi_comm, mpi_ierr)    


      end subroutine MAKE_SEISMIC_MOMENT_POINTSOURCES_NHFAULT

