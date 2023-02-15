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


!> @brief ...Writing VTK file to visualise in Paraview
!! @author Srihari Sangaraju
!> @date July, 2021 
!> @version 1.0

!> @param[in] loc_n_num. Global node number of 'i'th local node is loc_n_num(i)
!> @param[in] nn_loc. No. of nodes in Local/Current Partition
!> @param[in] nmat_nhe No. of Blocks specified with NHE case
!> @param[in] nhe_mat Tag/Labels of Blocks where NHE has to be implemented
!> @param[in] xs_loc x-coordinate of spectral nodes
!> @param[in] ys_loc y-coordinate of spectral nodes
!> @param[in] zs_loc z-coordinate of spectral nodes
! propname = exaple 'MASS', 'FORCEZ', 'DAMPMATRIX'
    subroutine   WRITE_VTK_MECH_PROP(nn_loc, cs_nnz_loc, cs_loc, &
                                  nmat, sdeg, prop_mat, tag_mat, &
                                  nmat_nlp, tag_mat_nlp, & 
                                  xx_loc, yy_loc, zz_loc, mpi_id, nn, vtk_numbering_map)

      implicit none

      character*70 :: file_name, prop_name, temp_char

      integer*4 :: nel_loc, nn_loc, cs_nnz_loc, mpi_id, nmat, nmat_nlp, nn
      integer*4 :: ie, inode, im, im_nlp, i, j, k
      integer*4 :: unit_mpi
    !   integer*4 :: ic1, ic2, ic3, ic4, ic5, ic6, ic7, ic8
      integer*4, dimension(nmat) :: sdeg, tag_mat
      integer*4, dimension(nn*nn*nn) :: vtk_numbering_map, node_numbering_vtkwrite, loc_nod_indx
      integer*4, dimension(nmat_nlp) :: tag_mat_nlp
      integer*4, dimension(0:cs_nnz_loc) :: cs_loc
      integer*4, dimension(nn_loc) :: nlp_flag
    !   integer*4, dimension(:), allocatable :: nlp_flag_el
      
      real*8, dimension(nmat,4) :: prop_mat
      real*8, dimension(nn_loc) :: rho, lambda, mu
      real*8, dimension(nn_loc) :: xx_loc, yy_loc, zz_loc
      
      nel_loc = cs_loc(0) - 1

    !   allocate(nlp_flag_el(nel_loc))
      rho = 0.d0; lambda = 0.d0; mu = 0.d0; nlp_flag = 0;

      if (mpi_id.eq.0) write(*,'(A)')
      if (mpi_id.eq.0) write(*,'(A)')'------Writing VTK file - SCALAR ----------' 
        
      prop_name = 'MECH_PROP_NLP'
      write(file_name,'(a,i5.5,a)') trim(prop_name),mpi_id,'.vtk'
      unit_mpi = 2500 + mpi_id                                 
      
      !----------------------------------------------------------------------
      open(unit_mpi,file=file_name)
      write(unit_mpi,'(a)') '# vtk DataFile Version 3.1'
      write(unit_mpi,'(a)') 'material model VTK file'
      write(unit_mpi,'(a)') 'ASCII'
      write(unit_mpi,'(a)') 'DATASET UNSTRUCTURED_GRID'
      write(unit_mpi, '(a,i12,a)') 'POINTS ', nn_loc, ' float'

      ! Node Coordinates
      do inode=1,nn_loc
          write(unit_mpi,'(3e20.12)') xx_loc(inode),yy_loc(inode),zz_loc(inode)
      enddo
      write(unit_mpi,*) ''

      ! Connectivity (note: node indices for vtk start at 0)
      write(unit_mpi,'(a,i12,i12)') "CELLS ",nel_loc,nel_loc*(nn*nn*nn + 1)

      write(temp_char,*)'(i12,',nn*nn*nn,'i12)'
      do ie=1,nel_loc
        im = cs_loc(cs_loc(ie -1) + 0 )

        ! nn = sdeg(im) +1
        do i=1,(nn*nn*nn)
            loc_nod_indx(i) = cs_loc(cs_loc(ie -1) + i) - 1
        enddo
        do i=1,(nn*nn*nn)
            node_numbering_vtkwrite(i) = loc_nod_indx(vtk_numbering_map(i))
        enddo

        ! Writing Element-node connectivity
        ! write(unit_mpi,'(9i12)') 8, ic2, ic6, ic7, ic3, &       ! Bottom Surface (front-left node, and then anticlosckwise)
        !                             ic1, ic5, ic8, ic4      ! Top Surface
       
        write(unit_mpi,temp_char) nn*nn*nn, (node_numbering_vtkwrite(i), i=1,(nn*nn*nn))

        ! dum_int = 0
        ! do im_nlp=1,nmat_nlp
        !     if (tag_mat_nlp(im_nlp) .eq. tag_mat(im)) dum_int = 1
        ! enddo

        do k=1,nn
            do j=1,nn
                do i=1,nn
                    inode = cs_loc( cs_loc(ie-1) + nn*nn*(k-1) + nn*(j-1) + i )
                    rho(inode) = prop_mat(im,1); lambda(inode) = prop_mat(im,2);  
                    mu(inode) = prop_mat(im,3); !nlp_flag(inode) = nlp_flag_el(ie);
                enddo
            enddo
        enddo

      enddo
      write(unit_mpi,*) ''

      ! vtkCellType: hexahedrons (ID = 12 for linear 8 noded hex)
      ! ID = 29 for HEX27 (triquadratic Hex)
      ! 72 for lagrange-higher order hex
      ! Reference : https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
      ! https://visit-sphinx-github-user-manual.readthedocs.io/en/v3.3.0/data_into_visit/VTKFormat.html
      write(unit_mpi,'(a,i12)') "CELL_TYPES ",nel_loc
      write(unit_mpi,'(6i12)') (72,ie=1,nel_loc)
      write(unit_mpi,*) ''

      ! Writing node data------------------------------------------------------------------
      write(unit_mpi,'(a,i12)') "POINT_DATA ",nn_loc

      ! Density
      write(unit_mpi,'(a)') 'SCALARS DENSITY float'
      write(unit_mpi,'(a)') "LOOKUP_TABLE default"
      do inode = 1,nn_loc
          write(unit_mpi,*) rho(inode)
      enddo
      write(unit_mpi,*) ''

      ! lambda
      write(unit_mpi,'(a)') 'SCALARS Lambda float'
      write(unit_mpi,'(a)') "LOOKUP_TABLE default"
      do inode = 1,nn_loc
          write(unit_mpi,*) lambda(inode)
      enddo
      write(unit_mpi,*) ''

      ! mu
      write(unit_mpi,'(a)') 'SCALARS Mu float'
      write(unit_mpi,'(a)') "LOOKUP_TABLE default"
      do inode = 1,nn_loc
          write(unit_mpi,*) mu(inode)
      enddo
      write(unit_mpi,*) ''

    !   ! nonlinear flag
    !   write(unit_mpi,'(a)') 'SCALARS NLP_FLAG integer'
    !   write(unit_mpi,'(a)') "LOOKUP_TABLE default"
    !   do inode = 1,nn_loc
    !       write(unit_mpi,*) nlp_flag(inode)
    !   enddo
    !   write(unit_mpi,*) ''
      !-----------------------------------------------------------------------------------

      ! Writing Element data------------------------------------------------------------------
    !   write(unit_mpi,'(a,i12)') "CELL_DATA ",nel_loc

    !   ! nonlinear flag
    !   write(unit_mpi,'(a)') 'SCALARS NLP_FLAG_EL integer'
    !   write(unit_mpi,'(a)') "LOOKUP_TABLE default"
    !   do ie = 1,nel_loc
    !       write(unit_mpi,*) nlp_flag_el(ie)
    !   enddo
    !   write(unit_mpi,*) ''
      !-----------------------------------------------------------------------------------


      close(unit_mpi)   
      !------------------------------------------------------------------

      if (mpi_id.eq.0) write(*,'(A)')'Completed.' 

    end subroutine WRITE_VTK_MECH_PROP





    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Maps the conversion -> speed node numbering to VTK nodenumbering for VTK_LAGRANGE_HEXAHEDRAN of order n
    ! Reference: https://www.kitware.com/main/wp-content/uploads/2018/09/Source_Issue_43.pdf
    ! https://gitlab.kitware.com/vtk/vtk/-/issues/17746
    subroutine VTK_NODE_NUMBERING_MAP(nn, vtk_numbering)  
  
        integer*4, intent(in) :: nn
        integer*4, dimension(nn*nn*nn), intent(out) :: vtk_numbering

        integer*4 :: i, j, k, i1, i2, i3, ncount
       
        vtk_numbering = 0;

        if (nn.le.1) call exit()

        ! 8 Vertices (first 4 belong to bottom face), next 4 belong to top face
        k= 1; j=nn; i=nn;  vtk_numbering(1) = nn*nn*(k-1) + nn*(j-1) + i;
        k=nn; j=nn; i=nn;  vtk_numbering(2) = nn*nn*(k-1) + nn*(j-1) + i;
        k=nn; j=nn; i= 1;  vtk_numbering(3) = nn*nn*(k-1) + nn*(j-1) + i;
        k= 1; j=nn; i= 1;  vtk_numbering(4) = nn*nn*(k-1) + nn*(j-1) + i;

        k= 1; j= 1; i=nn;  vtk_numbering(5) = nn*nn*(k-1) + nn*(j-1) + i;
        k=nn; j= 1; i=nn;  vtk_numbering(6) = nn*nn*(k-1) + nn*(j-1) + i;
        k=nn; j= 1; i= 1;  vtk_numbering(7) = nn*nn*(k-1) + nn*(j-1) + i;
        k= 1; j= 1; i= 1;  vtk_numbering(8) = nn*nn*(k-1) + nn*(j-1) + i;

        if (nn.le.2) return

        ! 12 Edges
        do i1 = 2,(nn-1)
            !Bottom face
            ! Edge 1 - Bottom - X1
            k=i1; j=nn; i=nn;  vtk_numbering(8 + i1 -1) = nn*nn*(k-1) + nn*(j-1) + i;

            ! Edge 2 - Bottom - Y2
            k=nn; j=nn; i=nn-i1+1;  vtk_numbering(8 + (nn-2) + i1 -1) = nn*nn*(k-1) + nn*(j-1) + i;

            ! Edge 3 - Bottom - X2
            k=i1; j=nn; i= 1;  vtk_numbering(8 + 2*(nn-2) + i1 -1) = nn*nn*(k-1) + nn*(j-1) + i;

            ! Edge 4 - Bottom - Y2
            k= 1; j=nn; i=nn-i1+1;  vtk_numbering(8 + 3*(nn-2) + i1-1) = nn*nn*(k-1) + nn*(j-1) + i;


            !Top Face
            ! Edge 5 - Top - X1
            k=i1; j= 1; i=nn;  vtk_numbering(8 + 4*(nn-2) + i1-1) = nn*nn*(k-1) + nn*(j-1) + i;

            ! Edge 6 - Top - Y2
            k=nn; j= 1; i=nn-i1+1;  vtk_numbering(8 + 5*(nn-2) + i1-1) = nn*nn*(k-1) + nn*(j-1) + i;

            ! Edge 7 - Top - X2
            k=i1; j= 1; i= 1;  vtk_numbering(8 + 6*(nn-2) + i1-1) = nn*nn*(k-1) + nn*(j-1) + i;

            ! Edge 8 - Top - Y2
            k= 1; j= 1; i=nn-i1+1;  vtk_numbering(8 + 7*(nn-2) + i1-1) = nn*nn*(k-1) + nn*(j-1) + i;


            !Vertical Edges
            ! Edge 9 - front-left
            k= 1; j=nn-i1+1; i=nn;  vtk_numbering(8 + 8*(nn-2) + i1-1) = nn*nn*(k-1) + nn*(j-1) + i;

            ! Edge 10 - front-right
            k=nn; j=nn-i1+1; i=nn;  vtk_numbering(8 + 9*(nn-2) + i1-1) = nn*nn*(k-1) + nn*(j-1) + i;

            ! Edge 11 - back-left
            k= 1; j=nn-i1+1; i= 1;  vtk_numbering(8 + 10*(nn-2) + i1-1) = nn*nn*(k-1) + nn*(j-1) + i;

            ! Edge 12 - back-right
            k=nn; j=nn-i1+1; i= 1;  vtk_numbering(8 + 11*(nn-2) + i1-1) = nn*nn*(k-1) + nn*(j-1) + i;
        enddo


        !Other Exterior nodes present on 6 faces
        ncount = 0;
        do i1 = 2,(nn-1)
            do i2 = 2,(nn-1)
                ncount = ncount+1;

                !1st Face - side left (Xmin)
                k= 1; j=nn-i1+1; i=nn-i2+1;  vtk_numbering(8 + 12*(nn-2) + ncount) = nn*nn*(k-1) + nn*(j-1) + i;

                !2nd Face - side left (Xmax)
                k=nn; j=nn-i1+1; i=nn-i2+1;  vtk_numbering(8 + 12*(nn-2) + 1*(nn-2)*(nn-2) + ncount) = nn*nn*(k-1) + nn*(j-1) + i;

                !3rd Face - front (Ymin)
                k=i2; j=nn-i1+1; i=nn;  vtk_numbering(8 + 12*(nn-2) + 2*(nn-2)*(nn-2) + ncount) = nn*nn*(k-1) + nn*(j-1) + i;

                !4th Face - Back (Ymax)
                k=i2; j=nn-i1+1; i= 1;  vtk_numbering(8 + 12*(nn-2) + 3*(nn-2)*(nn-2) + ncount) = nn*nn*(k-1) + nn*(j-1) + i;

                !5th Face - Bottom (Zmin)
                k=i2; j=nn; i=nn-i1+1;  vtk_numbering(8 + 12*(nn-2) + 4*(nn-2)*(nn-2) + ncount) = nn*nn*(k-1) + nn*(j-1) + i;

                !6th Face - Top (Zmax)
                k=i2; j=1; i=nn-i1+1;  vtk_numbering(8 + 12*(nn-2) + 5*(nn-2)*(nn-2) + ncount) = nn*nn*(k-1) + nn*(j-1) + i;
            enddo
        enddo

        ! Volume nodes
        ncount = 0;
        do i1=2,(nn-1)
            do i2=2,(nn-1)
                do i3=2,(nn-1)
                    ncount=ncount+1;
                    k=i3; j=nn-i1+1; i=nn-i2+1;  vtk_numbering(8 + 12*(nn-2) + 6*(nn-2)*(nn-2) + ncount) = nn*nn*(k-1) + nn*(j-1) + i;
                enddo
            enddo
        enddo
      
    end subroutine VTK_NODE_NUMBERING_MAP  
