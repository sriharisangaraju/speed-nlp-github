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
    subroutine  WRITE_VTK_TIMEDATA(nn_loc, cs_nnz_loc, cs_loc, &
                                  nmat, sdeg, prop_mat, tag_mat, &
                                  nmat_nlp, tag_mat_nlp, & 
                                  xx_loc, yy_loc, zz_loc, mpi_id )

      implicit none

      character*70 :: file_name, prop_name

      integer*4 :: nel_loc, nn_loc, cs_nnz_loc, mpi_id, nmat, nmat_nlp
      integer*4 :: ie, inode, im, nn, im_nlp, i, j, k
      integer*4 :: unit_mpi
      integer*4 :: ic1, ic2, ic3, ic4, ic5, ic6, ic7, ic8
      integer*4, dimension(nmat) :: sdeg, tag_mat
      integer*4, dimension(nmat_nlp) :: tag_mat_nlp
      integer*4, dimension(0:cs_nnz_loc) :: cs_loc
      integer*4, dimension(nn_loc) :: nlp_flag
      integer*4, dimension(:), allocatable :: nlp_flag_el
      
      real*8, dimension(nmat,4) :: prop_mat
      real*8, dimension(nn_loc) :: rho, lambda, mu
      real*8, dimension(:), allocatable :: rho_el, lambda_el, mu_el
      real*8, dimension(nn_loc) :: xx_loc, yy_loc, zz_loc
      
      nel_loc = cs_loc(0) - 1

      allocate(rho_el(nel_loc), lambda_el(nel_loc), mu_el(nel_loc), nlp_flag_el(nel_loc))
      rho_el = 0.d0; lambda_el = 0.d0; mu_el = 0.d0; nlp_flag_el = 0;
      rho = 0.d0; lambda = 0.d0; mu = 0.d0; nlp_flag = 0;

      if (mpi_id.eq.0) write(*,'(A)')
      if (mpi_id.eq.0) write(*,'(A)')'------Writing VTK file - SCALAR ----------' 
        
      prop_name = 'MECH_PROP_NLP_test'
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
      write(unit_mpi,'(a,i12,i12)') "CELLS ",nel_loc,nel_loc*9

      do ie=1,nel_loc
        im = cs_loc(cs_loc(ie -1) + 0 )
        nn = sdeg(im) +1

        ic1 = cs_loc(cs_loc(ie -1) + nn*nn*(1 -1) + nn*(1 -1) + 1) - 1
        ic2 = cs_loc(cs_loc(ie -1) + nn*nn*(1 -1) + nn*(1 -1) + nn) - 1
        ic3 = cs_loc(cs_loc(ie -1) + nn*nn*(1 -1) + nn*(nn -1) + nn) - 1
        ic4 = cs_loc(cs_loc(ie -1) + nn*nn*(1 -1) + nn*(nn -1) + 1) - 1
        ic5 = cs_loc(cs_loc(ie -1) + nn*nn*(nn -1) + nn*(1 -1) + 1) - 1
        ic6 = cs_loc(cs_loc(ie -1) + nn*nn*(nn -1) + nn*(1 -1) + nn) - 1
        ic7 = cs_loc(cs_loc(ie -1) + nn*nn*(nn -1) + nn*(nn -1) + nn) - 1
        ic8 = cs_loc(cs_loc(ie -1) + nn*nn*(nn -1) + nn*(nn -1) + 1) - 1

        ! HEX 8 
        ! Node Numbering Order for VTK
        !
        !   7----6
        !   /|   /|
        ! 4----5 |
        ! | 3--|-2
        ! |/   |/
        ! 0----1
        !
        ! !When using Python scripts to export mesh from cubit (HEX8)
        ! write(unit_mpi,'(9i12)') 8, ic2, ic6, ic7, ic3, ic1, ic5, ic8, ic4

        ! !When using NCDUMPS to export mesh from cubit (HEX8)
        ! write(unit_mpi,'(9i12)') 8, ic1, ic2, ic3, ic4, ic5, ic6, ic7, ic8

        ! HEX27
        !When using Python scripts to export mesh from cubit (HEX27)
        write(unit_mpi,'(9i12)') 8, ic2, ic6, ic7, ic3, &       ! Bottom Surface (front-left node, and then anticlosckwise)
                                    ic1, ic5, ic8, ic4       ! Top Surface

        do im_nlp=1,nmat_nlp
            if (tag_mat_nlp(im_nlp) .eq. tag_mat(im)) nlp_flag_el(ie) = 1
        enddo
        rho_el(ie) = prop_mat(im,1);
        lambda_el(ie) = prop_mat(im,2);
        mu_el(ie) = prop_mat(im,3);

        if (ie.eq.1) write(unit_mpi,'(a)') 'Node ids Start' 
        do k=1,nn
            do j=1,nn
                do i=1,nn
                    inode = cs_loc( cs_loc(ie-1) + nn*nn*(k-1) + nn*(j-1) + i )
                    rho(inode) = prop_mat(im,1); lambda(inode) = prop_mat(im,2);  
                    mu(inode) = prop_mat(im,3); nlp_flag(inode) = nlp_flag_el(ie);
                    if (ie.eq.1) write(unit_mpi,'(i12)') inode
                enddo
            enddo
        enddo
        if (ie.eq.1) write(unit_mpi,'(a)') 'Node ids End' 

      enddo
      write(unit_mpi,*) ''

      ! vtkCellType: hexahedrons (ID = 12 for linear 8 noded hex
      ! ID = 29 for HEX27 (triquadratic Hex)
      ! 72 for lagrange-higher order hex
      ! Reference : https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
      ! https://visit-sphinx-github-user-manual.readthedocs.io/en/v3.3.0/data_into_visit/VTKFormat.html
      write(unit_mpi,'(a,i12)') "CELL_TYPES ",nel_loc
      write(unit_mpi,'(6i12)') (12,ie=1,nel_loc)
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

      ! nonlinear flag
      write(unit_mpi,'(a)') 'SCALARS NLP_FLAG integer'
      write(unit_mpi,'(a)') "LOOKUP_TABLE default"
      do inode = 1,nn_loc
          write(unit_mpi,*) nlp_flag(inode)
      enddo
      write(unit_mpi,*) ''
      !-----------------------------------------------------------------------------------

      ! Writing Element data------------------------------------------------------------------
      write(unit_mpi,'(a,i12)') "CELL_DATA ",nel_loc

      ! Density
      write(unit_mpi,'(a)') 'SCALARS DENSITY_EL float'
      write(unit_mpi,'(a)') "LOOKUP_TABLE default"
      do ie = 1,nel_loc
          write(unit_mpi,*) rho_el(ie)
      enddo
      write(unit_mpi,*) ''

      ! lambda
      write(unit_mpi,'(a)') 'SCALARS LAMBDA_EL float'
      write(unit_mpi,'(a)') "LOOKUP_TABLE default"
      do ie = 1,nel_loc
          write(unit_mpi,*) lambda_el(ie)
      enddo
      write(unit_mpi,*) ''

      ! mu
      write(unit_mpi,'(a)') 'SCALARS Mu_el float'
      write(unit_mpi,'(a)') "LOOKUP_TABLE default"
      do ie = 1,nel_loc
          write(unit_mpi,*) mu_el(ie)
      enddo
      write(unit_mpi,*) ''

      ! nonlinear flag
      write(unit_mpi,'(a)') 'SCALARS NLP_FLAG_EL integer'
      write(unit_mpi,'(a)') "LOOKUP_TABLE default"
      do ie = 1,nel_loc
          write(unit_mpi,*) nlp_flag_el(ie)
      enddo
      write(unit_mpi,*) ''
      !-----------------------------------------------------------------------------------


      close(unit_mpi)   
      !------------------------------------------------------------------

      if (mpi_id.eq.0) write(*,'(A)')'Completed.' 

    end subroutine WRITE_VTK_TIMEDATA
