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
    subroutine  WRITE_VTU_TIMEDATA(nn_loc, cs_nnz_loc, cs_loc, &
                                  nmat, sdeg, prop_mat, tag_mat, &
                                  xx_loc, yy_loc, zz_loc, mpi_id, nmpi, &
                                  nn, vtk_numbering_map, &
                                  its, strain, stress, disp )

      implicit none

      character*70 :: file_name_vtu, temp_char, file_name_pvtu, num_str1, num_str2
      character*200 :: buffer

      integer*4 :: nel_loc, nn_loc, cs_nnz_loc, mpi_id, nmat, nn, its, nmpi
      integer*4 :: ie, inode, im, i
      integer :: i1
      integer*4 :: unit_mpi, pvtu_unit
      integer*4, dimension(nmat) :: sdeg, tag_mat
      integer*4, dimension(nn*nn*nn) :: vtk_numbering_map, node_numbering_vtkwrite, loc_nod_indx
      integer*4, dimension(0:cs_nnz_loc) :: cs_loc
      
      real*8, dimension(nmat,4) :: prop_mat
      real*8, dimension(nn_loc) :: xx_loc, yy_loc, zz_loc
      real*8, dimension(6*nn_loc) :: strain, stress 
      real*8, dimension(3*nn_loc) :: disp

      
      
      nel_loc = cs_loc(0) - 1

      !Coverting strain/stress to satisfy real*4 (Float32 in vtu) range
      do inode=1,6*nn_loc
        if (strain(inode).lt. 1.0e-29) strain(inode) = 1.0e-29;
        if (strain(inode).gt. 1.0e+29) strain(inode) = 1.0e+29;
      enddo

      

      write(file_name_vtu,'(a,i8.8,a,i5.5,a)') './VTKOUT/SNAPSHOT_',its,'.',mpi_id,'.vtu'
      !-------------- VTU FILE -------------------------------------------------------------------------------------
      write(temp_char,*)'(i12,',nn*nn*nn,'i12)'
      
      unit_mpi = 2500 + mpi_id             
      open(unit_mpi,file=file_name_vtu)

        ! write header
        buffer='<?xml version="1.0" encoding="utf-8"?>'
        write(unit_mpi,'(a)')trim(buffer)
        buffer='<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">';
        write(unit_mpi,'(a)')trim(buffer)
        buffer='<UnstructuredGrid>'
        write(unit_mpi,'(a)')trim(buffer)
        write(num_str1,*)nn_loc; write(num_str2,*)nel_loc;
        buffer='<Piece NumberOfPoints="'//trim(adjustl(num_str1))//'" NumberOfCells="'//trim(adjustl(num_str2))//'">'
        write(unit_mpi,'(a)')trim(buffer)
        
        ! Write POINTS Info (Node Coordinates)--------------------------------
        buffer='<Points>'
        write(unit_mpi,'(a)')trim(buffer)
        buffer='<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
        write(unit_mpi,'(a)')trim(buffer)
        do inode=1,nn_loc
            write(unit_mpi,'(3e20.12)') xx_loc(inode),yy_loc(inode),zz_loc(inode)
        enddo
        buffer='</DataArray>'
        write(unit_mpi,'(a)')trim(buffer)
        buffer='</Points>'
        write(unit_mpi,'(a)')trim(buffer)
        !--------------------------------------------------------------------

        ! Write CELL info-------------------------------------------------------
        buffer='<Cells>'
        write(unit_mpi,'(a)')trim(buffer)

        buffer='<DataArray type="Int32" Name="connectivity" format="ascii">'
        write(unit_mpi,'(a)')trim(buffer)
        do ie=1,nel_loc
            im = cs_loc(cs_loc(ie -1) + 0 )
            ! nn = sdeg(im) +1
            do i=1,(nn*nn*nn)
                loc_nod_indx(i) = cs_loc(cs_loc(ie -1) + i) - 1
            enddo
            do i=1,(nn*nn*nn)
                node_numbering_vtkwrite(i) = loc_nod_indx(vtk_numbering_map(i))
            enddo
            write(unit_mpi,temp_char) (node_numbering_vtkwrite(i), i=1,(nn*nn*nn))
        enddo
        buffer='</DataArray>'
        write(unit_mpi,'(a)')trim(buffer)

        buffer='<DataArray type="Int32" Name="offsets" format="ascii">'
        write(unit_mpi,'(a)')trim(buffer)
        write(unit_mpi,'(6i12)') (ie*(nn*nn*nn), ie=1,nel_loc)
        buffer='</DataArray>'
        write(unit_mpi,'(a)')trim(buffer)

        buffer='<DataArray type="Int32" Name="types" format="ascii">'
        write(unit_mpi,'(a)')trim(buffer)
        ! vtkCellType: hexahedrons (ID = 12 for linear 8 noded hex)
        ! ID = 29 for HEX27 (triquadratic Hex)
        ! 72 for lagrange-higher order hex
        ! Reference : https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
        ! https://visit-sphinx-github-user-manual.readthedocs.io/en/v3.3.0/data_into_visit/VTKFormat.html
        write(unit_mpi,'(6i12)') (72,ie=1,nel_loc)
        buffer='</DataArray>'
        write(unit_mpi,'(a)')trim(buffer)

        buffer='</Cells>'
        write(unit_mpi,'(a)')trim(buffer)
        !--------------------------------------------------------------------

        ! Write POINTS data --------------------------------------------
        buffer='<PointData>'
        write(unit_mpi,'(a)')trim(buffer)

        !Strain
        buffer='<DataArray type="Float32" NumberOfComponents="6" Name="Strain" format="ascii">'
        write(unit_mpi,'(a)')trim(buffer)
        do inode = 1,nn_loc
            im = 6*(inode-1)
            write(unit_mpi,'(6e20.12)') strain(im+1), strain(im+2), strain(im+3), strain(im+4), strain(im+5), strain(im+6)
        enddo
        buffer='</DataArray>'
        write(unit_mpi,'(a)')trim(buffer)

        !displacement
        buffer='<DataArray type="Float32" NumberOfComponents="3" Name="Disp" format="ascii">'
        write(unit_mpi,'(a)')trim(buffer)
        do inode = 1,nn_loc
            im = 3*(inode-1)
            write(unit_mpi,'(3e20.12)') disp(im+1), disp(im+2), disp(im+3)
        enddo
        buffer='</DataArray>'
        write(unit_mpi,'(a)')trim(buffer)


        buffer='</PointData>'
        write(unit_mpi,'(a)')trim(buffer)
        !--------------------------------------------------------------------


        ! Write CELLS data --------------------------------------------
        ! buffer='<CellData>'
        ! write(unit_mpi,'(a)')trim(buffer)
        ! buffer='</CellData>'
        ! write(unit_mpi,'(a)')trim(buffer)
        !--------------------------------------------------------------------

        buffer='</Piece>'
        write(unit_mpi,'(a)')trim(buffer)
        buffer='</UnstructuredGrid>'
        write(unit_mpi,'(a)')trim(buffer)
        buffer='</VTKFile>'
        write(unit_mpi,'(a)')trim(buffer)
        close(unit_mpi)


      close(unit_mpi)   
      !-----------------------------------------------------------------------------------------------------------------




      
      !-------------- PVTU FILE -------------------------------------------------------------------------------------
      if (mpi_id.eq.0) then
            write(file_name_pvtu,'(a,i8.8,a)') './VTKOUT/SNAPSHOT_',its,'.pvtu'
            pvtu_unit = 2499          
            open(pvtu_unit,file=file_name_pvtu)

            buffer='<?xml version="1.0" encoding="utf-8"?>'
            write(pvtu_unit,'(a)')trim(buffer)
            buffer='<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">';
            write(pvtu_unit,'(a)')trim(buffer)
            buffer='<PUnstructuredGrid GhostLevel="0">'
            write(pvtu_unit,'(a)')trim(buffer)
            

            buffer='<PPoints>'
            write(pvtu_unit,'(a)')trim(buffer)
            buffer='<PDataArray type="Float32" NumberOfComponents="3" format="ascii"/>'
            write(pvtu_unit,'(a)')trim(buffer)
            buffer='</PPoints>'
            write(pvtu_unit,'(a)')trim(buffer)


            buffer='<PCells>'
            write(pvtu_unit,'(a)')trim(buffer)
            buffer='<PDataArray type="Int32" Name="connectivity" NumberOfComponents="1" format="ascii"/>'
            write(pvtu_unit,'(a)')trim(buffer)
            buffer='<PDataArray type="Int32" Name="offsets" NumberOfComponents="1" format="ascii"/>'
            write(pvtu_unit,'(a)')trim(buffer)
            buffer='<PDataArray type="Int32" Name="types" NumberOfComponents="1" format="ascii"/>'
            write(pvtu_unit,'(a)')trim(buffer)
            buffer='</PCells>'
            write(pvtu_unit,'(a)')trim(buffer)


            buffer='<PPointData>'
            write(pvtu_unit,'(a)')trim(buffer)
            buffer='<PDataArray type="Float32" NumberOfComponents="6" Name="Strain" format="ascii"/>'
            write(pvtu_unit,'(a)')trim(buffer)
            buffer='<PDataArray type="Float32" NumberOfComponents="3" Name="Disp" format="ascii"/>'
            write(pvtu_unit,'(a)')trim(buffer)
            buffer='</PPointData>'
            write(pvtu_unit,'(a)')trim(buffer)


            ! buffer='<PCellData>'
            ! write(pvtu_unit,'(a)')trim(buffer)
            ! buffer='</PCellData>'
            ! write(pvtu_unit,'(a)')trim(buffer)

            do i=0,(nmpi-1)
                write(file_name_vtu,'(a,i8.8,a,i5.5,a)') './SNAPSHOT_',its,'.',i,'.vtu'
                buffer='<Piece Source="'//trim(file_name_vtu)//'"/>'
                write(pvtu_unit,'(a)')trim(buffer)
            enddo

            buffer='</PUnstructuredGrid>'
            write(pvtu_unit,'(a)')trim(buffer)
            buffer='</VTKFile>'
            write(pvtu_unit,'(a)')trim(buffer)
            close(pvtu_unit);
      endif

    end subroutine WRITE_VTU_TIMEDATA

    subroutine WRITE_PVD_TIMEDATA(nmpi, nts, ncount, dt)  

        integer*4, intent(in) :: nmpi, nts, ncount
        integer*4 :: unit_mpi

        real*8 :: dt

        character*70 :: filename_pvtu, filename_pvd, str_time
        character*200 :: buffer

        write(filename_pvd,'(a,i8.8,a)') './VTKOUT/SNAPSHOT.pvd'
        unit_mpi = 2499          
        open(unit_mpi,file=filename_pvd)

        buffer='<?xml version="1.0" encoding="utf-8"?>'
        write(unit_mpi,'(a)')trim(buffer)
        buffer='<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
        write(unit_mpi,'(a)')trim(buffer)
        buffer='<Collection>'
        write(unit_mpi,'(a)')trim(buffer)

        do i=0,nts,ncount
            write(filename_pvtu,'(a,i8.8,a)') './SNAPSHOT_',i,'.pvtu'
            write(str_time,'(f16.6)') i*dt 
            buffer='<DataSet timestep="'//trim(adjustl(str_time))//'" part="001" file="'//trim(filename_pvtu)//'"/>'
            write(unit_mpi,'(a)')trim(buffer)
        enddo

        buffer='</Collection>'
        write(unit_mpi,'(a)')trim(buffer)
        buffer='</VTKFile>'
        write(unit_mpi,'(a)')trim(buffer)
        close(unit_mpi);

    end subroutine WRITE_PVD_TIMEDATA