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

!> @brief Writes output results for Restart.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in,out] file_name directory where saving files
!> @param[in] count index for snapshot
!> @param[in] proc mpi process id
!> @param[in] nv number of values to print
!> @param[in] vec values to print
!> @param[in] xx,yy,zz coordinates of the grid nodes to print
!> @param[in] local to global numbering for nodes

      subroutine WRITE_FILEOUT_GRID_SV(file_name,file_xyz,count, proc, nv,ncol, &
                                        vec,xx,yy,zz,loc_n_num,tstart)
      

      
      character*70 :: file_name, file_xyz

      integer*4 :: count,proc,nv, ncol
      
      real*8 :: tstart
      real*8, dimension(nv,ncol) :: vec
      
      real*8, dimension(nv/6) :: xx,yy,zz
      integer*4, dimension(nv/6) :: loc_n_num
      
      character*70 :: out_file, xyz_file
      integer*4 :: i,lname, lnamexyz
      
      lname = len_trim(file_name)
      lnamexyz = len_trim(file_xyz)
      
      out_file = file_name(1:lname) // '000000_000000.out'
      xyz_file = file_xyz(1:lname) // '000000.out'
      
      if (proc .lt. 10) then
         write(out_file(lname+6:lname+6),'(i1)') proc
         if(tstart .eq. 0.d0) write(xyz_file(lname+6:lname+6),'(i1)') proc
      else if (proc .lt. 100) then
         write(out_file(lname+5:lname+6),'(i2)') proc
         if(tstart .eq. 0.d0) write(xyz_file(lname+5:lname+6),'(i2)') proc
      else if (proc .lt. 1000) then
         write(out_file(lname+4:lname+6),'(i3)') proc     
         if(tstart .eq. 0.d0) write(xyz_file(lname+4:lname+6),'(i3)') proc     
      else if (proc .lt. 10000) then
         write(out_file(lname+3:lname+6),'(i4)') proc      
         if(tstart .eq. 0.d0) write(xyz_file(lname+3:lname+6),'(i4)') proc      
      else if (proc .lt. 100000) then
         write(out_file(lname+2:lname+6),'(i5)') proc
         if(tstart .eq. 0.d0) write(xyz_file(lname+2:lname+6),'(i5)') proc
      else
         write(out_file(lname+1:lname+6),'(i6)') proc
         if(tstart .eq. 0.d0) write(xyz_file(lname+1:lname+6),'(i6)') proc
      endif
      
      if (count .lt. 10) then
         write(out_file(lname+13:lname+13),'(i1)') count
      else if (count .lt. 100) then
         write(out_file(lname+12:lname+13),'(i2)') count
      else if (count .lt. 1000) then
         write(out_file(lname+11:lname+13),'(i3)') count    
      else if (count .lt. 10000) then
         write(out_file(lname+10:lname+13),'(i4)') count    
      else if (count .lt. 100000) then
         write(out_file(lname+9:lname+13),'(i5)') count    
      else
         write(out_file(lname+8:lname+13),'(i6)') count
      endif
               
      open(20+proc, file=out_file)
      do i = 1,nv
         write(20+proc,*)  vec(i,:)
      enddo
      close(20+proc)


!         write(20+proc,*) loc_n_num(i), xx(i), yy(i), zz(i), &
!                                        vec(3*(i-1)+1), vec(3*(i-1)+2),vec(3*(i-1)+3)
      open(20+proc, file=xyz_file)
      do i = 1,nv/6
         write(20+proc,*)  loc_n_num(i), xx(i), yy(i), zz(i)
      enddo
      close(20+proc)

      
      
      return
      
      end subroutine WRITE_FILEOUT_GRID_SV

