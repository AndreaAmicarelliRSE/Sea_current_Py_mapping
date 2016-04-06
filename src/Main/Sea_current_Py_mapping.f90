PROGRAM Sea_current_Py_mapping
!-------------------------------------------------------------------------------
! "Sea current Py mapping v.1.0" 
! Copyright 2016 (RSE SpA)
! "Sea current Py mapping v.1.0" authors and email contact are provided on 
! the documentation file.
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Description: "Sea current Py mapping v.1.0" is a minor tool, which reads a
!              formatted ".nc" file from the free dataset of MetOcean-INGV at 
!              www.myocean.eu on daily mean current synoptic-scale velocity (of 
!              a given year) and provides the yearly average of the specific 
!              power flow of marine currents at the same scale 
!              (P=500(radq(u_x^2+u_y^2))^3. The output is both in binary and 
!              formatted ".vtk" file format.
!-------------------------------------------------------------------------------  
!------------------------
! Modules
!------------------------ 
use netcdf
!------------------------
! Declarations
!------------------------
implicit none
! It reads 3D data, a 72 x 677 x 253 matrix each day
integer,parameter :: ndepth=72,nlat=253,nlon=677,ntime=1
! These will be the netCDF IDs for the file and data variable
integer :: ncid,varid
integer :: i,j,k,IOstatus,ndata,read_coord,nfile
real :: vomecrty_max,vomecrty_min,vozocrtx_max,vozocrtx_min
character(100) FILE_NAME
integer n_power(nlon,nlat,ndepth)
real power(nlon,nlat,ndepth)
! Dynamic matrices for reading from ".nc" files
real,dimension(:,:,:,:),allocatable :: vozocrtx_in,vomecrty_in
real,dimension(:),allocatable :: lon,lat,depth
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
allocate(vozocrtx_in(nlon,nlat,ndepth,ntime))
allocate(vomecrty_in(nlon,nlat,ndepth,ntime))
allocate(lon(nlon))
allocate(lat(nlat))
allocate(depth(ndepth))
!------------------------
! Initializations
!------------------------
vomecrty_max = -1000.
vomecrty_min = 1000.
vozocrtx_max = -1000.
vozocrtx_min = 1000.
do i=1,nlon
   do j=1,nlat
      do k=1,ndepth
         power(i,j,k) = 0.
         n_power(i,j,k) = 0
      enddo
   enddo
enddo
read_coord = 1
nfile = 0
!------------------------
! Statements
!------------------------
write(*,*)                                                                     &
"Sea current Py mapping v.1.0 is a minor tool, which reads a formatted ",      &
".nc file from the free dataset of MetOcean-INGV at www.myocean.eu on daily ", &
"mean current synoptic-scale velocity (of a given year) and provides the ",    &
"yearly average of the specific power flow of marine currents ",               &
"(P=500(radq(u_x^2+u_y^2))^3. The output is both in binary and formatted ",    &
".vtk file format."
! Reading and pre-processing (file loop start)
write(*,*) "   Start reading nc files and pre-processing data"
! Start reading file name list
write(*,*) "   Start reading file name list "
open(1,file='file_name_list.txt')
do
   read (1,'(a)',IOSTAT=IOstatus) FILE_NAME
   if (IOstatus/=0) exit
! Start reading nc file
   write(*,*) "   Start reading nc file ",FILE_NAME
   nfile = nfile + 1
! It opens the ".nc" file and tells netCDF it wants read-only access
   call check(nf90_open(FILE_NAME,NF90_NOWRITE,ncid))
! Get the "varid" of "u", based on its name.
   call check(nf90_inq_varid(ncid,"vozocrtx",varid))
! To read "u"
   call check(nf90_get_var(ncid,varid,vozocrtx_in))
! To get the "varid" of "v", based on its name.
   call check(nf90_inq_varid(ncid,"vomecrty",varid))
! To read "v"
   call check(nf90_get_var(ncid,varid,vomecrty_in))
   if (read_coord==1) then
! To get the "varid" of "lonu", based on its name.
! For MyOcean dataset
      call check(nf90_inq_varid(ncid,"lonu",varid))
! For Copernicus dataset
!      call check(nf90_inq_varid(ncid,"lon",varid))
! To read "lon"
      call check(nf90_get_var(ncid,varid,lon))
! To get the "varid" of "latu", based on its name.
! For MyOcean dataset
      call check(nf90_inq_varid(ncid,"latu",varid))
! For Copernicus dataset
!      call check(nf90_inq_varid(ncid,"lat",varid))
! To read "lat"
      call check(nf90_get_var(ncid,varid,lat))
! To get the "varid" of "depth", based on its name.
      call check(nf90_inq_varid(ncid,"depth",varid))
! To read "depth"
      call check(nf90_get_var(ncid,varid,depth))
      read_coord = 0
   endif
! To close the file, freeing all resources.
   call check(nf90_close(ncid))
   write(*,*) "   End reading nc file ",FILE_NAME
! End reading file
! Start pre-processing file
   write(*,*) "   Start pre-processing nc file ",FILE_NAME
! First elaboration loop: treatment of the missing values, data check (max and 
! min), to compute 24h specific power flow, contribution to compute the yearly 
! average.
   do k=1,ndepth
      do j=1,nlat
         do i=1,nlon  
            if (vozocrtx_in(i,j,k,1)>1000) then
               vozocrtx_in(i,j,k,1) = -999.
               else
                  if (vozocrtx_in(i,j,k,1)>vozocrtx_max) vozocrtx_max =        &
                     vozocrtx_in(i,j,k,1)
                  if (vozocrtx_in(i,j,k,1)<vozocrtx_min) vozocrtx_min =        &
                     vozocrtx_in(i,j,k,1)
            endif
            if (vomecrty_in(i,j,k,1)>1000) then
               vomecrty_in(i,j,k,1) = -999.
               else
                  if (vomecrty_in(i,j,k,1)>vomecrty_max) vomecrty_max =        &
                     vomecrty_in(i,j,k,1)
                  if (vomecrty_in(i,j,k,1)<vomecrty_min) vomecrty_min =        &
                     vomecrty_in(i,j,k,1)
            endif
            if ((vozocrtx_in(i,j,k,1)/=-999.).and.                             &
               (vomecrty_in(i,j,k,1)/=-999.)) then
               power(i,j,k) = power(i,j,k) + 500. * sqrt(vozocrtx_in(i,j,k,1)  &
                              ** 2 + vomecrty_in(i,j,k,1) ** 2) ** 3
               n_power(i,j,k) = n_power(i,j,k) + 1
            endif
         enddo
      enddo
   enddo
   write(*,*) "   End pre-processing nc file ",FILE_NAME
! End pre-processing file
enddo
close(1)
write(*,*) "   End reading file name list "
write(*,*) "   u_max(m/s) ",vozocrtx_max," u_min(m/s) ",vozocrtx_min
write(*,*) "   v_max(m/s) ",vomecrty_max," v_min(m/s) ",vomecrty_min
write(*,*) "End reading nc files and pre-processing data"
! End reading and pre-processing data
! Yearly specific power flow computation
write(*,*) "   Start yearly power computation "
do k=1,ndepth
   do j=1,nlat
      do i=1,nlon 
         if (n_power(i,j,k)>=(nfile/2)) then
            power(i,j,k) = power(i,j,k) / n_power(i,j,k)
            else
               power(i,j,k) = -999.
         endif
      enddo
   enddo
enddo
write(*,*) "End yearly power computation "
! To write the 3D yearly power ".dat" binary files
write(*,*) "   Start writing .dat binary output"
open(2,file='lon.dat',form='unformatted')
write(2) lon
close(2)
open(2,file='lat.dat',form='unformatted')
write(2) lat
close(2)
open(2,file='depth.dat',form='unformatted')
write(2) depth
close(2)
open(2,file='yearly_power.dat',form='unformatted')
write(2) power
close(2)
write(*,*) "End writing .dat binary output"
! test start: ASCII format 
! Write 3D yearly power .txt file
!  print *,"Start writing .txt output"
!  open (2,file='yearly_power.txt')
!  write (2,"(1X,a,(F12.4,1X),a,(F12.4,1X))") '         lon         lat       depth         P(W) vomecrty_max(m/s) ',vomecrty_max,' vomecrty_min(m/s) ',vomecrty_min
!  do i=1,nlon
!     do j=1,nlat
!        do k=1,ndepth
!           write (2,"(1X,4(F12.4,1X))") lon(i),lat(j),depth(k),power(i,j,k)
!        end do
!     end do
!  end do
!  close (2)
!  print *,"End writing .txt output"
! test end: ASCII format 
! Write 3D yearly power on a ".vtk" file 
write(*,*) "   Start writing .vtk output"
open(3,file='yearly_power.vtk')
write(3,'(a)') '# vtk DataFile Version 3.0'
write(3,'(a)') 'YEARLY POWER'
write(3,'(a)') 'ASCII'
write(3,'(a)') 'DATASET RECTILINEAR_GRID'
write(3,'(a,3(i5,1X))') 'DIMENSIONS ',nlon,nlat,ndepth
ndata = nlat * nlon * ndepth
write(3,'(a,i10,a)') 'X_COORDINATES ',nlon,' float '
do i=1,nlon
   write(3,'(F12.4)',ADVANCE='NO') lon(i)
enddo
write(3,'(/a,i10,a)') 'Y_COORDINATES ',nlat,' float '
do i=1,nlat
   write(3,'(F12.4)',ADVANCE='NO') lat(i)
enddo
write(3,'(/a,i10,a)') 'Z_COORDINATES ',ndepth,' float '
do i=1,ndepth
   write(3,'(F12.4)',ADVANCE='NO') depth(i)
enddo
write(3,'(/a,i10)') 'POINT_DATA ',ndata
write(3,'(a)') 'SCALARS A float '
write(3,'(a)') 'LOOKUP_TABLE default ' 
! Start test: structured grid
!  write (3,'(a)') 'DATASET STRUCTURED_GRID'
!  write (3,'(a,3(i5,1X))') 'DIMENSIONS ',nlon,nlat,ndepth
!  ndata = nlat*nlon*ndepth
!  write (3,'(a,i10,a)') 'POINTS ', ndata, ' float '
!  do k=1,ndepth
!     do j=1,nlat
!        do i=1,nlon
!           write (3,'(3(F12.4))') lon(i),lat(j),depth(k)
!        end do
!     end do
!  end do  
!  write (3,'(a,i10)') 'POINT_DATA ', ndata
!  write (3,'(a)') 'SCALARS A float '
!  write (3,'(a)') 'LOOKUP_TABLE default '
! End test: structured grid
do k=1,ndepth
   do j=1,nlat
      do i=1,nlon 
         write(3,'(F12.4)') power(i,j,k)
      enddo
   enddo
enddo  
close(3) 
write(*,*) "End writing .vtk output" 
write(*,*) "Sea current Py mapping v.1.0 has terminated. "
!------------------------
! Deallocations
!------------------------
deallocate(vozocrtx_in)
deallocate(vomecrty_in)
deallocate(lon)
deallocate(lat)
deallocate(depth)
contains
   subroutine check(error_status)
      integer,intent(in) :: error_status
      if(error_status/=nf90_noerr) then
         write(*,*) trim(nf90_strerror(error_status))
         stop
      endif
   end subroutine check
end program Sea_current_Py_mapping

