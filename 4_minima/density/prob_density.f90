program prob_density
implicit none
real*8,allocatable :: grid_min(:),grid_max(:),grid_width(:),s(:,:)
real*8,allocatable :: prob(:,:)
integer :: nd,ndata,norm,i,j,k,m,ios
integer,allocatable :: ngrid(:),igrid(:)
real*8:: etc, sum

open(1,file='input')
open(2,file='coordinate_file')
open(20,file='read_data',status='replace')
open(50,file='output',status='replace')

read(1,*) nd

allocate(grid_min(nd),grid_max(nd),grid_width(nd))

do i=1,nd
   read(1,*) grid_min(i),grid_max(i),grid_width(i)
end do

ndata=0

do 
  read(2,*,iostat=ios)
  if(ios.ne.0) exit
  ndata=ndata+1
end do
print *, 'ndata=', ndata

rewind(2)

allocate(s(nd,ndata))

! to print the data                                          ! nd columns and ndata rows
do m=1,ndata
  read(2,*) etc,(s(k,m),k=1,nd)
  write(20,*)(s(k,m),k=1,nd)
end do
                                     !computing probability density
allocate(ngrid(nd))

! to print the data
!do k=1,nd
!ngrid(k)=nint((grid_max(k)-grid_min(k))/grid_width(k))       ! number of grids along nd dimensions
!print *, ngrid(k)
!end do


ngrid(1:nd)=nint( ((grid_max(1:nd))- grid_min(1:nd))/grid_width(1:nd) )

   allocate( prob(ngrid(1),ngrid(2)) )

prob=0.d0

norm=0.d0


allocate(igrid(nd))

print *,'grids_x=',ngrid(1), 'grid_y=', ngrid(2)

       do m=1,ndata                                                             ! number of students

!      igrid(1:nd) = nint( (s(1:nd,idata)-grid_min(1:nd))/grid_width(1:nd) )
       igrid(1:nd) = nint( (s(1:nd,m)-grid_min(1:nd))/grid_width(1:nd) + 0.5 )                                                                                

! lets say we have 23 as a value, it should be in the 3 grid point if 0-100 distribution is there and the gris size is 10. but without 0.5 factor it is in 2
! grid.


       prob(igrid(1),igrid(2)) = prob(igrid(1),igrid(2)) + 1.d0                  !     prob(igrid(1),igrid(2),1)=prob(igrid(1),igrid(2),1)+1.d0

       norm=norm+1.d0

   end do


sum =0.d0
do i=1,ngrid(1)
   do j=1,ngrid(2)
   write(50,*)((i*grid_width(1))+grid_min(1)), ((j*grid_width(2))+grid_min(2)), prob(i,j)/norm
   sum = sum + prob(i,j)/norm
   end do
end do
print *, 'sum=', sum
end program prob_density
