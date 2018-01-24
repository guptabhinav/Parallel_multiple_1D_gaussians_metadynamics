PROGRAM fes
IMPLICIT NONE

REAL*8 :: deltaT, t              ! inputs

REAL*8, ALLOCATABLE :: hx(:), hy(:), sx(:),sy(:), xtau(:), ytau(:)
INTEGER, ALLOCATABLE :: bins(:)

REAL*8 :: dsx, dsy, alpha, disp2_x, disp2_y, x, y, vb
INTEGER ::  mtd_step, lx, ly, dummy1, i, j, p, m, ios

OPEN(1,FILE="input")
OPEN(2,FILE="metadynamics_file")
OPEN(10,FILE="data_input",STATUS="REPLACE")
OPEN(20,FILE="output",STATUS="REPLACE")

mtd_step=0
do
read(2,*,iostat=ios)
if (ios.ne.0) exit
mtd_step=mtd_step+1
end do

REWIND(2)

print *, 'mtd_step=', mtd_step

ALLOCATE(xtau(mtd_step))
ALLOCATE(ytau(mtd_step))
ALLOCATE(sx(3))
ALLOCATE(sy(3))
ALLOCATE(hx(mtd_step))
ALLOCATE(hy(mtd_step))
ALLOCATE(bins(2))

READ(1,*) sx(1:3)                       ! xmin, xmax, xbinsize          
READ(1,*) sy(1:3)                       ! ymin, ymax, ybinsize
READ(1,*) deltaT, t                     ! temp_factor, temperature


alpha= (t+deltaT)/deltaT
print*, 'alpha=', alpha


DO i=1,mtd_step
  READ(2,*)   dummy1,hx(i),hy(i),dsx,dsy,xtau(i),ytau(i)
  write(10,*) dummy1,hx(i),hy(i),dsx,dsy,xtau(i),ytau(i)
END DO


bins(1) = nint( (sx(2)-sx(1))/sx(3) )
bins(2) = nint( (sy(2)-sy(1))/sy(3) )

if(bins(1) .le.0) stop 'error in grid definitions'
if(bins(2) .le.0) stop 'error in grid definitions'

WRITE(*,*) "bin size =",bins(1:2)

print *, 'working'

DO lx=1,bins(1)                                                                      ! sx(*) number of grids along x
   x=sx(1)+sx(3)*float(lx-1)                      
       DO ly=1,bins(2) 
         y=sy(1)+sy(3)*float(ly-1)                                                      ! sy(*) number of grids along y
         vb=0.0                                                                         !    IF x(l,2).le.x(l,2)+ABS(bin_size) THEN
           DO i=1,mtd_step
           disp2_x=0.5*(x-xtau(i))**2
           disp2_y=0.5*(y-ytau(i))**2   ! exponential term for bias
           vb=vb+hx(i)*dexp(-disp2_x/dsx**2)+hy(i)*dexp(-disp2_y/dsy**2)
           END DO
       WRITE(20,'(3f16.6)') x, y, -vb*alpha                      ! alpha is wt_mtd factor
       END DO
write(20,*)
END DO

CLOSE(1)
CLOSE(2)
CLOSE(10)
CLOSE(20)
DEALLOCATE(xtau)
DEALLOCATE(ytau)
DEALLOCATE(sx)
DEALLOCATE(sy)
DEALLOCATE(hx)
DEALLOCATE(hy)


END PROGRAM
