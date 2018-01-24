PROGRAM fes
IMPLICIT NONE

REAL*8 :: deltaT, t, x, y, vb

REAL*8, ALLOCATABLE :: ht(:), sx(:), sy(:), xtau(:), ytau(:)
INTEGER, ALLOCATABLE :: bins(:)

REAL*8 :: ds, dummy, alpha, disp2

INTEGER ::  mtd_step, lx, ly, k, i, j, p, m, ios

OPEN(1,FILE="input")
OPEN(2,FILE="metadynamics_file")
OPEN(10,FILE="data_input",STATUS="REPLACE")
OPEN(20,FILE="gaussian",STATUS="REPLACE")

READ(1,*) sx(1:3)                       ! xmin, xmax, xgrids           
READ(1,*) sy(1:3)                       ! ymin, ymax, ygrids
READ(1,*) deltaT, t                     ! temp_factor, temperature


mtd_step=0
do
read(2,*,iostat=ios)
if (ios.ne.0) exit
mtd_step=mtd_step+1
end do

print *, 'mtd_step=', mtd_step

ALLOCATE(xtau(mtd_step))
ALLOCATE(ytau(mtd_step))
ALLOCATE(sx(mtd_step))
ALLOCATE(sy(mtd_step))
ALLOCATE(ht(mtd_step))

!ALLOCATE(ht(mtd_step))

alpha= (t+deltaT)/deltaT
print*, 'alpha=', alpha

REWIND(2)

DO i=1,mtd_step
  READ(2,*) k, ht(i), dummy, ds, xtau(i), ytau(i)
END DO


bins(1) =nint( (sx(2)-sx(1))/sx(3) )
bins(2) =nint( (sy(2)-sy(1))/sy(3) )

if(bins(1) .le.0) stop 'error in grid definitions'
if(bins(2) .le.0) stop 'error in grid definitions'

WRITE(*,*) "bin size =",bins(1:2)

!PRINT*,bins,mtd_step,smin
DO lx=1,bins(1)
   x=sx(1)+sx(3)*float(lx-1)
   DO ly=1,bins(2) 
     y=sy(1)+sy(3)*float(ly-1)
     vb=0.0                                                                         !    IF x(l,2).le.x(l,2)+ABS(bin_size) THEN
     DO i=1,mtd_step
        disp2=((x-xtau(i))**2+(y-ytau(i))**2)*0.5
        vb=vb+ht(i)*dexp(-disp2/ds**2)
     END DO
     WRITE(20,'(3f16.6)') x, y, -vb*alpha
   END DO
   WRITE(20,*) 
END DO

CLOSE(1)
CLOSE(2)
CLOSE(10)
CLOSE(20)
DEALLOCATE(xtau)
DEALLOCATE(ytau)
DEALLOCATE(sx)
DEALLOCATE(sy)
DEALLOCATE(ht)


END PROGRAM
