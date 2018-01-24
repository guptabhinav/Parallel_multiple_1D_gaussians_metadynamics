PROGRAM mtd
IMPLICIT NONE
REAL*8 :: x,vel_x,force_x,hx0,hy0,dsx,dsy,dd,deltaT
REAL*8 :: y,vel_y,force_y,vel,vbias
REAL*8, ALLOCATABLE :: xtau(:), ytau(:), hx(:), hy(:)
REAL*8, PARAMETER :: m=100000.0d0, k=0.001d0, kb=3.16e-6            ! au
REAL*8 :: system_t,bath_temp,au_t,dt,t,ket,pet,tot_energy,a,b,vel_sum,f,au_tbath,noise
INTEGER :: nsteps,i,j,d,natoms,mtd_step,c,n,inter_meta,mtd_time_step
LOGICAL :: displaced
a=1.0d0

OPEN(10,FILE="coordinate_file",STATUS="REPLACE")
OPEN(20,FILE="energy_file",STATUS="REPLACE")
OPEN(40,FILE="metadynamics_file",STATUS="REPLACE")
OPEN(50,FILE="noise_file",STATUS="REPLACE")
OPEN(60,FILE="hills",STATUS="REPLACE")
OPEN(1,FILE="input")
READ(1,*) natoms, system_t
READ(1,*) bath_temp, dt
READ(1,*) nsteps,f
READ(1,*) hx0,hy0,dsx,dsy,deltaT

deltaT=deltaT*kb
print *, 'deltaT=',deltaT
mtd_step=0

ALLOCATE(xtau(nsteps))
ALLOCATE(ytau(nsteps))
ALLOCATE(hx(nsteps))
ALLOCATE(hy(nsteps))
hx=0.d0
hy=0.d0

PRINT*,'height of gaussain hx and hy =',hx0,hy0,'width of gaussian =',dsx,dsy
  
  DO i=1,natoms
  READ(1,*) x,y
  END DO


au_t = system_t*kb

!...setting initial velocities...................................................................
  CALL temp(bath_temp,kb,au_tbath)
  CALL rand_vel(vel,natoms)
vel_x = vel*sqrt(m/(2.d0*3.14*au_tbath))                                   ! velocities at bath temperature
PRINT *, "vel_x", vel_x  
CALL rand_vel(vel,natoms)
vel_y = vel*sqrt(m/(2.d0*3.14*au_tbath))
PRINT *, "vel_x", vel_y  

 call kin_ene(ket,vel_x,vel_y,m,kb,t,natoms)                               ! calling the temperature in kelvin
print *, 'temperature of randomly generated velocity=', t
 vel_x=sqrt(bath_temp/t)*vel_x
 vel_y=sqrt(bath_temp/t)*vel_y
print *, 'velocity after scaling vx=', vel_x, 'and vy=', vel_y
!.....................................................................................................
PRINT*, 'bath Tempearture in kelvin=',bath_temp, 'and in atomic units=', au_tbath
  
  CALL fion_x(x,force_x,m,a,k,kb,natoms,f,au_tbath,dt,noise,mtd_step,xtau,hx,dsx)
  CALL fion_y(y,force_y,m,b,k,kb,natoms,f,au_tbath,dt,noise,mtd_step,ytau,hy,dsy)

c = 0
inter_meta=0
mtd_time_step=100

  DO i=1,nsteps
  inter_meta=inter_meta+1
  
  CALL new_x(x,vel_x,force_x,natoms,dt)
  CALL new_y(y,vel_y,force_y,natoms,dt)
  CALL new_vel_x(vel_x,force_x,dt,natoms,f,m)   
  CALL new_vel_y(vel_y,force_y,dt,natoms,f,m)   
 
  CALL fion_x(x,force_x,m,a,k,kb,natoms,f,au_tbath,dt,noise,mtd_step,xtau,hx,dsx)
  CALL fion_y(y,force_y,m,b,k,kb,natoms,f,au_tbath,dt,noise,mtd_step,ytau,hy,dsy)
  CALL new_vel_x(vel_x,force_x,dt,natoms,f,m)
  CALL new_vel_y(vel_y,force_y,dt,natoms,f,m)

  CALL kin_ene(ket,vel_x,vel_y,m,kb,t,natoms)

!  call check_displacement(x,xtau(mtd_step),1.5d0*ds,dd,displaced)
!  call check_displacement(y,ytau(mtd_step),1.5d0*ds,dd,displaced)
 displaced = .false.

     IF(((i-1).eq.0).or.(inter_meta.ge.mtd_time_step).or.displaced) THEN       ! displaced is logical here. 
     mtd_step=mtd_step+1

         call eval_vbias(x,xtau,mtd_step,dsx,hx,vbias)
          hx(mtd_step)=hx0*dexp(-vbias/deltaT)             ! deltaT is in atomic units.

         call eval_vbias(y,ytau,mtd_step,dsy,hy,vbias)
          hy(mtd_step)=hy0*dexp(-vbias/deltaT)             ! deltaT=kb*deltaT is in atomic units.

    inter_meta=0
    ytau(mtd_step)=y
    xtau(mtd_step)=x
    write(40,'(i10,6f16.6)') i,hx(mtd_step),hy(mtd_step),dsx,dsy,xtau(mtd_step),ytau(mtd_step)
  END IF
   
  CALL pot_ene(pet,x,y,a,b,natoms,k,mtd_step,xtau,ytau,i,hx,dsx,hy,dsy)
tot_energy = ket + pet

   IF (MOD(c,100) == 0) THEN
       DO j=1,natoms
       WRITE(10,*) i,x,y,vel_x,vel_y        
       WRITE(20,*) i,x,y,ket,pet,tot_energy,t     
       WRITE(30,*) i,x,y,force_x,force_y     
       END DO   
   END IF
c = c+1
END DO

DEALLOCATE(xtau)
DEALLOCATE(ytau)
DEALLOCATE(hx)
DEALLOCATE(hy)
CONTAINS

SUBROUTINE temp(bath_t,kb,au_tbath)
implicit none
REAL*8 :: bath_t,kb,au_tbath 
au_tbath = bath_t*kb
END SUBROUTINE

SUBROUTINE new_x(x,vel_x,force_x,natoms,dt)
implicit none
REAL*8 :: x,vel_x,force_x
REAL*8 :: dt
INTEGER :: natoms,i
!  DO i=V1,natoms
      x= x + dt*vel_x + 0.5*dt**2*force_x
!PRINT*,i
!  END DO
END SUBROUTINE

SUBROUTINE new_y(y,vel_y,force_y,natoms,dt)
implicit none
REAL*8 :: y,vel_y,force_y
REAL*8 :: dt
INTEGER :: natoms,i
!  DO i=1,natoms
      y= y + dt*vel_y + 0.5*dt**2*force_y
!PRINT*,i
!  END DO
END SUBROUTINE

SUBROUTINE fion_x(x,force_x,m,a,k,kb,natoms,f,au_tbath,dt,noise,mtd_step,xtau,h,ds)
implicit none
REAL*8 :: x,force_x,xtau(*),h(*)
REAL*8 :: a,k,m,f,au_tbath,gauss,ds
REAL*8 :: noise,d1,dt,kb,biased_force
REAL*8, PARAMETER :: pi=3.14d0, v0=0.001d0
INTEGER :: i,l,natoms,mtd_step
 CALL RANDOM_NUMBER(d1)
noise = dsqrt(2*f*au_tbath/dt)*2.0d0*(d1-0.5)                       ! kb is not multiplied because the temperature is already in atomic units.
    DO i=1,natoms
    force_x = -v0*(4.0d0*x/(m*a**4))*(x**2-a**2) + (1.0/m)*noise        ! the noise term is not coming from surface potential. 
         DO l=1,mtd_step-1
         biased_force = (h(l)/(ds**2))*(x - xtau(l))*exp(-((x-xtau(l))**2/(2*(ds**2))))
         force_x =  force_x + biased_force/m
         END DO
    END DO
END SUBROUTINE

SUBROUTINE fion_y(y,force_y,m,b,k,kb,natoms,f,au_tbath,dt,noise,mtd_step,ytau,h,ds)
implicit none
REAL*8 :: y,force_y,ytau(*),h(*)
REAL*8 :: b,k,m,f,au_tbath,gauss,ds
REAL*8 :: noise,d1,dt,kb,biased_force
REAL*8, PARAMETER :: pi=3.14d0,v0=0.001d0
INTEGER :: i,l,natoms,mtd_step
 CALL RANDOM_NUMBER(d1)
noise = dsqrt(2*f*au_tbath/dt)*2.0d0*(d1-0.5)                ! kb is not multiplied because the temperature is already in atomic units f= kg/dt units 
    DO i=1,natoms
    force_y = -k*y/m + (1.0/m)*noise        ! the stochastic term is not coming from surface potential. 
         DO l=1,mtd_step-1
         biased_force = (h(l)/(ds**2))*(y - ytau(l))*exp(-((y-ytau(l))**2)/(2*(ds**2)))
         force_y =  force_y + biased_force/m
         END DO
    END DO
END SUBROUTINE

SUBROUTINE new_vel_x(vel_x,force_x,dt,natoms,f,m)   
implicit none
REAL*8 :: vel_x,force_x,m
REAL*8 :: dt,f
INTEGER :: natoms,i
   DO i = 1,natoms
    vel_x =  vel_x*(1-0.5*f*dt/m) + 0.5*dt*force_x                       ! this velocity friction is due to stochastic force. kg/dt units 
  END DO
END SUBROUTINE
   
SUBROUTINE new_vel_y(vel_y,force_y,dt,natoms,f,m)   
implicit none
REAL*8 :: vel_y,force_y,m
REAL*8 :: dt,f
INTEGER :: natoms,i
   DO i = 1,natoms
    vel_y =  vel_y*(1-0.5*f*dt/m) + 0.5*dt*force_y                       ! this velocity friction is due to stochastic force. kg/dt units
  END DO
END SUBROUTINE
  
SUBROUTINE rand_vel(vel,natoms)
implicit none
REAL*8 :: vel
REAL*8 :: a1,a2,pi
INTEGER :: i,j,natoms
pi=3.14d0
     DO i=1,natoms
            CALL RANDOM_NUMBER(a1)
            CALL RANDOM_NUMBER(a2)
            vel = (-2*log(a1))**0.5*cos(2*pi*a2)
     END DO
END SUBROUTINE

SUBROUTINE pot_ene(pet,x,y,a,b,natoms,k,mtd_step,xtau,ytau,i,hx,dsx,hy,dsy)
implicit none
REAL*8 :: x,y, xtau(*),ytau(*),hx(*),hy(*)
REAL*8 :: pet,a,k,bias,b,v0,dsx,dsy
INTEGER :: natoms,i,l,mtd_step,p
bias = 0.0d0
v0=0.001d0
      DO p=1,natoms
      pet=  v0*((x**2 - a**2)/(a**2))**2+ 0.5d0*k*y**2
      
         DO l=1,mtd_step-1
         bias = bias + hx(l)*exp(-(x-xtau(l))**2/2*dsx**2)+ hy(l)*exp(-(y-ytau(l))**2/(2*dsy**2))
         END DO
      pet = pet + bias
      END DO
END SUBROUTINE

SUBROUTINE kin_ene(ket,vel_x,vel_y,m,kb,t,natoms)
implicit none
REAL*8 :: vel_x,vel_y
REAL*8 :: ket,m,t,kb
INTEGER :: natoms,i
   DO i=1,natoms
   ket =  0.5d0*m*(vel_x**2+vel_y**2) 
   END DO
t = ket/kb
END SUBROUTINE

  subroutine check_displacement(x,xtau,dis_cut,dd,displaced)
  implicit none
  real*8 :: x, xtau,dis_cut,dd
  logical:: displaced
  displaced=.false.
  dd=dabs(x-xtau)
  if(dd.ge.dis_cut)displaced=.true.
  end subroutine 

subroutine eval_vbias(x,xtau,mtd_step,ds,h,vbias)
implicit none
real*8 :: x,ds,vbias,xtau(*),h(*)
integer :: mtd_step,l
vbias=0.d0
         DO l=1,mtd_step-1
         vbias = vbias + h(l)*exp(-((x-xtau(l))**2/2*ds**2))
         END DO
end subroutine


END PROGRAM
