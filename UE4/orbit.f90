PROGRAM ORBIT
    USE constants
    IMPLICIT NONE

!   Define position and velocity Vector/Arrays
    REAL(KIND=p), DIMENSION(2) :: r, v, dummy

!   Time and time step
    REAL(KIND=p) :: dt
    INTEGER :: N, MAXN

!   Intermediate quantities for Runge-Kutta (RK4)
    REAL(KIND=p), DIMENSION(2) :: k

!   Define simulation length and time step size first
	dt=0.01
	MAXN=365*100
	
    WRITE(*,*)'#Welcome to Orbit calculator!'
    WRITE(*,*)'#Time step is set to:       ',dt,' days.'
    WRITE(*,*)'#Number of steps is set to: ',MAXN,'.'

!   Define initial paramters
!    CALL INIT_EARTH(r,v)
!   
!    OPEN(unit=11,file="EarthOrbit_Euler_001.dat")
!    OPEN(unit=12,file="EarthOrbit_Energy_Euler_001.dat")
!    WRITE(*,*)' '
!    WRITE(*,*)'#Calculating EarthOrbit using first order Euler method.'
!    WRITE(*,*)'#Results will be written to EarthOrbit_Euler.dat and EarthOrbit_Energy_Euler.dat files.'
!    DO N=1,MAXN
!	
!	dummy = r
!	r = r + dt * v
!	v = v + dt * FORCE_SUN(Ms,dummy)
!
!   WRITE(11,*) N*dt, r
!    WRITE(12,*) N*dt, ENERGY(Me,r,v)
	
!    ENDDO

!    CLOSE(11)
!    CLOSE(12)


!  Define initial paramters for earth again
   CALL INIT_EARTH(r,v)

   OPEN(unit=11,file="EarthOrbit_RK4_001.dat")
   OPEN(unit=12,file="EarthOrbit_Energy_RK4_001.dat")
   WRITE(*,*)' '
   WRITE(*,*)'#Calculating EarthOrbit using RK4 method.'
   WRITE(*,*)'#Results will be written to EarthOrbit_RK4.dat and EarthOrbit_Energy_RK4.dat files.'

   DO N=1,MAXN

!     Integrate r and v using RK4
		dummy = r
		k = v
		r = r + (dt/6) * k
		k = v + 0.5*dt*k
		r = r + (dt/3) * k
		k = v + 0.5*dt*k
		r = r + (dt/3) * k
		k = v + dt * k
		r = r + (dt/6) * k	

		k = FORCE_SUN(Ms,dummy)
		v = v + (dt/6) * k
		k = FORCE_SUN(Ms,dummy+0.5*dt*k)
		v = v + (dt/3) * k
		k = FORCE_SUN(Ms,dummy+0.5*dt*k)
		v = v + (dt/3) * k
		k = FORCE_SUN(Ms,dummy+dt*k)
		v = v + (dt/6) * k
		
		k0_r=r
		k0_v=v

		k1_r=dt*k0_v
		k1_v=dt*FORCE_SUN(Me,k0_r)

		k2_r=dt*(k0_v+0.5*k1_v)
		k2_v=dt*FORCE_SUN(Me,k0_r+0.5*k1_r)

		k3_r=dt*(k0_v+0.5*k2_v)
		k3_v=dt*FORCE_SUN(Me,k0_r+0.5*k2_r)

		k4_r=dt*(k0_v+k3_v)
		k4_v=dt*FORCE_SUN(Me,k0_r+k3_r)

		r = k0_r + 1./6*(k1_r+2*k2_r+2*k3_r+k4_r)
		v = k0_v + 1./6*(k1_v+2*k2_v+2*k3_v+k4_v)/Me

	  IF ((abs(r(2))<10E-3).and.(r(1)<0)) THEN
		write(*,*) 'Successss      ', r(1) 
	  END IF

      WRITE(11,*) N*dt, r
      WRITE(12,*) N*dt, ENERGY(Me,r,v)
   ENDDO
   
   CLOSE(11)
   CLOSE(12)



CONTAINS

! This function returns the gravitational force of the sun acting on a mass m at coordinate  r  
FUNCTION FORCE_SUN(m,r)
    USE constants
    IMPLICIT NONE
    REAL(KIND=p), DIMENSION(2) :: FORCE_SUN
    REAL(KIND=p), DIMENSION(2) :: r      ! position of mass
    REAL(KIND=p) :: m                    ! and mass m

	FORCE_SUN = -G*m*1/(sqrt(r(1)**2+r(2)**2))**3 * r	

    RETURN
END

! This function returns the kinetic and potential energy of a mass m at coordinate  r  in the grav. field of the sun
FUNCTION ENERGY(m,r,v)
   USE constants
   IMPLICIT NONE
   REAL(KIND=p) :: ENERGY
   REAL(KIND=p) :: m
   REAL(KIND=p), DIMENSION(2) :: r,v

!   GESAMTENERGIE MUSS HIER RICHTIG BERECHNET WERDEN.
	ENERGY = (v(1)**2+v(2)**2)*0.5*m - G*m*Ms/sqrt(r(1)**2+r(2)**2)

   RETURN
END

! This routine defines the initial r and v values for the earth.
SUBROUTINE INIT_EARTH(r,v)
   IMPLICIT NONE
   REAL(KIND=p), DIMENSION(2) :: r,v

!   ANFANGSWERTE MUESSEN HIER EINGESETZT WERDEN
   r = (/ 1.0 , 0.0 /)
   v = (/ 0.0 , -0.017326 /)

END


END PROGRAM



