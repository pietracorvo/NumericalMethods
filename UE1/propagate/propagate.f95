program propagate_psi

	implicit none

	integer, parameter :: ind = 800
	complex :: psi(ind), psi_new(ind), dumc(ind), dumf(ind)
	real :: x(ind), V(ind) = 0	

	real :: dx=0.1, dt=0.01 		! iteration paramter	
	integer :: cycles=20000			!      -- " --	
	real :: k=1						! constant
	complex :: i=(0,1)				! imaginary unit
	integer :: j, l					! loopindices

	! Build spatial vector x(j) and potential V(x(j))
	do j=1,ind,1 
		x(j) = (j-int(0.5*ind))*dx
		! V(j) = -exp(-3*abs(x(j))/(abs(x(j))+1))    ! Yukawa potential, uncomment to use it ...		 	
	end do

	! Build starting psi for t=0
	psi = exp(-0.5*x**2)*exp(i*k*x)
 
!	Open output file
	open(20,file='data_psi.dat')


!	iteration:
!	----------
!	We now want to solve the equation
					
!					Hl * Psi_new = Hr * Psi_old

!	with Hr = Hl^*, with the matrices of the Haimltonians of the Crank-Nicholsons
!	algorithm.

!	To obtain the evolved state vector psi_new, first the matrix product on the 
!	right hand side and then the solution of the resulting system of linear 
!	equations has to be performed.  

	do j=1,cycles,1
!		Matrix product on right hand side (for tridiagonal Hr)
		psi_new(1:ind) = (1-0.5*i*dt/dx**2) * psi(1:ind) - (0.5*i*dt) * V(1:ind) * psi(1:ind)
		psi_new(2:ind) = psi_new(2:ind) + 0.25*i*dt/dx**2 * psi(1:ind-1)
		psi_new(1:ind-1) = psi_new(1:ind-1) + 0.25*i*dt/dx**2 * psi(2:ind)
		
!       Uncomment for optional output
		! if (real(psi_new(1))>0.1) write(*,*) j, real(psi_new(1)) 

!		Thomas algorithm to solve system of lineqs on left hand side (for tridig Hl)  
		dumc(1) = ( -0.25*i*dt/dx**2 ) / ( 1+i*0.5*dt/dx**2+(0.5*i*dt)*V(1)control )
		dumf(1) = psi_new(1) / ( -0.25*i*dt/dx**2 )		
		do l=2,ind,1
			dumc(l) = (-0.25*i*dt/dx**2)/((1+i*0.5*dt/dx**2+(0.5*i*dt)*V(l))-dumc(l-1)*(-0.25*i*dt/dx**2)) 
			dumf(l) = (psi_new(l)-dumf(l-1)*(-0.25*i*dt/dx**2))/((1+i*0.5*dt/dx**2+(0.5*i*dt)*V(1))-dumc(l-1)*(-0.25*i*dt/dx**2))		
		end do
		psi_new(ind) = dumf(ind)
		do l=ind-1,1,-1
			psi_new(l) = dumf(l)-dumc(l)*psi_new(l+1) 
		end do

!       Uncomment for optional output
		! if (real(psi_new(1))>0.1) write(*,*) j, real(psi_new(1))
		
		psi = psi_new		

!		Here a boundary condition is needed or bad thing will happen!
		psi(1) = psi(2) ! psi(1) = 0 ! psi(1) = psi(ind)
	
!		Writeout the result (on modulo steps)
		if (modulo(j,20)==0) then		
			do l=1,ind,1
				if (isnan(real(psi_new(l)))) then
					write(*,*) 'Nan on iteration step: ',j
					close(20)
					stop	! stop program if something bad happens 
				end if
				write(20,*) x(l), sqrt(real(psi_new(l))**2+aimag(psi_new(l))**2), V(l)
			end do
			write(20,*)
			write(20,*)		
		end if	

	end do

!	Close output file
	close(20)

end program
