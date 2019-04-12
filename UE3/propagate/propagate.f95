program propagate_psi

!	 compile with: gfortran propagate.f95 -lmathlib

	implicit none

	integer, parameter :: ind = 2**11
	real, parameter :: dx = 0.05

	complex :: psi(ind)
	real :: x(ind), k(ind), V(ind)=0, tmp(ind)	

	real :: dk = 2*3.142/(ind*dx), dt = 0.01, shift=10		
	integer :: cycles = 20000		
	real :: k0 = 1					
	complex :: i=(0,1)				
	integer :: j, l				


!	Open output files
	open(20,file='data_psi.dat')
	open(21,file='data_Fpsi.dat')	! that file writes FFT(psi)

! 	Build spatial vector x(j), potential V(x(j)) and k(j) ...
	do j=1,ind,1 
		x(j) = (j-int(0.5*ind))*dx
		k(j) = (j-int(0.5*ind))*dk
		V(j) = -exp(-3*abs(x(j)-shift)/(abs(x(j)-shift)+1))   ! Yukawa potential
	end do
!	... and flip k vector for FFT.
	tmp(1:ind/2) = k(1:ind/2)
	k(1:ind/2) = k(ind/2+1:ind)
	k(ind/2+1:ind) = tmp(1:ind/2)

!	Build starting psi for t=0
	psi = exp(-0.5*(x+shift)**2)*exp(i*k0*x)
	call writeout(x,psi,V,ind,20,0)


!	iteration:
!	----------
	do j=1,cycles,1

!		FFT
		call cfstft(11,psi)

!		Apply momentum operator
		psi = psi * exp(-0.5*i*dt*k**2)

!		Writeout FFT(psi) (on modulo steps)
		if (modulo(j,100)==0) call writeout(k,psi,V,ind,21,1)

!		Inverse FFT
		call cfstft(-11,psi)
		
!		Apply space operator
		psi = exp(-i*dt*V) * psi
	
!		Writeout psi (on modulo steps)
		if (modulo(j,100)==0) call writeout(x,psi,V,ind,20,0)		

!		Renormisation
		psi = psi / sqrt(sum(abs(psi)**2)*dx)

		write(*,*) 'Iteration step: ', j

	end do

!	Close output file
	close(20)
	close(21)

end program


subroutine writeout(x,psi,V,ind,channel,switch)

	implicit none

	complex, intent(in) :: 	psi(ind) 
	real, intent(in) 	:: 	x(ind), V(ind)
	integer, intent(in) ::	ind, channel, switch

	integer 			:: 	l 	
	real 				:: norm, dx = 0.05

	norm = 1 / sqrt(sum(abs(psi)**2)*dx)

	if(switch==1) then
!		For k sapce
		do l=ind/2+1,ind,1
			write(channel,*) x(l), norm*sqrt(abs(psi(l))**2)
		end do
		do l=1,ind/2,1
			write(channel,*) x(l), norm*sqrt(abs(psi(l))**2)
		end do
	else
!		For x space
		do l=1,ind,1
			write(channel,*) x(l), norm*sqrt(abs(psi(l))**2), V(l)
		end do	
	end if
	write(channel,*)
	write(channel,*)

end subroutine

