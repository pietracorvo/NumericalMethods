program wannier

!   NOTE: this code uses the LAPACK subroutine 'ssyev', which must be present at your system ...

	implicit none
	
	integer, parameter :: n=1500, steps=12
	real :: V(n)=0, H(n,n)=0, EW(n)=0, EV(n,n) = 0, work(3*n-1)

	real :: L=(19+2)*5, dx, E(steps)=(/ 0.0, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018, 0.02 /)
	integer :: i, j, k, info


!	Open file for saving the 20 EiegnValues and the first 20 EigenVectors
	open(unit=1, file="wannierEW.dat")
	open(unit=2, file="wannierEV.dat")

!	calculate steps
	dx = L/n

!	Loop over different fields E
	do k=1,steps,1	
!		Initialize potential	
		call pot(n,V,E(k))
!		Fill diagonal and offdiagonal elements of H	
		H(n,n) = 1/dx**2 + V(n)
		do j=1,n-1,1
			H(j,j) = 1/dx**2 + V(j)
			H(j+1,j) = -0.5/dx**2
			H(j,j+1) = -0.5/dx**2
		end do

!		Routine for calculation of EW and EV
		call ssyev('V', 'U', n, H, n, EW, work, 3*n-1, info)

		write(*,*) E(k), info
	
!		Writeout
		write(1,*) 'Electric field: ', E(k)
		do j=1,20,1
			write(1,*) EW(j)
		end do
		write(1,*)
		write(1,*)
		do j=1,n,1
			! just write the first 20 EV
			write(2,*) dx*j, V(j), (H(j,i)*H(j,i),i=1,20)
		end do
		write(2,*)
		write(2,*) 

	end do

	close(1)

end program


subroutine pot(n,V,E)
!	This subroutine builds the the vector with the potential V
	
	implicit none
	real, intent(in) :: E
	integer, intent(in) :: n
	real, intent(out) :: V(n)
	real :: dum(n)

	integer, parameter :: off=5	! off is number of potlengths vacuum per side

	real :: L=(19+2*off)*5, dx
	integer :: j, holecount, sw, loff

	holecount=0
	dum=0

	dx = L/n
	sw = int(n/(19+2*off))
	loff = int(0.5*(n-19*sw))	

!	Fill potential
	do j = loff,n,1
		if (holecount<10) then
			if (modulo(j-loff,2*sw)<=sw-1) then
				dum(j) = -1
			end if
			if (modulo(j-loff,2*sw)==sw-1) then
				holecount = holecount+1
			end if
			dum(j) = dum(j) + E*(j-loff)*dx 		
		end if
		!write(*,*) j, dum(j), holecount, E*(j-loff)*dx	
	end do
	V = dum

	return

end subroutine pot
