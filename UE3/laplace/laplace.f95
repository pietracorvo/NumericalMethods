program laplace

	implicit none
	
	integer, parameter 	:: 	n=500
	integer 			:: 	mask(n,n)=0, start, ende, counter
	real 				:: 	U(n,n)=0, Uold(n,n)=1, V0=0, V1=1000, tmp
	integer 			:: 	i, j, k=0

!	Open output file
	open(unit=10, file='laplace.dat')

!	Start and end index for inner potential
	start=int(n/3)
	ende=int(2*n/3)

!	Fill mask for potential
	mask(1,1:n) = -1
	mask(n,1:n) = -1
	mask(1:n,1) = -1
	mask(1:n,n) = -1
	mask(start:ende,start:ende) = +1

!	Fill potetnial U with boundary conditions
	do i=1,n,1
		do j=1,n,1
			if (mask(i,j)==1) then
				U(i,j) = V1
			elseif (mask(i,j)==-1) then
				U(i,j) = V0
			end if				
		end do
	end do

!	Iterate (with exit for changes smaller 1E-00)
	do while (sum(U-Uold)>1)

		write(*,*) sum(U-Uold)
		Uold = U

!		Just iterate over chip without borders
		do i=2,n-1,1
			do j=2,n-1,1

!				If value at i,j is not a fixed potential value
				if (mask(i,j)==0) then
					U(i,j) = 0.2 * ( U(i,j) + U(i-1,j) + U(i+1,j) + U(i,j-1) + U(i,j+1) )	
				end if						

			end do
		end do

		k = k+1
		write(*,*) 'Iteration', k

	end do

	do i=1,n,1
		write(10,*) (U(i,j),j=1,n,1)
	end do

!	Close output file
	close(10)


end program
