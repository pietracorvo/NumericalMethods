program calculation_binomial

	implicit none

	integer :: n, k
	integer :: binomialkoeff

	write (*,*) 'Bitte eingeben n: '
	read (*,*) n

	write (*,*) 'Bitte eingeben k: '
	read (*,*) k

	write (*,*) binomialkoeff(n,k)

end program

integer function binomialkoeff(n,k)

	implicit none 

	integer, intent(in) :: n, k
	integer :: fact

	binomialkoeff = fact(n)/(fact(k)*fact(n-k))
	
	return
end function

recursive function fact(m) result(x)
	implicit none
	integer, intent(in) :: m
	integer :: x

	if (m<=0) then
		x = 1
	else	
		x = m * fact(m-1)
	end if

	return
end function
