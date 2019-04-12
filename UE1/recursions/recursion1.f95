program calculate_neumann

implicit none

real :: neumann
real :: rho = 0.5
integer :: l = 5

write(*,*) 'Das ist der Wert:', neumann(rho,l)

end program


real function neumann(rho,l)

	implicit none
	real, intent(in) :: rho
	real :: tmp, tmp0, tmp1
	integer, intent(in) :: l
	integer :: i

	tmp0 = -cos(rho)/rho
	tmp1 = -sin(rho)/rho - cos(rho)/rho**2

	do i=2,1,l
		tmp = tmp1*real(2*i+1)/rho-tmp0			
		tmp0 = tmp1
		tmp1 = tmp
	end do

	neumann = tmp
	return

end function neumann



