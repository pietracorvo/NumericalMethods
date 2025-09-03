Program nmr

	implicit none

  external CalcE
  external CalcE_short

	real CalcE
	real CalcE_short

  integer Nimax,Njmax, dum1, dum2		!dimension of the image
  parameter (Nimax=254,Njmax=333)
	integer MR(Nimax,Njmax), mag(Nimax,Njmax), mag_old(Nimax,Njmax) !/* Picture/sigma array

  integer NN                   			!/* Maximal lattice size, # neighbors*/               
  parameter (NN=4)
	integer Inn(NN)                   !/*  Nearest neighbor array I  */
  parameter (Inn=(/1,-1,0,0/))      !/* to CHANGE!!!
  integer Jnn(NN)                   !/*  Nearest neighbor array J  */
  parameter (Jnn=(/0,0,1,-1/))      !/* to CHANGE!!!

	real*8 JJ
  parameter (JJ=1.d0)								!/* Coupling constant, magnetic field */

	integer MW(5)											!mean values and standard deviation
!                  BG WM  GM  CSF  SB
  parameter (MW=(/30,426,602,1223,167/))
  integer SIGMA(5)
  parameter (SIGMA=(/30,59,102,307,69/))

  integer i, j, k, l					      !/* Loop counters */
	integer newmag

  real*8 Energy,deltaE, Energy_old, dEnergy             !/* Total lattice energy */
  
  real*8 T, Tstart                          !/* temperature (in units of J/k_B) */
  integer sweeps                    !/* number of measurement sweeps */
  real*8 r, error, ran


  
open(unit=22, file='SimMRimage.dat')
open(unit=32, file='NMR.dat')

  !read in image
	do i=1,Nimax
     do j=1,Njmax
       read(22,'(3I5)') dum1, dum2, MR(i,j)
     enddo
     read(22,*)
  enddo

  Energy=0 

 ! /***************************
 !  * Initialization          * 
 !  ***************************/

  ! write (*,*) "Temperatur ?"
  !read(*,*) T
  Tstart = 6
	T = Tstart
  write (*,*) 'Temperatur: ', T

  !write (*,*) "# sweeps ?"
  !read(*,*) sweeps
  sweeps = 1E2
  write (*,*) '# sweeps ', sweeps
 
  call random_init(mag,MR,MW,Nimax,Njmax)

	do l = 1,118

  Energy_old = CalcE(mag,MR,Nimax,Njmax,Jnn,Inn,NN,JJ,MW,SIGMA)
 

! /***************************
! /*  sweeps */
! /***************************

  do k=1,sweeps,1
    
    do i=1,Nimax,1
      do j=1,Njmax,1

				newmag=INT(rand()*5.d0)+1
        deltaE = CalcE_short(i,j,newmag,mag,MR,Nimax,Njmax,Jnn,Inn,NN,JJ,MW,SIGMA)
				r = exp(-deltaE/T)
				!write(*,*) deltaE, r 
        if (rand() .lt. (r/(r+1.d0))) then
          mag(i,j) = newmag
          Energy = Energy + deltaE
        end if

      end do
    end do

  end do 


! /***************************
! /*  end sweeps */ 
! /***************************

	error = sum(abs(mag_old-mag))
	T = Tstart - 0.05*l

	write(*,*) T, error

	end do


  call  outputmag(mag,Nimax,Njmax)
  close(32)

end Program nmr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine random_init(mag,MR,MW,Nimax,Njmax) 

	implicit none
  integer Nimax,Njmax
  integer mag(Nimax,Njmax), MR(Nimax,Njmax), MW(5), vgl(4)

	integer i,j	

	vgl = 0.5*(MW(1:4)+MW(2:5))
  
  do i=1,Nimax
    do j=1,Njmax
			if ((MR(i,j) .ge. 0) .and. (MR(i,j) .lt. vgl(1))) then
	    	mag(i,j) = 1
      elseif ((MR(i,j) .ge. vgl(1)) .and. (MR(i,j) .lt. vgl(2))) then
      	mag(i,j) = 2
			elseif ((MR(i,j) .ge. vgl(2)) .and. (MR(i,j) .lt. vgl(3))) then
				mag(i,j) = 3
			elseif ((MR(i,j) .ge. vgl(3)) .and. (MR(i,j) .lt. vgl(4))) then
				mag(i,j) = 4
			elseif (MR(i,j) .ge. vgl(4)) then
				mag(i,j) = 5
			else
				write(*,*) "Uuuuups:  ", MR(i,j)
      endif
    enddo
  enddo

  return
end subroutine random_init


subroutine outputmag(mag,Nimax,Njmax)
	implicit none
  integer Nimax, Njmax
  integer mag(Nimax,Njmax)
 
  integer i,j

  do i=1,Nimax
    do j=1,Njmax
        write(32,'(3I5)') i,j,mag(i,j)
    enddo
    write (32,*)
  enddo

  return 
end subroutine outputmag


real function CalcE(mag,MR,Nimax,Njmax,Jnn,Inn,NN,JJ,MW,SIGMA)

	implicit none

  integer Nimax,Njmax, NN
  integer mag(Nimax,Njmax),MR(Nimax,Njmax),Jnn(NN),Inn(NN),MW(5),SIGMA(5)
  
  real*8 JJ
  integer i,j,k,Inew,Jnew
  real*8 Energy

  !/* Determine the initial energy */
  Energy = 0.0
  
  do i=1,Nimax
    do j=1,Njmax
        
			!/* Loop over nearest neighbors */
			do k=1, NN  
      	Inew = i + Inn(k)       
      	Jnew = j + Jnn(k)
    
	    	!/* Check periodic boundary conditions */
	    	if (Inew .le. 0) then
          Inew = Nimax
        else 
        	if(Inew .gt. Nimax ) then 
            Inew = 1
          endif
        endif
        if (Jnew .le. 0) then
        	Jnew = Njmax
        else 
        	if(Jnew .gt. Njmax) then 
            Jnew = 1
          endif
        endif
	    
	    	!/* Update the energy */
			  if (mag(i,j).eq.mag(Inew,Jnew)) then
					Energy = Energy - JJ
				endif 
      enddo
		!/*Calculate Energy */
		Energy = Energy + 0.5*(MR(i,j) - MW(mag(i,j)))**2/(SIGMA(mag(i,j)))**2 + 2*log(real(SIGMA(mag(i,j)))) 

    enddo
  enddo

  CalcE=Energy
  return
end function CalcE


real function CalcE_short(i,j,newmag,mag,MR,Nimax,Njmax,Jnn,Inn,NN,JJ,MW,SIGMA)

	implicit none

  integer Nimax,Njmax, NN
  integer mag(Nimax,Njmax),MR(Nimax,Njmax),Jnn(NN),Inn(NN),MW(5),SIGMA(5)
  
  real*8 JJ
  integer i,j,k,Inew,Jnew,newmag
  real*8 dEnergy

  !/* Determine the initial energy */
  dEnergy = 0.0
        
	!/* Loop over nearest neighbors */
	do k=1, NN  
    Inew = i + Inn(k)       
    Jnew = j + Jnn(k)
    
	  !/* Check periodic boundary conditions */
	  if (Inew .le. 0) then
      Inew = Nimax
    else 
     	if(Inew .gt. Nimax ) then 
        Inew = 1
      endif
    endif
    if (Jnew .le. 0) then
     	Jnew = Njmax
    else 
     	if(Jnew .gt. Njmax) then 
        Jnew = 1
      endif
    endif
	    
	  !/* Update the energy */
		if ((newmag .eq. mag(Inew,Jnew)) .and. (mag(i,j) .ne. mag(Inew,Jnew))) then
			dEnergy = dEnergy - JJ
		elseif ((newmag .ne. mag(Inew,Jnew)) .and. (mag(i,j) .eq. mag(Inew,Jnew))) then 
			dEnergy = dEnergy + JJ
		endif 
  enddo

  !/*Calculate Energy */
	dEnergy = dEnergy & 
						+ 0.5*((MR(i,j)-MW(newmag))/SIGMA(newmag))**2 &
						- 0.5*((MR(i,j)-MW(mag(i,j)))/SIGMA(mag(i,j)))**2 &
						+ log(real(SIGMA(newmag))) &
						- log(real(SIGMA(mag(i,j))))

  CalcE_short = dEnergy
  return

end function CalcE_short

