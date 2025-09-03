
Program Ising

  external CalcE
  external CalcE_short

  integer NMax,NN                   !/* Maximal lattice size, # neighbors*/
  real*8 JJ, HH                     !/* Coupling constant, magnetic field */
  parameter (NMax=50,NN=4)
  parameter (JJ=1.d0,HH=0.d0)


  integer mag(Nmax,Nmax)            !/* 2D Ising Lattice */
  integer i, j, k, a,b,iT                   !/* Loop counters */
  integer s,d                       !/* Lattice spin variables */  
  real*8 Energy                     !/* Total lattice energy */
  integer Inn(NN)                   !/*  Nearest neighbor array I  */
  parameter (Inn=(/1,-1,0,0/))      !/* to CHANGE!!!
  integer Jnn(NN)                   !/*  Nearest neighbor array J  */
  parameter (Jnn=(/0,0,1,-1/))      !/* to CHANGE!!!
  integer Inew, Jnew                !/* Nearest neighbot indices */ 
  real*8 Etemp, deltaE              !/* Temp energy variables for MC moves */ 
  integer accept                    !/* Number of accepted moves */   
  integer move                      !/* Number of moves total */   

  real*8 T                          !/* temperature (in units of J/k_B) */
  real*8 beta
  integer sweeps                    !/* number of measurement sweeps */
  integer warm                      !/* number of warm-up sweeps */
  integer L                         !/* lattice dimension */
  real*8 r, dum1, dum2


  real*8 :: M_abs=0, M_square=0, M_four=0, Binder=0 

open(unit=22, file='ergebnis.dat')

do iT=0,30

  M_abs=0
  M_square=0
  M_four=0
  Binder=0
  Energy=0 

  do i=1,Nmax
    do j=1,Nmax
      mag(i,j)=0
    end do
  end do

 ! /***************************
 !  * Initialization          * 
 !  ***************************/

 ! write (*,*) "Temperatur ?"
  !read(*,*) T
  T=1+iT*0.1
  beta = 1.d0/T
  write (*,*) 'Temperatur: ', T

  !write (*,*) "Gittergroesse (L x L) ?"
  !read(*,*) L
  L = 20
  write (*,*) 'Gittergroesse (L x L), L = ', L

  !write (*,*) "# sweeps ?"
  !read(*,*) sweeps
  sweeps = 10E4
  write (*,*) '# sweeps ', sweeps
  
  !write (*,*) "# warm up sweeps?"
  !read(*,*) warm
  warm = 0.1*sweeps
  write (*,*) "# warm up sweeps ", warm


  call random_init(mag,L,Nmax)  
  Energy = CalcE(mag,L,Jnn,Inn,NN,JJ,HH,Nmax)

  !call  outputmag(mag,L,Nmax)
 
! /***************************
! /* warum up sweeps */
! /***************************


 do k=1,warm,1
    
    do i=1,L,1
      do j=1,L,1

        deltaE = CalcE_short(i,j,mag,L,Jnn,Inn,NN,JJ,HH,Nmax)
        r = exp(-beta*deltaE)

        if (rand().lt.r/(r+1)) then
          mag(i,j) = mag(i,j)*(-1)
          Energy = Energy + deltaE
        end if

      end do
    end do

  end do  

 
! /***************************
! /* end warum up sweeps */ 
! /***************************

! /***************************
! /*  sweeps */
! /***************************

 do k=1,sweeps,1
    
    do i=1,L,1
      do j=1,L,1

        deltaE = CalcE_short(i,j,mag,L,Jnn,Inn,NN,JJ,HH,Nmax)
        r = exp(-beta*deltaE)

        if (rand().lt.r/(r+1.d0)) then
          mag(i,j) = mag(i,j)*(-1)
          Energy = Energy + deltaE
        end if

      end do
    end do

    dum=0
    do a=1,L
      do b=1,L  
      !dum = mag(a,b)+dum
      end do
    end do
  

    M_abs = M_abs + real(abs(sum(mag)))/(sweeps*L**2)
  !  write(*,*) sum(mag)
    M_square = M_square + (real(sum(mag))/L**2)**2/sweeps
    M_four = M_four + (real(sum(mag))/L**2)**4/sweeps
    Binder = 1 - M_four/(3*M_square**2)

  end do 


! /***************************
! /*  end sweeps */ 
! /***************************


  !call  outputmag(mag,L,Nmax)
  write(*,*) 'Average energy',Energy
  write(*,*) M_abs, M_square, M_four, Binder

 
  write(22,*) L, T, M_abs, M_square, M_four, Binder
  

end do
 
 close(22)

end Program Ising







subroutine random_init(mag,L,Nmax)
  integer Nmax
  integer mag(Nmax,Nmax), L
 
  integer i,j
  
  do i=1,L
     do j=1,L
       	if (rand().gt.0.5d0) then
	   mag(i,j) = 1
        else
           mag(i,j) = -1
        endif
     enddo
  enddo

  return
end subroutine random_init



subroutine outputmag(mag,L,Nmax)
  integer Nmax
  integer mag(Nmax,Nmax), L
 
  integer i,j
  write (*,*) 'Spin snapshot'

  do i=1,L
     do j=1,L
        if(mag(i,j).eq. 1) then
	    write(*,'(a)',ADVANCE='NO') '-> '
        else
           write(*,'(a)',ADVANCE='NO') '<- '
        endif
      
     enddo
     !     Newline to complete
     write (*,*) 
     write (*,*)
  enddo
  return 
end subroutine outputmag




real function CalcE(mag,L,Jnn,Inn,NN,JJ,HH,Nmax)
  integer Nmax
  integer mag(Nmax,Nmax),Jnn(NN),Inn(NN),L,NN
  
  real*8 JJ, HH
  integer i,j,k,Inew,Jnew
  real*8 Energy

  !/* Determine the initial energy */
  Energy = 0.0
  
  do i=1,L
     do j=1,L
        
	!/* Loop over nearest neighbors */
	do k=1, NN  
      Inew = i + Inn(k)       
      Jnew = j + Jnn(k)
    
	    !/* Check periodic boundary conditions */
	    if (Inew .le. 0) then
	      Inew = L
      else 
        if(Inew .gt. L) then 
          Inew = 1
        endif
      endif    
      if (Jnew .le. 0) then
	      Jnew = L
      else 
        if(Jnew .gt. L) then 
          Jnew = 1
        endif
      endif
	    
	    !/* Update the energy */
	    Energy = Energy-JJ * mag(i,j) * mag(Inew,Jnew)
        enddo
	!/*Calculate the contribution from the field H */
	Energy = Energy - 2.d0*HH*mag(i,j);
     enddo
  enddo

   !/* Account for double counting */
   Energy = Energy/2.d0;

   CalcE=Energy
   return
 end function CalcE


real function CalcE_short(i_index,j_index,mag,L,Jnn,Inn,NN,JJ,HH,Nmax)

  integer Nmax, i_index, j_index
  integer mag(Nmax,Nmax),Jnn(NN),Inn(NN),L,NN

  real*8 JJ, HH
  integer k,Inew,Jnew
  real*8 dEnergy

  !/* Determine the initial energy */
  dEnergy = 0.0
        
	!/* Loop over nearest neighbors */
	do k=1, NN  
      Inew = i_index + Inn(k)       
      Jnew = j_index + Jnn(k)
    
	    !/* Check periodic boundary conditions */
	    if (Inew .le. 0) then
	      Inew = L
      else 
        if(Inew .gt. L) then 
          Inew = 1
        endif
      endif    
      if (Jnew .le. 0) then
	      Jnew = L
      else 
        if(Jnew .gt. L) then 
          Jnew = 1
        endif
      endif
	    
	    !/* Update the energy */
	    dEnergy = dEnergy + 2 * JJ * mag(i_index,j_index) * mag(Inew,Jnew) 
  enddo
	!/*Calculate the contribution from the field H */
	dEnergy = dEnergy + 2.d0*HH*mag(i_index,j_index);

  CalcE_short = dEnergy
  return

end function CalcE_short


