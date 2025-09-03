! solves the He atom in the Hartree--Fock approximation
! on an equidistant radial grid.
! note: the approximation for the triplet is not the Hartree--Fock approximation
program helium
  implicit none
  ! grid parameters
  real(8), parameter :: r_max = 20
  integer, parameter :: samples = 262144
  ! desired energy precision
  real(8), parameter :: prec = 1e-7_8
  ! lowest possible eigenenergy of Schroedinger equation searched for
  real(8), parameter :: epsilon_min = -20
  ! nuclear charge
  integer, parameter :: Z = 2

  ! radial grid for solving the Schroedinger and the Poisson equation
  real(8), dimension(:), allocatable :: r

  ! orbitals u(r,n) = r*psi_n(r) where n is the orbital index (1 or 2)
  ! for singlet n=1 for triplet n=1,2
  real(8), dimension(:,:), allocatable :: u
  ! induced electrostatic potential q(r,n) = r*phi_n(r)
  real(8), dimension(:,:), allocatable :: q


!====================================================================
!
!                      MAIN PROGRAM
!
!====================================================================

  call setup_regular_grid()
  call compute_singlet()
!  call compute_triplet()

contains
  subroutine compute_singlet()
    real(8) :: e, last_total_energy, total_energy
		integer :: j

    write (0,*) 'Helium singlet'
    ! start self consistent loop with hydrogenic orbitals
    q(:,:) = 0
    call se_solve(1, prec/128, q(:,1), e, u(:,1))
    call pe_solve(u(:,1), q(:,1))
    total_energy = 2*e + repulsion_energy(u(:,1), q(:,1))
    write (0,*) e, total_energy

    ! start self consistent loop
    last_total_energy = 1
    total_energy = 0
    do while (abs(last_total_energy - total_energy) > prec)
      last_total_energy = total_energy
      ! solve updated Schroedinger equation with higher precision
      call se_solve(1, prec/128, q(:,1), e, u(:,1))

      ! compute the total energy
      total_energy = 2*e - repulsion_energy(u(:,1), q(:,1))
      write (0,*) total_energy

      ! solve updated poisson equation
      call pe_solve(u(:,1), q(:,1))
    end do
!		Writeout psi
		open(unit=21,file='psi.dat')
		do j=1,samples,1
			write(21,*) r(j), u(j,1)*r(j)
		end do
		close(21)
  end subroutine

! allocate and setup the grid
  subroutine setup_regular_grid()
    implicit none
    integer :: i
    allocate(r(samples))
    allocate(u(samples,2), q(samples,2))
    do i= 1, samples
      r(i) = r_max*(i-1)/samples
    end do
  end subroutine

! calculate the energy by integrating the given q_1(r)/r*u_2^2(r)
  real(8) function repulsion_energy(u,q)
    real(8), dimension(:), intent(in) :: u  ! electron's wavefunction
    real(8), dimension(:), intent(in) :: q  ! other electron's induced potential
    real(8) :: e, step
    integer :: i, samples

    samples = size(r)
    ! note: assumes regular grid
    step = r(2) - r(1)
    e = 0
    do i = 2, samples
      e = e + step*q(i)/r(i)*u(i)**2
    end do

!   Here we return the computed value for the repulsion energy
    repulsion_energy = e
  end function

  ! solves the Poisson equation for the given wave function u(r) = r*psi(r)
  ! q''(r) = -u^2(r) / r
  ! returning the cummulant charge q(r) = r*phi(r)
  ! with the boundary conditions q(0) = 0, q(r_max) = 1.
  subroutine pe_solve(u, q)
    implicit none
    real(8), dimension(:), intent(in) :: u   ! electron's wavefunction
    real(8), dimension(:), intent(out) :: q  ! same electron's induced potential
    integer :: i, samples
    real(8) :: alpha, step

    samples = size(r)
    ! note: assumes regular grid
    step = r(2) - r(1)

!====================================================================
!
! Special boundary conditions!
! The potential can not be determined uniquely by solving the Poisson equation.
! A constant and linear term has to be fixed
!
! q(r) + alpha_0 + r * alpha
!
!
!====================================================================


!====================================================================
! Outwards integration from first known.
! Boundary condition q(r=0) = 0: no charge within radius 0
! This determines alpha_0
!====================================================================
    q(1) = 0
!====================================================================
! linear component alpha of q(r) set to zero. it is later determined from
! r_max boundary condition matching after integration
!====================================================================
    q(2) = 0
    ! treat i=2 specially: lim_r->0 u(r)^2/r = 0
    q(3) = 2*q(2) - q(1) - step**2/12 * ( &
      u(3)**2/r(3) + 10*u(2)**2/r(2) + 0 &
    )

!   Implement Numerow-Method for solving the Poisson equation
    ! outwards integration
    do i = 3, samples-1, 1
			q(i+1) = 2*q(i)-q(i-1)- step**2/12*((u(i+1)**2)/r(i+1)+10*(u(i)**2)/r(i)+(u(i-1)**2)/r(i-1))
    end do

!   match boundary condition at r_max: full charge 1 within r_max
!   to determine linear component alpha of q(r)
    alpha = (1 - q(samples)) / r(samples)

!   update q(r) according to integration constant
    do i = 1, samples
      q(i) = q(i) + alpha * r(i)
    end do

  end subroutine

  ! solves a hydrogen like Schroedinger equation in the given state with the
  ! main quantum number n.
  ! the grid has to be specified in r and the other electron's potential
  ! is given by means of the cummulant charge q(r)=r*phi(r) on the same grid.
  ! the resulting radial part of the wave function is written in u=r*psi(r).
  subroutine se_solve(n, prec, q, e, u)
    implicit none
    integer, intent(in) :: n                ! main quantum number
    real(8), intent(in) :: prec             ! desired energy precision
    real(8), dimension(:), intent(in) :: q  ! other electron's potential
    real(8), intent(out) :: e               ! eigenvalue found
    real(8), dimension(:), intent(out) :: u ! wave function on the grid

    real(8) :: e_min, e_max
    integer :: i, nodes

    ! only search for bound states within [epsilon_min,0]
    ! Initially we set e_min to -20 Ha
    e_min = epsilon_min
    e_max = 0

    do while (abs(e_max-e_min) > prec)
      e = (e_min+e_max)/2
      call se_integrate(q, e, u)

      ! count nodes
      nodes = 0
      do i = 1, size(r)-1
        if (u(i)*u(i+1) < 0) nodes = nodes + 1
      end do

      ! continue search in the upper or lower half
      if (nodes > n-1) then
        e_max = e
      else
        e_min = e
      endif
    end do

    call se_normalize(u)
  end subroutine

  subroutine se_integrate(q, e, u)
    implicit none
    real(8), dimension(:), intent(in) :: q  ! other electron's potential
    real(8), intent(in) :: e                ! eigenvalue tried
    real(8), dimension(:), intent(out) :: u ! wave function on the grid
    integer :: i, samples
    real(8) :: step
    real(8) :: gip1, gi, gim1   ! g_(i+1), g_i and g_(i-1) for Numerov

    samples = size(r)
    ! note: assumes regular grid
    step = r(2) - r(1)

!   Implement Numerow-Method for solving the SE equation
!    ! inward integration
!    ! from known boundary condition u(r) = r*exp(-(Z-1)*r) for large r

    u(samples) = r(samples)*exp(-(Z-1)*r(samples))
    u(samples-1) = r(samples-1)*exp(-(Z-1)*r(samples-1))

!   compute g_(i+1) and g_i for these values
    gip1 = 2*(e + Z/r(samples) - q(samples)/r(samples))
    gi = 2*(e + Z/r(samples-1) - q(samples-1)/r(samples-1))

!    integrate inward until i=3
    do i = samples-1, 3, -1

!   compute g_(i-1) 
      gim1 = 2*(e + Z/r(i-1) - q(i-1)/r(i-1))  

!   compute u(i-1)
      u(i-1) = (2*u(i)*(1-step**2*5/12*gi) - u(i+1)*(1+step**2/12*gip1)) / (1+step**2/12*gim1) 

!   update g_(i+1) and g_i for next iteration
       gip1 = gi
       gi = gim1

    end do

    ! treat i=2 specially since g(r) becomes infinite at r->0
    ! u(1) only used for node counting, so just take sign of numerator
    u(1) = sign( &
      1.0_8, &
      2*u(2)*(1 - 5*step**2/12*gi) - u(3)*(1 + step**2/12*gip1) &
    )
  end subroutine

  subroutine se_normalize(u)
    implicit none
    real(8), dimension(:), intent(inout) :: u ! wave function to normalize
    real(8) :: norm
    integer :: i, samples
    real(8) :: step

    samples = size(r)
    ! note: assumes regular grid
    step = r(2) - r(1)
    norm = 0

!   Compute Norm and normalize u(r)
!   u(1) does not contribute
!   all other points count fully
    do i = 2, samples
      norm = norm + u(i)**2 * step
    end do

!   normalize:
    u = u/sqrt(norm)


  end subroutine

!====================================================================
!
! BONUS
!
!====================================================================
  subroutine compute_triplet()
    real(8) :: e1, e2, last_total_energy, total_energy

    write (0,*) 'Helium triplet'
!====================================================================
! TODO (BONUS) : 
! start self consistent loop with hydrogenic orbitals
    q(:,:) = 0

    ! start self consistent loop

    ! solve updated Schroedinger equation

    ! compute the total energy

    ! solve updated poisson equation
!====================================================================
  end subroutine


end program
