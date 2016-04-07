module evolution_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use rec_mod
  implicit none

  ! Accuracy parameters
  real(dp),     parameter, private :: a_init   = 1.d-8
  real(dp),     parameter, private :: k_min    = 0.1d0 * H_0 / c
  real(dp),     parameter, private :: k_max    = 1.d3  * H_0 / c
  integer(i4b), parameter          :: n_k      = 100
  integer(i4b), parameter, private :: lmax_int = 6

  ! Perturbation quantities
  real(dp), allocatable, dimension(:,:,:) :: Theta
  real(dp), allocatable, dimension(:,:)   :: delta
  real(dp), allocatable, dimension(:,:)   :: delta_b
  real(dp), allocatable, dimension(:,:)   :: Phi
  real(dp), allocatable, dimension(:,:)   :: Psi
  real(dp), allocatable, dimension(:,:)   :: v
  real(dp), allocatable, dimension(:,:)   :: v_b
  real(dp), allocatable, dimension(:,:)   :: dPhi
  real(dp), allocatable, dimension(:,:)   :: dPsi
  real(dp), allocatable, dimension(:,:)   :: dv_b
  real(dp), allocatable, dimension(:,:,:) :: dTheta

  ! Fourier mode list
  real(dp), allocatable, dimension(:) :: ks

  ! Book-keeping variables
  real(dp),     private :: k_current
  integer(i4b), private :: npar = 6+lmax_int

  !With or without polarization
  logical(lgt) :: polarize = False
contains


  ! NB!!! New routine for 4th milestone only; disregard until then!!!
  subroutine get_hires_source_function(k, x, S)
    implicit none

    real(dp), pointer, dimension(:),   intent(out) :: k, x
    real(dp), pointer, dimension(:,:), intent(out) :: S

    integer(i4b) :: i, j
    real(dp)     :: g, dg, ddg, tau, dt, ddt, H_p, dH_p, ddHH_p, Pi, dPi, ddPi
    real(dp), allocatable, dimension(:,:) :: S_lores

    ! Task: Output a pre-computed 2D array (over k and x) for the 
    !       source function, S(k,x). Remember to set up (and allocate) output 
    !       k and x arrays too. 
    !
    ! Substeps:
    !   1) First compute the source function over the existing k and x
    !      grids
    !   2) Then spline this function with a 2D spline
    !   3) Finally, resample the source function on a high-resolution uniform
    !      5000 x 5000 grid and return this, together with corresponding
    !      high-resolution k and x arrays

  end subroutine get_hires_source_function




  ! Routine for initializing and solving the Boltzmann and Einstein equations
  subroutine initialize_perturbation_eqns
    implicit none

    integer(i4b) :: l, i

    ! Task: Initialize k-grid, ks; quadratic between k_min and k_max

    ! Allocate arrays for perturbation quantities
    allocate(Theta(0:n_t, 0:lmax_int, n_k))
    allocate(delta(0:n_t, n_k))
    allocate(delta_b(0:n_t, n_k))
    allocate(v(0:n_t, n_k))
    allocate(v_b(0:n_t, n_k))
    allocate(Phi(0:n_t, n_k))
    allocate(Psi(0:n_t, n_k))
    allocate(dPhi(0:n_t, n_k))
    allocate(dPsi(0:n_t, n_k))
    allocate(dv_b(0:n_t, n_k))
    allocate(dTheta(0:n_t, 0:lmax_int, n_k))

    ! Task: Set up initial conditions for the Boltzmann and Einstein equations
    Phi(0,:)     = 1.d0
    delta(0,:)   = 1.5d0*Phi(0,:)
    delta_b(0,:) = delta(0,:)
       
    do i = 1, n_k
        v(0,i)       = c*ks(i)/(2.d0*get_H_p(x_rec(i))*Phi(0,:)
        v_b(0,i)     = v(0,i)
        Theta(0,0,i) = 0.5d0*Phi(0,:)
        Theta(0,1,i) = -c*ks(i)/(6.d0*get_H_p(x_rec(i)))*Phi(0,:)

        if(polarize==True)
            Theta(0,2,i) = -8.d0*c*ks(i)/(15.d0*get_H_p(x_rec(i))*get_dtau(x_rec(i)))*Theta(0,1,i) !with polarization
        else
            Theta(0,2,i) = -20.d0*c*ks(i)/(45.d0*get_H_p(x_rec(i))*get_dtau(x_rec(i)))*Theta(0,1,i) !without polarization
        end if
        do l = 3, lmax_int
            Theta(0,l,i) = -l/(2.d0*l+1.d0)*c*ks(i)/(get_H_p(x_rec(i))*get_dtau(x_rec(i)))*Theta(0,l-1,i)
        end do
    end do

  end subroutine initialize_perturbation_eqns

  subroutine integrate_perturbation_eqns
    implicit none

    integer(i4b) :: i, j, k, l
    real(dp)     :: x1, x2, x_init
    real(dp)     :: eps, hmin, h1, x_tc, H_p, dt, t1, t2

    real(dp), allocatable, dimension(:) :: y, y_tight_coupling, dydx

    x_init = log(a_init)
    eps    = 1.d-8
    hmin   = 0.d0

    allocate(y(npar))
    allocate(dydx(npar))
    allocate(y_tight_coupling(7))

    ! Propagate each k-mode independently
    do k = 1, n_k

       k_current = ks(k)  ! Store k_current as a global module variable
       h1        = 1.d-5

       ! Initialize equation set for tight coupling
       y_tight_coupling(1) = delta(0,k)
       y_tight_coupling(2) = delta_b(0,k)
       y_tight_coupling(3) = v(0,k)
       y_tight_coupling(4) = v_b(0,k)
       y_tight_coupling(5) = Phi(0,k)
       y_tight_coupling(6) = Theta(0,0,k)
       y_tight_coupling(7) = Theta(0,1,k)
       
       ! Find the time to which tight coupling is assumed, 
       ! and integrate equations to that time
       x_tc = get_tight_coupling_time(k_current)

       ! Task: Integrate from x_init until the end of tight coupling, using
       !       the tight coupling equations

       q = (-((1.d0-2.d0*R)*get-dtau(x_t(i)) + (1+R)*get_ddtau(x_t(i)))*(3.d0*theta(i,1,k_current)) &
           +v_b(i,k_current)) -c*k_current/get_H_p(x_t(i))*Psi(i,k) +(1.d0-get_dH_p(x_t(i)) &
           /get_H_p(x_t(i)))*c*k_current/get_H_p(x_t(i))*(-Theta(i,0,k)+2.d0*Theta(i,2,k))-c*k_current &
           /get_H_p(x_t(i))*dTheta(i,0,k))/((1.d0+R)*get_dtau(x_t(i))+get_dH_p(x_t(i))/get_H_p(x_t(i)) -1.d0)

       dv_b(i,k) = 1.d0/(1.d0+R)*(-v_b(i,k)-c*k_current/get_H_p(x_t(i))*Psi(i,k) &
                   +R*(q+c*k_current/get_H_p(x_t(i))*(-Theta(i,0,k)+2.d0*Theta(i,2,k)) &
                   -c*k_current/get_H_p(x_t(i))*Psi(i,k)))
       d_Theta(i,1,k) = 1.d0/3.d0*(q-dv_b(i,k))
       ! Task: Set up variables for integration from the end of tight coupling 
       ! until today
       y(1:7) = 
       y(8)   = 
       do l = 3, lmax_int
          y(6+l) = 
       end do
!
       
       do i = 1, n_t
          ! Task: Integrate equations from tight coupling to today

          ! Task: Store variables at time step i in global variables
          delta(i,k)   = 
          delta_b(i,k) = 
          v(i,k)       = 
          v_b(i,k)     = 
          Phi(i,k)     = 
          do l = 0, lmax_int
             Theta(i,l,k) = 
          end do
          Psi(i,k)     =  - Phi(i,k) - 12.d0*H_0**2/(c*k_current*a_t(i))**2*Omega_r*Theta(i,2,k)

          ! Task: Store derivatives that are required for C_l estimation
          dPhi(i,k)     = Psi(i,k) -c**2*k_current**2/(3.d0*get_H_p(x_t(i))**2)*Phi(i,k) +H_0**2/(2.d0*get_H_p(x_t(i))) &
                          *(Omega_m/a_t(i)*delta(i,k) +Omega_b/a_t(i)*delta_b(i,k_current) + 4.d0*Omega_r/a_t(i)**2 &
                          *Theta(i,0,k_current) )
          dv_b(i,k)     = -v_b(i,k) -c*k_current/get_H_p(x_t(i))*Psi(i,k) +get_dtau(x_t(i))*R*(3.d0*Theta(i,1,k)+ v_b(i,k))
          dTheta(i,:,k) = 
          dPsi(i,k)     = 
       end do

    end do

    deallocate(y_tight_coupling)
    deallocate(y)
    deallocate(dydx)

  end subroutine integrate_perturbation_eqns


  ! Task: Complete the following routine, such that it returns the time at which
  !       tight coupling ends. In this project, we define this as either when
  !       dtau < 10 or c*k/(H_p*dt) > 0.1 or x > x(start of recombination)
  function get_tight_coupling_time(k)
    implicit none

    real(dp), intent(in)  :: k
    real(dp)              :: get_tight_coupling_time

  end function get_tight_coupling_time

end module evolution_mod
