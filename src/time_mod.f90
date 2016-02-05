module time_mod
  use healpix_types
  use params
  use spline_1D_mod
  use ode_solver
  implicit none

  integer(i4b)                           :: n_t                ! Number of x-values
  real(dp),    allocatable, dimension(:) :: x_t                ! Grid of relevant x-values
  real(dp),    allocatable, dimension(:) :: a_t                ! Grid of relevant a-values

  integer(i4b)                           :: n_eta              ! Number of eta grid poins
  real(dp),    allocatable, dimension(:) :: x_eta              ! Grid points for eta
  real(dp),    allocatable, dimension(:) :: eta, eta2          ! Eta and eta'' at each grid point


contains

  subroutine initialize_time_mod
    implicit none

    integer(i4b) :: i, n, n1, n2
    real(dp)     :: z_start_rec, z_end_rec, z_0, x_start_rec, x_end_rec, x_0, dx, x_eta1, x_eta2, a_init

    ! Define two epochs, 1) during and 2) after recombination.
    n1          = 200                       ! Number of grid points during recombination
    n2          = 300                       ! Number of grid points after recombination
    n_t         = n1 + n2                   ! Total number of grid points

    z_start_rec = 1630.4d0                  ! Redshift of start of recombination
    z_end_rec   = 614.2d0                   ! Redshift of end of recombination
    z_0         = 0.d0                      ! Redshift today

    x_start_rec = -log(1.d0 + z_start_rec)  ! x of start of recombination
    x_end_rec   = -log(1.d0 + z_end_rec)    ! x of end of recombination
    x_0         = 0.d0                      ! x today
    
    n_eta       = 1000                      ! Number of eta grid points (for spline)
    a_init      = 1.d-10                    ! Start value of a for eta evaluation
    x_eta1      = log(a_init)               ! Start value of x for eta evaluation
    x_eta2      = 0.d0                      ! End value of x for eta evaluation

    ! Task: Fill in x and a grids
    allocate(x_t(n_t))

    do i = 0,n1 ! Fill interval during recombination
        x_t(i+1)= x_start_rec + i*(x_end_rec-x_start_rec)/n1
    end do
    do i = 2,n2 !Fill from end of recomb to today
        x_t(n1+i) = x_end_rec + i*(x_0-x_end_rec)/n2
    end do

    !write(*,*) x_t !print x_t to terminal



    allocate(a_t(n_t))
    a_t = exp(x_t) !fill a grid using the x grid

    !write(*,*) a_t !print a_t to terminal




    ! Task: 1) Compute the conformal time at each eta time step
    !       2) Spline the resulting function, using the provided "spline" routine in spline_1D_mod.f90
    allocate(x_eta(n_eta+1))
    do i = 1,n_eta+1
        x_eta(i) = x_eta1 + (x_eta2-x_eta1)*(i-1)/n_eta
    end do
    !write(*,*) x_eta !print x_eta to terminal





    allocate(eta(n_eta))
    do i=1,n_eta
	eta(i+1) = eta(i)
!        call odeint(ystart, x1, x2, eps, h1, hmin, derivs, rkqs, output)
!        call derivs(x,y,dydx)
!        call odeint(eta(i+1), x_eta(i+1), x_eta(i+2), 1d-10, 1d30, hmin, derivs, rkqs, output)
    end do




    allocate(eta2(n_eta))
    
  end subroutine initialize_time_mod



  !Stuff needed to make odeint work
  subroutine derivs(x, y, dydx)
       use healpix_types
         implicit none
         real(dp),               intent(in)  :: x
         real(dp), dimension(:), intent(in)  :: y
         real(dp), dimension(:), intent(out) :: dydx
         !dydx = 0
  end subroutine derivs

  subroutine output(x, y)
         use healpix_types
         implicit none
         real(dp),               intent(in)  :: x
         real(dp), dimension(:), intent(in)  :: y
  end subroutine output


  ! Task: Write a function that computes H at given x
  function get_H(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H
    real(dp) 		 :: a
    a = exp(x)

    get_H = H_0*sqrt((Omega_b+Omega_m)*a**-3 + (Omega_r+Omega_nu)*a**-4 + Omega_lambda)
  end function get_H

  ! Task: Write a function that computes H' = a*H  at given x
  function get_H_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H_p
    real(dp) 		 :: a
    a = exp(x)

    get_H_p = a*get_H(x)
  end function get_H_p

  ! Task: Write a function that computes dH'/dx at given x
  function get_dH_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dH_p
    get_dH_p = 0
  end function get_dH_p

  ! Task: Write a function that computes eta(x), using the previously precomputed splined function
  function get_eta(x_in)
    implicit none

    real(dp), intent(in) :: x_in
    real(dp)             :: get_eta
    get_eta = 0
  end function get_eta

end module time_mod
