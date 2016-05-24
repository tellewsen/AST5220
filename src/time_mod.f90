module time_mod
    use healpix_types
    use params
    use spline_1D_mod
    use ode_solver
    implicit none


    integer(i4b)                           :: n1, n2,n3
    integer(i4b)                           :: n_t                ! Number of x-values
    integer(i4b)                           :: n_eta              ! Number of eta grid poins
    real(dp),    allocatable, dimension(:) :: x_t                ! Grid of relevant x-values
    real(dp),    allocatable, dimension(:) :: a_t                ! Grid of relevant a-values
    real(dp),    allocatable, dimension(:) :: eta_t              ! Grid of relevant eta-values
 
    real(dp),    allocatable, dimension(:) :: Omega_mx           !Relative densities
    real(dp),    allocatable, dimension(:) :: Omega_bx
    real(dp),    allocatable, dimension(:) :: Omega_rx 
    real(dp),    allocatable, dimension(:) :: Omega_nux 
    real(dp),    allocatable, dimension(:) :: Omega_lambdax 


    real(dp),    allocatable, dimension(:) :: H                  !Hubble constant as func of x
    real(dp)				   :: x_init,x_start_rec, x_end_rec
    real(dp)                               :: a_init,a_end


    real(dp),    allocatable, dimension(:) :: x_eta              ! Grid points for eta
    real(dp),    allocatable, dimension(:) :: a_eta              ! Grid points for eta
    real(dp),    allocatable, dimension(:) :: z_eta 		 ! Grid points for eta
    real(dp),    allocatable, dimension(:) :: eta, eta2          ! Eta and eta'' at each grid point

contains
  subroutine initialize_time_mod
    implicit none

    integer(i4b) :: i
    real(dp)     :: z_start_rec, z_end_rec, z_0, x_0
    real(dp)     :: dx, x_eta1, x_eta2,h1,eta_init,rho_crit0,rho_crit
    real(dp)     :: eps,hmin,yp1,ypn

    real(dp)                               :: rho_m0             !Matter density today
    real(dp) 				   :: rho_b0             !Baryon density today
    real(dp)				   :: rho_r0             !Radiation density today
    real(dp) 				   :: rho_nu0            !Neutrino density today
    real(dp)				   :: rho_lambda0        !Vacuum energy density today

    real(dp),    allocatable, dimension(:) :: rho_m              !Matter density
    real(dp),    allocatable, dimension(:) :: rho_b              !Baryon density
    real(dp),    allocatable, dimension(:) :: rho_r              !Radiation density
    real(dp),    allocatable, dimension(:) :: rho_nu             !Neutrino density
    real(dp),    allocatable, dimension(:) :: rho_lambda         !Vacuum energy density


    ! Define two epochs, 1) during and 2) after recombination.
    n1          = 300                       !Number of grid points before recombination
    n2          = 200                       !Number of grid points during recombination
    n3          = 300                       !Number of grid points after recombination
    n_t         = n1 + n2 +n3               !Total number of grid points

    z_start_rec = 1630.4d0                  ! Redshift of start of recombination
    z_end_rec   = 614.2d0                   ! Redshift of end of recombination
    z_0         = 0.d0                      ! Redshift today

    a_init      = 1.d-8                    ! Start value of a for eta evaluation
    a_end       = 1.d0

    x_start_rec = -log(1.d0 + z_start_rec)  ! x of start of recombination
    x_end_rec   = -log(1.d0 + z_end_rec)    ! x of end of recombination
    x_0         = 0.d0                      ! x today
    
    n_eta       = 1000                      ! Number of eta grid points (for spline)

    x_init      = log(a_init)               ! Start value of x for eta evaluation
    x_0         = 0.d0                      ! End value of x for eta evaluation
    eta_init    = a_init/(H_0*sqrt(Omega_r+Omega_nu))

    eps = 1.d-10
    hmin = 0.d0

    ! Task: Fill in x and a grids ( These will be used in later milestones)
    allocate(x_t(n_t))
    allocate(a_t(n_t))

    !Fill in x,a,z (rec) grids
    do i = 1,n1
        x_t(i)       = x_init + (i-1)*(x_start_rec-x_init)/(n1-1)
    end do
    do i = 1,n2 ! Fill interval during recombination
        x_t(n1+i)    = x_start_rec + i*(x_end_rec-x_start_rec)/(n2)
    end do
    do i = 1,n3 !Fill from end of recomb to today
        x_t(n1+n2+i) = x_end_rec + i*(x_0-x_end_rec)/(n3)
    end do
    !write(*,*) x_t !print x_t to terminal

    a_t = exp(x_t) !fill the a grid using the x grid
    !write(*,*) a_t !print a_t to terminal


    !Allocate and fill a,x, and z arrays
    allocate(a_eta(n_eta))
    allocate(x_eta(n_eta))
    allocate(z_eta(n_eta))

    x_eta(1) = x_init
    do i = 1,(n_eta-1)
        x_eta(i+1) = x_init + i*(x_0-x_init)/(n_eta-1)
    end do
    a_eta = exp(x_eta)
    z_eta = 1.d0/a_eta -1.d0
    
    !write(*,*) z_eta
    !write(*,*) size(z_eta)
    !print *, "x"
    !write(*,*) x_eta(1)
    !write(*,*) x_eta(-1)
    !print *, "a"
    !write(*,*) a_eta(1)
    !write(*,*) a_eta(-1)
    !print *, "z"
    !write(*,*) z_eta(1)
    !write(*,*) z_eta(-1)   

    !Calculate the various densities for each scale factor
    rho_crit0   = 3.d0*H_0**2.d0/(8.d0*pi*G_grav)
    rho_m0  	= Omega_m     *rho_crit0
    rho_b0  	= Omega_b     *rho_crit0
    rho_r0  	= Omega_r     *rho_crit0
    rho_nu0 	= Omega_nu    *rho_crit0
    rho_lambda0 = Omega_lambda*rho_crit0

    allocate(rho_m(n_eta))
    allocate(rho_b(n_eta))
    allocate(rho_r(n_eta))
    allocate(rho_nu(n_eta))
    allocate(rho_lambda(n_eta))

    allocate(Omega_mx(n_eta))
    allocate(Omega_bx(n_eta))
    allocate(Omega_rx(n_eta))
    allocate(Omega_nux(n_eta))
    allocate(Omega_lambdax(n_eta))
    allocate(H(n_eta))


    do i=1,n_eta
    H(i) = get_H(x_eta(i))
    Omega_mx(i) 	= Omega_m 	*(H_0/H(i))**2	*a_eta(i)**-3.d0
    Omega_bx(i) 	= Omega_b 	*(H_0/H(i))**2	*a_eta(i)**-3.d0
    Omega_rx(i) 	= Omega_r 	*(H_0/H(i))**2	*a_eta(i)**-4.d0
    Omega_nux(i) 	= Omega_nu 	*(H_0/H(i))**2	*a_eta(i)**-4.d0
    Omega_lambdax(i) 	= Omega_lambda	*(H_0/H(i))**2
    end do
    !End of density calculations



    allocate(eta(n_eta))
    eta(1) = eta_init !Start value of eta 

    h1 = abs(1.d-2*(a_eta(1)-a_eta(2))) !Defines the steplength
    !allocate(dydx(1))


    do i =2,n_eta
        eta(i) =eta(i-1)
        call odeint(eta(i:i),a_eta(i-1) ,a_eta(i), eps, h1, hmin, eta_derivs, bsstep, output) 
    end do
    !write(*,*) eta !check that eta gives reasonable values


    !Spline eta and place the second derivative of
    !this function in eta2
    allocate(eta2(n_eta))
    yp1 = 1.d30
    ypn = 1.d30
    call spline(a_eta, eta, yp1, ypn, eta2)
    

    allocate(eta_t(n_t))
    do i=1,n_t
        eta_t(i) = get_eta(x_t(i))
    end do

  end subroutine initialize_time_mod



  !Begin Stuff needed to make odeint work
  subroutine eta_derivs(a, eta, dydx) !Define the derivative d/da(eta)
    use healpix_types
        implicit none
        real(dp),               intent(in)  :: a
        real(dp), dimension(:), intent(in)  :: eta
        real(dp), dimension(:), intent(out) :: dydx
	real(dp)                            :: x
        x = log(a)
        dydx = c/(a*get_H_p(x))
  end subroutine eta_derivs

  subroutine output(x, y)
        use healpix_types
        implicit none
        real(dp),               intent(in)  :: x
        real(dp), dimension(:), intent(in)  :: y
  end subroutine output
  !End Stuff needed to make odeint work



  ! Task: Write a function that computes H at given x
  function get_H(x)
    implicit none
      real(dp), intent(in) :: x
      real(dp)             :: get_H
      real(dp) 		   :: a
      a = exp(x)
      get_H = H_0*sqrt((Omega_b+Omega_m)*a**-3.d0 + (Omega_r+Omega_nu)*a**-4.d0 + Omega_lambda)
  end function get_H

  ! Task: Write a function that computes H' = a*H  at given x
  function get_H_p(x)
    implicit none

      real(dp), intent(in) :: x
      real(dp)             :: get_H_p
      real(dp) 		   :: a
      a = exp(x)
      get_H_p = a*get_H(x)
  end function get_H_p

  ! Task: Write a function that computes dH'/dx at given x
  function get_dH_p(x)
      implicit none

      real(dp), intent(in) :: x
      real(dp)             :: get_dH_p
      real(dp)             :: exp2x,expx
      exp2x  = exp(2.d0*x)
      expx   = exp(x)
      get_dH_p = H_0/2.d0/sqrt((Omega_m+Omega_b)/expx+(Omega_r+Omega_nu)/exp2x &
                 + Omega_lambda*exp2x) * (-(Omega_m+Omega_b)/expx - 2.d0*(Omega_r+Omega_nu)/exp2x &
                 + 2.d0*Omega_lambda*exp2x)
  end function get_dH_p


  ! Task: Write a function that computes eta(x), using the previously precomputed splined function
  function get_eta(x_in)
    implicit none

      real(dp), intent(in) :: x_in
      real(dp)             :: get_eta
      real(dp) 		   :: a_in
      a_in = exp(x_in)
      get_eta = splint(a_eta, eta, eta2, a_in)
  end function get_eta

end module time_mod
