module rec_mod
    use healpix_types
    use params
    use time_mod
    use ode_solver
    use spline_1D_mod
    implicit none

    integer(i4b)                                 :: n                        ! Number of grid points
    real(dp), allocatable, dimension(:)          :: x_rec,a_rec,z_rec        ! Grid
    real(dp), allocatable, dimension(:)          :: tau, dtau, ddtau         ! Splined tau and derivatives
    real(dp), allocatable, dimension(:)          :: d4tau,logd4tau
    real(dp), allocatable, dimension(:)          :: n_e, n_e2,logn_e,logn_e2 ! Splined (log of) electron density, n_e
    real(dp), allocatable, dimension(:)          :: g, g2, g22               ! Splined visibility function
    real(dp), allocatable, dimension(:)          :: H_rec,X_e  	             ! Variables for H and X_e
    integer(i4b)                                 :: j,k                        !Used for choosing right indexes
    real(dp), allocatable, dimension(:)          :: x_test,n_etest,z_test,a_test  !Used for testing the splines
    real(dp), allocatable, dimension(:)          :: tau_test,dtau_test,ddtau_test  !Used for testing the splines
    real(dp)                                     :: x_0
    real(dp)                                     :: x_test_start
    real(dp)                                     :: x_test_end
contains

    subroutine initialize_rec_mod
    implicit none
    
    integer(i4b) :: i,n1,n2
    real(dp)     :: saha_limit, y, T_b, n_b, dydx, xmin, xmax, dx, f, n_e0, X_e0, &
		    X_econst, phi2,alpha2,beta,beta2,n1s,lambda_alpha,C_r
    real(dp)     :: eps,hmin,yp1,ypn,h1,h2
    real(dp)     :: z_start_rec, z_end_rec, z_0, x_start_rec, x_end_rec
    logical(lgt) :: use_saha


    x_test_start = -7.0d0
    x_test_end   = -1.0d0
    saha_limit   = 0.99d0       ! Switch from Saha to Peebles when X_e < 0.99

    !ODE int variables
    eps  = 1.d-10
    hmin = 0.d0
    yp1  = 1.d30
    ypn  = 1.d30

    !Grid sizes 
    n1          = 200                       ! Number of grid points during recombination
    n2          = 300                       ! Number of grid points after recombination
    n           = n1 + n2                   ! Total number of grid points
    z_start_rec = 1630.4d0                  ! Redshift of start of recombination
    z_end_rec   = 614.2d0                   ! Redshift of end of recombination
    z_0         = 0.d0                      ! Redshift today
    x_start_rec = -log(1.d0 + z_start_rec)  ! x of start of recombination
    x_end_rec   = -log(1.d0 + z_end_rec)    ! x of end o

    allocate(x_rec(n))
    allocate(a_rec(n))
    allocate(z_rec(n))
    allocate(H_rec(n))
    allocate(X_e(n))
    allocate(tau(n))
    allocate(dtau(n))
    allocate(ddtau(n))
    allocate(d4tau(n))
    allocate(n_e(n))
    allocate(n_e2(n))
    allocate(logn_e(n))
    allocate(logn_e2(n))
    allocate(g(n))
    allocate(g2(n))
    allocate(g22(n))
    allocate(x_test(n))
    allocate(z_test(n))
    allocate(a_test(n))
    allocate(n_etest(n))
    allocate(tau_test(n))
    allocate(dtau_test(n))
    allocate(ddtau_test(n))



    !fill test 
    do i=1,n
        x_test(i) = x_test_start + i*(x_test_end-x_test_start)/n
    end do
    z_test = 1.d0/exp(x_test) -1.d0

    !Fill in x,a,z (rec) grids
    do i = 1,n1 ! Fill interval during recombination
        x_rec(i)  = x_start_rec + (i-1)*(x_end_rec-x_start_rec)/(n1-1)
    end do
    do i = 1,n2 !Fill from end of recomb to today
        x_rec(n1+i) = x_end_rec + i*(x_0-x_end_rec)/(n2)
    end do
    write(*,*) x_rec
    a_rec = exp(x_rec)
    z_rec = 1.d0/a_rec -1.d0

    do i = 1,n
        H_rec(i) = get_H(x_rec(i))
    end do


    h1 = abs(1.d-2*(x_rec(1)-x_rec(2)))     !Defines the steplength to 100th of length between     
    h2 = abs(1.d-2*(x_rec(n-1)-x_rec(n-2))) !neighbouring x values, for both intervals

    !Since we have two different steplengths in our x array 
    !we choose the steplength that is smallest of the two parts
    if (h2<h1) then
    h1=h2
    end if
    

    !Compute X_e and n_e at all grid times
    use_saha = .true.
    do j = 1, n
        n_b = Omega_b*rho_c/(m_H*a_rec(j)**3)	
        if (use_saha) then
            ! Use the Saha equation
	    T_b = T_0/a_rec(j)
            X_econst = 1.d0/n_b*(m_e*k_b*T_b/(2.d0*hbar**2*pi))**1.5d0*exp(-epsilon_0/(k_b*T_b))
            X_e(j) = (-X_econst + sqrt(X_econst**2 +4.d0*X_econst))/2.d0

        if (X_e(j) < saha_limit) use_saha = .false.
        else
            ! Use the Peebles equation
            X_e(j) =X_e(j-1)
            call odeint(X_e(j:j),x_rec(j-1) ,x_rec(j), eps, h1, hmin, dX_edx, bsstep, output1) 
        end if
	n_e(j) = X_e(j)*n_b !Calculate electron density
    end do

    !Compute splined (log of) electron density function
    logn_e =log(n_e) !Turn n_e into its logarithm
    call spline(x_rec, logn_e, yp1, ypn,logn_e2)
    !call splint_deriv()
    !Test spline for x values between those used for spline
    do i=1,n  
        n_etest(i) = get_n_e(x_test(i))
    end do



    !Compute optical depth,and first deriv at all grid points
    tau(n) = 0 !Optical depth today is 0
    do k=n-1,1,-1
        tau(k) = tau(k+1)
        call odeint(tau(k:k),x_rec(k+1),x_rec(k),eps,h1,hmin,dtaudx,bsstep,output1)
    end do

    !Compute splined (log of) optical depth
    call spline(x_rec, tau, yp1, ypn,ddtau)
    !write(*,*) ddtau
    
    !Test the get_tau function
    do i=1,n
        tau_test(i) = get_tau(x_test(i))
    end do


    !Compute firste derivative of optical depth
    do i=1,n
        dtau(i) = -n_e(i)*sigma_T*c/H_rec(i)
    end do


    !Test the get_dtau function
    do i=1,n
        dtau_test(i) = get_dtau(x_test(i))
    end do

    !Compute splined second derivative of optical depth
    call spline(x_rec,ddtau,yp1,ypn,d4tau)
    !Test get_ddtau function
    do i=1,n
        ddtau_test(i) = get_ddtau(x_test(i))
    end do



    ! Task: Compute splined visibility function
    ! Task: Compute splined second derivative of visibility function


  end subroutine initialize_rec_mod

  !Begin Stuff needed to make odeint work
  subroutine dX_edx(x_rec, X_e, dydx) 
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x_rec
    real(dp), dimension(:), intent(in)  :: X_e
    real(dp), dimension(:), intent(out) :: dydx
    real(dp) :: T_b,n_b,phi2,alpha2,beta,beta2,n1s,lambda_alpha,C_r
    real(dp) :: Xe
    Xe = X_e(1)
    T_b          = T_0/a_rec(j)
    n_b          = Omega_b*rho_c/m_H/a_rec(j)**3
    phi2         = 0.448d0*log(epsilon_0/(k_b*T_b))
    alpha2       = 64.d0*pi/sqrt(27.d0*pi)*(alpha/m_e)**2*sqrt(epsilon_0/(k_b*T_b))*phi2 *hbar**2/c
    beta         = alpha2 *(m_e*k_b*T_b/(2.d0*hbar**2*pi))**1.5d0*exp(-epsilon_0/(k_b*T_b))
    !This part is needed since the exponent
    !in beta2 becomes so large that the computer 
    !sets it to infinity. However beta goes to zero before that 
    !so it should be 0 even if the exponent is enormous.
    if(T_b <= 169.d0) then
        beta2    = 0.d0
    else
        beta2    = beta*exp((3.d0*epsilon_0)/(4.d0*k_b*T_b))
    end if

    n1s          = (1.d0-Xe)*n_b
    lambda_alpha = H_rec(j)*(3.d0*epsilon_0)**3/((8.d0*pi)**2*n1s) /(c*hbar)**3

    C_r          = (lambda_2s1s +lambda_alpha)/(lambda_2s1s+lambda_alpha+beta2)
    dydx         = C_r/H_rec(j)*(beta*(1.d0-Xe)-n_b*alpha2*Xe**2)

    !Print values for testing
    !write(*,*) 'j =',j
    !write(*,*) T_b,n_b,X_econst
    !write(*,*) phi2,alpha2,beta
    !write(*,*) beta2,n1s,lambda_alpha
    !write(*,*) C_r,dydx,Xe    
    !write(*,*) beta,beta2,C_r
  end subroutine dX_edx

  subroutine dtaudx(x_rec,tau, dydx) 
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x_rec
    real(dp), dimension(:), intent(in)  :: tau
    real(dp), dimension(:), intent(out) :: dydx
    dydx         = -n_e(k)*sigma_T/H_rec(k)*c
  end subroutine dtaudx

  subroutine output1(x, y)
         use healpix_types
         implicit none
         real(dp),               intent(in)  :: x
         real(dp), dimension(:), intent(in)  :: y
  end subroutine output1
  !End Stuff needed to make odeint work


  !Complete routine for computing n_e at arbitrary x, using precomputed information
  function get_n_e(x_in)
    implicit none

    real(dp), intent(in) :: x_in
    real(dp)             :: get_n_e
    get_n_e = splint(x_rec, logn_e, logn_e2, x_in)
    !Return the actual n_e instead of log(n_e)
    get_n_e = exp(get_n_e)
  end function get_n_e

  ! Task: Complete routine for computing tau at arbitrary x, using precomputed information
  function get_tau(x_in)
    implicit none

    real(dp), intent(in) :: x_in
    real(dp)             :: get_tau
    get_tau  = splint(x_rec,tau,ddtau,x_in)
  end function get_tau

  ! Task: Complete routine for computing the derivative of tau at arbitrary x, using precomputed information
  function get_dtau(x_in)
    implicit none

    real(dp), intent(in) :: x_in
    real(dp)             :: get_dtau
    get_dtau = splint_deriv(x_rec,tau,ddtau,x_in)
  end function get_dtau

  ! Task: Complete routine for computing the second derivative of tau at arbitrary x, 
  ! using precomputed information
  function get_ddtau(x_in)
    implicit none

    real(dp), intent(in) :: x_in
    real(dp)             :: get_ddtau
    get_ddtau = splint(x_rec,ddtau,d4tau,x_in)
  end function get_ddtau

  ! Task: Complete routine for computing the visibility function, g, at arbitray x
  function get_g(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_g
get_g = 0
  end function get_g

  ! Task: Complete routine for computing the derivative of the visibility function, g, at arbitray x
  function get_dg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dg
get_dg=0
  end function get_dg

  ! Task: Complete routine for computing the second derivative of the visibility function, g, at arbitray x
  function get_ddg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddg
get_ddg=0
  end function get_ddg


end module rec_mod
