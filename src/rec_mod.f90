module rec_mod
    use healpix_types
    use params
    use time_mod
    use ode_solver
    use spline_1D_mod
    implicit none

    real(dp), allocatable, dimension(:)   :: x_rec,a_rec,z_rec             !Grid
    real(dp), allocatable, dimension(:)   :: tau_rec, dtau_rec, ddtau_rec  !Splined tau and derivatives
    real(dp), allocatable, dimension(:)   :: d4tau
    real(dp), allocatable, dimension(:)   :: n_e, n_e2,logn_e,logn_e2      !Splined (log of) n_e
    real(dp), allocatable, dimension(:)   :: g,dg,ddg,d4g                  !Splined visibility function
    real(dp), allocatable, dimension(:)   :: g_test,dg_test,ddg_test       !Splined visibility function
    real(dp), allocatable, dimension(:)   :: H_rec,X_e                     !Variables for H and X_e
    real(dp), allocatable, dimension(:)   :: x_test,n_etest,z_test,a_test  !Used for testing the splines
    real(dp), allocatable, dimension(:)   :: tau_test,dtau_test,ddtau_test !Used for testing the splines
    real(dp)                              :: x_0
    real(dp)                              :: x_test_start
    real(dp)                              :: x_test_end
    integer(i4b),private :: i
    integer(i4b) :: n
    !Variables used or for splines in this file, as well as rec_mod and cl_mod
    real(dp)     :: eps,hmin,yp1,ypn,h1,h2,h3

contains

  subroutine initialize_rec_mod
    implicit none

    real(dp)     :: saha_limit, y, T_b, n_b, dydx, xmin, xmax, dx, f, n_e0, X_e0, &
		    X_econst, phi2,alpha2,beta,beta2,n1s,lambda_alpha,C_r
    real(dp)     :: z_start_rec, z_end_rec, z_0, x_before_rec,x_start_rec, x_end_rec
    logical(lgt) :: use_saha

    n = n1+n2+n3
    x_test_start = -10.0d0
    x_test_end   = -1.0d0
    saha_limit   = 0.99d0       ! Switch from Saha to Peebles when X_e < 0.99

    !ODE int variables
    eps  = 1.d-10
    hmin = 0.d0
    !Spline variables
    yp1  = 1.d30
    ypn  = 1.d30

    !Grid sizes 

    z_start_rec = 1630.4d0                  ! Redshift of start of recombination
    z_end_rec   = 614.2d0                   ! Redshift of end of recombination
    z_0         = 0.d0                      ! Redshift today

    x_init      = log(a_init)!-10.5d0      ! x at start of array
    !write(*,*) x_init
    x_start_rec = -log(1.d0 + z_start_rec)  ! x of start of recombination
    x_end_rec   = -log(1.d0 + z_end_rec)    ! x of end of recombination

    !Allocate all the arrays which will be used
    allocate(x_rec(n),a_rec(n),z_rec(n),H_rec(n),X_e(n) &
            ,tau_rec(n),dtau_rec(n),ddtau_rec(n),d4tau(n),n_e(n) &
            ,n_e2(n),logn_e(n),logn_e2(n),g(n),dg(n) &
            ,ddg(n),d4g(n),x_test(n),z_test(n),a_test(n) &
            ,n_etest(n),tau_test(n),dtau_test(n)   &
            ,ddtau_test(n),g_test(n),dg_test(n),ddg_test(n))

    !x_rec = x_t
    !fill test 
    !do i=1,n
    !    x_test(i) = x_test_start + i*(x_test_end-x_test_start)/n
    !end do
    !z_test = 1.d0/exp(x_test) -1.d0

    !Fill in x,a,z (rec) grids
    do i = 1,n1
        x_rec(i)       = x_init + (i-1)*(x_start_rec-x_init)/(n1-1)
    end do
    do i = 1,n2+1 ! Fill interval during recombination
        x_rec(n1+i)    = x_start_rec + i*(x_end_rec-x_start_rec)/(n2)
    end do
    do i = 1,n3 !Fill from end of recomb to today
        x_rec(n1+n2+i) = x_end_rec + i*(x_0-x_end_rec)/(n3)
    end do

    !write(*,*) 'x_rec(1),x_rec(100)'
    !write(*,*) x_rec(1),x_rec(100)
    a_rec = exp(x_rec)
    z_rec = 1.d0/a_rec -1.d0
    
    do i = 1,n
        H_rec(i) = get_H(x_rec(i))
    end do
        
    h1 = abs(1.d-3*(x_rec(1)   -x_rec(2)))  !Defines the steplength to 100th of length between     
    h2 = abs(1.d-3*(x_rec(n1+1)-x_rec(n1))) !neighbouring x values, for all three intervals
    h3 = abs(1.d-3*(x_rec(n2+1)-x_rec(n2))) 
    !Since we have three different steplengths in our x array 
    !we choose the steplength that is smallest of the three parts
    if (h3<h2) then
    h2 = h3
    end if
    if (h2<h1) then
    h1=h2
    end if

    !Compute X_e and n_e at all grid times
    use_saha = .true.
    do i = 1, n
        n_b = Omega_b*rho_c/(m_H*a_rec(i)**3)	
        if (use_saha) then
            ! Use the Saha equation
            T_b = T_0/a_rec(i)
            X_econst = ((m_e*k_b*T_b)/(2.d0*pi*hbar**2))**1.5d0*exp(-epsilon_0/(k_b*T_b))/n_b
            X_e(i) = (-X_econst + sqrt(X_econst**2 +4.d0*X_econst))/2.d0

        if (X_e(i) < saha_limit) use_saha = .false.
        else
            ! Use the Peebles equation
            X_e(i) =X_e(i-1)
            call odeint(X_e(i:i),x_rec(i-1) ,x_rec(i), eps, h1, hmin, dX_edx, bsstep, output) 
        end if
        n_e(i) = X_e(i)*n_b !Calculate electron density
        !write(*,*) use_saha,x_rec(i), X_e(i)
    end do

    !Compute splined (log of) electron density function
    logn_e =log(n_e)
    call spline(x_rec, logn_e, yp1, ypn,logn_e2)


    !Test spline for x values between those used for spline
    do i=1,n  
        n_etest(i) = get_n_e(x_test(i))
    end do

    !Compute optical depth,and first deriv at all grid points
    tau_rec(n) = 0.d0 !Optical depth today is 0
    do i=n-1,1,-1
        tau_rec(i) = tau_rec(i+1)
        call odeint(tau_rec(i:i),x_rec(i+1),x_rec(i),eps,h1,hmin,dtaudx,bsstep,output)
    end do
    !write(*,*)'sum(tau_rec)'
    !write(*,*) sum(tau_rec)

    !Compute splined optical depth,and second derivative
    call spline(x_rec, tau_rec, yp1, ypn,ddtau_rec)
    call spline(x_rec,ddtau_rec,yp1,ypn,d4tau)
    !write(*,*) ddtau_rec
 
    !Compute first derivative of optical depth
    do i=1,n
        dtau_rec(i) = -n_e(i)*sigma_T*c/H_rec(i)
    end do
   
    !Test the get_tau,get_dtau,get_ddtau function
    do i=1,n
        tau_test(i) = get_tau(x_test(i))
        dtau_test(i) = get_dtau(x_test(i))
        ddtau_test(i) = get_ddtau(x_test(i))
    end do

    
    !Compute g for values in x_rec
    do i=1,n
        g(i) = -dtau_rec(i)*exp(-tau_rec(i))
    end do
    !write(*,*)'sum(g)'
    !write(*,*) sum(g)
    !stop
    !Compute splined visibility function, and second derivative
    call spline(x_rec,g,yp1,ypn,ddg)
    call spline(x_rec,ddg,yp1,ypn,d4g)
    !Test get_g and get_dg
    do i=1,n
        dg(i)       = get_dg(x_rec(i))
        g_test(i)   = get_g(x_test(i))
        dg_test(i)  = get_dg(x_test(i))
        ddg_test(i) = get_ddg(x_test(i))
    end do


  end subroutine initialize_rec_mod

  !Begin Stuff needed to make odeint work
  subroutine dX_edx(x, X_e, dydx) 
      use healpix_types
      implicit none
      real(dp),               intent(in)  :: x
      real(dp), dimension(:), intent(in)  :: X_e
      real(dp), dimension(:), intent(out) :: dydx
      real(dp) :: T_b,n_b,phi2,alpha2,beta,beta2,n1s,lambda_alpha,C_r
      real(dp) :: Xe
      real(dp) :: a
      real(dp) :: H 
      Xe = X_e(1)
      a  = exp(x)
      H  = get_H(x)
      T_b          = T_0/a
      n_b          = Omega_b*rho_c/(m_H*a**3)

      phi2         = 0.448d0*log(epsilon_0/(k_b*T_b))
      alpha2       = 64.d0*pi/sqrt(27.d0*pi)*(alpha/m_e)**2*sqrt(epsilon_0/(k_b*T_b))*phi2 *hbar**2/c
      beta         = alpha2 *((m_e*k_b*T_b)/(2.d0*pi*hbar**2))**1.5*exp(-epsilon_0/(k_b*T_b))

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
      lambda_alpha = H*(3.d0*epsilon_0)**3/((8.d0*pi)**2*n1s) /(c*hbar)**3
  


      C_r          = (lambda_2s1s +lambda_alpha)/(lambda_2s1s+lambda_alpha+beta2)
      dydx         = C_r/H*(beta*(1.d0-Xe)- n_b*alpha2*Xe**2)

      !Print values for testing
      !write(*,*) 'i =',i
      !write(*,*) T_b,n_b,X_econst
      !write(*,*) phi2,alpha2,beta
      !write(*,*) beta2,n1s,lambda_alpha
      !write(*,*) C_r,dydx,Xe    
      !write(*,*) beta,beta2,C_r
  end subroutine dX_edx

  subroutine dtaudx(x,tau, dydx) 
      use healpix_types
      implicit none
      real(dp),               intent(in)  :: x
      real(dp), dimension(:), intent(in)  :: tau
      real(dp), dimension(:), intent(out) :: dydx
      real(dp)                            :: n_e
      real(dp)                            :: H
      n_e  = get_n_e(x)
      H    = get_H(x)
      dydx = -n_e*sigma_T/H*c
  end subroutine dtaudx

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

  !Routine for computing tau at arbitrary x, using precomputed information
  function get_tau(x_in)
      implicit none
      real(dp), intent(in) :: x_in
      real(dp)             :: get_tau
      get_tau  = splint(x_rec,tau_rec,ddtau_rec,x_in)
  end function get_tau

  !Routine for computing the derivative of tau at arbitrary x, using precomputed information
  function get_dtau(x_in)
      implicit none
      real(dp), intent(in) :: x_in
      real(dp)             :: get_dtau
      real(dp)             :: n_e,a,H_p
      H_p = get_H_p(x_in)
      a = exp(x_in)
      n_e = get_n_e(x_in)
      get_dtau = -n_e*sigma_T*a*c/H_p!splint_deriv(x_rec,tau_rec,ddtau_rec,x_in)
  end function get_dtau

  !Routine for computing the second derivative of tau at arbitrary x, 
  !using precomputed information
  function get_ddtau(x_in)
      implicit none
      real(dp), intent(in) :: x_in
      real(dp)             :: get_ddtau
      get_ddtau = splint(x_rec,ddtau_rec,d4tau,x_in)
  end function get_ddtau

  !Routine for computing the visibility function, g, at arbitray x
  function get_g(x_in)
      implicit none
      real(dp), intent(in) :: x_in
      real(dp)             :: get_g
      real(dp)             :: dtau
      real(dp)             :: tau
      dtau  =  get_dtau(x_in)
      tau   =  get_tau(x_in)
      get_g = -dtau*exp(-tau)
  end function get_g

  !Routine for computing the derivative of the visibility function, g, at arbitrary x
  function get_dg(x_in)
      implicit none
      real(dp), intent(in) :: x_in
      real(dp)             :: get_dg
      get_dg= splint_deriv(x_rec,g,ddg,x_in)
  end function get_dg

  !Routine for computing the second derivative of the visibility function, g, at arbitrary x
  function get_ddg(x_in)
      implicit none
      real(dp), intent(in) :: x_in
      real(dp)             :: get_ddg
      get_ddg = splint(x_rec,ddg,d4g,x_in)
  end function get_ddg

end module rec_mod
