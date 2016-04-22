module evolution_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use rec_mod
  implicit none

  ! Accuracy parameters
  !real(dp),     parameter, private :: a_init   = 1.d-8
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
  !logical(lgt) :: polarize = False
contains


  ! NB!!! New routine for 4th milestone only; disregard until then!!!
  !subroutine get_hires_source_function(k, x, S)
  !  implicit none

  !  real(dp), pointer, dimension(:),   intent(out) :: k, x
  !  real(dp), pointer, dimension(:,:), intent(out) :: S

  !   integer(i4b) :: i, j
  !  real(dp)     :: g, dg, ddg, tau, dt, ddt, H_p, dH_p, ddHH_p, Pi, dPi, ddPi
  !  real(dp), allocatable, dimension(:,:) :: S_lores

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

 ! end subroutine get_hires_source_function




  ! Routine for initializing and solving the Boltzmann and Einstein equations
  subroutine initialize_perturbation_eqns
    implicit none

    integer(i4b) :: i,l, k
    real(dp) :: k_min = 0.1d0*H_0
    real(dp) :: k_max = 1.0d3*H_0

    !Initialize k-grid, ks; quadratic between k_min and k_max
    allocate(ks(n_k))
    do i=1,n_k
        ks(i) = k_min +(k_max -k_min)*(i/1.0d2)**2
    end do

    !Allocate arrays for perturbation quantities
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
       
    do k = 1, n_k
        v(0,k)       = c*ks(k)/(2.d0*get_H_p(x_t(1)))*Phi(0,1)
        v_b(0,k)     = v(0,k)
        Theta(0,0,k) = 0.5d0*Phi(0,1)
        Theta(0,1,k) = -c*ks(k)/(6.d0*get_H_p(x_t(1)))*Phi(0,1)

        !if(polarize==True)
        !    Theta(0,2,i) = -8.d0*c*ks(i)/(15.d0*get_H_p(x_t(i))*get_dtau(x_t(i)))*Theta(0,1,i) !with polarization
        !else
        Theta(0,2,k) = -20.d0*c*ks(k)/(45.d0*get_H_p(x_t(1))*get_dtau(x_t(1)))*Theta(0,1,k) !without polarization
        !end if

        do l = 3, lmax_int
            Theta(0,l,k) = -l/(2.d0*l+1.d0)*c*ks(k)/(get_H_p(x_t(1))*get_dtau(x_t(1)))*Theta(0,l-1,k	)
        end do

    end do
  end subroutine initialize_perturbation_eqns

  subroutine integrate_perturbation_eqns
    implicit none

    integer(i4b) :: i, j, k, l
    real(dp)     :: x1, x2, x_init
    real(dp)     :: eps, hmin, h1, x_tc, i_tc,H_p, dt, t1, t2
    real(dp)     :: R,d_v,d_v_b,d_Theta0,d_Theta_1,Theta_2,q

    real(dp), allocatable, dimension(:) :: y, y_tight_coupling, dydx

    eps    = 1.d-8
    hmin   = 0.d0

    allocate(y(npar))
    allocate(dydx(npar))
    allocate(y_tight_coupling(7))

    ! Propagate each k-mode independently
    do k = 1, n_k
       write(*,*) 'Current k', k
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
       write(*,*) 'over x_tc'
       x_tc = get_tight_coupling_time(k_current)
       write(*,*) 'x_tc =',x_tc
       write(*,*) 'under x_tc'
       ! Task: Integrate from x_init until the end of tight coupling, using
       !       the tight coupling equations
       write(*,*) 'Start of tight coupling'
       do i=2,n_t
           write(*,*) 'i=',i
           if (x_t(i)< x_tc) then 
           !Solve next step
           call odeint(y_tight_coupling,x_t(i-1),x_t(i),eps,h1,hmin,derivs_tc, bsstep, output3)
           !Save variables
           delta(i,k)   = y_tight_coupling(1)
           delta_b(i,k) = y_tight_coupling(2)
           v(i,k)       = y_tight_coupling(3)
           v_b(i,k)     = y_tight_coupling(4)
           Phi(i,k)     = y_tight_coupling(5)
           Theta(i,0,k) = y_tight_coupling(6)
           Theta(i,1,k) = y_tight_coupling(7)
           Theta(0,2,k) = -20.d0*c*k_current/(45.d0*get_H_p(x_t(i))*get_dtau(x_t(i)))*Theta(i,1,k) !without polarization
           do l = 3, lmax_int
              Theta(i,l,k) = -l/(2.d0*l+1.d0)*c*k_current/(get_H_p(x_t(i))*get_dtau(x_t(i)))*Theta(0,l-1,k)
           end do

           else
               i_tc = i
               exit
           end if
       end do
       write(*,*) 'End of tight coupling'

       ! Task: Set up variables for integration from the end of tight coupling 
       ! until today
       y(1:7) = y_tight_coupling(1:7)
       y(8)   = Theta(0,2,k)
       do l = 3, lmax_int
          y(6+l) = Theta(0,l,k)
       end do


       !Continue after tight coupling
       write(*,*) 'start of rec'       
       do i = i_tc, n_t

          ! Task: Integrate equations from tight coupling to today
          write(*,*) 'running odeint with i =', i
          call odeint(y, x_t(i-1) ,x_t(i), eps, h1, hmin, derivs, bsstep, output3)

          ! Task: Store variables at time step i in global variables
          delta(i,k)   = y(1)
          delta_b(i,k) = y(2)
          v(i,k)       = y(3)
          v_b(i,k)     = y(4)
          Phi(i,k)     = y(5)
          
          do l = 0, lmax_int
             Theta(i,l,k) = y(6+l)
          end do
          Psi(i,k)     =  - Phi(i,k) - 12.d0*H_0**2/(c*k_current*a_t(i))**2*Omega_r*Theta(i,2,k)

          ! Task: Store derivatives that are required for C_l estimation
          dPhi(i,k)     = Psi(i,k) -c**2*k_current**2/(3.d0*get_H_p(x_t(i))**2)*Phi(i,k) +H_0**2/(2.d0*get_H_p(x_t(i))) &
                          *(Omega_m/a_t(i)*delta(i,k) +Omega_b/a_t(i)*delta_b(i,k_current) + 4.d0*Omega_r/a_t(i)**2 &
                          *Theta(i,0,k_current))

          dv_b(i,k)     = -v_b(i,k) -c*k_current/get_H_p(x_t(i))*Psi(i,k) +get_dtau(x_t(i))*R*(3.d0*Theta(i,1,k)+ v_b(i,k))


          dTheta(i,0,k) = -c*k_current/get_H_p(x_t(i))*Theta(i,1,k) -dPhi(i,k)
          dTheta(i,1,k) = c*k_current/(3.d0*get_H_p(x_t(i)))*Theta(i,0,k) - &
                          2.d0*c*k_current/(3.d0*get_H_p(x_t(i)))*Theta(i,2,k)+&
                          c*k_current/(3.d0*get_H_p(x_t(i)))*Psi(i,k) + &
                          get_dtau(x_t(i))*(Theta(i,1,k)+ 1.d0/3.d0*v_b(i,k))
          dTheta(i,2,k) = l*c*k_current/((2.d0*l+1.d0)*get_H_p(x_t(i)))*Theta(i,l-1,k) -&
                          (l+1.d0)*c*k_current/((2.d0*l+1.d0)*get_H_p(x_t(i)))*Theta(i,l+1,k)+&
                          get_dtau(x_t(i))*(Theta(i,l,k) - 1.d0/10.d0*Theta(i,l,k))
          do l=3,lmax_int-1
              dTheta(i,l,k) = l*c*k_current/((2.d0*l+1.d0)*get_H_p(x_t(i)))*Theta(i,l-1,k) -&
                              (l+1.d0)*c*k_current/((2.d0*l+1.d0)*get_H_p(x_t(i)))*Theta(i,l+1,k)+&
                              get_dtau(x_t(i))*Theta(i,l,k)
          end do
          dTheta(i,lmax_int,k) = c*k_current/get_H_p(x_t(i))*Theta(i,l-1,k) -&
                                 c*(l+1.d0)/(get_H_p(x_t(i))*get_eta(x_t(i)))*&
                                 Theta(i,l,k) + get_dtau(x_t(i))*Theta(i,l,k)
          !dPsi(i,k)     = 
       end do
       write(*,*) 'today'
    end do
    deallocate(y_tight_coupling)
    deallocate(y)
    deallocate(dydx)

  end subroutine integrate_perturbation_eqns

  subroutine derivs_tc(x,y, dydx)
      use healpix_types
      implicit none
      real(dp),               intent(in)  :: x
      real(dp), dimension(:), intent(in)  :: y
      real(dp), dimension(:), intent(out) :: dydx

      real(dp) :: a,H_p,dtau,dH_p
      real(dp) :: delta,d_delta
      real(dp) :: delta_b,d_delta_b
      real(dp) :: v,d_v
      real(dp) :: v_b,d_v_b
      real(dp) :: Phi,d_Phi
      real(dp) :: Theta_0,d_Theta_0
      real(dp) :: Theta_1,d_Theta_1
      real(dp) :: Theta_2
      real(dp) :: Psi,q,R

      H_p     = get_H_p(x)
      dH_p    = get_dH_p(x)
      dtau    = get_dtau(x)
      a       = exp(x)
      delta   = y(1)
      delta_b = y(2)
      v       = y(3)
      v_b     = y(4)
      Phi     = y(5)
      Theta_0 = y(6)
      Theta_1 = y(7)

      Theta_2   = -20.d0*c*k_current/(45.d0*H_p*dtau)*Theta_1
      R         = 4.d0*Omega_r/(3.d0*Omega_b*a)
      Psi       = -Phi - 12.d0*H_0**2/(c*k_current*a)**2*Omega_r*Theta_2

      d_Phi     = Psi -c**2*k_current**2/(3.d0*H_p**2)*Phi +H_0**2/(2.d0*H_p) &
                  *(Omega_m/a*delta +Omega_b/a*delta_b + 4.d0*Omega_r/a**2 &
                  *Theta_0)

      d_Theta_0 = -c*k_current/H_p*Theta_1 -d_Phi

      d_delta   = c*k_current/H_p*v  -3d0*d_Phi
      d_delta_b = c*k_current/H_p*v_b-3d0*d_Phi
      d_v       = -v -c*k_current/H_p*Psi

      q         = (-((1.d0-2.d0*R)*dtau + (1.d0+R)*get_ddtau(x))* &
                  (3.d0*Theta_1 +v_b) -c*k_current/H_p* &
                  Psi +(1.d0-dH_p /H_p)*c*k_current/ &
                  H_p*(-Theta_0+2.d0*Theta_2)-c*k_current / &
                  H_p*d_Theta_0)/((1.d0+R)*dtau+ &
                  dH_p/H_p -1.d0)
      d_v_b     = 1.d0/(1.d0+R)*(-v_b-c*k_current/H_p*Psi &
                  +R*(q+c*k_current/H_p*(-Theta_0+2.d0*Theta_2) &
                  -c*k_current/H_p*Psi))
      d_Theta_1 = 1.d0/3.d0*(q-d_v_b)
      dydx(1) = d_delta
      dydx(2) = d_delta_b
      dydx(3) = d_v
      dydx(4) = d_v_b
      dydx(5) = d_Phi
      dydx(6) = d_Theta_0
      dydx(7) = d_Theta_1      
  end subroutine derivs_tc

  subroutine derivs(x,y, dydx) 
      use healpix_types
      implicit none
      real(dp),               intent(in)  :: x
      real(dp), dimension(:), intent(in)  :: y
      real(dp), dimension(:), intent(out) :: dydx

      real(dp) :: a,H_p,dH_p
      real(dp) :: Psi,dtau
      real(dp) :: Phi,d_Phi
      real(dp) :: delta,d_delta
      real(dp) :: delta_b,d_delta_b
      real(dp) :: v,d_v
      real(dp) :: v_b,d_v_b
      real(dp) :: Theta_0,d_Theta_0
      real(dp) :: Theta_1,d_Theta_1
      real(dp) :: Theta_2
      real(dp) :: R,q


      H_p     = get_H_p(x)
      dH_p    = get_dH_p(x)
      dtau    = get_dtau(x)
      a       = exp(x)
      delta   = y(1)
      delta_b = y(2)
      v       = y(3)
      v_b     = y(4)
      Phi     = y(5)
      Theta_0 = y(6)
      Theta_1 = y(7)
      Theta_2 = y(8)

      R         = 4.d0*Omega_R/(3.d0*Omega_b*a)

      Psi       = -Phi - 12.d0*H_0**2/(c*k_current*a)**2*Omega_r*Theta_2
      d_Phi     = Psi -c**2*k_current**2/(3.d0*H_p**2)*Phi +H_0**2/(2.d0*H_p) &
                  *(Omega_m/a*delta +Omega_b/a*delta_b + 4.d0*Omega_r/a**2 &
                  *Theta_0)
      d_delta   = c*k_current/H_p*v-3d0*d_Phi
      d_delta_b = c*k_current/H_p*v_b-3d0*d_Phi
      d_v       = -v -c*k_current/H_p*Psi
      d_v_b     = -v_b -c*k_current/H_p*Psi + dtau*R*(3.d0*Theta_1+v_b)


      d_Theta_0 = -    c*k_current/H_p*Theta_1 -d_Phi
      d_Theta_1 =      c*k_current/(3.d0*H_p)*Theta_0 - &
                  2.d0*c*k_current/(3.d0*H_p)*Theta_2+&
                       c*k_current/(3.d0*H_p)*Psi + &
                  dtau*(Theta_1+ 1.d0/3.d0*v_b)

      dydx(1) = d_delta
      dydx(2) = d_delta_b
      dydx(3) = d_v
      dydx(4) = d_v_b
      dydx(5) = d_Phi
      dydx(6) = d_Theta_0
      dydx(7) = d_Theta_1
  end subroutine derivs

  subroutine output3(x, y)
      use healpix_types
      implicit none
      real(dp),               intent(in)  :: x
      real(dp), dimension(:), intent(in)  :: y
  end subroutine output3



  ! Task: Complete the following routine, such that it returns the time at which
  !       tight coupling ends. In this project, we define this as either when
  !       dtau < 10 or c*k/(H_p*dt) > 0.1 or x > x(start of recombination)
  function get_tight_coupling_time(k)
    implicit none

    real(dp), intent(in)  :: k
    real(dp)              :: get_tight_coupling_time
    integer(i4b)          :: i,n
    real(dp)              :: x,a
    n =1d6
    do i=0,n
        x = x_init +i*(0.d0-x_init)/n
        !write(*,*) x,x_start_rec
        if (x < x_start_rec .and. abs(c*k/(get_H_p(x)*get_dtau(x))) <= 0.1d0 .and. abs(get_dtau(x)) > 10.d0) then 
            get_tight_coupling_time = x
        end if
    end do
  end function get_tight_coupling_time

end module evolution_mod
