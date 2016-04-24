module evolution_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use rec_mod
  implicit none

  !Use j,k,l as global variable
  integer(i4b) :: j,k,l

  ! Accuracy parameters
  real(dp),     parameter, private :: k_min    = 0.1d0 * H_0 / c
  real(dp),     parameter, private :: k_max    = 1.d3  * H_0 / c
  integer(i4b), parameter          :: n_k      = 1!00
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

  real(dp), allocatable, dimension(:) :: dtau
  real(dp), allocatable, dimension(:) :: ddtau
  real(dp), allocatable, dimension(:) :: H_p
  real(dp), allocatable, dimension(:) :: dH_p

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
    integer(i4b) :: i
    real(dp)     :: k_min = 0.1d0*H_0
    real(dp)     :: k_max = 1.0d3*H_0

    !Initialize k-grid, ks; quadratic between k_min and k_max
    allocate(ks(n_k))
    do k=1,n_k
        ks(k) = k_min +(k_max -k_min)*(k/100.d0)**2
    end do

    !Allocate arrays for perturbation quantities
    allocate(Theta(1:n_t, 0:lmax_int, n_k))
    allocate(dTheta(1:n_t, 0:lmax_int, n_k))
    allocate(delta(1:n_t, n_k))
    allocate(delta_b(1:n_t, n_k))
    allocate(v(1:n_t, n_k))
    allocate(v_b(1:n_t, n_k))
    allocate(Phi(1:n_t, n_k))
    allocate(Psi(1:n_t, n_k))
    allocate(dPhi(1:n_t, n_k))
    allocate(dPsi(1:n_t, n_k))
    allocate(dv_b(1:n_t, n_k))

    !Allocate arrays for precomputed variables
    allocate(dtau(n_t),H_p(n_t),dH_p(n_t))
    allocate(ddtau(n_t))
    !Precompute useful variables
    do i=1,n_t
       dtau(i)  = get_dtau(x_t(i))
       ddtau(i) = get_ddtau(x_t(i))
       H_p(i)   = get_H_p(x_t(i))
       dH_p(i)  = get_dH_p(x_t(i))
    end do


    ! Task: Set up initial conditions for the Boltzmann and Einstein equations
    Phi(1,:)     = 1.d0
    delta(1,:)   = 1.5d0*Phi(1,:)
    delta_b(1,:) = delta(1,:)

    do k = 1, n_k
        v(1,k)       = c*ks(k)/(2.d0*get_H_p(x_t(1)))*Phi(1,k)
        v_b(1,k)     = v(1,k)
        Theta(1,0,k) = 0.5d0*Phi(1,k)
        Theta(1,1,k) = -c*ks(k)/(6.d0*get_H_p(x_t(1)))*Phi(1,k)
        Theta(1,2,k) = -20.d0*c*ks(k)/(45.d0*get_H_p(x_t(1))*get_dtau(x_t(1)))*Theta(1,1,k) !without polarization

        do l = 3, lmax_int
            Theta(1,l,k) = -l/(2.d0*l+1.d0)*c*ks(k)/(get_H_p(x_t(1))*get_dtau(x_t(1)))*Theta(1,l-1,k)
        end do

    end do
  end subroutine initialize_perturbation_eqns

  subroutine integrate_perturbation_eqns
    implicit none
    real(dp)     :: x1, x2, x_init
    real(dp)     :: eps, hmin, h1, x_tc, j_tc, dt, t1, t2
    real(dp)     :: R,d_v,d_v_b,q
    real(dp), allocatable, dimension(:) :: y, y_tight_coupling, dydx
    eps    = 1.d-8
    hmin   = 0.d0
    h1     = 1.d-5
    allocate(y(npar))
    allocate(dydx(npar))
    allocate(y_tight_coupling(7))

    ! Propagate each k-mode independently
    do k = 1, n_k
       write(*,*) 'Current k', k
       k_current = ks(k)  ! Store k_current as a global module variable

       ! Initialize equation set for tight coupling
       y_tight_coupling(1) = delta(1,k)
       y_tight_coupling(2) = delta_b(1,k)
       y_tight_coupling(3) = v(1,k)
       y_tight_coupling(4) = v_b(1,k)
       y_tight_coupling(5) = Phi(1,k)
       y_tight_coupling(6) = Theta(1,0,k)
       y_tight_coupling(7) = Theta(1,1,k)
       
       ! Find the time to which tight coupling is assumed, 
       ! and integrate equations to that time
       x_tc = get_tight_coupling_time(k_current)
       write(*,*) 'x_tc =',x_tc
       write(*,*) 'under x_tc'

       ! Task: Integrate from x_init until the end of tight coupling, using
       !       the tight coupling equations
       write(*,*) 'Start of tight coupling'
       write (*,'(*(2X, ES14.6))') Psi(1,k) ,Phi(1,k) ,delta(1,k),delta_b(1,k),v(1,k),v_b(1,k) ,Theta(1,0,k) ,Theta(1,1,k)
       !write (*,'(*(2X, ES14.6))') dPsi(1,k),dPhi(1,k),dv_b(1,k),dTheta(1,0,k),dTheta(1,1,k)

       do j=2,n_t
           !write(*,*) 'i=',i
           if (x_t(j)< x_tc) then 
               !Solve next step
               call odeint(y_tight_coupling,x_t(j-1),x_t(j),eps,h1,hmin,derivs_tc, bsstep, output3)
               !Save variables
               delta(j,k)   = y_tight_coupling(1)
               delta_b(j,k) = y_tight_coupling(2)
               v(j,k)       = y_tight_coupling(3)
               v_b(j,k)     = y_tight_coupling(4)
               Phi(j,k)     = y_tight_coupling(5)
               Psi(j,k)     = -Phi(j,k) - 12.d0*H_0**2/(c*k_current*a_t(j))**2*Omega_r*Theta(j,2,k)
               Theta(j,0,k) = y_tight_coupling(6)
               Theta(j,1,k) = y_tight_coupling(7)
               Theta(j,2,k) = -(20.d0*c*k_current)/(45.d0*H_p(j)*dtau(j))*Theta(j,1,k)
               do l = 3, lmax_int
                  Theta(j,l,k) = -l/(2.d0*l+1.d0)*c*k_current/(get_H_p(x_t(j))*get_dtau(x_t(j)))*Theta(j,l-1,k)
               end do

               ! Task: Store derivatives that are required for C_l estimation
               dPsi(j,k)    = -dPhi(j,k) - 12.d0*H_0**2/(c*k_current*a_t(j))**2*Omega_r*(-2.d0*Theta(j,2,k)+dTheta(j,2,k))

               !write (*,'(*(2X, ES14.6))') Psi(j,k) ,Phi(j,k) ,delta(j,k),delta_b(j,k),v(j,k),v_b(j,k) ,Theta(j,0,k) ,Theta(j,1,k)
               write (*,'(*(2X, ES14.6))') dPsi(j,k),dPhi(j,k),dv_b(j,k),dTheta(j,0,k),dTheta(j,1,k),dTheta(j,2,k)
           else
               j_tc = j
               exit
           end if
       end do
       !write(*,*) 'End of tight coupling'

       ! Task: Set up variables for integration from the end of tight coupling 
       ! until today
       !y(1:7) = y_tight_coupling(1:7)
       !y(8)   = Theta(0,2,k)
       !do l = 3, lmax_int
       !   y(6+l) = Theta(0,l,k)
       !end do


       !Continue after tight coupling
       !write(*,*) 'start of rec'       
       !do j = i_tc, n_t

          ! Task: Integrate equations from tight coupling to today
!          write(*,*) 'running odeint with j =', j
 !         call odeint(y, x_t(i-1) ,x_t(i), eps, h1, hmin, derivs, bsstep, output3)

          ! Task: Store variables at time step i in global variables
!          delta(i,k)   = y(1)
 !         delta_b(i,k) = y(2)
  !        v(i,k)       = y(3)
   !       v_b(i,k)     = y(4)
    !      Phi(i,k)     = y(5)
          
!          do l = 0, lmax_int
 !            Theta(i,l,k) = y(6+l)
  !        end do
   !       Psi(i,k)     =  - Phi(i,k) - 12.d0*H_0**2/(c*k_current*a_t(i))**2*Omega_r*Theta(i,2,k)

          ! Task: Store derivatives that are required for C_l estimation
!          dPhi(i,k)     = Psi(i,k) -c**2*k_current**2/(3.d0*get_H_p(x_t(i))**2)*Phi(i,k) +H_0**2/(2.d0*get_H_p(x_t(i))) &
 !                         *(Omega_m/a_t(i)*delta(i,k) +Omega_b/a_t(i)*delta_b(i,k_current) + 4.d0*Omega_r/a_t(i)**2 &
  !                        *Theta(i,0,k_current))

 !         dv_b(i,k)     = -v_b(i,k) -c*k_current/get_H_p(x_t(i))*Psi(i,k) +get_dtau(x_t(i))*R*(3.d0*Theta(i,1,k)+ v_b(i,k))


!          dTheta(i,0,k) = -c*k_current/get_H_p(x_t(i))*Theta(i,1,k) -dPhi(i,k)
 !         dTheta(i,1,k) = c*k_current/(3.d0*get_H_p(x_t(i)))*Theta(i,0,k) - &
  !                        2.d0*c*k_current/(3.d0*get_H_p(x_t(i)))*Theta(i,2,k)+&
   !                       c*k_current/(3.d0*get_H_p(x_t(i)))*Psi(i,k) + &
    !                      get_dtau(x_t(i))*(Theta(i,1,k)+ 1.d0/3.d0*v_b(i,k))
     !     dTheta(i,2,k) = l*c*k_current/((2.d0*l+1.d0)*get_H_p(x_t(i)))*Theta(i,l-1,k) -&
      !                    (l+1.d0)*c*k_current/((2.d0*l+1.d0)*get_H_p(x_t(i)))*Theta(i,l+1,k)+&
       !                   get_dtau(x_t(i))*(Theta(i,l,k) - 1.d0/10.d0*Theta(i,l,k))
!          do l=3,lmax_int-1
 !             dTheta(i,l,k) = l*c*k_current/((2.d0*l+1.d0)*get_H_p(x_t(i)))*Theta(i,l-1,k) -&
  !                            (l+1.d0)*c*k_current/((2.d0*l+1.d0)*get_H_p(x_t(i)))*Theta(i,l+1,k)+&
   !                           get_dtau(x_t(i))*Theta(i,l,k)
    !      end do
     !     dTheta(i,lmax_int,k) = c*k_current/get_H_p(x_t(i))*Theta(i,l-1,k) -&
      !                           c*(l+1.d0)/(get_H_p(x_t(i))*get_eta(x_t(i)))*&
       !                          Theta(i,l,k) + get_dtau(x_t(i))*Theta(i,l,k)
          !dPsi(i,k)     = 
!       end do
 !      write(*,*) 'today'
    end do
    deallocate(y_tight_coupling)
    deallocate(y)
    deallocate(dydx)

  end subroutine integrate_perturbation_eqns

  subroutine derivs_tc(x,y_tc, dydx)
      use healpix_types
      implicit none
      real(dp),               intent(in)  :: x
      real(dp), dimension(:), intent(in)  :: y_tc
      real(dp), dimension(:), intent(out) :: dydx

      real(dp) :: d_delta
      real(dp) :: d_delta_b
      real(dp) :: d_v
      real(dp) :: q,R

      delta(j,k)   = y_tc(1)
      delta_b(j,k) = y_tc(2)
      v(j,k)       = y_tc(3)
      v_b(j,k)     = y_tc(4)
      Phi(j,k)     = y_tc(5)
      Theta(j,0,k) = y_tc(6)
      Theta(j,1,k) = y_tc(7)
      Theta(j,2,k) = -20.d0*c*k_current/(45.d0*H_p(j)*dtau(j))*Theta(j,1,k)

      do l=3,lmax_int
          Theta(j,l,k) = -l/(2.d0*l+1.d0)*c*k_current/(H_p(j)*dtau(j))*Theta(j,l-1,k)
      end do

      R         = 4.d0*Omega_r/(3.d0*Omega_b*a_t(j))

      Psi(j,k)  = -Phi(j,k) - 12.d0*H_0**2/(c*k_current*a_t(j))**2*Omega_r*Theta(j,2,k)

      dPhi(j,k) = Psi(j,k) -(c*k_current)**2/(3.d0*H_p(j)**2)*Phi(j,k) +H_0**2/(2.d0*H_p(j)**2) &
                  *(Omega_m/a_t(j)*delta(j,k) +Omega_b/a_t(j)*delta_b(j,k) + 4.d0*Omega_r/a_t(j)**2 &
                  *Theta(j,0,k))

      dTheta(j,0,k) = -c*k_current/H_p(j)*Theta(j,1,k) - dPhi(j,k)

      d_delta   = c*k_current/H_p(j)*v(j,k)   - 3.d0*dPhi(j,k)

      d_delta_b = c*k_current/H_p(j)*v_b(j,k) - 3.d0*dPhi(j,k)

      d_v       = -v(j,k) -c*k_current/H_p(j)*Psi(j,k)

      q         = (-((1.d0-2.d0*R)*dtau(j) + &
                  (1.d0+R)*ddtau(j))*(3.d0*Theta(j,1,k)+v_b(j,k)) - &
                  c*k_current/H_p(j)*Psi(j,k) +&
                  (1.d0-dH_p(j)/H_p(j))*c*k_current/H_p(j)*(-Theta(j,0,k)+2.d0*Theta(j,2,k))-&
                  c*k_current/H_p(j)*dTheta(j,0,k))/((1.d0+R)*dtau(j)+dH_p(j)/H_p(j) -1.d0)

      dv_b(j,k) = 1.d0/(1.d0+R)*(-v_b(j,k)-c*k_current/H_p(j)*Psi(j,k)+&
                  R*(q+c*k_current/H_p(j)*(-Theta(j,0,k)+2.d0*Theta(j,2,k))-&
                  c*k_current/H_p(j)*Psi(j,k)))

      dTheta(j,1,k) = 1.d0/3.d0*(q-dv_b(j,k))

      dydx(1) = d_delta
      dydx(2) = d_delta_b
      dydx(3) = d_v
      dydx(4) = dv_b(j,k)
      dydx(5) = dPhi(j,k)
      dydx(6) = dTheta(j,0,k)
      dydx(7) = dTheta(j,1,k)      

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
    real(dp)              :: x
    n =1d4
    do i=0,n
        x = x_init +i*(0.d0-x_init)/n
        !write(*,*) x,x_start_rec
        if (x < x_start_rec .and. abs(c*k/(get_H_p(x)*get_dtau(x))) <= 0.1d0 .and. abs(get_dtau(x)) > 10.d0) then 
            get_tight_coupling_time = x
        end if
    end do
  end function get_tight_coupling_time

end module evolution_mod
