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

  real(dp), allocatable, dimension(:) :: dtau
  real(dp), allocatable, dimension(:) :: ddtau
  real(dp), allocatable, dimension(:) :: H_p
  real(dp), allocatable, dimension(:) :: dH_p
  real(dp), allocatable, dimension(:),private :: eta_precomp

  ! Fourier mode list
  real(dp), allocatable, dimension(:) :: ks

  ! Book-keeping variables
  real(dp),     private :: k_current,ck_current,ckH_p
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
    real(dp)     :: k_min = 0.1d0*H_0/c
    real(dp)     :: k_max = 1000.d0*H_0/c

    !Initialize k-grid, ks; quadratic between k_min and k_max
    allocate(ks(n_k))
    do k=1,n_k
        ks(k) = k_min +(k_max -k_min)*((k-1)/100.d0)**2
    end do

    !Allocate arrays for perturbation quantities
    allocate(delta(1:n_t, n_k))
    allocate(delta_b(1:n_t, n_k))
    allocate(v(1:n_t, n_k))
    allocate(v_b(1:n_t, n_k))
    allocate(Phi(1:n_t, n_k))
    allocate(Theta(1:n_t, 0:lmax_int, n_k))
    allocate(Psi(1:n_t, n_k))

    allocate(dPhi(1:n_t, n_k))
    allocate(dPsi(1:n_t, n_k))
    allocate(dv_b(1:n_t, n_k))
    allocate(dTheta(1:n_t, 0:lmax_int, n_k))

    !Allocate arrays for precomputed variables
    allocate(dtau(n_t),H_p(n_t),dH_p(n_t))
    allocate(ddtau(n_t),eta_precomp(n_t))

    !Precompute useful variables
    do i=1,n_t
       dtau(i)  = get_dtau(x_t(i))
       ddtau(i) = get_ddtau(x_t(i))
       H_p(i)   = get_H_p(x_t(i))
       dH_p(i)  = get_dH_p(x_t(i))
       eta_precomp(i)   = get_eta(x_t(i))
    end do
    write(*,'(*(2X, ES14.6))') H_p(1),dH_p(1),ddtau(1),dtau(1),ks(1)
    ! Task: Set up initial conditions for the Boltzmann and Einstein equations
    Phi(1,:)     = 1.d0
    delta(1,:)   = 1.5d0*Phi(1,:)
    delta_b(1,:) = delta(1,:)
    Theta(1,0,:) = 0.5d0*Phi(1,:)
    do k = 1, n_k
        v(1,k)       = c*ks(k)/(2.d0*H_p(1))*Phi(1,k)
        v_b(1,k)     = v(1,k)
        Theta(1,1,k) = -c*ks(k)/(6.d0*H_p(1))*Phi(1,k)
        Theta(1,2,k) = -20.d0*c*ks(k)/(45.d0*H_p(1)*dtau(1))*Theta(1,1,k) !without polarization
        do l = 3, lmax_int
            Theta(1,l,k) = -l/(2.d0*l+1.d0)*c*ks(k)/(H_p(1)*dtau(1))*Theta(1,l-1,k)
        end do
        Psi(1,k)     = -Phi(1,k) - 12.d0*H_0**2/(ks(k)*c*a_t(1))**2*Omega_r*Theta(1,2,k)
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

    dydx(:) = 0

    ! Propagate each k-mode independently
    do k = 1, n_k
       write(*,*) 'Current k', k
       k_current = ks(k)  ! Store k_current as a global module variable
       ck_current = c*ks(k) !store c*k

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
       !write(*,*) 'x_tc =',x_tc
       !write(*,*) 'under x_tc'

       ! Task: Integrate from x_init until the end of tight coupling, using
       !       the tight coupling equations
       !write(*,*) 'Start of tight coupling'
       !write (*,'(*(2X, ES14.6))') delta(1,k), delta_b(1,k), &
       !v(1,k), v_b(1,k), Phi(1,k), Theta(1,0,k), Theta(1,1,k),Psi(1,k)
       !write (*,'(*(2X, ES14.6))') x_t(1),dv_b(1,k),dPsi(1,k),dPhi(1,k),dTheta(1,0,k),dTheta(1,1,k),dTheta(1,2,k)

       do j=2,n_t
           if (x_t(j)< x_tc) then 
               !precompute some variables
               ckH_p = ck_current/H_p(j)

               !Solve next step
               call odeint(y_tight_coupling,x_t(j-1),x_t(j),eps,h1,hmin,derivs_tc, bsstep, output3)

               !Save variables
               delta(j,k)   = y_tight_coupling(1)
               delta_b(j,k) = y_tight_coupling(2)
               v(j,k)       = y_tight_coupling(3)
               v_b(j,k)     = y_tight_coupling(4)
               Phi(j,k)     = y_tight_coupling(5)
               Theta(j,0,k) = y_tight_coupling(6)
               Theta(j,1,k) = y_tight_coupling(7)
               Theta(j,2,k) = -(20.d0*ckH_p)/(45.d0*dtau(j))*Theta(j,1,k)
               do l = 3, lmax_int
                  Theta(j,l,k) = -l/(2.d0*l+1.d0)*ckH_p/dtau(j)*Theta(j,l-1,k)
               end do	
               Psi(j,k)      = -Phi(j,k) - 12.d0*H_0**2/(ck_current*a_t(j))**2*Omega_r*Theta(j,2,k)

               !Store derivatives that are required for C_l estimation
               dPhi(j,k)     = Psi(j,k) - (ckH_p)**2/3.d0*Phi(j,k) + (H_0**2/H_p(j))**2/2.d0 &
                              *(Omega_m/a_t(j)*delta(j,k) + Omega_b/a_t(j)*delta_b(j,k) &
                              + 4.d0*Omega_r/a_t(j)**2 *Theta(j,0,k))

               dPsi(j,k)     = -dPhi(j,k) - 12.d0*H_0**2/(ck_current*a_t(j))**2*Omega_r*(-2.d0*Theta(j,2,k)+dTheta(j,2,k))

               dTheta(j,0,k) = -ckH_p*Theta(j,1,k) - dPhi(j,k)

               R             = 4.d0*Omega_r/(3.d0*Omega_b*a_t(j))

               q             = (-((1.d0-2.d0*R)*dtau(j) + &
                              (1.d0+R)*ddtau(j))*(3.d0*Theta(j,1,k)+v_b(j,k)) - &
                              ckH_p*Psi(j,k) +&
                              (1.d0-dH_p(j)/H_p(j))*ckH_p*(-Theta(j,0,k) + 2.d0*Theta(j,2,k))-&
                              ckH_p*dTheta(j,0,k))/((1.d0+R)*dtau(j)+dH_p(j)/H_p(j) -1.d0)
 
               dv_b(j,k)     = 1.d0/(1.d0+R)*(-v_b(j,k)-ckH_p*Psi(j,k)+&
                              R*(q+ckH_p*(2.d0*Theta(j,2,k)-Theta(j,0,k))-&
                              ckH_p*Psi(j,k)))

               dTheta(j,1,k) = 1.d0/3.d0*(q-dv_b(j,k))
               dTheta(j,2,k) = 0
               do l = 3, lmax_int
                   dTheta(j,l,k) = 0 
               end do	

               !write (*,'(*(2X, ES14.6))') delta(j,k), delta_b(j,k), v(j,k), v_b(j,k), Phi(j,k),  Theta(j,0,k), Theta(j,1,k),Psi(j,k)
               !write (*,'(*(2X, ES14.6))') x_t(j),dPsi(j,k),dPhi(j,k),dv_b(j,k),dTheta(j,0,k),dTheta(j,1,k),dTheta(j,2,k)


           else
               j_tc = j
               exit
           end if
       end do
       !write(*,*) 'End of tight coupling'

       ! Task: Set up variables for integration from the end of tight coupling 
       ! until today
       y(1:7) = y_tight_coupling(1:7)
       y(8)   = Theta(1,2,k)
       do l = 3, lmax_int
          y(6+l) = Theta(1,l,k)
       end do


       !Continue after tight coupling
       !write(*,*) 'start of rec'       
       do j = j_tc, n_t

          !Precompute some variables
          ckH_p = ck_current/H_p(j)

          !Integrate equations from tight coupling to today
          !write(*,*) 'running odeint with j =', j
          call odeint(y, x_t(j-1) ,x_t(j), eps, h1, hmin, derivs, bsstep, output3)

          ! Task: Store variables at time step i in global variables
          delta(j,k)   = y(1)
          delta_b(j,k) = y(2)
          v(j,k)       = y(3)
          v_b(j,k)     = y(4)
          Phi(j,k)     = y(5)
          
          do l = 0, lmax_int
             Theta(j,l,k) = y(6+l)
          end do
          Psi(j,k)     =  - Phi(j,k) - 12.d0*H_0**2/(ck_current*a_t(j))**2*Omega_r*Theta(j,2,k)

          ! Task: Store derivatives that are required for C_l estimation
          dPhi(j,k)     = Psi(j,k) -c**2*k_current**2/(3.d0*H_p(j)**2)*Phi(j,k) +H_0**2/(2.d0*H_p(j)) &
                          *(Omega_m/a_t(j)*delta(j,k) +Omega_b/a_t(j)*delta_b(j,k) + 4.d0*Omega_r/a_t(j)**2 &
                          *Theta(j,0,k))

          dv_b(j,k)     = -v_b(j,k) -ckH_p*Psi(j,k) +dtau(j)*R*(3.d0*Theta(j,1,k)+ v_b(j,k))

          dTheta(j,0,k) = -ckH_p*Theta(j,1,k) -dPhi(j,k)
          dTheta(j,1,k) = ckH_p/3.d0*Theta(j,0,k) - &
                          2.d0*ckH_p/3.d0*Theta(j,2,k)+&
                          ckH_p/3.d0*Psi(j,k) + &
                          dtau(j)*(Theta(j,1,k)+ 1.d0/3.d0*v_b(j,k))
          dTheta(j,2,k) = 2.d0*ckH_p/5.d0*Theta(j,1,k) -&
                          3.d0*ckH_p/5.d0*Theta(j,3,k)+&
                          dtau(j)*0.9d0*Theta(j,2,k)
          do l=3,lmax_int-1
              dTheta(j,l,k) = l*ckH_p/(2.d0*l+1.d0)*Theta(j,l-1,k) -&
                              (l+1.d0)*ckH_p/(2.d0*l+1.d0)*Theta(j,l+1,k)+&
                              dtau(j)*Theta(j,l,k)
          end do
          dTheta(j,lmax_int,k) = ckH_p*Theta(j,l-1,k) -&
                                 c*(l+1.d0)/(H_p(j)*eta_precomp(j))*&
                                 Theta(j,l,k) + dtau(j)*Theta(j,l,k)
          dPsi(j,k)     = -dPhi(j,k) - 12.d0*H_0**2/(ck_current*a_t(j))**2*Omega_r*(-2.d0*Theta(j,2,k)+dTheta(j,2,k))
       end do
       !write(*,*) 'today'
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

      real(dp) :: delta,delta_b,v,v_b,Phi,Theta0,Theta1,Theta2
      real(dp) :: Psi,dPhi,dTheta0,dv_b,dTheta1

      delta   = y_tc(1)
      delta_b = y_tc(2)
      v       = y_tc(3)
      v_b     = y_tc(4)
      Phi     = y_tc(5)
      Theta0  = y_tc(6)
      Theta1  = y_tc(7)

      Theta2    = -20.d0*ckH_p/(45.d0*dtau(j))*Theta1

      R         = (4.d0*Omega_r)/(3.d0*Omega_b*a_t(j))

      Psi       = -Phi - 12.d0*(H_0/ck_current/a_t(j))**2.d0*Omega_r*Theta2

      dPhi      = Psi - ckH_p**2/3.d0*Phi + (H_0/H_p(j))**2/2.d0*(Omega_m/a_t(j)*delta + Omega_b/a_t(j)*delta_b + 4.d0*Omega_r/a_t(j)**2*Theta0)

      dTheta0   = -ckH_p*Theta1 - dPhi

      d_delta   = ckH_p*v   - 3.d0*dPhi

      d_delta_b = ckH_p*v_b - 3.d0*dPhi

      d_v       = -v -ckH_p*Psi

      q         = ( -((1.d0-2.d0*R)*dtau(j) + (1.d0+R)*ddtau(j)) *(3.d0*Theta1+v_b) - ckH_p*Psi +(1.d0-dH_p(j)/H_p(j))*ckH_p*(-Theta0 + 2.d0*Theta2) - ckH_p*dTheta0) / ((1.d0+R)*dtau(j)+dH_p(j)/H_p(j) -1.d0)

      dv_b      = (1.d0/(1.d0+R)) *(-v_b - ckH_p*Psi + R*(q+ckH_p*(-Theta0 + 2.d0*Theta2)-ckH_p*Psi))

      dTheta1   = (1.d0/3.d0)*(q-dv_b)

      dydx(1) = d_delta
      dydx(2) = d_delta_b
      dydx(3) = d_v
      dydx(4) = dv_b
      dydx(5) = dPhi
      dydx(6) = dTheta0
      dydx(7) = dTheta1
      !write(*,*) 'dydx(1) =',dydx(1)
      !write(*,*) 'dydx(2) =',dydx(2)
  end subroutine derivs_tc

  subroutine derivs(x,y, dydx) 
      use healpix_types
      implicit none
      real(dp),               intent(in)  :: x
      real(dp), dimension(:), intent(in)  :: y
      real(dp), dimension(:), intent(out) :: dydx

      real(dp) :: d_delta
      real(dp) :: d_delta_b
      real(dp) :: d_v
      real(dp) :: q,R
      integer(i4b) :: i
      real(dp) :: delta,delta_b,v,v_b,Phi,Theta0,Theta1,Theta2,Theta3,Theta4,Theta5,Theta6
      real(dp) :: Psi,dPhi,dTheta0,dv_b,dTheta1,dTheta2

      delta   = y(1)
      delta_b = y(2)
      v       = y(3)
      v_b     = y(4)
      Phi     = y(5)
      Theta0  = y(6)
      Theta1  = y(7)
      Theta2  = y(8)
      Theta3  = y(9)
      Theta4  = y(10)
      Theta5  = y(11)
      Theta6  = y(12)

      R         = (4.d0*Omega_r)/(3.d0*Omega_b*a_t(j))

      Psi       = -Phi - 12.d0*(H_0/ck_current/a_t(j))**2.d0*Omega_r*Theta2

      dPhi      = Psi - ckH_p**2/3.d0*Phi + (H_0/H_p(j))**2/2.d0*(Omega_m/a_t(j)*delta + Omega_b/a_t(j)*delta_b + 4.d0*Omega_r/a_t(j)**2*Theta0)

      dTheta0   = -ckH_p*Theta1 - dPhi

      d_delta   = ckH_p*v   - 3.d0*dPhi

      d_delta_b = ckH_p*v_b - 3.d0*dPhi

      d_v       = -v -ckH_p*Psi

      dv_b      = -v_b -ckH_p*Psi +dtau(j)*R*(3.d0*Theta1+v_b)

      dTheta1   = ckH_p/3.d0*Theta0 -2.d0/3.d0*ckH_p*Theta2 +ckH_p/3.d0*Psi +dtau(j)*(Theta1+v_b/3.d0)
      dTheta2   = l/(2.d0*l+1)*ckH_p*Theta1 - (l+1.d0)/(2.d0*l+1.d0)*ckH_p*Theta3+dtau(j)*0.9d0*Theta2

      do i=3,lmax_int-1
          dydx(6+i) = l/(2.d0*l+1)*ckH_p*y(5+i) - (l+1.d0)/(2.d0*l+1.d0)*ckH_p*y(7+i) +dtau(j)*y(6+i)
      end do

      dydx(12) = ckH_p*Theta5 -c*(l+1.d0)/H_p(j)/eta_precomp(j)*Theta6 +dtau(j)*Theta6

      dydx(1) = d_delta
      dydx(2) = d_delta_b
      dydx(3) = d_v
      dydx(4) = dv_b
      dydx(5) = dPhi
      dydx(6) = dTheta0
      dydx(7) = dTheta1
      dydx(8) = dTheta2
      !write(*,*) 'dydx(1) =',dydx(1)
      !write(*,*) 'dydx(2) =',dydx(2)
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
