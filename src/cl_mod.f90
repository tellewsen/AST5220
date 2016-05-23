module cl_mod
  use healpix_types
  use evolution_mod
  use sphbess_mod
  implicit none

  real(dp), allocatable, dimension(:,:) :: S, S2
  real(dp), allocatable, dimension(:)   :: x_hires, k_hires
  real(dp),     allocatable, dimension(:)       :: z_spline
  real(dp),     allocatable, dimension(:,:)     :: j_l,j_l2
  real(dp),     allocatable, dimension(:,:,:)   :: integrand
contains

  ! Driver routine for (finally!) computing the CMB power spectrum
  subroutine compute_cls
    implicit none

    integer(i4b) :: i, j, l, l_num, x_num, n_spline
    real(dp)     :: dx, S_func, j_func, z, eta, eta0, x0, x_min, x_max, d, e
    integer(i4b), allocatable, dimension(:)       :: ls,cl_hires

    real(dp),     allocatable, dimension(:)       :: x_arg, cls, cls2, ls_dp,l_hires
    !real(dp),     pointer,     dimension(:)       :: k, x
    !real(dp),     pointer,     dimension(:,:,:,:) :: S_coeff
    real(dp),     allocatable, dimension(:,:)     :: Theta_l
    real(dp),     allocatable, dimension(:)       :: j_l_spline, j_l_spline2

    real(dp)                                      :: a,b,h,C_lint

    real(dp)           :: t1, t2, integral
    logical(lgt)       :: exist
    character(len=128) :: filename
    real(dp), allocatable, dimension(:) :: y, y2

    ! Set up which l's to compute
    l_num = 5!44
    allocate(ls(l_num))

    ls = (/2,3,4,5,6/)

    !ls = (/ 2, 3, 4, 6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, &
    !     & 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, &
    !     & 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 /)

    ! Task: Get source function from evolution_mod
    allocate(S(n_x_highres,n_k_highres))
    allocate(x_hires(n_x_highres),k_hires(n_k_highres))


    write(*,*) 'Compute hires Source function'
    call get_hires_source_function(x_hires, k_hires, S)
    
    !Test source func
    write(*,*) 'S'
    write(*,*) S(1,1),S(n_x_highres,n_k_highres)


    !Test of x and k grid.
    !write(*,*) 'x_hires'
    !write(*,*) x_hires(1),x_hires(n_x_highres)
    !write(*,*) 'k_hires'
    !write(*,*) k_hires(1),k_hires(n_k_highres)

    


    ! Task: Initialize spherical Bessel functions for each l; use 5400 sampled points between 
    !       z = 0 and 3500. Each function must be properly splined
    n_spline = 5400
    allocate(z_spline(n_spline))

    do i=1,n_spline
        z_spline(i) = 0.d0 + (i-1)*(3500.d0-0.d0)/(n_spline-1.d0)
    end do 

    !Test z_spline
    write(*,*) 'z_spline'
    write(*,*) z_spline(1), z_spline(n_spline)


    allocate(j_l(n_spline,l_num))
    allocate(j_l2(n_spline,l_num))
    !allocate(z_spline(n_spline))    ! Note: z is *not* redshift, but simply the dummy argument of j_l(z)


    !Calculate spherical bessel functions for select ls
    write(*,*) 'Compute spherical Bessel functions'
    do i =1,n_spline
        do l=1,l_num
            if(z_spline(i) > 2.d0) then
                call sphbes(ls(l),z_spline(i), j_l(i,l))
            endif
        end do
    end do

    ! Hint: It may be useful for speed to store the splined objects on disk in an unformatted
    !       Fortran (= binary) file, so that these only has to be computed once. Then, if your
    !       cache file exists, read from that; if not, generate the j_l's on the fly.

    !Spline across z for each l
    write(*,*) 'splining bessel'
    do l=1,l_num
          call spline(z_spline, j_l(:,l), yp1, ypn, j_l2(:,l))
    end do

    

    ! Overall task: Compute the C_l's for each given l

    !Compute the transfer function, Theta_l(k)
    ! For this I use trapezoidal intergration. Better methods should be implemented
    ! for better precision.
    write(*,*) 'Starting integration of Theta_l'
    allocate(Theta_l(l_num,n_k_highres))
    allocate(integrand(l_num,n_x_highres,n_k_highres))
    allocate(cls(l_num))
    allocate(cls2(l_num))

    do l =1,l_num
        write(*,*)'l =',l
        do k=1,n_k_highres
            !write(*,*)'k = ',k
            !trapezoidal method start
            a = x_hires(1)
            b = x_hires(n_x_highres)
            h = (b-a)/n_x_highres
            !write(*,*)'still working'

            
            do i=1,n_x_highres
                !write(*,*) 'i =',i
                integrand(l,i,k) = S(i,k)*j_lfunc(l,k_hires(k),x_hires(i))
            end do
            Theta_l(l,k) = 0.5d0*(integrand(l,1,k)+integrand(l,n_x_highres,k))

            do i=2,n_x_highres-1
                Theta_l(l,k) = Theta_l(l-1,k) +integrand(l,i,k)
            end do
        end do
        !trapezoidal method end

        !Integrate P(k) * (Theta_l^2 / k) over k to find un-normalized C_l's
        !trapezoidal method start
        !write(*,*)'doing c_l integration'
        a = k_hires(1)
        b = k_hires(n_k_highres)
        h = (b-a)/n_k_highres
        C_lint = 0.5d0*( (a/H_0)**(n_s-1.d0)*Theta_l(l,1)**2/a + (b/H_0)**(n_s-1.d0)*Theta_l(l,n_k_highres)**2/b)
        do k=2,n_k_highres-1
                !write(*,*) 'k=',k
                C_lint = C_lint + (k_hires(k)/H_0)**(n_s-1.d0)*Theta_l(l,k)**2/k_hires(k)
        end do

        !Store C_l in an array. Optionally output to file
        cls(l) = h*C_lint *ls(l)*(ls(l)+1)/(2.d0*pi)
        !trapezoidal method end
    end do


    ! Task: Spline C_l's found above, and output smooth C_l curve for each integer l

    !This is needed to make the spline funciton happy(it demands double precision)

    write(*,*) 'converting ls to double precision'
    allocate(ls_dp(l_num))
    do l=1,l_num
        ls_dp(l) = ls(l)
    end do

    write(*,*)'splining cls'
    call spline(ls_dp, cls, yp1, ypn, cls2)

    write(*,*)'done splining cls'

    allocate(l_hires(int(maxval(ls))))
    allocate(cl_hires(int(maxval(ls))))

    write(*,*) 'making l_hires'
    do l=1,int(maxval(ls))
        l_hires(l) = l
    end do

    !Find Cls for all ls, also those not in the original list

    write(*,*)'saving splined cls'
    do l=1,int(maxval(ls))
        cl_hires(l) = splint(ls_dp, cls, cls2, l_hires(l))
    end do

  end subroutine compute_cls
  
  function j_lfunc(l,k,x)
      implicit none
      integer(i4b),intent(in) ::l
      real(dp), intent(in)    :: x,k
      real(dp)                :: j_lfunc
      !write(*,*)'inside j_lfunc'
      j_lfunc = splint(z_spline,j_l(:,l),j_l2(:,l),-k*get_eta(x))
      !write(*,*)'j_lfunc calculated'
  end function j_lfunc


end module cl_mod
