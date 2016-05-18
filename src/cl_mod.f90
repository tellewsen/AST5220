module cl_mod
  use healpix_types
  use evolution_mod
  use sphbess_mod
  implicit none

  real(dp), allocatable, dimension(:,:) :: S, S2
  real(dp), allocatable, dimension(:)   :: x_hires, k_hires
contains

  ! Driver routine for (finally!) computing the CMB power spectrum
  subroutine compute_cls
    implicit none

    integer(i4b) :: i, j, l, l_num, x_num, n_spline
    real(dp)     :: dx, S_func, j_func, z, eta, eta0, x0, x_min, x_max, d, e
    integer(i4b), allocatable, dimension(:)       :: ls
    real(dp),     allocatable, dimension(:)       :: integrand,zs
    real(dp),     allocatable, dimension(:,:)     :: j_l, j_l2
    real(dp),     allocatable, dimension(:)       :: x_arg, int_arg, cls, cls2, ls_dp
    !real(dp),     pointer,     dimension(:)       :: k, x
    !real(dp),     pointer,     dimension(:,:,:,:) :: S_coeff
    real(dp),     allocatable, dimension(:,:)     :: Theta_l
    real(dp),     allocatable, dimension(:)       :: z_spline, j_l_spline, j_l_spline2

    real(dp)           :: t1, t2, integral
    logical(lgt)       :: exist
    character(len=128) :: filename
    real(dp), allocatable, dimension(:) :: y, y2

    ! Set up which l's to compute
    l_num = 44
    allocate(ls(l_num))
    ls = (/ 2, 3, 4, 6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, &
         & 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, &
         & 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 /)

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

    !Test zs
    !write(*,*) 'zs'
    !write(*,*) zs(1), zs(n_spline)


    allocate(j_l(n_spline,l_num))
    allocate(j_l2(n_spline,l_num))
    !allocate(z_spline(n_spline))    ! Note: z is *not* redshift, but simply the dummy argument of j_l(z)


    !Calculate spherical bessel functions for select ls
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
    do l=1,l_num
          call spline(z_spline, j_l(:,l), yp1, ypn, j_l2(:,l))
    end do


    ! Overall task: Compute the C_l's for each given l
       !Compute the transfer function, Theta_l(k)
       ! For this I use trapezoidal intergration. Better methods should be implemented
       ! for better precision.
    allocate(Theta_l(l_num,n_k_highres))
    do l = 1, l_num
       Theta_l()

       do i=1,n_intergration
           Theta_l()
       end do

       ! Task: Integrate P(k) * (Theta_l^2 / k) over k to find un-normalized C_l's


       ! Task: Store C_l in an array. Optionally output to file

    end do


    ! Task: Spline C_l's found above, and output smooth C_l curve for each integer l


  end subroutine compute_cls
  
end module cl_mod
