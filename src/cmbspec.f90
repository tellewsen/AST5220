program cmbspec
    use healpix_types
    use params
    use rec_mod
    use time_mod
  
    implicit none

    integer, parameter :: out_unit1=10
    integer, parameter :: out_unit2=20
    integer, parameter :: out_unit3=30
    integer, parameter :: out_unit4=40
    integer, parameter :: out_unit5=50
    integer, parameter :: out_unit6=60
    integer, parameter :: out_unit7=70
    integer, parameter :: out_unit8=80
    integer :: i

    !Initialize time and rec
    call initialize_time_mod
    call initialize_rec_mod

    !Output to file desired quantities here

    !This part opens all the files we want to write to
    open (unit=out_unit1,file="xaeta.dat"         ,action="write",status="replace")
    open (unit=out_unit2,file="omega_mbr.dat"     ,action="write",status="replace")
    open (unit=out_unit3,file="omega_nulambda.dat",action="write",status="replace")
    open (unit=out_unit4,file="Hz.dat"            ,action="write",status="replace")
    open (unit=out_unit5,file="eta_t.dat"         ,action="write",status="replace")
    open (unit=out_unit6,file="X_e.dat"           ,action="write",status="replace")
    open (unit=out_unit7,file="n_e.dat"           ,action="write",status="replace")
    open (unit=out_unit8,file="n_etest.dat"       ,action="write",status="replace")
    !Write to all the files
    do i=1,n_eta+1
        write (out_unit1,*) x_eta(i) , a_eta(i),eta(i)/Mpc
        write (out_unit2,*) Omega_mx(i) , Omega_bx(i),Omega_rx(i)
        write (out_unit3,*) Omega_nux(i) , Omega_lambdax(i)
        write (out_unit4,*) H(i)*Mpc/1.d3 , z_eta(i)
    end do
    do i=1,n_t
        write (out_unit5,*) x_t(i),eta_t(i)/Mpc
    end do
    do i=1,n
        write (out_unit6,*) x_rec(i),z_rec(i),X_e(i)
	write (out_unit7,*) n_e(i)
        write (out_unit8,*) x_test(i),z_test(i),n_etest(i)
    end do
    !Close all the files
    close (out_unit1)
    close (out_unit2)
    close (out_unit3)
    close (out_unit4)
    close (out_unit5)
    close (out_unit6)

end program cmbspec

