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
    integer, parameter :: out_unit9=90
    integer, parameter :: out_unit10=100
    integer, parameter :: out_unit11=110
    integer, parameter :: out_unit12=120
    integer :: i

    !Initialize time and rec
    call initialize_time_mod
    call initialize_rec_mod

    !Output to file desired quantities here

    !This part opens all the files we want to write to
!    open (unit=out_unit1 ,file="xaeta.dat"         ,action="write",status="replace")
!    open (unit=out_unit2 ,file="omega_mbr.dat"     ,action="write",status="replace")
!    open (unit=out_unit3 ,file="omega_nulambda.dat",action="write",status="replace")
!    open (unit=out_unit4 ,file="Hz.dat"            ,action="write",status="replace")
!    open (unit=out_unit5 ,file="eta_t.dat"         ,action="write",status="replace")
    open (unit=out_unit6 ,file="X_e.dat"           ,action="write",status="replace")
    open (unit=out_unit7 ,file="n_e.dat"           ,action="write",status="replace")
 !   open (unit=out_unit8 ,file="n_etest.dat"       ,action="write",status="replace")
    open (unit=out_unit9 ,file="tau2.dat"          ,action="write",status="replace")
  !  open (unit=out_unit10,file="tau3.dat"          ,action="write",status="replace")
    !open (unit=out_unit11,file="g.dat"             ,action="write",status="replace")
    !open (unit=out_unit12,file="g2.dat"            ,action="write",status="replace")

    !Write to all the files
  !  do i=1,n_eta
   !     write (out_unit1 ,*) x_eta(i) , a_eta(i),eta(i)/Mpc
    !    write (out_unit2 ,*) Omega_mx(i) , Omega_bx(i),Omega_rx(i)
    !    write (out_unit3 ,*) Omega_nux(i) , Omega_lambdax(i)
    !    write (out_unit4 ,*) H(i)*Mpc/1.d3 , z_eta(i)
    !end do
    !do i=1,n_t
    !    write (out_unit5 ,*) x_t(i),eta_t(i)/Mpc
    !end do
    do i=1,n
        write (out_unit6 ,*) x_rec(i),z_rec(i),X_e(i)
	write (out_unit7 ,*) n_e(i),tau(i),dtau(i)
 !       write (out_unit8 ,*) x_test(i),z_test(i),n_etest(i)
        write (out_unit9 ,*) ddtau(i),tau_test(i),dtau_test(i)
!        write (out_unit9 ,*) tau2(i)!
   !     write (out_unit10,*) ddtau_test(i),g(i),g_test(i)
    !    write (out_unit11,*) dg(i),dg_test(i),ddg(i)
     !   write (out_unit12,*) ddg_test(i)
    end do

    !Close all the files
    !close (out_unit1)
!    close (out_unit2)
 !   close (out_unit3)
  !  close (out_unit4)
   ! close (out_unit5)
    close (out_unit6)
    close (out_unit7)
 !   close (out_unit8)
    close (out_unit9)
   ! close (out_unit10)
    !close (out_unit11)
    !close (out_unit12)
end program cmbspec

