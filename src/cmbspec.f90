program cmbspec
    use healpix_types
    use params
    use rec_mod
    use time_mod
    use evolution_mod
    use cl_mod

    implicit none

    !Output unit
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
    integer, parameter :: out_unit13=130
    integer, parameter :: out_unit14=140
    integer, parameter :: out_unit15=150
    integer, parameter :: out_unit16=160
    integer, parameter :: out_unit17=170
    integer, parameter :: out_unit18=180
    integer, parameter :: out_unit19=190
    integer, parameter :: out_unit20=200
    integer, parameter :: out_unit21=210
    integer, parameter :: out_unit22=220
    integer, parameter :: out_unit23=230
    integer, parameter :: out_unit24=240

    integer :: i

    !Time keeping variables
    real(dp) :: start_time,end_time

    !Start timer
    call cpu_time(start_time)

    !Initialize time_mod and save to file
    write(*,*) 'initialize_time_mod'
    call initialize_time_mod

    !Opens all the files we want to write to
    open (unit=out_unit1 ,file="xaeta.dat"         ,action="write",status="replace")
    open (unit=out_unit2 ,file="omega_mbr.dat"     ,action="write",status="replace")
    open (unit=out_unit3 ,file="omega_nulambda.dat",action="write",status="replace")
    open (unit=out_unit4 ,file="Hz.dat"            ,action="write",status="replace")
    open (unit=out_unit5 ,file="eta_t.dat"         ,action="write",status="replace")

    !write to file
    do i=1,n_eta
        write (out_unit1 ,*) x_eta(i) , a_eta(i),eta(i)/Mpc
        write (out_unit2 ,*) Omega_mx(i) , Omega_bx(i),Omega_rx(i)
        write (out_unit3 ,*) Omega_nux(i) , Omega_lambdax(i)
        write (out_unit4 ,*) H(i)*Mpc/1.d3 , z_eta(i)
    end do
    
    do i=1,n_t
        write (out_unit5 ,*) x_t(i),eta_t(i)/Mpc
    end do

    !Close the files
    close (out_unit1)
    close (out_unit2)
    close (out_unit3)
    close (out_unit4)
    close (out_unit5)

    !Initialize and save values from rec_mod
    write(*,*) 'initialize_rec_mod'
    call initialize_rec_mod
    !Open all the files
    open (unit=out_unit6 ,file="X_e.dat"           ,action="write",status="replace")
    open (unit=out_unit7 ,file="n_e.dat"           ,action="write",status="replace")
    open (unit=out_unit8 ,file="n_etest.dat"       ,action="write",status="replace")
    open (unit=out_unit9 ,file="tau2.dat"          ,action="write",status="replace")
    open (unit=out_unit10,file="tau3.dat"          ,action="write",status="replace")
    open (unit=out_unit11,file="g.dat"             ,action="write",status="replace")
    open (unit=out_unit12,file="g2.dat"            ,action="write",status="replace")

    !Write to all the files
    do i=1,n
        write (out_unit6 ,*) x_rec(i),z_rec(i),X_e(i)
	write (out_unit7 ,*) n_e(i),tau_rec(i),dtau_rec(i)
        write (out_unit8 ,*) x_test(i),z_test(i),n_etest(i)
        write (out_unit9 ,*) ddtau_rec(i),tau_test(i),dtau_test(i)
        write (out_unit10,*) ddtau_test(i),g(i),g_test(i)
        write (out_unit11,*) dg(i),dg_test(i),ddg(i)
        write (out_unit12,*) ddg_test(i)
    end do

    !Close all the files
    close (out_unit6)
    close (out_unit7)
    close (out_unit8)
    close (out_unit9)
    close (out_unit10)
    close (out_unit11)
    close (out_unit12)

    !Intialize and save to file for evolution_mod
    write(*,*) 'initialize_perturbation_eqns'
    call initialize_perturbation_eqns
    write(*,*) 'integrate_perturbation_eqns'
    call integrate_perturbation_eqns

    write(*,*) 'Saving perturbations to file'
    open (unit=out_unit13,file="delta.dat",action="write",status="replace")
    open (unit=out_unit14,file="delta_b.dat",action="write",status="replace")
    open (unit=out_unit15,file="v.dat",action="write",status="replace")
    open (unit=out_unit16,file="v_b.dat",action="write",status="replace")
    open (unit=out_unit17,file="Phi.dat",action="write",status="replace")
    open (unit=out_unit18,file="Psi.dat",action="write",status="replace")
    open (unit=out_unit19,file="Theta0.dat",action="write",status="replace")
    open (unit=out_unit20,file="dPhi.dat",action="write",status="replace")
    open (unit=out_unit21,file="dPsi.dat",action="write",status="replace")

    do i=1,n_t            
        write (out_unit13,'(*(2X, ES14.6))') delta(i,1),delta(i,5),delta(i,10),delta(i,40),delta(i,60),delta(i,100)
        write (out_unit14,'(*(2X, ES14.6))') delta_b(i,1),delta_b(i,5),delta_b(i,10),delta_b(i,40),delta_b(i,60),delta_b(i,100)
        write (out_unit15,'(*(2X, ES14.6))') v(i,1),v(i,5),v(i,10),v(i,40),v(i,60),v(i,100)
        write (out_unit16,'(*(2X, ES14.6))') v_b(i,1),v_b(i,5),v_b(i,10),v_b(i,40),v_b(i,60),v_b(i,100)
        write (out_unit17,'(*(2X, ES14.6))') Phi(i,1),Phi(i,5),Phi(i,10),Phi(i,40),Phi(i,60),Phi(i,100)
        write (out_unit18,'(*(2X, ES14.6))') Psi(i,1),Psi(i,5),Psi(i,10),Psi(i,40),Psi(i,60),Psi(i,100)
        write (out_unit19,'(*(2X, ES14.6))') Theta(i,0,1),Theta(i,0,5),Theta(i,0,10),Theta(i,0,40),Theta(i,0,60),Theta(i,0,100)
        write (out_unit20,'(*(2X, ES14.6))') dPhi(i,1),dPhi(i,5),dPhi(i,10),dPhi(i,40),dPhi(i,60),dPhi(i,100)
        write (out_unit21,'(*(2X, ES14.6))') dPsi(i,1),dPsi(i,5),dPsi(i,10),dPsi(i,40),dPsi(i,60),dPsi(i,100)
    end do

    open (unit=123,file="dTheta2.dat",action="write",status="replace")
    open (unit=124,file="Theta2.dat",action="write",status="replace")
    do i=1,n_t            
        write (123,'(*(2X, ES14.6))') dTheta(i,2,1),dTheta(i,2,5),dTheta(i,2,10),dTheta(i,2,40),dTheta(i,2,60),dTheta(i,2,100)
        write (124,'(*(2X, ES14.6))') Theta(i,2,1),Theta(i,2,5),Theta(i,2,10),Theta(i,2,40),Theta(i,2,60),Theta(i,2,100)
    end do
    close (123)
    close (124)



    close (out_unit13)
    close (out_unit14)
    close (out_unit15)
    close (out_unit16)
    close (out_unit17)
    close (out_unit18)
    close (out_unit19)
    close (out_unit20)
    close (out_unit21)



    !Milestone 4 starts here

    !Initialize and compute the C_l
    write(*,*) 'initialize cl_mod'
    call compute_cls

    !Writing source func to file
    write(*,*) 'Writing source func to file'
    open (unit=out_unit22, file="Source.dat", action="write", status="replace")
    do i = 1,n_x_highres
        write (out_unit22,'(*(2X, ES14.6E3))') x_hires(i),S(i,50),S(i,250),S(i,500),S(i,2000),S(i,3000),S(i,5000)
    end do
    close (out_unit22)

    !write low res source for testing
    open (unit=out_unit23, file="S_low.dat", action="write", status="replace")
    do i=1,n_t
        write (out_unit23,'(*(2X, ES14.6E3))') S_lores(i,1),S_lores(i,5),S_lores(i,10),S_lores(i,40),S_lores(i,60),S_lores(i,100)   
    end do
    close (out_unit23)

    !write Theta_l for six different ks
    open (unit=12, file="Theta_l.dat", action="write", status="replace")
    open (unit=13, file="ls.dat", action="write", status="replace")
    do l=1,44
        write (12,'(*(2X, ES14.6E3))') Theta_l(l,50),Theta_l(l,250),Theta_l(l,500),Theta_l(l,2000),Theta_l(l,3000),Theta_l(l,5000)
        write (13,*) ls(l)
    end do
    close (12)
    close (13)

    !write cls and ls to file
    open (unit=out_unit24, file="C_l.dat", action="write", status="replace")
    do i=1,1200
        write (out_unit24,'(*(2X, ES14.6E3))') l_hires(i),cl_hires(i)
    end do
    close (out_unit24)
     

    !Print time used
    call cpu_time(end_time)
    print'("Time used = ",f7.2," seconds.")',end_time-start_time


end program cmbspec



