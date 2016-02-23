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
  integer :: i

  ! Initialize time grids
  call initialize_time_mod
  call initialize_rec_mod

  ! Output to file desired quantities here

  !This part opens a file results.txt and prints x,a,eta for each step on one line
  !at a time.
  open (unit=out_unit1,file="xaeta.dat",action="write",status="replace")
  open (unit=out_unit2,file="omega_mbr.dat",action="write",status="replace")
  open (unit=out_unit3,file="omega_nulambda.dat",action="write",status="replace")
  open (unit=out_unit4,file="Hz.dat",action="write",status="replace")
  open (unit=out_unit5,file="eta_t.dat",action="write",status="replace")

  do i=1,n_eta+1
      write (out_unit1,*) x_eta(i) , a_eta(i),eta(i)/Mpc
      write (out_unit2,*) Omega_mx(i) , Omega_bx(i),Omega_rx(i)
      write (out_unit3,*) Omega_nux(i) , Omega_lambdax(i)
      write (out_unit4,*) H(i)*Mpc/1.d3 , z_eta(i)
  end do
  do i=1,n_t
      write (out_unit5,*) x_t(i),eta_t(i)/Mpc
  end do
  close (out_unit1)
  close (out_unit2)
  close (out_unit3)
  close (out_unit4)
  close (out_unit5)

end program cmbspec

