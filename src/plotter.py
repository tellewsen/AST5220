import matplotlib.pyplot as plt
from numpy import *
from matplotlib import rc 
rc('font',**{'family':'serif'}) # This is for Latex writing

#Import datafiles
#Relative densities
nu,lamda                 = loadtxt("omega_nulambda.dat",unpack=True)
m,b,r                    = loadtxt("omega_mbr.dat",unpack=True)
total                    = nu+lamda+m+b+r

#x,scalefactor,conformal time
x,a,eta                  = loadtxt("xaeta.dat",unpack=True)
#Hubbleparam,redshift
H,z                      = loadtxt("Hz.dat",unpack=True)

x_rec,z_rec,X_e          = loadtxt("X_e.dat",unpack=True)
n_e,tau,dtau             = loadtxt("n_e.dat",unpack=True)
x_test,z_test,n_etest    = loadtxt("n_etest.dat",unpack=True)

ddtau,tau_test,dtau_test = loadtxt("tau2.dat",unpack=True)
ddtau_test,g,g_test      = loadtxt("tau3.dat",unpack=True)
dg,dg_test,ddg           = loadtxt("g.dat",unpack=True)
ddg_test                 = loadtxt("g2.dat",unpack=True)


#Test x,test eta
x_t,eta_t                = loadtxt("eta_t.dat",unpack=True)

"""
plt.figure(0)
plt.plot(x,m,label = r'$\Omega_m$')
plt.plot(x,b,label = r'$\Omega_b$')
plt.plot(x,r,label = r'$\Omega_r$')
plt.plot(x,nu,label = r'$\Omega_\nu$')
plt.plot(x,lamda,label = r'$\Omega_\lambda$')
plt.plot(x,total,label = r'Total')
plt.ylim([0,1.1])
plt.xlim([min(x),max(x)])
plt.legend(loc='best')
plt.xlabel(r'x')
plt.ylabel(r'$\Omega_x$')


plt.figure(1)
plt.plot(x,eta,label = r'$\eta(x)$')
plt.plot(x_t,eta_t,label=r'$\eta_2(x)$')
plt.yscale('log')
plt.xlim([min(x),max(x)])
plt.xlabel(r'x')
plt.ylabel(r'$\eta$(x)[Mpc]')
plt.legend(loc='best')

plt.figure(2)
plt.plot(x,H,label = r'H(x)')
plt.yscale('log')
plt.xlim([min(x),max(x)])
plt.xlabel(r'x')
plt.ylabel(r'H(x)[km s$^{-1}$Mpc$^{-1}$]')


plt.figure(3)
plt.plot(z,H,label = r'H(z)')
plt.yscale('log')
plt.xscale('log')
plt.xlim([z[0],z[-1]])
plt.xlabel(r'z')
plt.ylabel(r'H(z)[km s$^{-1}$Mpc$^{-1}$]')


plt.figure(4)
plt.plot(z_rec,X_e,label = r'$X_e$')
plt.yscale('log')
#plt.xscale('log')
plt.xlim([1800,0])
plt.xlabel(r'z')
plt.ylabel(r'$X_e$')


plt.figure(5)
plt.plot(z_rec ,n_e,'-',label = r'$n_{e}$')
plt.plot(z_test,n_etest,'-',label = r'Splined')
plt.yscale('log')
#plt.xscale('log')
plt.xlim([z_rec[0],z_rec[-1]])
plt.xlabel(r'z')
plt.ylabel(r'$n_e$[m$^{-3}$]')
plt.legend()

plt.figure(6)
plt.plot(x_rec , tau,            '-',label = r'$\tau(x)$')
#plt.plot(x_test, tau_test,       '-' ,label = r'$\tau_{test}(x)$')
plt.plot(x_rec , abs(dtau),      '--',label = r'$|\tau^\prime(x)|$')
#plt.plot(x_test, abs(dtau_test), '-' ,label = r'$|\tau^{\prime}_{test}(x)|$')
plt.plot(x_rec , abs(ddtau),     '-.',label = r'$|\tau^{\prime\prime}(x)|$')
#plt.plot(x_test, abs(ddtau_test),'-' ,label = r'$|\tau^{\prime\prime}_{test}(x)|$')
plt.yscale('log')
#plt.xscale('log')
plt.xlim([x_rec[0],x_rec[-1]])
plt.xlabel(r'x')
plt.ylabel(r'$\tau$,$|\tau^\prime|$,$|\tau^{\prime\prime}|$')
plt.legend()


plt.figure(7)
plt.plot(x_rec , g,            '-',label = r'$\tilde g(x)$')
#plt.plot(x_test, g_test,       '-' ,label = r'$\tilde g_{test}(x)$')
plt.plot(x_rec , dg/10.,       '--',label = r'$\tilde g^\prime(x)$')
#plt.plot(x_test, dg_test/10.,  '-' ,label = r'$\tilde g^\prime_{test}(x)$')
plt.plot(x_rec , ddg/300.,     '-.',label = r'$\tilde g^{\prime\prime}(x)$')
#plt.plot(x_test, ddg_test/300.,'-' ,label = r'$\tilde g^{\prime\prime}_{test}(x)$')
#plt.yscale('log')
#plt.xscale('log')
#plt.xlim([-7.5,-6])
plt.xlabel(r'x')
plt.ylabel(r'$\tilde g$,$\tilde g^\prime/10$,$\tilde g^{\prime\prime}/300$')
plt.legend()
plt.show()
"""

#Milestone 3 starts here
"""
delta  = loadtxt("delta.dat",unpack=True)
deltab = loadtxt("delta_b.dat",unpack=True)
Phi    = loadtxt("Phi.dat",unpack=True)
dPhi   = loadtxt("dPhi.dat",unpack=True)
Psi    = loadtxt("Psi.dat",unpack=True)
dPsi   = loadtxt("dPsi.dat",unpack=True)
Theta0 = loadtxt("Theta0.dat",unpack=True)
v      = loadtxt("v.dat",unpack=True)
vb     = loadtxt("v_b.dat",unpack=True)
"""
"""
plt.figure(9)
plt.plot(x_t,delta[0],label = r'$\delta_{1}$')
plt.plot(x_t,delta[1],label = r'$\delta_{5}$')
plt.plot(x_t,delta[2],label = r'$\delta_{10}$')
plt.plot(x_t,delta[3],label = r'$\delta_{40}$')
plt.plot(x_t,delta[4],label = r'$\delta_{60}$')
plt.plot(x_t,delta[5],label = r'$\delta_{100}$')
plt.xlim([min(x_t),max(x_t)])
plt.legend(loc='best')
plt.yscale('symlog')
plt.xlabel(r'x')
plt.ylabel(r'$\delta_{k}$')


plt.figure(10)
plt.plot(x_t,deltab[0],label = r'$\delta_{b,1}$')
plt.plot(x_t,deltab[1],label = r'$\delta_{b,5}$')
plt.plot(x_t,deltab[2],label = r'$\delta_{b,10}$')
plt.plot(x_t,deltab[3],label = r'$\delta_{b,40}$')
plt.plot(x_t,deltab[4],label = r'$\delta_{b,60}$')
plt.plot(x_t,deltab[5],label = r'$\delta_{b,100}$')
plt.xlim([min(x_t),max(x_t)])
plt.legend(loc='best')
plt.yscale('symlog')
plt.xlabel(r'x')
plt.ylabel(r'$\delta_{b,k}$')

plt.figure(11)
plt.plot(x_t,vb[0],label = r'$v_{b,1}$')
plt.plot(x_t,vb[1],label = r'$v_{b,5}$')
plt.plot(x_t,vb[2],label = r'$v_{b,10}$')
plt.plot(x_t,vb[3],label = r'$v_{b,40}$')
plt.plot(x_t,vb[4],label = r'$v_{b,60}$')
plt.plot(x_t,vb[5],label = r'$v_{b,100}$')
plt.xlim([min(x_t),max(x_t)])
plt.legend(loc='best')
plt.yscale('symlog')
plt.xlabel(r'x')
plt.ylabel(r'$v_{b,k}$')

plt.figure(12)
plt.plot(x_t,v[0],label = r'$v_{1}$')
plt.plot(x_t,v[1],label = r'$v_{5}$')
plt.plot(x_t,v[2],label = r'$v_{10}$')
plt.plot(x_t,v[3],label = r'$v_{40}$')
plt.plot(x_t,v[4],label = r'$v_{60}$')
plt.plot(x_t,v[5],label = r'$v_{100}$')
plt.xlim([min(x_t),max(x_t)])
plt.legend(loc='best')
plt.yscale('symlog')
plt.xlabel(r'x')
plt.ylabel(r'$v_{k}$')

plt.figure(15)
plt.plot(x_t,Theta0[0],label = r'$\Theta_{0,1}$')
plt.plot(x_t,Theta0[1],label = r'$\Theta_{0,5}$')
plt.plot(x_t,Theta0[2],label = r'$\Theta_{0,10}$')
plt.plot(x_t,Theta0[3],label = r'$\Theta_{0,40}$')
plt.plot(x_t,Theta0[4],label = r'$\Theta_{0,60}$')
plt.plot(x_t,Theta0[5],label = r'$\Theta_{0,100}$')
plt.xlim([min(x_t),max(x_t)])
plt.legend(loc='best')
plt.yscale('symlog')
plt.xlabel(r'x')
plt.ylabel(r'$\Theta_{0,k}$')
"""
"""
dTheta2 = loadtxt("dTheta2.dat",unpack=True)
Theta2 = loadtxt("Theta2.dat",unpack=True)
plt.figure(22)
plt.plot(x_t,dTheta2[0],label = r'$d\Theta_{2,1}/dx$')
plt.plot(x_t,dTheta2[1],label = r'$d\Theta_{2,5}/dx$')
plt.plot(x_t,dTheta2[2],label = r'$d\Theta_{2,10}/dx$')
plt.plot(x_t,dTheta2[3],label = r'$d\Theta_{2,40}/dx$')
plt.plot(x_t,dTheta2[4],label = r'$d\Theta_{2,60}/dx$')
plt.plot(x_t,dTheta2[5],label = r'$d\Theta_{2,100}/dx$')
plt.xlim([min(x_t),max(x_t)])
plt.legend(loc='best')
#plt.yscale('symlog')
plt.xlabel(r'x')
plt.ylabel(r'$d\Theta_{2,k}/dx$')

plt.figure(23)
plt.plot(x_t,Theta2[0],label = r'$\Theta_{2,1}$')
plt.plot(x_t,Theta2[1],label = r'$\Theta_{2,5}$')
plt.plot(x_t,Theta2[2],label = r'$\Theta_{2,10}$')
plt.plot(x_t,Theta2[3],label = r'$\Theta_{2,40}$')
plt.plot(x_t,Theta2[4],label = r'$\Theta_{2,60}$')
plt.plot(x_t,Theta2[5],label = r'$\Theta_{2,100}$')
plt.xlim([min(x_t),max(x_t)])
plt.legend(loc='best')
#plt.yscale('symlog')
plt.xlabel(r'x')
plt.ylabel(r'$\Theta_{2,k}$')
plt.show()
"""
"""
plt.figure(13)
plt.plot(x_t,Phi[0],label = r'$\Phi_{1}$')
plt.plot(x_t,Phi[1],label = r'$\Phi_{5}$')
plt.plot(x_t,Phi[2],label = r'$\Phi_{10}$')
plt.plot(x_t,Phi[3],label = r'$\Phi_{40}$')
plt.plot(x_t,Phi[4],label = r'$\Phi_{60}$')
plt.plot(x_t,Phi[5],label = r'$\Phi_{100}$')
plt.xlim([min(x_t),max(x_t)])
plt.legend(loc='best')
#plt.yscale('symlog')
plt.xlabel(r'x')
plt.ylabel(r'$\phi_{k}$')

plt.figure(14)
plt.plot(x_t,Psi[0],label = r'$\Psi_{1}$')
plt.plot(x_t,Psi[1],label = r'$\Psi_{5}$')
plt.plot(x_t,Psi[2],label = r'$\Psi_{10}$')
plt.plot(x_t,Psi[3],label = r'$\Psi_{40}$')
plt.plot(x_t,Psi[4],label = r'$\Psi_{60}$')
plt.plot(x_t,Psi[5],label = r'$\Psi_{100}$')
plt.xlim([min(x_t),max(x_t)])
plt.legend(loc='best')
#plt.yscale('symlog')
plt.xlabel(r'x')
plt.ylabel(r'$\Psi_{k}$')

plt.figure(16)
plt.plot(x_t,dPsi[0],label = r'$\Psi^\prime_{1}$')
plt.plot(x_t,dPsi[1],label = r'$\Psi^\prime_{5}$')
plt.plot(x_t,dPsi[2],label = r'$\Psi^\prime_{10}$')
plt.plot(x_t,dPsi[3],label = r'$\Psi^\prime_{40}$')
plt.plot(x_t,dPsi[4],label = r'$\Psi^\prime_{60}$')
plt.plot(x_t,dPsi[5],label = r'$\Psi^\prime_{100}$')
plt.xlim([min(x_t),max(x_t)])
plt.legend(loc='best')
#plt.yscale('symlog')
plt.xlabel(r'x')
plt.ylabel(r'$\Psi^\prime_{k}$')

plt.figure(17)
plt.plot(x_t,dPhi[0],label = r'$\Phi^\prime_{1}$')
plt.plot(x_t,dPhi[1],label = r'$\Phi^\prime_{5}$')
plt.plot(x_t,dPhi[2],label = r'$\Phi^\prime_{10}$')
plt.plot(x_t,dPhi[3],label = r'$\Phi^\prime_{40}$')
plt.plot(x_t,dPhi[4],label = r'$\Phi^\prime_{60}$')
plt.plot(x_t,dPhi[5],label = r'$\Phi^\prime_{100}$')
plt.xlim([min(x_t),max(x_t)])
plt.legend(loc='best')
#plt.yscale('symlog')
plt.xlabel(r'x')
plt.ylabel(r'$\Phi^\prime_{k}$')
"""

#Milestone 4 
"""
x,S1,S2,S3,S4,S5,S6     = loadtxt("Source.dat",unpack=True)
S_lores1,S_lores2,S_lores3,S_lores4,S_lores5,S_lores6  =  loadtxt("S_low.dat",unpack=True)
"""
"""
#Low res for testing
plt.figure(18)
plt.title('THIS IS THE LOW RES SOURCE FOR TESTING')
#plt.plot(x_t,S_lores1,label = r'$S_{1}$')
#plt.plot(x_t,S_lores2,label = r'$S_{5}$')
#plt.plot(x_t,S_lores3,label = r'$S_{10}$')
plt.plot(x_t,S_lores4,label = r'$S_{40}$')
#plt.plot(x_t,S_lores5,label = r'$S_{60}$')
#plt.plot(x_t,S_lores6,label = r'$S_{100}$')
plt.xlim([min(x),max(x)])
plt.legend(loc='best')
plt.yscale('symlog')
plt.xlabel(r'x')
plt.ylabel(r'$S_k$')
"""

"""
#High res for further use
plt.figure(19)
#plt.plot(x,S1,label = r'$S_{50}$')
#plt.plot(x,S2,label = r'$S_{250}$')
#plt.plot(x,S3,label = r'$S_{500}$')
plt.plot(x,S4,label = r'$S_{2000}$')
#plt.plot(x,S5,label = r'$S_{3000}$')
#plt.plot(x,S6,label = r'$S_{5000}$')
plt.xlim([min(x),max(x)])
plt.legend(loc='best')
plt.yscale('symlog')
plt.xlabel(r'x')
plt.ylabel(r'$S_k$')

bessel = loadtxt("besseltest.dat",unpack=True)
plt.figure(20)
plt.plot(x,bessel)#,label=r'$$')
plt.xlim([min(x),max(x)])
#plt.legend(loc='best')
#plt.yscale('symlog')
plt.xlabel(r'x')
plt.ylabel(r'$j_l[k(\eta_0-\eta(x))]$')
"""
"""
Sj_l = loadtxt("Sj_l.dat",unpack=True)
plt.figure(21)
plt.plot(x,Sj_l/1e-3)#,label=r'$$')
plt.xlim([-8,max(x)])
#plt.legend(loc='best')
#plt.yscale('symlog')
plt.xlabel(r'x')
plt.ylabel(r'$\tilde{S}(k,x)j_l[k(\eta_0-\eta(x))]/10^{-3}$')

"""
"""
Theta_l = loadtxt("Theta_l.dat",unpack=True)
ls = loadtxt("ls.dat",unpack=True)
"""
#plot trasnfer function
"""
plt.figure(22)
plt.plot(ls,Theta_l[0],label = r'$\Theta_l(50)$')
plt.plot(ls,Theta_l[1],label = r'$\Theta_l(250)}$')
plt.plot(ls,Theta_l[2],label = r'$\Theta_l(500)}$')
plt.plot(ls,Theta_l[3],label = r'$\Theta_l(2000)}$')
plt.plot(ls,Theta_l[4],label = r'$\Theta_l(3000)}$')
plt.plot(ls,Theta_l[5],label = r'$\Theta_l(5000)}$')
plt.xlim([min(ls),max(ls)])
plt.legend(loc='best')
#plt.yscale('symlog')
plt.xlabel(r'l')
plt.ylabel(r'$\Theta_l(k)$')
"""

#Plot integrand in cmb spec
"""
integrand1 = loadtxt("integrand1.dat",unpack=True)
integrand2 = loadtxt("integrand2.dat",unpack=True)
integrand3 = loadtxt("integrand3.dat",unpack=True)
integrand4 = loadtxt("integrand4.dat",unpack=True)
integrand5 = loadtxt("integrand5.dat",unpack=True)
integrand6 = loadtxt("integrand6.dat",unpack=True)
ks = integrand1[0]
integrand1 = integrand1[1]

plt.figure(23)
plt.plot(ks,integrand1,label = 'l=2')
plt.plot(ks,integrand2,label = 'l=50')
plt.plot(ks,integrand3,label = 'l=200')
plt.plot(ks,integrand4,label = 'l=500')
plt.plot(ks,integrand5,label = 'l=800')
plt.plot(ks,integrand6,label = 'l=1200')
plt.xlim([min(ls),500])
plt.legend(loc='best')
#plt.yscale('symlog')
plt.xlabel(r'$ck/H_0$')
plt.ylabel(r'$l(l+1)\Theta_l^2(k)H_0/ck$')
plt.show()
"""

#Plot Power spectrum
l,Cl = loadtxt("C_l_default.dat",unpack=True,skiprows=1)
l,Clh066 = loadtxt("C_l_h066.dat",unpack=True,skiprows=1)
l,Clh074 = loadtxt("C_l_h074.dat",unpack=True,skiprows=1)
#load planck data
planck1  = loadtxt("COM_PowerSpect_CMB-TT-loL-full_R2.02.txt",unpack=True,skiprows=3)
planck2  = loadtxt("COM_PowerSpect_CMB-TT-hiL-full_R2.02.txt",unpack=True,skiprows=3)

planck_l1 = planck1[0]
planck_l2 = planck2[0]

C_llow = planck1[1]
C_lhi  = planck2[1]

error1 = planck1[2]
error2 = planck2[2]

Clplanck = hstack([C_llow,C_lhi])
planck_l = hstack([planck_l1,planck_l2])
error = hstack([error1,error2])

#planck3  = loadtxt("COM_PowerSpect_CMB-TT-hiL-binned_R2.02.txt",unpack=True,skiprows=3)
#planck_l = planck3[0]
#Clplanck = planck3[3]
#error    = planck3[4]




#normalize Cl to planck
Cl = Cl/max(Cl)*5775
Clh066 = Clh066/max(Clh066)*5775
Clh074 = Clh074/max(Clh074)*5775

#plot


plt.figure(123)
plt.errorbar(planck_l,Clplanck,yerr=error,label='Planck data',color='grey')#,fmt='-o')
plt.plot(l,Cl,label='Simulated')
#plt.plot(planck_l,Clplanck,label='Planck data')
plt.xlim([min(l),max(l)])
plt.legend(loc='best')
#plt.yscale('symlog')
plt.xlabel(r'l')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.show()

plt.figure(24)
plt.errorbar(planck_l,Clplanck,yerr=error,label='Planck data',color='grey')#,fmt='-o')
plt.plot(l,Cl,label='Default')
plt.plot(l,Clh066,label='h=0.66')
plt.plot(l,Clh074,label='h=0.74')
#plt.plot(planck_l,Clplanck,label='Planck data')
plt.xlim([min(l),max(l)])
plt.legend(loc='best')
#plt.yscale('symlog')
plt.xlabel(r'l')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.show()

l,Clm0200 = loadtxt("C_l_m0200.dat",unpack=True,skiprows=1)
l,Clm0248 = loadtxt("C_l_m0248.dat",unpack=True,skiprows=1)

Clm0200 = Clm0200/max(Clm0200)*5775
Clm0248 = Clm0248/max(Clm0248)*5775

plt.figure(24)
plt.errorbar(planck_l,Clplanck,yerr=error,label='Planck data',color='grey')#,fmt='-o')
plt.plot(l,Cl,label='Default ')
plt.plot(l,Clm0200,label=r'$\Omega_m=0.200$')
plt.plot(l,Clm0248,label=r'$\Omega_m=0.248$')
#plt.plot(planck_l,Clplanck,label='Planck data')
plt.xlim([min(l),max(l)])
plt.legend(loc='best')
#plt.yscale('symlog')
plt.xlabel(r'l')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.show()

l,Clb0042 = loadtxt("C_l_b0042.dat",unpack=True,skiprows=1)
l,Clb0050 = loadtxt("C_l_b0050.dat",unpack=True,skiprows=1)

Clb0042 = Clb0042/max(Clb0042)*5775
Clb0050 = Clb0050/max(Clb0050)*5775



plt.figure(24)
plt.errorbar(planck_l,Clplanck,yerr=error,label='Planck data',color='grey')#,fmt='-o')
plt.plot(l,Cl,label='Default ')
plt.plot(l,Clb0042,label=r'$\Omega_b=0.042$')
plt.plot(l,Clb0050,label=r'$\Omega_b=0.050$')
#plt.plot(planck_l,Clplanck,label='Planck data')
plt.xlim([min(l),max(l)])
plt.legend(loc='best')
#plt.yscale('symlog')
plt.xlabel(r'l')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.show()


