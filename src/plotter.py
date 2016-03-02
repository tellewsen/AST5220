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
#Test x,test eta
x_t,eta_t                = loadtxt("eta_t.dat",unpack=True)


x_rec,z_rec,X_e          = loadtxt("X_e.dat",unpack=True)
n_e,tau,dtau             = loadtxt("n_e.dat",unpack=True)
x_test,z_test,n_etest    = loadtxt("n_etest.dat",unpack=True)
ddtau,tau_test,dtau_test = loadtxt("tau2.dat",unpack=True)
ddtau_test               = loadtxt("tau3.dat",unpack=True)
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
plt.xlim([z_rec[0],z_rec[-1]])
plt.xlabel(r'z')
plt.ylabel(r'$X_e$')
plt.show()
"""
"""
plt.figure(5)
plt.plot(z_rec ,n_e,'--',label = r'$n_{e}$')
plt.plot(z_test,n_etest,'-',label = r'Splined')
plt.yscale('log')
#plt.xscale('log')
plt.xlim([z_rec[0],z_rec[-1]])
plt.xlabel(r'z')
plt.ylabel(r'$n_e$[m$^{-3}$]')
plt.legend()
plt.show()
"""
plt.figure(6)
plt.plot(x_rec , tau,            '--',label = r'$\tau(x)$')
plt.plot(x_test, tau_test,       '-' ,label = r'$\tau_{test}(x)$')
plt.plot(x_rec , abs(dtau),      '--',label = r'$|\tau^\prime(x)|$')
plt.plot(x_test, abs(dtau_test), '-' ,label = r'$|\tau^{\prime}_{test}(x)|$')
plt.plot(x_rec , abs(ddtau),     '--' ,label = r'$|\tau^{\prime\prime}(x)|$')
plt.plot(x_test, abs(ddtau_test),'-' ,label = r'$|\tau^{\prime\prime}_{test}(x)|$')
plt.yscale('log')
#plt.xscale('log')
plt.xlim([x_rec[0],x_rec[-1]])
plt.xlabel(r'x')
plt.ylabel(r'$\tau$,$|\tau^\prime|$,$|\tau^{\prime\prime}|$')
plt.legend()
plt.show()

