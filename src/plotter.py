import matplotlib.pyplot as plt
from numpy import *
from matplotlib import rc 
rc('font',**{'family':'serif'}) # This is for Latex writing

#Import datafiles
#Relative densities
nu,lamda = loadtxt("omega_nulambda.dat",unpack=True)
m,b,r = loadtxt("omega_mbr.dat",unpack=True)
total=nu+lamda+m+b+r

#x,scalefactor,conformal time
x,a,eta = loadtxt("xaeta.dat",unpack=True)
#Hubbleparam,redshift
H,z = loadtxt("Hz.dat",unpack=True)
#Test x,test eta
x_t,eta_t = loadtxt("eta_t.dat",unpack=True)



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
plt.show()



