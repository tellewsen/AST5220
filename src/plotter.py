import matplotlib.pyplot as plt
from numpy import *
from matplotlib import rc 
rc('font',**{'family':'serif'}) # This is for Latex writing

nu,lamda = loadtxt("omega_nulambda.dat",unpack=True)
m,b,r = loadtxt("omega_mbr.dat",unpack=True)
x,a,eta = loadtxt("xaeta.dat",unpack=True)
H,z = loadtxt("Hz.dat",unpack=True)
x_t,eta2 = loadtxt("eta2.dat",unpack=True)
total=nu+lamda+m+b+r

plt.figure(0)
plt.plot(x,m,label = r'$\Omega_m$')
plt.plot(x,b,label = r'$\Omega_b$')
plt.plot(x,r,label = r'$\Omega_r$')
plt.plot(x,nu,label = r'$\Omega_\nu$')
plt.plot(x,lamda,label = r'$\Omega_\lambda$')
plt.plot(x,total,label = r'Total')
plt.ylim([0,1.1])
plt.xlim([min(x),max(x)])
#plt.yscale('log')
plt.legend(loc='best')
plt.title('Relative densities vs x')
plt.xlabel(r'x')
plt.ylabel(r'$\Omega_x$')

plt.figure(1)
plt.plot(x,eta,label = r'$\eta(x)$')
plt.plot(x_t,eta2,label=r'$\eta2(x)$')
#plt.ylim([0,1.1])
#plt.xlim([min(x),max(x)])
plt.yscale('log')
#plt.legend(loc='best')
plt.title('Conformal time vs x')
plt.xlabel(r'x')
plt.ylabel(r'$\eta(x)/Mpc$')
plt.legend(loc='best')

plt.figure(2)
plt.plot(x,H,label = r'H(x)')
#plt.ylim([0,1.1])
#plt.xlim([min(x),max(x)])
plt.yscale('log')
#plt.legend(loc='best')
plt.title('Hubble parameter vs x')
plt.xlabel(r'x')
plt.ylabel(r'H(x)')

plt.figure(3)
plt.plot(z,H,label = r'H(z)')
#plt.ylim([0,1.1])
#plt.xlim([min(x),max(x)])
#plt.yscale('log')
#plt.legend(loc='best')
plt.xlim([z[0],z[-1]])
plt.title('Hubble parameter vs z')
plt.xlabel(r'z')
plt.ylabel(r'H(z)')
plt.show()



