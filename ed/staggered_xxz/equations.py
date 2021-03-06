from staggered import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import scipy.special as sp

# Define the expression whose roots we want to find

def E(n=0,P=0):
    func = lambda pa : -2.*pa*np.cos(2.*pa) + np.sin(2.*pa) - (np.pi*h*(n-0.25))/(np.cos(P))

# Plot it

#pa = np.linspace(-np.pi, np.pi, 201)

#plt.plot(pa, func(pa))
#plt.xlabel("pa")
#plt.ylabel("expression value")
#plt.grid()
#plt.show()

# Use the numerical solver to find the roots

    pa_initial_guess = 0.
    pa_solution = fsolve(func, pa_initial_guess)

#print("The solution is pa = ",pa_solution[0])
#print("at which the value of the expression is",func(pa_solution)[0])
    return 2.-4.*epsilon*np.cos(P)*np.cos(2.*pa_solution[0])


def E_bessel(n=0,P=0,nmax=10):
    prec=10**(-2)
    step=0.02
    x=np.arange(-nmax,nmax,step)
    func = lambda x : sp.jv(x,h**(-1)*np.cos(P))
    z = func(x)*func(x+step)
    zeros=np.sort(fsolve(func,x[z<0]))[::-1]
    return 2. -4.*epsilon*h*zeros[n]
