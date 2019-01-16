import scipy.special as sp
from fullspectrum import *
import matplotlib.pyplot as plt

h=0.05
epsilon=0.001
#L=8
#if __name__=="__main__":
colors = { 1:'r', 2:'b', 3:'g', 4:'y' }
maxdev = 4
for i in range(6,13):
        E=fullspectrum(i)[0]
        plt.scatter(i*np.ones(len(E)),E,color='r')

L = np.linspace(6,12,101)
for i in range(4):
    exact_single = lambda l : 1.*(1.-2.*epsilon*np.cos(4.*i*np.pi/l))
    plt.plot(L,exact_single(L),'g-')
def exact(n,l,k):
    zn=-sp.ai_zeros(4)[0]
    return 2.-4.*epsilon*np.cos(2.*np.pi*k/l)+ 4.*epsilon*((0.5*np.cos(2.*np.pi*k/l))**(1./3.))*h**(2./3.)*zn[n]
for k in range(1,5):
    for n in range(0,4):
        plt.plot(L,np.array([exact(n,i,k) for i in L] ),colors[n+1])

plt.show()
