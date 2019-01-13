import scipy.special as sp
from staggered import *
import matplotlib.pyplot as plt

#h=0.001
#epsilon=0.001
#L=8
#if __name__=="__main__":
colors = { 1:'r', 2:'b', 3:'g', 4:'y' }
maxdev = 4
for dev in range(maxdev,0,-1):
    for k in range(L//2+1):
        Ek=np.sort(spectrum_totsz_k(Sz=dev-1,m=k))
      #  Ek=Ek[Ek<=2.5]
      #  Ek=Ek[Ek>=1.5]
        Ks = 2*np.pi*k/L*np.ones(len(Ek))
        plt.scatter(Ks,Ek,color=colors[dev])

k = np.linspace(0,0.5*np.pi,101)
#exact_single = lambda k : 1.*(1.-2.*epsilon*np.cos(2.*k))
#plt.plot(k,exact_single(k),'r-')
#Ek[k]=E[:]
#    plt.plot(spec)
zn=-sp.ai_zeros(6)[0]
#for i in range(len(zn)):
#    exact= lambda k : 2.-4.*epsilon*np.cos(k) + 4.*epsilon*((0.5*np.cos(k))**(1./3.))*(h**(2./3.))*zn[i]
#    plt.plot(k,exact(k))
#exact= lambda k : 2.-4.*epsilon*np.cos(k) + 4.*epsilon*((0.5*np.cos(k))**(1./3.))*(h**(2./3.))*zn[1]
#plt.plot(k,exact(k),'g-')
plt.show()
