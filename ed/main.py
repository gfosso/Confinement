import scipy.special as sp
from staggered import *
import matplotlib.pyplot as plt
from equations import *
#h=0.001
#epsilon=0.001
#L=8
#if __name__=="__main__":
colors = { 1:'r', 2:'b', 3:'g', 4:'y' }
maxdev = 1
for dev in range(maxdev,0,-1):
    for k in range(L//2+1):
        Ek=np.sort(spectrum_totsz_k(Sz=dev-1,m=k))
        Ek=Ek[Ek<=2.5]
        Ek=Ek[Ek>=1.5]
        Ks = 2*np.pi*k/L*np.ones(len(Ek))
        plt.scatter(Ks,Ek,color=colors[dev])

k = np.linspace(0,0.5*np.pi,101)
#exact_single = lambda k : 1.*(1.-2.*epsilon*np.cos(2.*k))
#plt.plot(k,exact_single(k),'r-')
zn=-sp.ai_zeros(6)[0]
for n in range(6):
    plt.plot(k,np.array([E(n=n,P=i)-epsilon*h*L for i in k]),'--')
    exact= lambda k : 2.-4.*epsilon*np.cos(k) + 4.*epsilon*((0.5*np.cos(k))**(1./3.))*(h**(2./3.))*zn[n] -epsilon*h*L
    plt.plot(k,exact(k))

#plt.plot(k,np.array([E(n=0,P=i) for i in k]),'g-')
plt.show()
