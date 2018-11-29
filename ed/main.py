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
        Ek=Ek[Ek<=2.]
        Ks = 2*np.pi*k/L*np.ones(len(Ek))
        plt.scatter(Ks,Ek,color=colors[dev])

k = np.linspace(0,np.pi,101)
exact_single = lambda k : 1.*(1.-2.*epsilon*np.cos(2.*k))
plt.plot(k,exact_single(k),'r-')
#Ek[k]=E[:]
#    plt.plot(spec)
plt.show()
