import scipy.special as sp
from staggered import *
import matplotlib.pyplot as plt

#h=0.001
#epsilon=0.001
#L=8
#if __name__=="__main__":
colors = { 1:'r', 2:'b', 3:'g', 4:'y' }
maxdev = 4
E0=-0.25*L
for dev in range(maxdev,0,-1):
    for k in range(L//2+1):
        Ek=np.sort(spectrum_totsz_k(Sz=L//2-dev,m=k))-E0
        Ks = 2*np.pi*k/L*np.ones(len(Ek))
        plt.scatter(Ks,Ek,color=colors[dev])

#Ek[k]=E[:]
#    plt.plot(spec)
plt.show()
