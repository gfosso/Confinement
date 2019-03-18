from mps import *
import sys
from scipy.linalg import expm
import matplotlib.pyplot as plt
from scipy.fftpack import fft

def main(argv):
    cristo=mps(bond_dimension=300,site_dimension=4)
    hilbertsize=4**2
    gpar=-10
    gort=0.5 #~h in the staggered

    HamA = np.diag([ -0.5*gpar*(4.*SzSz(conf,0,2)+1)- 0.5*gpar*(4.*SzSz(conf,1,3)+1)   for conf in range(hilbertsize)])
    HamA += np.diag([ -gort*2.*SzSz(conf,0,1) -gort*2.*SzSz(conf,2,3)    for conf in range(hilbertsize)])
    # off-diagonal part
    for conf in range(hilbertsize):
            value, newconf = Spinflip(conf,0,2)
            HamA[newconf,conf] -=2.*value     
            value, newconf = Spinflip(conf,1,3)
            HamA[newconf,conf] -=2.*value    
    HamB=HamA
#    HamA=epsilon*HamA
#    HamB=epsilon*HamB
    delta=0.01;T=100;L=int(T/delta)
    epsilon=np.abs(gpar**(-1))
    UA=np.reshape(expm(-complex(0,delta)*HamA),(4,4,4,4))
    UB=np.reshape(expm(-complex(0,delta)*HamB),(4,4,4,4))
#    print(U)
    cristo.product_state(neel=True)
    corr=[]
    corr.append(np.zeros(L))
    corr.append(np.zeros(L))
    corr.append(np.zeros(L))
    for step in range(L):
        corr[0][step]=4.*cristo.expectation_SzSz(1,connected=False)
        corr[1][step]=4.*cristo.expectation_SzSz(2,connected=False)
        corr[2][step]=4.*cristo.expectation_SzSz(3,connected=False)  
 #       corr[0][step]=4.*cristo.staggered_Sz()
        cristo.evol(UA,UB)
#    np.savetxt("cristo.out",corr)
#    corrf= fft(corr[2])
#    time=np.arange(0,T,delta)
#    plt.plot(time,corr[0])
    #plt.plot(time,corr[1])
    #plt.plot(time,corr[2])
    plt.plot(np.array([2*np.pi*i/T for i in range(1,L//2)]),np.abs(fft(corr[0])[1:L//2]))
    plt.plot(np.array([2*np.pi*i/T for i in range(1,L//2)]),np.abs(fft(corr[1])[1:L//2]))
    plt.plot(np.array([2*np.pi*i/T for i in range(1,L//2)]),np.abs(fft(corr[2])[1:L//2]))
    plt.yscale("log")
    plt.show()
if __name__=="__main__":
    main(sys.argv[1:])
