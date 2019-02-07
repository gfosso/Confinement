from mps import *
import sys
from scipy.linalg import expm
import matplotlib.pyplot as plt
from scipy.fftpack import fft

def main(argv):
    cristo=mps(bond_dimension=300)
    hilbertsize=2**2
    epsilon=0.001
    gpar=-epsilon**(-1)
    h=5

    HamA = np.diag([ -0.5*gpar*(4.*SzSz(conf,0,1)+1)-4.*h*Sz(conf,0)+4.*h*Sz(conf,1)   for conf in range(hilbertsize)])
    HamB = np.diag([ -0.5*gpar*(4.*SzSz(conf,0,1)+1)+4.*h*Sz(conf,0)-4.*h*Sz(conf,1)   for conf in range(hilbertsize)])
    # off-diagonal part
    for conf in range(hilbertsize):
            value, newconf = Spinflip(conf,0,1)
            HamA[newconf,conf] -=2.*value     
            HamB[newconf,conf] -=2.*value     
#    HamA=epsilon*HamA
#    HamB=epsilon*HamB
    delta=0.0001;T=100;L=int(T/delta)
    UA=np.reshape(expm(-complex(0,delta)*HamA),(2,2,2,2))
    UB=np.reshape(expm(-complex(0,delta)*HamB),(2,2,2,2))
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
        cristo.evol(UA,UB)
#    np.savetxt("cristo.out",corr)
#    corrf= fft(corr[2])
    #time=np.arange(0,T,delta)
    #plt.plot(time,corr[0])
    #plt.plot(time,corr[1])
    #plt.plot(time,corr[2])
    plt.plot(np.array([2*np.pi*i/T for i in range(1,L//2)]),np.abs(fft(corr[0])[1:L//2]))
    plt.plot(np.array([2*np.pi*i/T for i in range(1,L//2)]),np.abs(fft(corr[1])[1:L//2]))
    plt.plot(np.array([2*np.pi*i/T for i in range(1,L//2)]),np.abs(fft(corr[2])[1:L//2]))
    plt.yscale("log")
    plt.show()
if __name__=="__main__":
    main(sys.argv[1:])
