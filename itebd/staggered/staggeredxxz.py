from mps import *
import sys
from scipy.linalg import expm
import matplotlib.pyplot as plt

def main(argv):
    cristo=mps()
    hilbertsize=2**2
    gpar=2
    
    Ham = np.diag([ +gpar*0.5*SzSz(conf,0,1)   for conf in range(hilbertsize)])
    # off-diagonal part
    for conf in range(hilbertsize):
            value, newconf = Spinflip(conf,0,1)
            Ham[newconf,conf] +=value     

    delta=0.1;T=10;L=int(T/delta)
    U=np.reshape(expm(-complex(0,delta)*Ham),(2,2,2,2))
    cristo.product_state(neel=True)
    corr=[]
    corr.append(np.zeros(L))
    corr.append(np.zeros(L))
    corr.append(np.zeros(L))
    for step in range(L):
        corr[0][step]=4.*cristo.expectation_SzSz(1,connected=False)
        corr[1][step]=4.*cristo.expectation_SzSz(2,connected=False)
        corr[2][step]=4.*cristo.expectation_SzSz(3,connected=False)  
        cristo.evol(U)
#    np.savetxt("cristo.out",corr)

    plt.plot(corr[0])
    plt.plot(corr[1])
    plt.plot(corr[2])
    plt.show()
if __name__=="__main__":
    main(sys.argv[1:])
