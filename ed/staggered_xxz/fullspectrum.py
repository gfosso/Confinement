import numpy as np
import scipy.special as sp


def fullspectrum(L):

    #hilbertsize
    hilbertsize=2**L
    epsilon=0.01
    delta=-epsilon**(-1)
    h=0.00
    def binconf(c): return np.binary_repr(c,L)

    def readsite(conf,i): return (conf&(1<<i))>>i

    def SzSz(conf,i,j):
        si = readsite(conf,i)-1/2
        sj = readsite(conf,j)-1/2
        return si*sj

    def Sz(conf,i):
        return readsite(conf,i)-0.5


#for conf in range(hilbertsize):
#    print(IsingHam(conf))

    def Sx(conf,i):
        return 0.5, conf^(1<<i)


#for conf in range(hilbertsize):
#    print(binconf(conf),binconf(Sx(conf,0)[1]))


    def Spinflip(conf,i,j):
        "SxSx + SySy = 0.5(S+S- + S-S+) term"
        if readsite(conf,i) != readsite(conf,j):
            return 0.5, conf^((1<<i)^(1<<j))
        else: return 0.0, conf


    def XXZHam(conf):
        return sum([ -0.5*delta*(4.0*SzSz(conf,i,(i+1)%L) + 1.0  ) -2.*h*((-1)**i)*Sz(conf,i) for i in range(L) ])


    # diagonal part
    Ham = np.diag([XXZHam(conf) for conf in range(hilbertsize)])
    # off-diagonal part
    for conf in range(hilbertsize):
        for i in range(L):
            value, newconf = Spinflip(conf,i,(i+1)%L)
            Ham[newconf,conf] -= 2.*value
    
    return np.linalg.eigh((np.abs(delta)**(-1))*Ham)
#    return np.sort(np.linalg.eigvals(np.abs(delta)**(-1)*Ham))
