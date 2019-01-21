import numpy as np
import scipy.special as sp
import scipy.linalg as spl

def fullspectrum(L):

    #hilbertsize
    hilbertsize=2**(2*L)
    gort=0.5
    epsilon=0.001
    gpar=-epsilon**(-1)
    def binconf(c): return np.binary_repr(c,L)

    def readsite(conf,i): return (conf&(1<<i))>>i

    def SzSz(conf,i,j):
        si = readsite(conf,i)-1/2
        sj = readsite(conf,j)-1/2
        return si*sj

    def Sz(conf,i):
        return readsite(conf,i)-0.5

    def Sx(conf,i):
        return 0.5, conf^(1<<i)

    def Spinflip(conf,i,j):
        "SxSx + SySy = 0.5(S+S- + S-S+) term"
        if readsite(conf,i) != readsite(conf,j):
            return 0.5, conf^((1<<i)^(1<<j))
        else: return 0.0, conf

    #longitudinal term chain 1
    def XXZHam1(conf):
        return sum([-0.5*gpar*(4.*SzSz(conf,i,(i+1)%L)+1.)  for i in range(L) ])

    #longitudinal term chain 2
    def XXZHam2(conf):
        return sum([-0.5*gpar*(4.*SzSz(conf,i+L,(i+1)%L+L)+1.)  for i in range(L) ])

    #coupling between the chains
    def XXZHamort(conf):
        return sum([- 4.*gort*SzSz(conf,i,i+L)  for i in range(L) ])


    #Hamiltonian
    # diagonal part
    Ham = np.diag([XXZHam1(conf) for conf in range(hilbertsize)])
    Ham += np.diag([XXZHam2(conf) for conf in range(hilbertsize)])
    Ham += np.diag([XXZHamort(conf) for conf in range(hilbertsize)])
    # off-diagonal part
    for conf in range(hilbertsize):
        for i in range(L):
            value, newconf = Spinflip(conf,i,(i+1)%L)
            Ham[newconf,conf] -=2.*value     
        for i in range(L):
            value, newconf = Spinflip(conf,i+L,(i+1)%L+L)
            Ham[newconf,conf] -=2.*value     
        
        
    
#    return spl.eigh((np.abs(gpar)**(-1))*Ham)
    return np.sort(np.real(np.linalg.eigvals(epsilon*Ham)))
#    return np.sort(np.linalg.eigvals(np.abs(delta)**(-1)*Ham))
