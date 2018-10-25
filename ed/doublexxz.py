import numpy as np

#size
L=6
#hilbertsize
hilbertsize=2**(2*L)
gort=0.0
gpar=-4.0

def binconf(c): return np.binary_repr(c,L)

def readsite(conf,i): return (conf&(1<<i))>>i

def SzSz(conf,i,j):
    si = readsite(conf,i)-1/2
    sj = readsite(conf,j)-1/2
    return si*sj

def Sz(conf,i):
    return readsite(conf,i)-0.5

#longitudinal term chain 1
def XXZHam1(conf):
    return sum([- gpar*SzSz(conf,i,(i+1)%L)  for i in range(L) ])

#longitudinal term chain 2
def XXZHam2(conf):
    return sum([- gpar*SzSz(conf,i+L,(i+1)%L+L)  for i in range(L) ])

#coupling between the chains
def XXZHamort(conf):
    return sum([- gort*SzSz(conf,i,i+L)  for i in range(L) ])

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


#for conf in range(hilbertsize):
#        print(binconf(conf),binconf(Spinflip(conf,0,1)[1]))




#Hamiltonian
# diagonal part
Ham = np.diag([XXZHam1(conf) for conf in range(hilbertsize)])
Ham += np.diag([XXZHam2(conf) for conf in range(hilbertsize)])
Ham += np.diag([XXZHamort(conf) for conf in range(hilbertsize)])
# off-diagonal part
for conf in range(hilbertsize):
    for i in range(L):
        value, newconf = Spinflip(conf,i,(i+1)%L)
        Ham[newconf,conf] -=value     
    for i in range(L):
        value, newconf = Spinflip(conf,i+L,(i+1)%L+L)
        Ham[newconf,conf] -=value     
        
        
print(Ham)
en ,pin= np.linalg.eigh(Ham)
print(en[0]/L)



