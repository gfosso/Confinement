import numpy as np
import scipy.special as sp
#size
L=12
#hilbertsize
hilbertsize=2**L
epsilon=0.001
delta=-epsilon**(-1)
h=0.05

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


def count(conf):
	return sum(readsite(conf,i) for i in range(L))
#for conf in range(hilbertsize):
#        print(binconf(conf),binconf(Spinflip(conf,0,1)[1]))


#translation operator 
def translate(conf):
    zero=readsite(conf,0)
    return (conf>>1)|(zero<<(L-1))

#total flip operator
def totalflip(conf):
        return conf^(hilbertsize-1)

#modified momentum operator
def modtrans(conf):
        return translate(totalflip(conf))


#gives the lowest integer representative and the periodicity
def lowestrepr(conf):
	conf0=conf
	conf1=conf
	for i in range(L):
		conf=modtrans(conf)
		if conf1>conf: conf1=conf
		elif conf0==conf: return conf1,i+1 
	return conf1,L

#gives the lowest integer representative and how many translations you need to reach that configuration	
def repr(conf):
	if lowestrepr(conf)[0]==conf: return conf,0
	lowest=lowestrepr(conf)[0]
	for i in range(L):
		conf=modtrans(conf)
		if conf==lowest:return lowest,i+1


def hilbertspace(Sz=0,m=0):
    #reduced hilbert space for symmetry sector tot_sz=0 
    c=[]
    #return True if the state doesn't exist, False if it's already in the configuration vector
    def checkstate(lw):
            for i in range(len(c)):
                    if c[i]==lw: return False
            return True
    for conf in range(hilbertsize):
    	if (count(conf) == L//2+Sz)&checkstate(lowestrepr(conf)):
        		lw=lowestrepr(conf)
        		c.append(lw)
    #reduced hilbert space for momentum states k=0
    ck=[]
    for i in range(len(c)):
        if m%(L/c[i][1]) == 0: ck.append(c[i])

    return ck

def hilbertspace_totsz(Sz=0):
    #reduced hilbert space for symmetry sector tot_sz=0 
    c=[]
    for conf in range(hilbertsize):
        if (count(conf) == L//2+Sz):
        		c.append(conf)
    return c

def XXZHam(conf):
    return sum([ -0.5*delta*(4.0*SzSz(conf,i,(i+1)%L) + 1.0  ) -2.*h*((-1)**i)*Sz(conf,i) for i in range(L) ])

def spectrum_totsz_k(Sz=0,m=0):
    # Hamiltonian in symmetry sector tot_sz=0 and k=0
    # diagonal part
    ck=hilbertspace(Sz,m)
    Ham = np.diag(np.array([XXZHam(ck[i][0]) for i in range(len(ck))],dtype=np.complex))
    # off-diagonal part
    for j in range(len(ck)):
        for i in range(L):
            value, newconf = Spinflip(ck[j][0],i,(i+1)%L)
            d=lowestrepr(newconf) 
            if m%(L/d[1]) == 0:
                Ham[ck.index(d),j] -= 2.*value*np.sqrt(ck[j][1]/d[1])*np.exp(-complex(0,(2.*np.pi*m)/L*repr(newconf)[1]))
    
    return np.sort(np.real(np.linalg.eigvals(np.abs(delta)**(-1)*Ham)))
       
def spectrum_totsz(Sz=0):
    # Hamiltonian in symmetry sector tot_sz=0 and k=0
    # diagonal part
    c=hilbertspace_totsz(Sz)
    Ham = np.diag([XXZHam(c[i]) for i in range(len(c))])
    # off-diagonal part
    for j in range(len(c)):
        for i in range(L):
            value, newconf = Spinflip(c[j],i,(i+1)%L)
            Ham[c.index(newconf),j] -= 2.*value
    
    return np.sort(np.linalg.eigvals(np.abs(delta)**(-1)*Ham))

def spectrum():
    # Hamiltonian in symmetry sector tot_sz=0 and k=0
    # diagonal part
    Ham = np.diag([XXZHam(conf) for conf in range(hilbertsize)])
    # off-diagonal part
    for conf in range(hilbertsize):
        for i in range(L):
            value, newconf = Spinflip(conf,i,(i+1)%L)
            Ham[newconf,conf] -= 2.*value
    
    return np.linalg.eigh((np.abs(delta)**(-1))*Ham)
#    return np.sort(np.linalg.eigvals(np.abs(delta)**(-1)*Ham))

#print(Ham)
#print(en[0]/L)

#E=[]

#for k in range(L//2+1):
#    E=np.append(E,spectrum(0,k))

#E=np.sort(E)

#print(spectrum(0,0))
