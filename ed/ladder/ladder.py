import numpy as np
#size
L=6
#hilbertsize
hilbertsize=2**(2*L)
gort=0.5
epsilon=0.001
gpar=-epsilon**(-1)

def binconf(c): return np.binary_repr(c,2*L)

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


def count(conf,c):
	if c==0:
		return sum(readsite(conf,i) for i in range(L))
	else if c==1:
		return sum(readsite(conf,i) for i in range(L,2*L))
def translate(conf):
	elle=readsite(conf,L)
	zero=readsite(conf,0)
	if(elle==0):
		return (conf>>1)|(zero<<(L-1))
	else:
		return((conf>>1)|(elle<<(2*L-1)))^((elle^zero)<<(L-1))

#gives the lowest integer representative and the periodicity
def lowestrepr(conf):
	conf0=conf
	conf1=conf
	for i in range(L):
		conf=translate(conf)
		if conf1>conf: conf1=conf
		elif conf0==conf: return conf1,i+1 
	return conf1,L

#gives the lowest integer representative and how many translations you need to reach that configuration	
def repr(conf):
	if lowestrepr(conf)[0]==conf: return conf,0
	lowest=lowestrepr(conf)[0]
	for i in range(L):
		conf=translate(conf)
		if conf==lowest:return lowest,i+1
#must be modified
def hilbertspace(Sz=0,m=0):
    #reduced hilbert space for symmetry sector tot_sz=0 
    c=[]
    #return True if the state doesn't exist, False if it's already in the configuration vector
    def checkstate(lw):
            for i in range(len(c)):
                    if c[i]==lw: return False
            return True

    for conf in range(hilbertsize):
    	if (count(conf) == L//2+Sz)&checkstate(lowestrepr(conf))&twokinks(conf):
        		lw=lowestrepr(conf)
        		c.append(lw)
    #reduced hilbert space for momentum states k=0
    ck=[]
    for i in range(len(c)):
        if m%(L/c[i][1]) == 0: ck.append(c[i])

    return ck

#Hamiltonian
#longitudinal term chain 1
def XXZHam1(conf):
	return sum([-0.5*gpar*(4.*SzSz(conf,i,(i+1)%L)+1.)  for i in range(L) ])

#longitudinal term chain 2
def XXZHam2(conf):
	return sum([-0.5*gpar*(4.*SzSz(conf,i+L,(i+1)%L+L)+1.)  for i in range(L) ])

#coupling between the chains
def XXZHamort(conf):
	return sum([- 4.*gort*SzSz(conf,i,i+L)  for i in range(L) ])
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
        
        
#print(Ham)
#en ,pin= np.linalg.eigh(Ham)
#print(en[0]/L)



