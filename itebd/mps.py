from ham import *


#expectation value of sz
def magnetization(s,B,d):
    sz=np.diag([Sz(conf,0) for conf in range(0,d)])
 #   sz=np.array([[0,1],[1,0]])
    mag=[]
    for i_bond in range(2):
        sB = np.tensordot(np.diag(s[np.mod(i_bond-1,2)]),B[i_bond],axes=(1,1))
        C=np.tensordot(sB,np.conj(sB),axes=([0,2],[0,2]))
        mag.append( np.tensordot(C,sz,axes=([0,1],[0,1])))
    return np.mean(np.real(mag))

#correlator between sz in different positions
def corrszsz(dist,s,B,d):
    sz=np.diag([Sz(conf,0) for conf in range(0,d)])
    corr=[]
    if dist ==0:
        sz2= np.tensordot(sz,sz,axes=(1,0))
        for i_bond in range(2):
            sB = np.tensordot(np.diag(s[np.mod(i_bond-1,2)]),B[i_bond],axes=(1,1))
            C=np.tensordot(sB,np.conj(sB),axes=([0,2],[0,2]))
            corr.append( np.tensordot(C,sz2,axes=([0,1],[0,1])) - np.tensordot(C,sz,axes=([0,1],[0,1]))*np.tensordot(C,sz,axes=([0,1],[0,1])))
        return np.mean(np.real(corr))

    if dist !=0:
        dist=np.abs(dist)
        for i_bond in range(2):
            sB = np.tensordot(np.diag(s[np.mod(i_bond-1,2)]),B[i_bond],axes=(1,1))
            C=np.tensordot(sB,np.conj(sB),axes=(0,0))
            R = np.tensordot(C,sz,axes=([0,2],[0,1]))
            mean1= np.trace(R)
            for i in range(dist-1):
                T=np.tensordot(R,B[np.mod(i_bond+1+i,2)],axes=(0,1))
                T=np.tensordot(T,np.conj(B[np.mod(i_bond+1+i,2)]),axes=(0,1))
                R=np.trace(T,axis1=0,axis2=2)
            C=np.tensordot(B[np.mod(i_bond+dist,2)],np.conj(B[np.mod(i_bond+dist,2)]),axes=(2,2))
            L=np.tensordot(C,sz,axes=([0,2],[0,1]))
            mean2=np.trace(L)
            corr.append( np.tensordot(R,L,axes=([0,1],[0,1])) - mean1*mean1)
        return np.mean(np.real(corr))


#time evolution
def evol(s,B,U,chi,d):
    for i_bond in [0,1]:
        ia = np.mod(i_bond-1,2); ib = np.mod(i_bond,2); ic = np.mod(i_bond+1,2)
        chia = B[ib].shape[1]; chic = B[ic].shape[2]
        # Construct theta matrix and time evolution #
        theta = np.tensordot(B[ib],B[ic],axes=(2,1)) # i a j b
        theta = np.tensordot(U,theta,axes=([2,3],[0,2])) # ip jp a b 
        theta = np.tensordot(np.diag(s[ia]),theta,axes=([1,2])) # a ip jp b 
        theta = np.reshape(np.transpose(theta,(1,0,2,3)),(d*chia,d*chic)) # ip a jp b
        # Schmidt decomposition #
        X, Y, Z = np.linalg.svd(theta,full_matrices=0)
        chi2 = np.min([np.sum(Y>10.**(-10)), chi])	
        piv = np.zeros(len(Y), np.bool)
        piv[(np.argsort(Y)[::-1])[:chi2]] = True
        Y = Y[piv]; invsq = np.sqrt(sum(Y**2))
        X = X[:,piv] 
        Z = Z[piv,:]
        # Obtain the new values for B and s #
        s[ib] = Y/invsq 
        X=np.reshape(X,(d,chia,chi2))
        X = np.transpose(np.tensordot(np.diag(s[ia]**(-1)),X,axes=(1,1)),(1,0,2))
        B[ib] = np.tensordot(X, np.diag(s[ib]),axes=(2,0))
        B[ic] = np.transpose(np.reshape(Z,(chi2,d,chic)),(1,0,2))

    return s,B

