from ham import *
import cupy as cp
import numpy as np


#expectation value of sz
def magnetization(s,B,d):
    sz=cp.diag([Sz(conf,0) for conf in range(0,d)])
 #   sz=cp.array([[0,1],[1,0]])
    mag=cp.array(0.,dtype=np.float32)
    for i_bond in range(2):
        sB = cp.tensordot(cp.diag(s[np.mod(i_bond-1,2)]),B[i_bond],axes=(1,1))
        C=cp.tensordot(sB,cp.conj(sB),axes=([0,2],[0,2]))
        mag+=cp.real( cp.tensordot(C,sz,axes=([0,1],[0,1])).get())
    return mag*0.5

#correlator between sz in different positions
def corrszsz(dist,s,B,d):
    sz=cp.diag([Sz(conf,0) for conf in range(0,d)])
    corr=cp.array(0.,dtype=np.float32)
    if dist ==0:
        sz2= cp.tensordot(sz,sz,axes=(1,0))
        for i_bond in range(2):
            sB = cp.tensordot(cp.diag(s[np.mod(i_bond-1,2)]),B[i_bond],axes=(1,1))
            C=cp.tensordot(sB,cp.conj(sB),axes=([0,2],[0,2]))
            corr += cp.real(cp.tensordot(C,sz2,axes=([0,1],[0,1])) - cp.tensordot(C,sz,axes=([0,1],[0,1]))*cp.tensordot(C,sz,axes=([0,1],[0,1])))
        return 0.5*corr

    if dist !=0:
        dist=cp.abs(dist)
        for i_bond in range(2):
            sB = cp.tensordot(cp.diag(s[np.mod(i_bond-1,2)]),B[i_bond],axes=(1,1))
            C=cp.tensordot(sB,cp.conj(sB),axes=(0,0))
            R = cp.tensordot(C,sz,axes=([0,2],[0,1]))
            mean1= cp.trace(R)
            for i in range(dist-1):
                T=cp.tensordot(R,B[np.mod(i_bond+1+i,2)],axes=(0,1))
                T=cp.tensordot(T,cp.conj(B[np.mod(i_bond+1+i,2)]),axes=(0,1))
                R=cp.trace(T,axis1=0,axis2=2)
            C=cp.tensordot(B[cp.mod(i_bond+dist,2)],cp.conj(B[cp.mod(i_bond+dist,2)]),axes=(2,2))
            L=cp.tensordot(R,C,axes=([0,1],[1,3]))
            corr += cp.real( cp.tensordot(L,sz,axes=([0,1],[0,1])) - mean1*mean1)
        return 0.5*corr


#time evolution
def evol(s,B,U,chi,d):
    for i_bond in [0,1]:
        ia = np.mod(i_bond-1,2); ib = np.mod(i_bond,2); ic = np.mod(i_bond+1,2)
        chia = B[ib].shape[1]; chic = B[ic].shape[2]
        # Construct theta matrix and time evolution #
        theta = cp.tensordot(B[ib],B[ic],axes=(2,1)) # i a j b
        theta = cp.tensordot(U,theta,axes=([2,3],[0,2])) # ip jp a b 
        theta = cp.tensordot(cp.diag(s[ia]),theta,axes=([1,2])) # a ip jp b 
        theta = cp.reshape(cp.transpose(theta,(1,0,2,3)),(d*chia,d*chic)) # ip a jp b
        # Schmidt decomposition #
        X, Y, Z = cp.linalg.svd(theta,full_matrices=0)
        chi2 = np.min([cp.sum(Y>10.**(-10)).get(), chi])	
        piv = cp.zeros(len(Y), cp.bool)
        piv[(cp.argsort(Y)[::-1])[:chi2]] = True
        Y = Y[piv]; invsq = cp.sqrt(sum(Y**2))
        X = X[:,piv] 
        Z = Z[piv,:]
        # Obtain the new values for B and s #
        s[ib] = Y/invsq
        X=cp.reshape(X,(d,chia,chi2))
        X = cp.transpose(cp.tensordot(cp.diag(s[ia]**(-1)),X,axes=(1,1)),(1,0,2))
        B[ib] = cp.tensordot(X, cp.diag(s[ib]),axes=(2,0))
        B[ic] = cp.transpose(cp.reshape(Z,(chi2,d,chic)),(1,0,2))

    return s,B

