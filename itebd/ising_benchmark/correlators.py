from ham import *



def corrszsz(dist,s,B):
    sz=np.diag([Sz(conf,0) for conf in range(0,2)])
    corr=[]
    if dist ==0:
        sz2= np.tensordot(sz,sz,axes=(1,0))
        for i_bond in range(2):
            sB = np.tensordot(np.diag(s[np.mod(i_bond-1,2)]),B[i_bond],axes=(1,1))
            C=np.tensordot(sB,np.conj(sB),axes=([0,2],[0,2]))
            corr.append( np.tensordot(C,sz2,axes=([0,1],[0,1])) - np.tensordot(C,sz,axes=([0,1],[0,1]))*np.tensordot(C,sz,axes=([0,1],[0,1])))
        return np.mean(corr)

    if dist !=0:
        dist=np.abs(dist)
        T=[] #transfer matrix
        for i_bond in range(2):
            T.append(np.tensordot(B[i_bond],np.conj(B[i_bond]),axes=(0,0)))
        for i_bond in range(2):
            sB = np.tensordot(np.diag(s[np.mod(i_bond-1,2)]),B[i_bond],axes=(1,1))
            C = np.tensordot(sB,sz,axes=(1,0))
            sB=np.conj(sB)
            R = np.tensordot(C,sB,axes=([0,2],[0,1]))
            mean1= np.trace(R)
            for i in range(dist-1):
                R=np.tensordot(R,T[np.mod(i_bond+1+i,2)],axes=([0,1],[0,2]))
            C= np.tensordot(B[np.mod(i_bond+dist,2)],sz,axes=(0,0))
            sB=np.conj(B[np.mod(i_bond+dist,2)])
            L=np.tensordot(C,sB,axes=([2,1],[0,2]))
            mean2=np.trace(L)
            corr.append( np.tensordot(R,L,axes=([0,1],[0,1])) -mean1*mean2)
        return np.mean(corr)
    #mag=[]
#for i_bond in range(2):
#	sB = np.tensordot(np.diag(s[np.mod(i_bond-1,2)]),B[i_bond],axes=(1,1))
#	C = np.tensordot(sB,sz,axes=(1,0))
#	sB=np.conj(sB)
#	 mag.append(np.squeeze(np.tensordot(sB,C,axes=([0,2,1],[0,1,2]))).item()) 
#print("sigmazeta =", np.mean(mag))

#E=[]
#for i_bond in range(2):
#    BB = np.tensordot(B[i_bond],B[np.mod(i_bond+1,2)],axes=(2,1))
#    sBB = np.tensordot(np.diag(s[np.mod(i_bond-1,2)]),BB,axes=(1,1))
#    C = np.tensordot(sBB,np.reshape(H_bond,[d,d,d,d]),axes=([1,2],[2,3]))
#    sBB=np.conj(sBB)
#    E.append(np.squeeze(np.tensordot(sBB,C,axes=([0,3,1,2],[0,1,2,3]))).item()) 
#print("E_iTEBD =",np.mean(E))

