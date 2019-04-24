l=np.where( (x[0]-min(x[0]))<5.5)[0]
states=np.array([x[0][l[i]]for i in range(len(l))])
#coeff=[]
#for i in range(len(l)):
#    coeff.append(np.real(x[1].T[l[i]][np.where(np.abs(x[1].T[l[i]])>0.1)]))

s=[]
for i in range(len(l)):
    index=np.where(np.abs(x[1].T[l[i]])>0.1)[0]
    s.append([(str(round(np.real(x[1].T[l[i]][j]),2)),binconf(h[j][0])[0:L],binconf(h[j][0])[L:2*L]) for j in index])



        
