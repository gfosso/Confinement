import matplotlib.pyplot as plt
import numpy as np

T=30
dist=10
delta=0.1
t=np.arange(0,T,delta)
distance= np.arange(-dist+1,dist,1)
data=np.loadtxt('cristo.out')
A=np.zeros((len(t),dist*2-1))
for i in range(0,len(data)):
    A[i,:]=np.append(np.flip(data[i,1:]),data[i,:])
plt.contourf(distance,t,A,20)
        #,levels=np.arange(-0.02,0.24,0.01))
        #cmap='RdGy')
plt.title("correlator with gort=0.1  gpar=3. and htras=0.5 to 0.25")
plt.colorbar()
plt.show()
