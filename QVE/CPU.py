from pylab import *
from QPU import *
import time
from scipy.linalg import expm,logm,polar,inv as expm,logm,polar,inv

def E(mlist,x):
    for k in range(len(mlist)):
        u2x2=mlist[k][0]
        i=mlist[k][1]
        j=mlist[k][2]
        d=x.shape[1]
        Ev=fastprod(u2x2,i,j,x,d)
    E=np.dot(x.T,Ev)
    
    return E
    
def grad(A,x):
    d=A.shape[0]
    df=np.empty([d,1])
    for i in range(len(df)):
        Fi=np.take(A,i,axis=0).reshape(1,d)            #fila i pel vector
        Ci=np.take(A,i,axis=1).reshape(1,d)                       #Columna i per vector[i]
        df[i]=np.dot(Fi,x)+np.dot(Ci,x)
    return df
    
def numgrad(E1,E0,x1,x0):
    d=len(x1)
    gradf=np.empty([d,1])
    for i in range(d):
        gradf[i]=(E1-E0)/(x1[i]-x0[i])
    return gradf

def norm(v):
    dummy=0
    for i in range(len(v)):
        dummy+=v[i]**2
    return sqrt(dummy)
    
def newAnsatz(ansatz,A):
    curr=ansatz
    gamma = 1e-2
    #diff=ansatz+1/precision
   # while norm(diff)>precision:         #Com tractar situacio "O-O"???
    prev=curr
    curr=curr-gamma*(grad(A,prev))
    curr=curr/norm(curr)
        #diff=abs(curr-prev)
        
    return curr
#np.set_printoptions(precision=6,suppress=True)
#N=2
    
#P=np.array([[0,1,1,0],[1,0,0,1],[1,0,0,1],[0,1,1,0]])*(-1)          #-1 per Ising
#x=np.array([[0],[1],[1],[1]])/(2)     #00+10
"""    
curr=x
gamma=1e-2
diff=x+1 
i=0   
while norm(diff)>precision:         #Com tractar situacio "O-O"???
    i+=1
    prev=curr
    curr=curr-gamma*(grad(P,prev))
    curr=curr/norm(curr)
    diff=abs(curr-prev)
    #print(-gamma*grad(P,prev))
    #print(curr,prev)
#print(E(P,x),E(P,curr),i)
#print(curr)
"""

