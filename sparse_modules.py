from pylab import *
import qgates as q
import time
from scipy.sparse.linalg import expm as expm
from scipy import sparse as sp
np.set_printoptions(precision=4,suppress=True)


def Hising(N,J,h):
    unitarylist=[]
    dummylist=[]
    Z=sp.csr_matrix(q.Z())
    X=sp.csr_matrix(q.X())
    Z2=sp.kron(Z,Z) 
    H1=sp.csr_matrix((2**N,2**N),dtype="complex")           #SigmaX
    H2=sp.csr_matrix((2**N,2**N),dtype="complex")           #SigmaZ
    
    
    #Test Hamiltonian@Ising
    
    for i in range(N):        
        H1=-J*sp.kron(sp.kron(sp.eye(2**i),X),sp.eye(2**(N-1-i)))
        if i <N-1:
            H2=-h*sp.kron(sp.kron(sp.eye(2**i),Z2),sp.eye(2**(N-1-i-1)))
            dummylist.append(H2)
        unitarylist.append(H1)
        

    
    #ising=J*H1+h*H2
    return unitarylist+dummylist

    
#def Tx2(unitary,d):
#    u2x2list=[]
#    #Bucle that decomposes unitary U into 2-level unitary matrices
#    for i in range(d):        #cicla columnes
#        column=unitary.getcol(i)        #Stores column i
#        indexes=np.transpose(column.nonzero())        #Stores nonzero elements from that row
#        while len(indexes)!=1 and i<d-2:                  #Unitary into 2-level unitaries
#            U2x2=sp.identity(d,dtype="complex",format="csr")
#            I,J=indexes[1][0], i             #Get the index of first nonzero in column i
#            a, b = unitary[J,J], unitary[I,J]
#            r=sqrt(abs(a)**2+abs(b)**2)
#            c=a.conj()/r
#            s=-b.conj()/r            
#            #print("-------",I,J)
#            U2x2[J,J], U2x2[J,I] =c,-s
#            U2x2[I,J], U2x2[I,I] =s.conj(),c.conj()
#            #print(U2x2.toarray())
#            u2x2list.append(U2x2.getH())
#            unitary=sp.csr_matrix.dot(U2x2,unitary)
#            unitary.real[abs(unitary.real)<1e-14]=0             
#            unitary.imag[abs(unitary.imag)<1e-14]=0
#
#            column=unitary.getcol(i)        #Stores column i
#            indexes=np.transpose(column.nonzero())
#
#        U2x2=unitary.getH() 
#        u2x2list.append(U2x2.getH())
#        #unitary=sp.csr_matrix.dot(U2x2,unitary)
#    return u2x2list

def grad(A,x):                                #Millorable per H1 i H2. Solucio temporal. Swap per H1
    mat=A+A.transpose()    
                                                        #Es 2A millor que A+A.T?
    df=sp.csr_matrix.dot(mat,x)
    
    return df

#def grad2(A,x):                                #MSwap per H1
#    for i in range(nswaps):
#        
#    
#    return df

def E(A,x): 
                                        
    eigval=sp.csr_matrix.dot(A,x)
    E=sp.csr_matrix.dot(x.transpose(),eigval)
    

    return E.toarray()[0]

N=15
d=2**N

J=1
h=1                                   #J=0 superposa solucions
hamiltonian=Hising(N,J,h)


empty=sp.csr_matrix((d,1))
x=np.random.random((d,1))          #To get a quick ansatz
x=x/norm(x)
epsilon=1e-6
err=1/epsilon
gamma=sqrt(2)
ans=sp.csr_matrix(x)
En=0
gradsum=0
Elist=[]
start=time.time()
while err>epsilon:
    prev=ans
    for k in range(len(hamiltonian)):
        gradsum=gradsum+grad(hamiltonian[k],prev)

    
    ans=ans-gamma*gradsum
    ans=ans/sp.linalg.norm(ans)
    err=sp.linalg.norm(ans-prev)
    gradsum=0

end=time.time()    
print(end-start)
print(ans.toarray())

