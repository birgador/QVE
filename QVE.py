from pylab import *
import qgates as q
import time

Zero = np.array([[1],[0]])
One = np.array([[0],[1]])
N=3
Z2=np.kron(q.Z(),q.Z())              
H1=0
H2=0
for i in range(N):
    H1=H1+np.kron(np.kron(np.eye(2**i),q.X()),np.eye(2**(N-1-i)))
n=1    
for i in range(N-1):    
    H2=H2+np.kron(np.kron(np.eye(2**i),Z2),np.eye(2**(N-1-i-1)))
    
H=H1+H2

P=np.array([[1,1,1,1],[1,1j,-1,-1j],[1,-1,1,-1],[1,-1j,-1,1j]])
def QVE(H):
    U=H
    d=U.shape[0]
    list
    for i in range(d):
        if i==1:
            break
        for j in range(i+1,d):
            if U[j,i]!=0:
                U_i=np.eye(d,dtype="complex")
                a, b = U[i,i], U[j,i]
                r=sqrt(abs(a)**2+abs(b)**2)
                c=a.conj()/r
                s=-b.conj()/r
                U_i[i,i], U_i[i,j] =c,-s
                U_i[j,i], U_i[j,j] =s.conj(),c.conj()
                
                U=np.dot(U_i,U)
                print(U_i)
                print("_________",i)
                
    return U
K=QVE(P)/2


"""while True:
        U_i=np.eye(d)
        for i in range(d):
            if U[i,0]!=0 and i!=0:
                a, b = U[i,i], U[i,0]
                r=sqrt(abs(a)**2+abs(b)**2)
                c=a.conj()/r
                s=-b.conj()/r
                U_i[0,0], U_i[0,i] =c,-s
                U_i[i,0], U_i[i,i] =s.conj(),c.conj()"""
