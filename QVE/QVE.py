"""
PSEUDOCODI
"""
from CPU import *
from QPU import *
import qgates as q
from scipy.linalg import expm as expm
import time

Zero = np.array([[1],[0]])
One = np.array([[0],[1]])
N=8
Z2=np.kron(q.Z(),q.Z()) 
H1=np.zeros(2**N,dtype="complex")*(-1)
H2=np.zeros(2**N,dtype="complex")*(-1)


#Test Hamiltonian@Ising
for i in range(N):
    H1=H1+np.kron(np.kron(np.eye(2**i),q.X()),np.eye(2**(N-1-i)))   
for i in range(N-1):
    H2=H2+np.kron(np.kron(np.eye(2**i),Z2),np.eye(2**(N-1-i-1)))
    

h=H1+H2
u=expm(1j*h)
d=u.shape[0]
ansatz=np.empty([d,1])
ansatz=ansatz/norm(ansatz)
start=time.time()
matlist=Tx2(u,d)          #matrix list with 2x2 level matrixes metode 1 per energia
precision=1e-9
diff=ansatz+1
#E0=0
#ansatz0=ansatz-ansatz

#Decompose H for quick E

while norm(diff)>precision:
    #En=E(matlist,ansatz)
    new=newAnsatz(ansatz,h)
    diff=abs(new-ansatz)
    ansatz=new
    #print(norm(diff))
end=time.time()
print(new)
print(end-start)