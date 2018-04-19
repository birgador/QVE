from pylab import *
from scipy import linalg

"""1-qubit gates"""
def Id():
    Id=np.eye(2)
    return Id

def Hadamard():
    H0 = np.array([[1, 1],
                   [1, -1]])
    return H0

def X():
    H0 = np.array([[0, 1],
                   [1, 0]])
    return H0

def Y():
    H0 = np.array([[0, -1j],
                   [1j, 0]])
    return H0

def Z():
    H0 = np.array([[1, 0],
                   [0, -1]])
    return H0

"""2-qubit gates"""

def CNOT():
    CNOT = np.array([[1, 0, 0, 0],
                     [0, 1, 0, 0],
                     [0, 0, 0, 1],
                     [0, 0, 1, 0]])
    return CNOT

def CU(U):
    a = np.array([[1, 0],
                 [0, 0]])
                 
    b = np.array([[0, 0],
                 [0, 1]])
    return np.kron(a,Id())+np.kron(b,U)

H1=np.kron(X(),Id())
H2=np.kron(Id(),X())
H3=np.kron(Z(),Z())
H3=0

H=H1+H2+H3

la,v = linalg.eig(H)
#print(v)