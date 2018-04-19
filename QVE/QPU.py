from pylab import *
import qgates as q
import time
from scipy.linalg import expm,logm,polar,inv as expm,logm,polar,inv


def dagger(M):
    M=M.T
    M=M.conj()
    return M
    
def fastprod(u2x2,i,j,targ,d):                  #Falta arreglar 1e-17
    vk=linspace(0,1,2,dtype="complex")          #Solucio cutre a float to com
    for k in range(d):
        vk[0],vk[1]=targ[i,k],targ[j,k]
        vr=np.dot(u2x2,vk)
        targ[i,k],targ[j,k]=vr[0],vr[1]
    return targ
    
def Tx2(H,d):
    U=H
    actionlist=[]
    U2x2=np.empty([2,2],dtype="complex")
    #Bucle that decomposes unitary U into 2-level unitary matrices
    for i in range(d):
        """if i==1:
            break"""
        for j in range(i+1,d):                  #Unitary into 2-level unitaries
            if U[j,i]!=0:
                if j==d-2 and i ==d-2:
                    U2x2[0,0], U2x2[0,1] =U[d-2][d-2].conj(),U[d-1][d-2].conj()
                    U2x2[1,0], U2x2[1,1] =U[d-2][d-1].conj(),U[d-1][d-1].conj()
                    actionlist.append((U2x2,i,j))
                #U_i=np.eye(d,dtype="complex")
                a, b = U[i,i], U[j,i]
                r=sqrt(abs(a)**2+abs(b)**2)
                c=a.conj()/r
                s=-b.conj()/r

                #U_i[i,i], U_i[i,j] =c,-s
                #U_i[j,i], U_i[j,j] =s.conj(),c.conj()
                U2x2[0,0], U2x2[0,1] =c,-s
                U2x2[1,0], U2x2[1,1] =s.conj(),c.conj()
                actionlist.append((dagger(U2x2),i,j))
                U=fastprod(U2x2,i,j,U,d)
                #print(U)
                #print("___")
    #actionlist.append(U)
    #print(np.dot(U_i,U))
                
    return actionlist
"""np.set_printoptions(precision=4,suppress=True)        #Ho poso perque m'ajuda al programar
print(end-start)
alist=2x2(P)"""


