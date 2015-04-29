from numpy import *
from ChargeQubits import *
from numpy.random import normal
import pickle
from mpi4py import MPI
from time import time
from qutip import *

comm = MPI.COMM_WORLD

full_params = pickle.load(open('fullParams.p','rb'))

my_param_list = []

for i,params in enumerate(full_params):
    if i%comm.size==comm.rank:
        my_param_list+=[params]


op_basis = [[qeye(2), sigmax(), sigmay(), sigmaz()]] * 2

U = 1./sqrt(2.)*tensor(qeye(2)-1.j*sigmax(),qeye(2))
U = spre(U) * spost(U.dag())
idealChi = qpt(U, op_basis)

def findChi(ep1,ep2,Delta1,Delta2,CouplingG):
    u = zeros([4 * 4, 4 * 4], dtype=complex)
    
    for n in range(4 * 4):
        psi0 = basis(4 * 4, n)
        rho0 = Qobj(vec2mat(psi0.full()))
        myQubits = Charge_Qubits(None,rho=rho0,masterEquation=True)
        myQubits.setHam(myQubits.twoQubitHam(ep1,ep2,Delta1,Delta2,CouplingG))
        myQubits.setExpectedHam(myQubits.twoQubitHam(0.,0.,Delta1,Delta2,CouplingG))
        
        myQubits.X1Rot(0.1,pi/2.)
        
        u[:, n] = mat2vec(myQubits.getRho().full()).T
    
    u = Qobj(u)
    
    return qpt(u, op_basis)



for params in my_param_list:
    (ep1,ep2,Delta1,Delta2,CouplingG) = params
    try:
       chi = findChi(ep1,ep2,Delta1,Delta2,CouplingG)
    except:
        f = open('error.log','a')
        f.write("Epsilon_1: %d , Epsilon_2: %d , Delta_1: %d , Delta_2: %d , CouplingG: %d failed to converge.\n" % params)
        f.close()
        continue
    fid = trace(dot(chi,idealChi))
    f = open('fidelity_data_'+str(comm.rank)+'.txt','a')
    f.write("Epsilon_1: %d , Epsilon_2: %d , Delta_1: %d , Delta_2: %d , CouplingG: %d , Fidelity: " % params + str(fid)+"\n")
    f.close()


