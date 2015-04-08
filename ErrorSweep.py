from numpy import *
from ChargeQubits import *
from numpy.random import normal
import pickle

psi = Qobj(array([[1.],[0.],[0.],[0.]])).unit()
#psi = Qobj(array([[0],[0],[1],[0]])).unit()

Delta1 = 40.
Delta2 = 40.
CouplingG = 60.

myQubits = Charge_Qubits(None,psi=psi)
myQubits.setHam(myQubits.twoQubitHam(0.,0.,Delta1,Delta2,CouplingG))

myQubits.X1Rot(0.1,pi/2.)

targetDensity = ket2dm(myQubits.getPsi())

data = []

#for i in range(1000):
for ep1 in linspace(-5.,5.,51):
    for ep2 in linspace(-5.,5.,51):
        myQubits = Charge_Qubits(None,psi=psi)
        #ep1 += normal(scale=0.05)
        #ep2 += normal(scale=0.05)
        myQubits.setHam(myQubits.twoQubitHam(ep1,ep2,Delta1,Delta2,CouplingG))
        myQubits.setExpectedHam(myQubits.twoQubitHam(0.,0.,Delta1,Delta2,CouplingG))
        
        myQubits.X1Rot(0.1,pi/2.)
        
        randDensity = ket2dm(myQubits.getPsi())
        
        print ep1,ep2,fidelity(targetDensity,randDensity)
        data += [[ep1,ep2,fidelity(targetDensity,randDensity)]]
        #data += [[ep1,ep2,myQubits.typeA]]

pickle.dump(array(data),open('scatter_grid_fixed_40.p','wb'))

