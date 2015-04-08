from numpy import *
from ChargeQubits import *
from numpy.random import normal
import pickle

psi = Qobj(array([[1.],[0.],[0.],[0.]])).unit()
#psi = Qobj(array([[0],[0],[1],[0]])).unit()

Delta1 = 30.
Delta2 = 40.
CouplingG = 60.

myQubits = Charge_Qubits(None,psi=psi)
myQubits.setHam(myQubits.twoQubitHam(0.,0.,Delta1,Delta2,CouplingG))

myQubits.X1Rot(0.1,pi/2.)

targetDensity = ket2dm(myQubits.getPsi())

data = []

for ep2 in linspace(-0.2,-0.19,2):
    myQubits = Charge_Qubits(None,psi=psi)
    ep1 = 0.1 #normal(scale=5.)
    myQubits.setHam(myQubits.twoQubitHam(ep1,ep2,Delta1,Delta2,CouplingG))
    myQubits.setExpectedHam(myQubits.twoQubitHam(0.,0.,Delta1,Delta2,CouplingG))
    
    myQubits.X1Rot(0.1,pi/2.)
    
    randDensity = ket2dm(myQubits.getPsi())
    
    print ep1,ep2,fidelity(targetDensity,randDensity)
    data += [[ep1,ep2,fidelity(targetDensity,randDensity)]]
    list = myQubits.stateList[0]
    timeList = myQubits.timeList[0]
    #for j in range(0,500,10):
    #    timeU = array([[exp(1.j*myQubits.energies[0]*timeList[j]/myQubits.hbar),0,0,0],[0,exp(1.j*myQubits.energies[1]*timeList[j]/myQubits.hbar),0,0],[0,0,exp(1.j*myQubits.energies[2]*timeList[j]/myQubits.hbar),0],[0,0,0,exp(1.j*myQubits.energies[3]*timeList[j]/myQubits.hbar)]])
    #    #timeU = eye(4,dtype=complex)
    #    myQubits.densityBarchart(ket2dm(list[j].transform(timeU)),'/home/frees/Dropbox/_UW/Detuning_shift/Fidelity_tests/pictures/typeA_Xgate_ep2_'+str(ep2)+'_'+str(j/10)+'.png','')

#pickle.dump(array(data),open('scatter_small_grid.p','wb'))

