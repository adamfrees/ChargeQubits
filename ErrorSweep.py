from numpy import *
from ChargeQubits import *
from numpy.random import normal
import pickle
from mpi4py import MPI
import time

psi = Qobj(array([[1.],[0.],[0.],[0.]])).unit()
#psi = Qobj(array([[0],[0],[1],[0]])).unit()

fullDeltas = [(i,j) for i in linspace(10.,100.,10) for j in linspace(10.,100.,10)]

unfinishedDeltas = []

for deltas in fullDeltas:
    delta1 = deltas[0]
    delta2 = deltas[1]
    try:
        f = open('data/DeltaRun/g_50_d1_'+str(int(delta1))+'_d2_'+str(int(delta2))+'.p','rb')
    except IOError:
        unfinishedDeltas+=[(delta1,delta2)]
        continue

comm = MPI.COMM_WORLD

myDeltas = []

for i,pair in enumerate(unfinishedDeltas):
    if i%comm.size==comm.rank:
        myDeltas+=[pair]

for deltas in myDeltas:
    startTime = time.time()
    Delta1 = deltas[0]
    Delta2 = deltas[1]
    CouplingG = 50.
    print "Core",comm.rank,"beginning parameters",deltas
    
    myQubits = Charge_Qubits(None,psi=psi)
    myQubits.setHam(myQubits.twoQubitHam(0.,0.,Delta1,Delta2,CouplingG))
    
    try:
        myQubits.X1Rot(0.1,pi/2.)
    except IndexError:
        f = open('data/DeltaRun/error.log','a')
        f.write("Delta1: %d and Delta2: %d failed to converge.\n" % deltas)
        f.close()
        continue
    
    targetDensity = ket2dm(myQubits.getPsi())
    
    data = []
    
    #for i in range(1000):
    for ep1 in linspace(-5.,5.,26):
        for ep2 in linspace(-5.,5.,26):
            myQubits = Charge_Qubits(None,psi=psi)
            #ep1 += normal(scale=0.05)
            #ep2 += normal(scale=0.05)
            myQubits.setHam(myQubits.twoQubitHam(ep1,ep2,Delta1,Delta2,CouplingG))
            myQubits.setExpectedHam(myQubits.twoQubitHam(0.,0.,Delta1,Delta2,CouplingG))
            
            myQubits.X1Rot(0.1,pi/2.)
            
            randDensity = ket2dm(myQubits.getPsi())
            
            #print ep1,ep2,fidelity(targetDensity,randDensity)
            data += [[ep1,ep2,fidelity(targetDensity,randDensity)]]
            #data += [[ep1,ep2,myQubits.typeA]]
    print 'Core',comm.rank,'finished a run in ',time.time()-startTime,'s.'
    pickle.dump(array(data),open('data/DeltaRun/g_50_d1_'+str(int(Delta1))+'_d2_'+str(int(Delta2))+'.p','wb'))

