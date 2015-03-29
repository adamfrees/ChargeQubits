from numpy import *
from ChargeQubits import *

psi = Qobj(array([[1.],[1.],[0],[0.]])).unit()
#psi = Qobj(array([[0],[0],[1],[0]])).unit()

myQubits = Charge_Qubits(None,psi=psi)

myQubits.setHam(myQubits.twoQubitHam(0.,0.,11.,5.,3.))

#print myQubits.omega1Eff,myQubits.omega2Eff

#print myQubits.getPsi()
myQubits.wait(1.e-13)
myQubits.CNOT2(0.2)
for i in range(5):
    list = myQubits.stateList[i]
    timeList = myQubits.timeList[i]
    for j in range(0,500,10):
        timeU = array([[exp(1.j*myQubits.energies[0]*timeList[j]/myQubits.hbar),0,0,0],[0,exp(1.j*myQubits.energies[1]*timeList[j]/myQubits.hbar),0,0],[0,0,exp(1.j*myQubits.energies[2]*timeList[j]/myQubits.hbar),0],[0,0,0,exp(1.j*myQubits.energies[3]*timeList[j]/myQubits.hbar)]])
        myQubits.densityBarchart(ket2dm(list[j].transform(timeU)),'pictures/'+str(i)+'_minus3_'+str(j/10)+'.png','(1,1,0,0) Drive '+str(i+1))
#print myQubits.totalTime
#print myQubits.getPsi()
#a = myQubits.getPsi()*exp(-1.j*angle(myQubits.getPsi()[0]))
#print a
#print '+om1',myQubits.getPsi()*exp(1.j*myQubits.omega1Eff*myQubits.totalTime)
#print '-om1',myQubits.getPsi()*exp(-1.j*myQubits.omega1Eff*myQubits.totalTime)
#print '+om2',myQubits.getPsi()*exp(1.j*myQubits.omega2Eff*myQubits.totalTime)
#print '-om2',myQubits.getPsi()*exp(-1.j*myQubits.omega2Eff*myQubits.totalTime)
#print 'om1-om2',myQubits.getPsi()*exp(-1.j*myQubits.omega2Eff*myQubits.totalTime)*exp(1.j*myQubits.omega1Eff*myQubits.totalTime)
#print 'om2-om1',myQubits.getPsi()*exp(1.j*myQubits.omega2Eff*myQubits.totalTime)*exp(-1.j*myQubits.omega1Eff*myQubits.totalTime)
#for i,x in enumerate(myQubits.getPsi()):
#    print angle(x*exp(-1.j*angle(myQubits.getPsi()[0])))
#print myQubits.omega1Eff*myQubits.totalTime % (2.*pi)
#print myQubits.omega2Eff*myQubits.totalTime % (2.*pi)
#print (myQubits.omega1Eff-myQubits.omega2Eff)*myQubits.totalTime % (2.*pi)
#print (myQubits.omega2Eff-myQubits.omega1Eff)*myQubits.totalTime % (2.*pi)
#print (myQubits.omega2Eff+myQubits.omega1Eff)*myQubits.totalTime % (2.*pi)

#for field in linspace(0.1,0.5,40):
#    psi = Qobj(array([[1],[1],[1],[1]])).unit()
#    myQubits = Charge_Qubits(None,psi=psi)
#    myQubits.setHam(myQubits.twoQubitHam(0.,0.,11.,5.,0.5))
#    myQubits.CNOT2(field)
#    print myQubits.totalTime, (myQubits.getPsi()[1]*exp(-1.j*angle(myQubits.getPsi()[0]))).real
    

#for i in range(5):
#    myQubits.X2Rot(0.1,pi)
#    print myQubits.getPsi()
#    myQubits.X1Rot(0.1,pi)
#    print myQubits.getPsi()

#for i in range(10):
#    myQubits.wait(2.*pi/myQubits.omega2Eff)
#    print myQubits.totalTime,myQubits.getPsi()*exp(1.j*myQubits.energies[0]*myQubits.totalTime/myQubits.hbar)


