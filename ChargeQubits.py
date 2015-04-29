from numpy import *
from scipy.linalg import eigh
from scipy.integrate import ode
from math import cos, pi
from numpy.linalg import norm
from qutip import *
import matplotlib.pyplot as plt
import matplotlib as mpl
from random import random
qutip.settings.qutip_graphics = "YES"
from qutip.visualization import *

class Charge_Qubits:
    '''
    Units:
    Energies are in eV
    Everything else is in SI - m, s, etc.
    Basis Convention - 00,01,10,11, with energies -lambda_1, -lambda_2, +lambda_2, +lambda_1 (at 0 detunings)
                       omega1Eff = lambda_1 + lambda_2, omega2Eff = lambda_1-lambda_2
    '''
    def __init__(self, ham, rho=None, psi=None,visual=False,masterEquation=False):
        self.hbar = 6.58211928e-16 #eV*s
        self.ham = ham
        self.expectedHam = None
        self.psi = psi
        self.rho = rho
        if psi!= None:
            self.rho = ket2dm(psi)
        self.energies = None
        self.U = None
        self.expectedEnergies = None
        self.expectedU = None
        self.omega1Eff = None
        self.omega2Eff = None
        self.masterEquation = masterEquation
        self.typeA = None
        self.stateList = []
        self.timeList = []
        if ham!=None:
            ham_eig = ham.eigenstates()
            self.energies = ham_eig[0]
            self.U = ham_eig[1]
            self.omega1Eff = (self.energies[3]+self.energies[2])/self.hbar
            self.omega2Eff = (self.energies[3]-self.energies[2])/self.hbar
            self.typeA = self.checkType() # Depending on the Deltas and g, there are two ways that the off-resonant pulses can affect the system.
        self.visual = visual
        if visual:
            global mlab
            from mayavi import mlab
        self.totalTime = 0
    
    def getRho(self):
        timeU = array([[exp(1.j*self.energies[0]*self.totalTime/self.hbar),0,0,0],[0,exp(1.j*self.energies[1]*self.totalTime/self.hbar),0,0],[0,0,exp(1.j*self.energies[2]*self.totalTime/self.hbar),0],[0,0,0,exp(1.j*self.energies[3]*self.totalTime/self.hbar)]])
        return self.rho.transform(timeU)
    
    def setRho(self,rho):
        self.rho = rho
    
    def getPsi(self):
        timeU = array([[exp(1.j*self.energies[0]*self.totalTime/self.hbar),0,0,0],[0,exp(1.j*self.energies[1]*self.totalTime/self.hbar),0,0],[0,0,exp(1.j*self.energies[2]*self.totalTime/self.hbar),0],[0,0,0,exp(1.j*self.energies[3]*self.totalTime/self.hbar)]])
        return self.psi.transform(timeU)
    
    def setPsi(self,psi,setRho=True):
        self.psi = psi
        if setRho:
            self.rho = ket2dm(psi)
    
    def checkType(self):
        return self.pertMatrix(0).transform(self.U)[0,2]*self.pertMatrix(0).transform(self.U)[1,3]>=0
    
    def getHam(self):
        return self.ham
    
    def setHam(self,ham):
        self.ham = ham
        ham_eig = ham.eigenstates()
        self.energies = ham_eig[0]
        self.U = ham_eig[1]
        self.expectedU = self.U
        self.omega1Eff = (self.energies[3]+self.energies[2])/self.hbar
        self.omega2Eff = (self.energies[3]-self.energies[2])/self.hbar
        self.typeA = self.checkType()
    
    def setExpectedHam(self,ham):
        self.expectedHam = ham
        ham_eig = ham.eigenstates()
        self.expectedEnergies = ham_eig[0]
        self.expectedU = ham_eig[1]
        self.U = self.expectedU
        self.omega1Eff = (self.expectedEnergies[3]+self.expectedEnergies[2])/self.hbar
        self.omega2Eff = (self.expectedEnergies[3]-self.expectedEnergies[2])/self.hbar
        self.typeA = self.pertMatrix(0).transform(self.expectedU)[0,2]*self.pertMatrix(0).transform(self.expectedU)[1,3]>=0
    
    def getEnergies(self):
        return self.energies
    
    def twoQubitHam(self,eps1,eps2,delta1,delta2,g):
        '''
        4x4 matrix that describes 2 capacitively coupled qubits
        '''
        H1 = matrix([[eps1/2.,0.,delta1,0.],[0.,eps1/2.,0.,delta1],[delta1,0.,-eps1/2.,0.],[0.,delta1,0.,-eps1/2.]])
        H2 = matrix([[eps2/2.,delta2,0.,0.],[delta2,-eps2/2.,0.,0.],[0.,0.,eps2/2.,delta2],[0.,0.,delta2,-eps2/2.]])
        cross = matrix([[g,0.,0.,0.],[0.,-g,0.,0.],[0.,0.,-g,0.],[0.,0.,0.,g]])
        return  Qobj(H1+H2+cross)
    
    def wait(self,time):
        steps = 100
        if self.masterEquation:
            result = mesolve(self.ham.transform(self.U),self.rho,linspace(0.,time/self.hbar,steps),[],[]) #time is divided by hbar because of mesolve convention.
            self.rho = result.states[steps-1]
        else:
            result = mesolve(self.ham.transform(self.U),self.psi,linspace(0.,time/self.hbar,steps),[],[]) #time is divided by hbar because of mesolve convention.
            self.psi = result.states[steps-1]
        self.totalTime += time
    
    def perturb(self,pert,time,freq,phase):
        steps = 500
        opts = Options(nsteps = 10000)
        args = {'w':  freq*self.hbar,'phi': phase}
        if self.masterEquation:
            result = mesolve([self.ham.transform(self.U),[pert,'sin(w*t+phi)']],self.rho,linspace(0.,time/self.hbar,steps),[],[],args,options=opts) #time is divided by hbar because of mesolve convention.
            self.rho = result.states[steps-1]
        else:
            result = mesolve([self.ham.transform(self.U),[pert,'sin(w*t+phi)']],self.psi,linspace(0.,time/self.hbar,steps),[],[],args,options=opts) #time is divided by hbar because of mesolve convention.
            self.psi = result.states[steps-1]
        self.stateList += [result.states]
        self.timeList += [linspace(0.,time,steps)]
        self.totalTime += time
    
    def pulse(self,attr,strength,time,freq,phase):
        '''
        attr:
        0: epsilon_1
        1: epsilon_2
        2: delta_1
        3: delta_2
        '''
        pert = 0.5*strength*self.pertMatrix(attr).transform(self.U)
        self.perturb(pert,time,freq,phase)
    
    def pertMatrix(self,attr):
        '''
        attr:
        0: epsilon_1
        1: epsilon_2
        2: delta_1
        3: delta_2
        '''
        if attr==0:
            pert = matrix([[1.,0.,0.,0.],[0.,1.,0.,0.],[0.,0.,-1.,0.],[0.,0.,0.,-1.]])
        elif attr==1:
            pert = matrix([[1.,0.,0.,0.],[0.,-1.,0.,0.],[0.,0.,1.,0.],[0.,0.,0.,-1.]])
        elif attr==2:
            pert = matrix([[0.,0.,1.,0.],[0.,0.,0.,1.],[1.,0.,0.,0.],[0.,1.,0.,0.]])
        elif attr==3:
            pert = matrix([[0.,1.,0.,0.],[1.,0.,0.,0.],[0.,0.,0.,1.],[0.,0.,1.,0.]])
        else:
            pert = None
            print "Error: invalid attribute"
        return Qobj(pert)
    
    def X1Rot(self,strength,theta):
        
        waitTime = 2.*pi/self.omega1Eff - (self.totalTime % (2.*pi/self.omega1Eff))
        self.wait(waitTime)
        if self.typeA:
            a2 = abs(0.25*strength*self.pertMatrix(0).transform(self.expectedU)[0,2])
            time = 0.5*self.hbar*theta/a2
            self.pulse(0,strength,time,self.omega1Eff,pi/2.)
        else:
            c2 = abs(0.25*strength*self.pertMatrix(1).transform(self.expectedU)[0,2])
            time = 0.5*self.hbar*theta/c2
            self.pulse(1,strength,time,self.omega1Eff,pi/2.)
        if not self.masterEquation: self.psi = self.psi*exp(1.j*theta/2.) # This is purely bookkeeping, psi obviously doesn't depend on global phase
    
    def Y1Rot(self,strength,theta):
        waitTime = 2.*pi/self.omega1Eff - (self.totalTime % (2.*pi/self.omega1Eff))
        self.wait(waitTime)
        if self.typeA:
            a2 = 0.25*strength*self.pertMatrix(0).transform(self.expectedU)[0,2]
            time = abs(0.5*self.hbar*theta/a2)
            self.pulse(0,strength,time,self.omega1Eff,0.)
        else:
            c2 = 0.25*strength*self.pertMatrix(1).transform(self.expectedU)[0,2]
            time = abs(0.5*self.hbar*theta/c2)
            self.pulse(1,strength,time,self.omega1Eff,0.)
        if not self.masterEquation: self.psi = self.psi*exp(1.j*theta/2.) # This is purely bookkeeping, psi obviously doesn't depend on global phase
    
    def X2Rot(self,strength,theta):
        waitTime = 2.*pi/self.omega2Eff - (self.totalTime % (2.*pi/self.omega2Eff))
        self.wait(waitTime)
        if self.typeA:
            c1 = 0.25*strength*self.pertMatrix(1).transform(self.expectedU)[0,1]
            time = abs(0.5*self.hbar*theta/c1)
            self.pulse(1,strength,time,self.omega2Eff,pi/2.)
        else:
            a1 = 0.25*strength*self.pertMatrix(0).transform(self.expectedU)[0,1]
            time = abs(0.5*self.hbar*theta/a1)
            self.pulse(0,strength,time,self.omega2Eff,pi/2.)
        if not self.masterEquation: self.psi = self.psi*exp(-1.j*theta/2.) # This is purely bookkeeping, psi obviously doesn't depend on global phase
    
    def Y2Rot(self,strength,theta):
        waitTime = 2.*pi/self.omega2Eff - (self.totalTime % (2.*pi/self.omega2Eff))
        self.wait(waitTime)
        if self.typeA:
            c1 = 0.25*strength*self.pertMatrix(1).transform(self.expectedU)[0,1]
            time = abs(0.5*self.hbar*theta/c1)
            self.pulse(1,strength,time,self.omega2Eff,0.)
        else:
            a1 = 0.25*strength*self.pertMatrix(0).transform(self.expectedU)[0,1]
            time = abs(0.5*self.hbar*theta/a1)
            self.pulse(0,strength,time,self.omega2Eff,0.)
        if not self.masterEquation: self.psi = self.psi*exp(-1.j*theta/2.) # This is purely bookkeeping, psi obviously doesn't depend on global phase
    
    def CNOT2(self,strength):
        waitTime = 2.*pi/self.omega2Eff - (self.totalTime % (2.*pi/self.omega2Eff))
        self.wait(waitTime)
        myQubits.totalTime =0.
        if self.typeA:
            c2 = 0.25*strength*self.pertMatrix(1).transform(self.expectedU)[0,2]
            #c2 = 0.25*self.pertMatrix(1).transform(self.U)[0,2]
            time = abs(0.25*self.hbar*pi/c2)
            self.pulse(1,strength,time,self.omega1Eff,pi/2.)
        if not self.masterEquation: self.psi = self.psi*exp(1.j*pi/4.) # This is purely bookkeeping, psi obviously doesn't depend on global phase
        self.X1Rot(strength,3.*pi/2.)
        self.Z2Rot(strength,3.*pi*2.)
        if not self.masterEquation: self.psi = -self.psi*1.j # This is purely bookkeeping, psi obviously doesn't depend on global phase
    
    def Z1Rot(self,strength,theta):
        self.Y1Rot(strength,pi/2.)
        self.X1Rot(strength,theta)
        self.Y1Rot(strength,3.*pi/2.)
    
    def Z2Rot(self,strength,theta):
        self.Y2Rot(strength,pi/2.)
        self.X2Rot(strength,theta)
        self.Y2Rot(strength,3.*pi/2.)
    
    def densityBarchart(self,rho,fname,title):
        xlabels = ["00","01","10","11"]
        #fig, ax = matrix_histogram(0.5*(self.rho+self.rho.dag()), xlabels, xlabels, limits=[-1,1])
        fig, ax = matrix_histogram_complex(rho, xlabels, xlabels, limits=[-1,1])
        ax.view_init(azim=-55, elev=45)
        ax.set_title(title)
        plt.savefig(fname)
            