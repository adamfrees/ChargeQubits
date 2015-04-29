from qutip import *
from numpy import *
import pickle
from time import time



delta1 = 10.
delta2 = 80.
g = 50.
eps1 = 0.
eps2 = 0.


epsilonPairs = [(i,j) for i in linspace(-10.,10.,5) for j in linspace(-10.,10.,5)]


strength = 0.1
theta = pi / 2.
phi =pi / 2.

# pre-calculate the necessary operators
sx1 = tensor(sigmax(),qeye(2)); sx2 = tensor(qeye(2),sigmax())
sz1 = tensor(sigmaz(),qeye(2)); sz2 = tensor(qeye(2),sigmaz())
szz = tensor(sigmaz(),sigmaz())
# collapse operators
c_ops = []

# setup time-dependent Hamiltonian (list-string format)
H0 = delta1 * sx1 + delta2 * sx2 + g * szz
ham_eig = H0.eigenstates()
#H0 = spre(H0) * spost(H0.dag())
H1 = [sz1, 'eps1 / 2.0 + A1 * sin(w * t + phi)']
H2 = [sz2, 'eps2 / 2.0 + A2 * sin(w * t + phi)']
H_td = [H0,H1,H2]



energies = ham_eig[0]
w = energies[3] + energies[2]
hamIsTypeA = sz1.transform(ham_eig[1])[0,2] * sz1.transform(ham_eig[1])[1,3]>=0

print energies[3] + energies[2]
print energies[3] - energies[2]
print 2.*energies[2]
print 2.*energies[3]

if hamIsTypeA:
    A1 = strength / 2.
    A2 = 0.
    T = theta / 2. / abs(0.25*strength*sz1.transform(ham_eig[1])[0,2])
else:
    A1 = 0.
    A2 = strength / 2.
    T = theta / 2. / abs(0.25*strength*sz2.transform(ham_eig[1])[0,2])
Hargs = {'w': w, 'eps1': eps1, 'A1': A1, 'eps2': eps2, 'A2': A2, 'phi': phi}

op_basis = [[qeye(2), sigmax(), sigmay(), sigmaz()]] * 2
op_label = [["i", "x", "y", "z"]] * 2
# ODE settings (for reusing list-str format Hamiltonian)
opts = Options(rhs_reuse=True)
opts.nsteps = 10000000
opts.max_step = 1e-4
# pre-generate RHS so we can use parfor
rhs_clear()
rhs_generate(H_td, c_ops, Hargs, name='lz_func')

# a task function for the for-loop parallelization:
# the m-index is parallelized in loop over the elements of p_mat[m,n]

Hargs = {'w': w, 'eps1': eps1, 'A1': A1, 'eps2': eps2, 'A2': A2, 'phi': phi}
start = time()

timeU = 1.j * T * H0.transform(ham_eig[1])
timeU = timeU.expm()

U = propagator(H_td, T, c_ops, Hargs, opts) #<- IMPORTANT LINE
print 'time: ', time()-start
U = U.transform(ham_eig[1])
U = U * timeU
print U
pickle.dump(U,open('u_test_4.p','wb'))

'''
def task(args):
    m, (eps1,eps2) = args
    Hargs = {'w': w, 'eps1': eps1, 'A1': A1, 'eps2': eps2, 'A2': A2, 'phi': phi}
    U = propagator(H_td, T, c_ops, Hargs, opts) #<- IMPORTANT LINE
    U = spre(U) * spost(U.dag())
    chi = qpt(U, op_basis)
    print m, 'complete.'
    return [m, chi]


chi_list = parfor(task, enumerate(epsilonPairs),num_cpus=12)

pickle.dump(chi_list,open('chi_test.p','wb'))'''
