from numpy import *
from qutip import *
import matplotlib.pyplot as plt
import pickle
U = pickle.load(open('propagator_test_off_high.p','rb'))
U = 1./sqrt(2.)*tensor(qeye(2)-1.j*sigmax(),qeye(2))
U = Qobj(U)
print U

U = spre(U) * spost(U.dag())

op_basis = [[qeye(2), sigmax(), sigmay(), sigmaz()]] * 2
op_label = [["i", "x", "y", "z"]] * 2
chi = qpt(U, op_basis)
#fig = qpt_plot_combined(chi, op_label, r'$X_1$ $pi/2$')
print chi
#plt.show()
