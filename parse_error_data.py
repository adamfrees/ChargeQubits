from numpy import *
import pickle

fullDeltas = [(i,j) for i in linspace(10.,100.,10) for j in linspace(10.,100.,10)]

delta_scatter = []
a = []

for deltas in fullDeltas:
    delta1 = deltas[0]
    delta2 = deltas[1]
    try:
        f = open('data/DeltaRun/g_50_d1_'+str(int(delta1))+'_d2_'+str(int(delta2))+'.p','rb')
    except IOError:
        a+=[(delta1,delta2)]
        continue
    errors = pickle.load(f)
    f.close()
    delta_scatter += [[delta1,delta2,mean(errors[:,2])]]
pickle.dump(array(delta_scatter),open('data/delta_scatter.p','wb'))
print a
print len(a)