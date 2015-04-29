from numpy import *
import pickle
from random import gauss

fullParams = []

for i in linspace(10.,100.,10):
    for j in linspace(10.,100.,10):
        for k in range(100):
            fullParams += [(gauss(0.,5.),gauss(0.,5.),i,j,50.)]

pickle.dump(fullParams,open('fullParams_2.p','wb'))