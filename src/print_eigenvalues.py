import numpy as np
from numpy import linalg as LA

fileD = open("D.txt", "r")
D = []
for line in fileD:
    D.append([float(x) for x in line.split()])

fileD.close()

D = np.array(D)

DtD = np.dot(D.T, D)

evals, V = LA.eigh(DtD)

idx = evals.argsort()[::-1]   
evals = evals[idx]
V = V[:,idx]

print len(evals)
print evals[0]/evals[len(evals)-1]
print evals