##### Bures-Wasserstein Metric #####

import numpy as np
from scipy.linalg import sqrtm

# Bures-Wasserstein Distance
def BW_dist(A,B):
    temp = np.real(sqrtm(np.dot(A,B)))
    output = np.real(np.sqrt((np.trace(A+B)-2*np.trace(temp))+0j)) # accommodate negative values
    return(output)

# Bures-Wasserstein Logarithm Map
def BW_log(A,B):
    temp1 = np.real(sqrtm(np.dot(A,B)))
    temp2 = np.real(sqrtm(np.dot(B,A)))
    output = temp1+temp2-2*A
    return(output)

# Bures-Wasserstein Exponential Map
def BW_exp(A,B):
    D, U = np.linalg.eigh(A)
    Xu = np.transpose(U).dot(B).dot(U)
    W = np.zeros((len(D), len(D)))
    for i in range(len(D)):
        for j in range(len(D)):
            W[i,j] = 1/(D[i]+D[j])
    output = A+B+U.dot(W*Xu).dot(np.diag(D)).dot(W*Xu).dot(np.transpose(U))
    return(output)