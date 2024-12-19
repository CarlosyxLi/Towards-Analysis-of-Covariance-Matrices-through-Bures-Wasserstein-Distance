##### Bures-Wasserstein Metric #####

import numpy as np
from scipy.linalg import sqrtm
from scipy.linalg import logm
from scipy.linalg import expm

# Affine Invariant Distance
def AI_dist(A,B):
    temp = np.dot(np.linalg.inv(A),B)
    lambdas, _ = np.linalg.eig(temp)
    output = np.real(np.sqrt(np.sum(np.log(np.real(lambdas))**2)+0j)) # accommodate negative values
    return(output)

# Affine Invariant Logarithm Map
def AI_log(A,B):
    temp = logm(sqrtm(np.linalg.inv(A)).dot(B).dot(sqrtm(np.linalg.inv(A))))
    output = sqrtm(A).dot(temp).dot(sqrtm(A))
    return(output)

# Affine Invariant Exponential Map
def AI_exp(A,B):
    temp = expm(sqrtm(np.linalg.inv(A)).dot(B).dot(sqrtm(np.linalg.inv(A))))
    output = sqrtm(A).dot(temp).dot(sqrtm(A))
    return(output)