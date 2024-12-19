##### 
# n - sample size
# p - dimension
# l - length
import numpy as np
from BW_metric import BW_dist
from BW_metric import BW_log
from BW_metric import BW_exp
from AI_metric import AI_dist
from AI_metric import AI_log
from AI_metric import AI_exp


##### Covariance Matrices #####
def x2cov(data, centering = True):
    # input: signal data [N,p,l]
    # output: covariance matrices (psd) [N,p,p]
    output = np.zeros((data.shape[0], data.shape[1], data.shape[1]), dtype = np.float32)
    
    for trial_idx in range(data.shape[0]):
        temp = data[trial_idx,:,:].T
        if centering == True:
            output[trial_idx,:,:] = np.cov(temp, rowvar=False)
        else:
            # ex = np.mean(temp, axis = 0
            output[trial_idx,:,:] = (temp).T.dot(temp)/temp.shape[0] # (temp-ex).T.dot(temp-ex)/temp.shape[0]
    
    return output

##### Correlation Matrices #####
def x2corr(data):
    # input: signal data [N,p,l]
    # output: covariance matrices (psd) [N,p,p]
    output = np.zeros((data.shape[0], data.shape[1], data.shape[1]), dtype = np.float32)
    
    for trial_idx in range(data.shape[0]):
        temp = data[trial_idx,:,:].T
        output[trial_idx,:,:] = np.corrcoef(temp, rowvar=False)
    
    return output

##### Matrix Regularization #####
def matrix_regularization(s, lambda_value):
    # input: correlation/covariance matrices [N,p,p]
    # output: regularized correlation/covariance matrices [N,p,p]
    output = np.zeros(s.shape, dtype = np.float32)
    
    for trial_idx in range(s.shape[0]):
        # R'=1/(1+λ)*(R+λI)
        output[trial_idx] = (s[trial_idx]+(lambda_value*np.identity(s.shape[1])))/(1+lambda_value)
        
    return output

##### BW Projection Mean #####
def bw_projection_mean(x, eps=0.001, verbose = False):
    # input: psd matrices [N,p,p]
    # output: mean matrix (psd) [p,p]
    mean_new = np.mean(x, axis = 0)
    dist_mean = 10 # initial value > eps
    k = 0
    
    while dist_mean > eps:
    # stop when distance between two successive means exceeds eps
        
        U, s, VT = np.linalg.svd(mean_new, full_matrices=True)
        mean_old = (np.transpose(VT).dot(np.diag(s)).dot(VT)+np.transpose(np.transpose(VT).dot(np.diag(s)).dot(VT)))/2

        mean_tangent = np.mean(np.array([BW_log(mean_old,x[i]) for i in range(x.shape[0])]), axis = 0)
        mean_new = np.real(BW_exp(mean_old, mean_tangent))
        dist_mean = BW_dist(mean_new, mean_old)
        
        # update after each iteration
        k += 1
        mean_old = mean_new
        
        # track convergence process
        if verbose == True:
            print("Iter", k, 
                  ", dist_mean=", np.round(dist_mean,7), 
                  sep="")
        
    return mean_new

##### AI Projection Mean #####
def ai_projection_mean(x, eps=0.001, verbose = False):
    # input: psd matrices [N,p,p]
    # output: mean matrix (psd) [p,p]
    mean_new = np.mean(x, axis = 0)
    dist_mean = 10 # initial value > eps
    smaller_step = True
    k = 0
    
    while dist_mean>eps and smaller_step: 
    # stop when distance between two successive means exceeds eps, or step size starts to get bigger
    
        U, s, VT = np.linalg.svd(mean_new, full_matrices=True)
        mean_old = (np.transpose(VT).dot(np.diag(s)).dot(VT)+np.transpose(np.transpose(VT).dot(np.diag(s)).dot(VT)))/2

        mean_tangent = np.mean(np.array([AI_log(mean_old,x[i]) for i in range(x.shape[0])]), axis = 0)
        mean_new = np.real(AI_exp(mean_old, mean_tangent))
        dist_mean = AI_dist(mean_new, mean_old)
        
        # stop the change in steps size become smaller than 0.0001
        smaller_step = True if k==0 else (last_dist-dist_mean)>0.0001

        # track convergence process
        if verbose == True:
            print("Iter", k, 
                  ", Smaller step? ", "NA" if k==0 else smaller_step, 
                  ", dist_mean=", np.round(dist_mean,7), 
                  sep="")

        # update after each iteration
        k += 1
        last_dist = dist_mean
        mean_old = mean_new if smaller_step else mean_old # only update the mean when step size is getting smaller

    return mean_old

##### Logarithm Mapping All Matrices #####
def log_all(x, tangent, method):
    # input: psd matrices [N,p,p]
    # output: tangent vectors [N,p,p]
    if method == "BW":
        output = np.array([BW_log(tangent, x[i]) for i in range(x.shape[0])])
    else:
        output = np.array([AI_log(tangent, x[i]) for i in range(x.shape[0])])
    return output

##### Vectorization #####
def vec_all(x):
    # input: symmetric matrices [N,p,p]
    # output: long vectors [N,p^2]
    output = np.array([x[i].flatten() for i in range(x.shape[0])])
    return output

##### Half-Vectorization #####
def vech_all(x):
    # input: symmetric matrices [N,p,p]
    # output: long vectors [N,p*(p+1)/2]
    upper_idx = np.triu_indices(x.shape[1])
    if x.shape[1] == x.shape[2]:
        output = np.array([x[i][upper_idx] for i in range(x.shape[0])])
    else:
        warnings.warn("Input is not in the form of square matrices.")
        output = None
    return output

##### Dissimilarity (Binary)#####
def dissimilarity_binary(label):
    # input: label [N]
    # output: dissimilarity matrix [N,N]
    output = np.zeros((len(label), len(label)), dtype = np.float32)
    
    for i in range(len(label)):
        for j in range(len(label)):
            if label[i] == label[j]:
                output[i,j] = 0
            else:
                output[i,j] = 1
    
    return output

#####  Dissimilarity (Decay Factor) #####
def dissimilarity_decay(x, label, t, method):
    # input: psd matrices [N,p,p] and label [N]
    # output: dissimilarity matrix [N,N]
    output = np.zeros((x.shape[0], x.shape[0]), dtype = np.float32)
    
    for i in range(x.shape[0]):
        for j in range(x.shape[0]):
            if method == "BW": # distance between two psd matrices as weight
                output[i,j] = BW_dist(x[i], x[j])
            else:
                output[i,j] = AI_dist(x[i], x[j])
            if label[i] == label[j]: # whether these two matrices are in the same class as a scalor
                output[i,j] = t*output[i,j]
    
    return output

##### Frobenius Distance #####
def F_dist(A,B):
    output = np.sqrt(np.trace(np.dot((A-B),np.transpose(A-B))))
    return output

#####