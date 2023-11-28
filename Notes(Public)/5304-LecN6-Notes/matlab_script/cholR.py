# Row-wise Cholesky factorization

import numpy as np
from numpy import linalg as la
import copy

def cholR (A):
#---------------------------------------------
# function [U] = cholR (A) 
# performs row-Cholesky factorization of A 
#---------------------------------------------
    n = np.size(A,0)
    U = copy.deepcopy(A) 
    for k in range(n):
        if (U[k,k]<= 0):
            raise TypeError('->nonpositive pivot!') 
        U[k,k:n] /= np.sqrt(U[k,k]) 
        for i in range(k+1,n):
            U[i,i:n] -= U[k,i]*U[k,i:n]
    U = np.triu(U)
    return(U)