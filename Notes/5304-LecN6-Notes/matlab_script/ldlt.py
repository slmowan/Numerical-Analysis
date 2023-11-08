# % This is LDLT factorization of A
# % Adapted from Gaussian Elimination. Main observation: The working matrix
# % A(k+1:n, k+1:n) in standard LU remains symmetric.
# % notice the difference with gauss.m

import numpy as np
from numpy import linalg as la
import copy
def ldlt (A):
#---------------------------------------------
#  performs LDLT factorization of a matrix  A 
#---------------------------------------------
    n = np.size(A,0)
    L = copy.deepcopy(A)
    for k in range(n-1): 
        if (L[k,k]==0):
            raise TypeError('->zero diagonal - stop') 
        for i in range(k+1,n): 
            piv =  L[k,i]/L[k,k]
            L[i,i:n]=L[i,i:n]-piv*L[k,i:n]
## end loop 
    d = np.diag(L)
    L  = np.tril( (1.0 / d) * L.T) 
    return  L, d