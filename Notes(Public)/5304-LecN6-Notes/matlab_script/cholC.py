# column-wise Cholesky factorization

import numpy as np
from numpy import linalg as la
import copy

def ColChol (A): 
# function [G] = ColChol (A) 
# Column Cholesky factorization of A
# returns a low. triang. matrix.
#
 n = np.size(A,0) 
 G = copy.deepcopy(A)
 for j in range(n):
##-------------------- perform updates
     for k in range(j):
         G[j:,j] -= G[j:,k]*G[j,k] 
##-------------------- get diag. entry
     G[j,j] = np.sqrt(G[j,j])
##-------------------- scale to get column
     G[j+1:,j] /= G[j,j]
 return(np.tril(G))