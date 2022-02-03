#!/usr/bin/env python3

import numpy as np
import math

# Function that performs Strassen matrix multiplication for big matrices that cannot be hold in memory
def strassen(matrix1, matrix2):

    # Slice matrices in 4 submatrices each
    matrices = [matrix1, matrix2]
    submatrices = []
    for m in matrices:
        rows = m.shape[0]
        cols = m.shape[1]
        m1 = m[0:math.ceil(rows/2), 0:math.ceil(cols/2)]
        m2 = m[0:math.ceil(rows/2), math.ceil(cols/2):]
        m3 = m[math.ceil(rows/2):, 0:math.ceil(cols/2)]
        m4 = m[math.ceil(rows/2):, math.ceil(cols/2):]
        submatrices.append([m1, m2, m3, m4])
    
    # Compute dot product between submatrices 
    d1 = np.add(submatrices[0][0].dot(submatrices[1][0]), submatrices[0][1].dot(submatrices[1][2]))
    d2 = np.add(submatrices[0][0].dot(submatrices[1][1]), submatrices[0][1].dot(submatrices[1][3]))
    d3 = np.add(submatrices[0][2].dot(submatrices[1][0]), submatrices[0][3].dot(submatrices[1][2]))
    d4 = np.add(submatrices[0][2].dot(submatrices[1][1]), submatrices[0][3].dot(submatrices[1][3]))

    # Concatenate matrices along columns and rows 
    d12 = np.concatenate((d1, d2), axis=1)
    d34 = np.concatenate((d3, d4), axis=1)
    strassen_multiplied = np.concatenate((d12, d34), axis=0)

    # Return result of multiplication 
    return strassen_multiplied

