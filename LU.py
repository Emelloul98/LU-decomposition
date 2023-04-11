import numpy as np


def lu_solving(a_matrix, b_vector):
    rows, cols = np.shape(a_matrix)
    u_matrix = np.zeros((rows, cols))
    l_matrix = np.zeros((rows, cols))
    i = 1
    matrix_stat = True
    y_vector = np.zeros(len(a_matrix))
    x_vector = np.zeros(len(a_matrix))
    # the first line of u is the same as A:
    for k in range(rows):
        u_matrix[0][k] = a_matrix[0][k]
        l_matrix[k][k] = 1
    # the first col of L:
    for j in range(1, cols):
        l_matrix[j][0] = a_matrix[j][0]/u_matrix[0][0]
   # going each time on one line of each matrix
    while i < rows:
        sigma_num = 0
        # if the status is true we will go to L and if false to U
        if matrix_stat:
            for j in range(1, i):
                for k in range(0, j):
                    sigma_num += l_matrix[i][k]*u_matrix[k][j]
                l_matrix[i][j] = (1/u_matrix[j][j])*(a_matrix[i][j]-sigma_num)
                sigma_num = 0
            matrix_stat = False
        else:
            for j in range(i, rows):
                for k in range(0, i):
                    sigma_num += l_matrix[i][k] * u_matrix[k][j]
                u_matrix[i][j] = a_matrix[i][j] - sigma_num
                sigma_num = 0
            matrix_stat = True
            i = i+1
    # the linear solving  from the beginning to end for Ly=b
    for k in range(0, rows):
        array_sum = 0
        for j in range(0, rows):
            if j == k:
                continue
            array_sum += l_matrix[k][j]*y_vector[j]
        y_vector[k] = (b_vector[k]-array_sum) / l_matrix[k][k]
    # the linear solving in reverse for Ux=y
    for k in range(rows-1, -1, -1):
        array_sum = 0
        for j in range(0, rows):
            if j == k:
                continue
            array_sum += u_matrix[k][j]*x_vector[j]
        x_vector[k] = (y_vector[k]-array_sum) / u_matrix[k][k]

    return x_vector

# made by emanuel malloul

