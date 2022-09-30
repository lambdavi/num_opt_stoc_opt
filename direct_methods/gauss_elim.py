import numpy as np

# Defining the Gaussian Elimination Method

def gaussianElim(a,b):
    # Defining the number of equations
    n = len(b) 

    # Elimination Phase
    for k in range(0, n-1):
        for i in range(k+1, n):
            temp = a[i,k]/a[k,k]
            a[i,k+1:n]=a[i,k+1:n] - temp*a[k, k+1:n]
            b[i] = b[i] - temp*b[k]
    
    # Back substitution phase
    for k in range(n-1, -1, -1):
        b[k] = (b[k] - np.dot(a[k, k+1:n], b[k+1:n]))/a[k,k]
    return b # we are reusing the vector b because it has the same dimentions of x


# Let's play with it
A = np.array([[1.0, 1.0, -1.0], [2.0, -1.0, 1.0], [-1.0, 2.0, 2.0]])
b = np.array([-2.0, 5.0, 1.0])
print(f'A: {A}')
print(f'b: {b}')
x = gaussianElim(A, b)
print(f'Solution of the linear system: {x}')
