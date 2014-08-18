from numpy import array, dot, arange
import matplotlib.pyplot as plt

'''
result = gaussElimination(A,b)
solves matrix equation [A]{x}={b}
'''
def gaussElimination(A, b):
    n = len(b)
    # Elimination phase
    for k in range(0, n - 1):
        for i in range(k + 1, n):
            if A[i, k] != 0.:
                lam = A[i, k] / A[k, k]
                A[i, k:n] = A[i, k:n] - lam * A[k, k:n]
                b[i] = b[i] - lam * b[k]
    # Back substitution
    for k in range(n - 1, -1, -1):
        b[k] = (b[k] - dot(A[k, k + 1:n], b[k + 1:n])) / A[k, k]
    return b


A = array([[300., -100., 0., 0., 0.],
          [-100., 200., -100., 0., 0.],
          [0., -100., 200., -100., 0.],
          [0., 0., -100., 200., -100.],
          [0., 0., 0., -100., 300.]])

b = array([200. * 100., 0., 0., 0., 200. * 500.])

res_num = gaussElimination(A, b)

# Plotting results
x_num = arange(0.05, 0.55, 0.1)
x = arange(0., 0.6, 0.1)
plt.plot(x, 800. * x + 100, 'k-', x_num, res_num, 'ks',)
plt.legend( ('Exact', 'Numerical'), loc='upper left')
plt.xlabel('$x$  in  $m$', fontsize=16)
plt.ylabel('$T$  in  $^\circ C$', fontsize=16)
plt.axis([0, 0.5, 0, 600])
plt.title('Source free heat conduction in an insulated rod', fontsize=16)
plt.show()