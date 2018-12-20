import power_diagram as pdm
import scipy.sparse.linalg
import scipy.sparse
import numpy as np
import matplotlib.pyplot as plt
import math

def partial_transport_solve(pd, prescribed_area, max_iter=10):
    iad = pd.der_measures()
    i = 0
    while True:
        M = scipy.sparse.csr_matrix( ( iad.der_data, iad.der_indices, iad.der_indptr ),
                                     [ pd.nb_diracs(), pd.nb_diracs() ] )
        V = prescribed_area - iad.measures

        if max( abs( V ) ) < prescribed_area * 1e-2:
            break
        assert(np.min(iad.measures) > 0)

        R = scipy.sparse.linalg.spsolve( M, V )
        alpha = 1
        while 1:
            pd.add_to_weights(alpha * R)
            iad = pd.der_measures()
            if np.min(iad.measures) > 0:
                break
            else: #try with a lower alpha
                pd.add_to_weights(-alpha*R)
                alpha /= 2
        if i==max_iter:
            print("Warning: Maximum number of iterations reached!")
            return
        i += 1

def proj_noncongested(X, mass):
    pd = pdm.PowerDiagram2d()
    pd.add_box_shape([-10, -10, 10, 10])
    pd.set_ball_cut( True )
    r = np.sqrt(mass/np.pi)
    for p in X:
        pd.add_dirac( [p[0], p[1]], r*r)
    partial_transport_solve(pd, mass)
    return pd

n = 20
t = np.linspace(-2,2,n)
x,y = np.meshgrid(t,t)
X = np.hstack((np.reshape(x,(n*n,1)),
               np.reshape(y,(n*n,1))))
R2 = X[:,0]**2 + X[:,1]**2
X = X[R2 <= 4]
N = X.shape[0]
#plt.scatter(X[:,0], X[:,1])
#plt.show()
rho0 = 0.5
mass = rho0 * math.pi * 4/N
eps = np.sqrt(mass) # = h
dt = 0.5*eps

for i in range( 100 ):
    print("ITERATION", i, " -----------------------------------------------------")
    pd = proj_noncongested(X, mass)
    #pd.vtk_output( "vtk2/nv-n=%d_%03d.vtk" % (N, i ) )
    B = pd.centroids()
    n = np.sqrt(X[:,0]**2 + X[:,1]**2+0.1)
    NX = np.hstack((np.reshape(X[:,0]/n,(N,1)),
                    np.reshape(X[:,1]/n,(N,1))))
    X += dt*((B-X)/eps - NX)
    



