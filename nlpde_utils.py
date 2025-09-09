#### Utility functions for the nonlocal advection-diffusion equations, bifurcations, eigenvalues etc.

import numpy as np
from scipy.linalg import eig
from scipy.optimize import fsolve, approx_fprime
from scipy.integrate import quad

def cheb(N):
    """ From Trefethen's book. Returns Chebyshev differentiation matrix and collocation points.
    This is on [-1, 1], so may need rescaling if you are using another interval. """
    x = np.cos(np.pi*np.linspace(0, N, N+1)/N)
    D = np.zeros((N+1, N+1))
    for i in range(N+1):
        for j in range(N+1):
            if i==j:
                if i==0:
                    D[i,j] = (2*N**2+1)/6
                elif i==N:
                    D[i,j] = -(2*N**2+1)/6
                else:
                    D[i,j] = -x[j]/2/(1-x[j]**2)
            else:
                ci = 1
                if i==0 or i==N:
                    ci=2
                cj = 1
                if j==0 or j==N:
                    cj=2
                D[i,j] = ci/cj * (-1)**(i+j)/(x[i]-x[j])

    return D, x


def clencurt(N):
    """ returns a vector with weights w and the interpolation points x, to determine the integral of a function.
    this is from Trefethen's book.
    So the dot product between w and f(x) gives a spectral approximation to the integral of f over the interval [-1, 1]
    Note: have to rescale for other intervals. """

    theta = np.pi*np.linspace(0, N, N+1)/N
    x = np.cos(theta)
    w = np.zeros(N+1)
    ii = np.arange(1, N)
    v = np.ones(N-1)
    if N%2==0:
        w[0]=1/(N**2-1)
        w[-1]=w[0]
        for k in range(1, N//2):
            v = v-2*np.cos(2*k*theta[ii])/(4*k**2-1)
        v = v - np.cos(N*theta[ii])/(N**2-1)
    else:
        w[0] = 1/N**2
        w[-1]=w[0]
        for k in range(1, (N-1)//2+1):
            v = v-2*np.cos(2*k*theta[ii])/(4*k**2-1)
    w[ii]=2*v/N
    return w,x

def pj(j,x):
    """ returns polynomial interpolating a function that is zero at all values of the vector x, but one at x[j]
    We use this to set up the matrix of the nonlocal operator."""
    pp= lambda xxx: np.prod([xxx-x[k] for k in range(len(x)) if k != j], )/np.prod([x[j]-x[k] for k in range(len(x)) if k!=j])
    return np.vectorize(pp)
    
    
def get_A_cheb(x, g, L, method='quad'):
    """ returns the matrix corresponding to the discretiation of the nonlocal operator.
    x is the vector of interpolation points (comes from the function cheb, for example). g is a function that corresponds to the interaction kernel (on positive values). x has length N+1
    The entry A_{ij} is given by integral from 0 to L of g(y-x[i])*pj(y), where pj is the polynomial interpolating a delta function at xj. 
    If method is quad, we compute these integrals numerically using scipy's quad. It can be quite slow for large N.
    If method is clencurt, the entry A_{ij} is approximated as g(x[j]-x[i])*wj, what it would be if the integral was computed using clenshaw curtis quadrature. In this case we can also directly work with matrices and we don't loop through every entry. 
    """
    N = len(x)-1

    if method=='quad':
        A = np.zeros((N+1,N+1))
        for i in range(N+1):
            for j in range(N+1):
                # set A_{ij} here
                A[i,j] = quad(lambda t: -g(x[i]-t)*pj(j,x)(t), 0, x[i])[0] + quad(lambda t: g(t-x[i])*pj(j,x)(t), x[i],L)[0]
    else:
        ### using matrices, hopefully faster
        XX = np.tile(x.reshape((1, -1)), (len(x),1)) ### each row is the X vector
        YY = XX - XX.T
        A = np.sign(YY) * g(abs(YY))
        # don't forget to multiply by the weights of the clenshaw curtis if needed
        w,_ = clencurt(N)
        w *=L/2
        A = A @ np.diag(w)
    return A