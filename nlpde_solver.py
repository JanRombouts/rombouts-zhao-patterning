import numpy as np
from scipy.integrate import quad

from scipy.fft import fft, ifft, fft2, ifft2

from scipy.linalg import circulant
from scipy.integrate import dblquad, solve_ivp, quad
from scipy.optimize import fsolve

########
#### Module containing the classes and methods to solve the nonlocal PDEs.
#### Also includes functions that are used in steady state solving, eigenvalue determination, etc. 
########


##### Some general functions used in different classes

def Phi1(x):
    # piecewise constant discretization
    return (x<=1) & (x>0)
def Phi2(x):
    return np.maximum(0, 1-abs(x-1/2))
    
    
def get_kernelweights_1d(G, L, N, R, phiN=1):
    """ creates the weight vector for one dimension.
    G is a function (the vector version of the kernel - in 1d this means sign change)
    R is the cutoff distance we use for the kernel
    L is the domain size
    N is domain discretization
    phiN is 1 or 2, 1 for pcw constant, 2 for pcw linear (see Gerisch 2010)
    """
    h = L/N
    M = int(R/h)+1 # to determine which coefficients will be nonzero
    if phiN==1:
        lmin = M-1
        lmax = M
        Phi = Phi1
    else:
        lmin=M
        lmax=M+1
        Phi=Phi2
        
    w = np.zeros(lmin + lmax + 1)
    for i in range(-lmin, lmax+1):
        if phiN==1: # integration domain...
            w[i+lmin] = quad(lambda r: G(r)*Phi(r/h+i), -h*i, (1-i)*h)[0] # use scipy's quadrature here
        elif phiN==2:
            w[i+lmin] = quad(lambda r: G(r)*Phi(r/h+i), h*(-i-1/2), (3/2-i)*h)[0]
            
    Vcol = np.zeros(N)
    Vcol[:lmax] = w[-lmax:]
    Vcol[-lmin-1:] = w[:lmin+1]

    return w, Vcol, [lmin, lmax]

#########################################################
##################### Solver classes
#########################################################

class Solver1DPeriodic(object):
    """ One population, periodic boundary conditions, linear packing function with carrying capacity K and constant diffusion constant
    Ie solve sthe equation u_t = -d_x ( -D u_x + chi0 p(u) u A(u))
    where A(u) is the integral for the cell-cell interactions. 
    """
    def __init__(self, **params):
        self.D = params['D']
        self.K = params['K']
        self.g = params['g'] # the interaction function (scaled), vectorized
        self.chi0 = params['chi0']
        self.sigma = params['sigma']
        def G_vec(x):
            return np.sign(x) * self.g(abs(x)/self.sigma)/self.sigma
        self.G_vec = G_vec
        self.p = lambda x: 1- x/self.K
        if 'd' in params.keys(): # pressure component
            self.d=params['d']
        else:
            self.d = 0
        
    def prepsolve(self, L, N, R, phiN=2):
        ## determine the weight vector
        self.w, self.Vcol, self.lminmax = get_kernelweights_1d(self.G_vec, L, N, R, phiN)
        self.L = L
        self.N = N
        self.R = R
        xx = np.linspace(0, self.L, self.N, endpoint=0) # no endpoint (periodicity)
        self.xx = xx
    
    def solve(self, T, dt, u0, save_every=1, method='scipy'):
        """
        Finite volume with computation of the integrals using FFTs (as in Gerisch 2010).
        The weight vector has to be predetermined using the prepsolve function.
        Timestepping is determined by 'method' argument; either standard forward euler, or if method is 'scipy', then use solve_ivp.
        Don't need to save at every time step. 
        """
        
        h = self.L/self.N
        
        Nsteps = int(T/dt)
        Nsave = Nsteps//save_every

        uu = np.zeros((Nsave, self.N)) # for keeping the results
        uu[0,:] = u0
        
        if method=='fwdeuler':  
            u = u0
            for i in range(1,Nsteps):
                # diffusive flux at half-points
                Fdiff = -self.D*(np.roll(u,-1) - u)/h
                # advective flux at half points, compute using the FFT
                uhalf = (np.roll(u, -1) + u)/2
                a = ifft(fft(self.Vcol)*fft(u)).real
                Fadv = a*uhalf*self.p(uhalf)*self.chi0
                ## pressure
                Fpr = -self.d*uhalf*(np.roll(u,-1) - u)/h
                ### Update
                u += dt/h*(-Fdiff + np.roll(Fdiff, 1) - Fadv + np.roll(Fadv, 1) - Fpr + np.roll(Fpr, 1))
                if i%save_every==0:
                    uu[i//save_every,:] = u
        elif method[:5]=='scipy':
            if len(method)>5: # the solver is also given (for example scipyBDF, scipyRK45,...)
                integration_method = method[5:]
            else:
                integration_method='RK45'
            # use solve_ivp
            def dudt(t,u):
                Fdiff = -self.D*(np.roll(u,-1) - u)/h
                # advective flux at half points, compute using the FFT
                uhalf = (np.roll(u, -1) + u)/2
                a = ifft(fft(self.Vcol)*fft(u)).real
                Fadv = a*uhalf*self.p(uhalf)*self.chi0
                
                return 1/h*(-Fdiff + np.roll(Fdiff, 1) - Fadv + np.roll(Fadv, 1))
            
            sol = solve_ivp(dudt, (0, T), t_eval = np.linspace(0, T, Nsave), y0=u0, method=integration_method)
            self.sol = sol # keep for inspection later possibly
            for i in range(sol.y.shape[1]):
                uu[i,:] = sol.y[:,i]
            
        
        self.uu = uu
        self.tt = np.linspace(0, T, uu.shape[0])

               
class Solver1DBoundary(object):
    """ One population, with boundary conditions, linear packing function with carrying capacity K and constant diffusion constant.
    The equation u_t = -d_x ( - (D+epsilon*u) u_x + u p(u) * (A(u) + b)) + f(u)
    where A(u) is the integral for the cell-cell interactions and b is another integral that indicates interaction with the boundary.
    Boundary is represented by a 'field' of boundary material v which is zero inside the domain and 1 outside. 
    There is a separate interaction function Gb for the boundary with its own range and strength. 

    This class allows solving of a more general set of equations than the one used in the paper. The epsilon-term indicates a pressure-like term that drives particles away from high-density regions. There is also a reaction function f which can be used to describe cell death/birth. 
    epsilon and f are zero in the paper. 
    """
    def __init__(self, **params):
        self.D = params['D']
        self.K = params['K']
        self.g = params['g'] # the interaction function (scaled), vectorized
        self.chi0 = params['chi0']
        self.sigma = params['sigma']
        if 'epsilon' not in params.keys(): # for backwards compatibility
            self.epsilon=0
        else:
            self.epsilon = params['epsilon']
        if 'reaction_term' not in params.keys():
            self.reaction_term = lambda u: 0
        else:
            self.reaction_term = params['reaction_term']

        def G_vec(x):
            return np.sign(x) * self.g(abs(x)/self.sigma)/self.sigma
        self.G_vec = G_vec
        self.p = lambda x: 1- x/self.K
        
        # boundary interaction 
        self.chib=params['chib']

        ## boundary asymmetry: boundary on one side is multiplied by 1 + boundary_asymmetry
        if 'boundary_asymmetry' in params.keys():
            self.boundary_asymmetry = params['boundary_asymmetry']
        else:
            self.boundary_asymmetry = 0.

        self.gb = params['gb'] # interaction function with boundary
        self.sigmab = params['sigmab'] # interaction length scale with boundary
        
        def Gb_vec(x):
            return np.sign(x)*self.gb(abs(x)/self.sigmab)/self.sigmab
        self.Gb_vec = Gb_vec
        
    def prepsolve(self, L, N, R, Rb, phiN=2, correct=False):
        """ determine the weight vectors. 
         We have w and wb, where wb is for the boundary interaction. 
         note, L is the actual domain size. We also define an extended domain size (basically L + 2R) 
         for the convolution
         if correct is True, use correct boundaries and centers for the finite volume discretization.
         if False, it is N points including 0 and L but this means that the boundaries are at -h/2 and L+h/2 """
        
        h = L/(N-1)
        if correct:
            h=L/N
        extra_space = np.max([R, Rb])
        Nb = int(extra_space/h)+1 # boundary points added on the left, we will add one less on the right (no endpoint- periodicity)
        # update extra_space to that we get exactly the right discretization
        extra_space = Nb*h
        
        self.w, self.Vcol, self.lminmax = get_kernelweights_1d(self.G_vec, L+2*extra_space, N+Nb+Nb-1, R, phiN)
        
        # also for the boundary
        self.wb, self.Vcolb, self.lminmaxb = get_kernelweights_1d(self.Gb_vec, L+2*extra_space, N+Nb+Nb-1, Rb, phiN)
        
        self.L = L
        self.N = N
        self.Nb = Nb
        self.R = R
        self.Rb = Rb
        
        self.extra_space = extra_space
        self.h=h
        
        xx = np.linspace(0, self.L, self.N, endpoint=1)
        xx_e = np.linspace(-self.extra_space, self.L + self.extra_space, self.N + 2*Nb-1, endpoint=0) # extended domain
        if correct:
            xx = np.linspace(h/2, self.L-h/2, self.N, endpoint=1)
            xx_e = np.linspace(-self.extra_space+h/2, self.L + self.extra_space-h/2, self.N + 2*Nb-1, endpoint=0) # extended domain
            
        # compute the function which will denote the boundary interactions.
        # i.e. the nonlocal part but with the boundary field
        # the boundary field: 1 outside of domain, zero inside
        v = np.ones_like(xx_e)
        v[Nb:-Nb+1] = 0
        self.v=v
        self.a_boundary = ifft(fft(self.Vcolb)*fft(v)).real

        ##### multiply the boundary 'velocity field' by the asymmetry
        ## linear gradient
        self.a_boundary = self.a_boundary * (1+ self.boundary_asymmetry/self.L*xx_e)
        
        self.xx=xx
        self.xxe=xx_e
        
    def solve(self, T, dt, u0, save_every=1,method='scipy', extra_args={}):
        """
        finite volume with computation of the integrals using FFTs (as in Gerisch 2010).
        The weight vector has to be predetermined using the prepsolve function.
        Timestepping is determined by 'method' argument; either standard forward euler, or if method is 'scipy', then use solve_ivp. extra_args can be a dictionary and it is passed as keywords to solve_ivp
        Don't need to save at every time step. 
        
        Special here: boundary conditions need to be enforced.
        """
        Nb = self.Nb # for ease of notation
        h = self.h
        
        Nsteps = int(T/dt)
        Nsave = Nsteps//save_every

        uu = np.zeros((Nsave, self.N)) # for keeping the results
        uu[0,:] = u0
        
        ## we will simulate using the extended domain, but save only the actual domain
        u0e = np.zeros(self.N + 2*self.Nb-1)
        u0e[Nb:-Nb+1] = u0
        self.u0e = u0e

        ### definition of the right-hand side
        def dudt(t,u):
            uhalf = (np.roll(u, -1) + u)/2

            # diffusive flux at half-points
            Fdiff = -(self.D+self.epsilon*uhalf)*(np.roll(u,-1) - u)/h

            # velocity due to cell-cell interactions at half points, compute using the FFT            
            a = ifft(fft(self.Vcol)*fft(u))


            # old behavior, use central differences
            Fadv = uhalf*self.p(uhalf)*(a*self.chi0 + self.a_boundary*self.chib)
            
            ### compensate for the boundary adhesion such that total flux across boundary is zero on both sides
            Fdiff[Nb-1] = -Fadv[Nb-1]
            Fdiff[-Nb] = -Fadv[-Nb]
            
            return 1/h*(-Fdiff + np.roll(Fdiff, 1) - Fadv + np.roll(Fadv, 1)) + self.reaction_term(u)*(1-self.v) #### Reaction term not on the outside! These values should stay zero
        
        if method=='fwdeuler':
            u=u0e.copy()
            dt = self.T/Nsteps
            for i in range(1,Nsteps):
                ### Update
                u += dt *dudt(u,i*dt)

                if i%save_every==0:
                    uu[i//save_every,:] = u[Nb:-Nb+1]
                    
        elif method[:5]=='scipy':
            if len(method)>5: # the solver is also given (for example scipyBDF, scipyRK45,...)
                integration_method = method[5:]
            else:
                integration_method='RK45'
            sol = solve_ivp(dudt, (0, T), t_eval = np.linspace(0, T, Nsave), y0=u0e, method=integration_method, **extra_args)
            self.sol = sol # keep for inspection later possibly
            for i in range(sol.y.shape[1]):
                uu[i,:] = sol.y[Nb:-Nb+1,i]

        self.uu = uu
        self.tt = np.linspace(0, T, uu.shape[0])
        
