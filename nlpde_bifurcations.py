### For the bifurcation analysis
### One dimensional advection-diffusion equation on bounded domain
### order of the parameters: 
import numpy as np
from scipy.linalg import eig
from scipy.optimize import fsolve, approx_fprime
from scipy.integrate import quad
from scipy.interpolate import BarycentricInterpolator

import nlpde_utils

class SteadyState():
    def __init__(self,rhs,X,pars,x):
        """ rhs is the RHS of the PDE. X is the vector of variables
        rhs should take argument (X, parameters)
        The x argument is the spatial points at which the profile is evaluated. 
        It is mainly for bookkeeping/plotting later.
        """
        self.X = X
        self.pars = np.array(pars)
        self.rhs = rhs
        self.mu = None # determined using the apply_measure function
        self.x = x
        
        # stability properties determined in 'get_stability'
        self.stable = None
        self.eigenvalues = None
        self.eigenvectors = None
        
    def correct(self):
        """ correct the u values for fixed parameters. Use current u values as guess."""
        self.X = fsolve(lambda x: self.rhs(x, self.pars), x0=self.X)
        
    def compute_stability(self,f_matrices, keep_n = 4, ignore_larger = 1e3, thr_zero = 1e-8):
        """ compute the stability of this point at the current parameter values.  
        The argument, f_matrices, is a function of X, pars and it returns two matrices A, B.
        The current method solves the generalized eigenvalue problem A v = omega B v.
        We then look at the real part of the eigenvalues. To avoid dealing with spurious eigenvalues there is the
        ignore_larger parameter. The argument thr_zero is used to determine when an eigenvalue is larger than zero. """
        A, B = f_matrices(self.X, self.pars)
        evals, evecs = eig(A, B)
        idx_toolarge = abs(evals)>ignore_larger
        evals = evals[~idx_toolarge]
        evecs = evecs[:, ~idx_toolarge]
        # sort
        idx_sorted = np.argsort(evals.real)[::-1] # from large to small
        evals=evals[idx_sorted]
        evecs = evecs[:,idx_sorted]
        self.eigenvalues = evals[:keep_n]
        self.eigenvectors = evecs[:, :keep_n]
        self.stable = evals[0]<thr_zero
        
    def compute_measure(self, measure_function):
        """ To obtain a scalar representation of this steady state. The measure function takes in
        two arguments: u and pars. """
        self.mu = measure_function(self.X, self.pars)
        
class BifurcationPoint():
    """ Similar to steady state but it satisfies extra equations (based on Seydel's book 2010). 
    We are not interested in stability here (the bifurcation point is a point where a switch of stability happens)."""
    def __init__(self, rhs, X, pars, x):
        """X contains u, vector of N+1 function values, and c, the dummy parameter for mass conservation. 
        Then h, a N+2-vector."""
        
        self.X = X
        self.pars = np.array(pars)
        self.rhs = rhs
        self.mu = None # determined using the apply_measure function
        self.x = x # only for bookkeeping
    
    def correct(self, par_ind):
        """ correct the u,c and h values + one parameter given by par_ind. Use current values as initial guess."""
        f_to_solve = lambda Y: self.rhs(Y[:-1], np.hstack((self.pars[:par_ind], Y[-1], self.pars[par_ind+1:])))
        Y = fsolve(f_to_solve, x0=np.hstack((self.X, self.pars[par_ind])))
        self.X = Y[:-1]
        self.pars[par_ind] = Y[-1]
        
    def compute_measure(self, measure_function):
        """ measure function takes in X and pars"""
        self.mu = measure_function(self.X, self.pars)
        
class Branch1Par():
    """ Basically a container for points that vary in one parameter.
    Contains the methods to perform the pseudo-arclength continuation. """
    def __init__(self, name, par_index, starting_point, starting_tv):
        """ argument name is just to identify the branch. par_index is which parameter will be varied along the branch. 
        Assume that the starting point is already corrected. """
        self.name = name
        self.rhs = starting_point.rhs
        self.par_index = par_index
        
        self.points = [starting_point]
        self.tangent_vectors = [starting_tv]
        
        # numerical parameters, can set them directly later
        self.h = 0.1 # step size to start with. 
        self.h_min = 1e-3
        self.h_max = 1
        
        self.par_min = 0
        self.par_max = 100
        
       
        # define adapted rhs, appending parameter to X
        def rhs1(Y):
            pp = starting_point.pars.copy()
            pp[self.par_index]=Y[-1]
            return self.rhs(Y[:-1], pp)
        self.rhs_1 = rhs1
        
    def correct(self, Y0, v0):
        """do the correction of the Y0 point with also the condition of orthogonality to vector v0.
        Also returns the number of iterations the solver needed. """
        def tosolve(y):
            return np.hstack((self.rhs_1(y), np.dot(v0, y-Y0)))
        Y1, infodict, _, _ =fsolve(tosolve, Y0, full_output=True, xtol=1e-14)
        return Y1, infodict['nfev']

    def continue_branch(self, Nsteps):
        """ use stepsize etc. from the numerical parameters 
        We use the notation Y for the vector containing X and the parameter. """
        for i in range(Nsteps):
            v_last = self.tangent_vectors[-1]
            P_last = self.points[-1]
            Y_last = np.hstack((P_last.X, P_last.pars[self.par_index]))
            Y_predicted = Y_last + v_last*self.h
            # correction
            Y_corrected, n_it = self.correct(Y_predicted, v_last)
            # turn this into a new point
            new_pars = P_last.pars.copy()
            new_pars[self.par_index] = Y_corrected[-1]
            new_P = SteadyState(self.rhs, Y_corrected[:-1], new_pars, x=P_last.x)
            
            ## compute current tangent vector with normalization based on previous tangent vector
            J = approx_fprime(Y_corrected, self.rhs_1)
            M = np.zeros((J.shape[0]+1, J.shape[1]))
            M[:-1,:] = J
            M[-1,:] = v_last
            b = np.zeros(M.shape[0])
            b[-1] = 1
            v_new = np.linalg.solve(M, b)
            v_new /= np.linalg.norm(v_new)
            ## append
            self.points.append(new_P)
            self.tangent_vectors.append(v_new)
            ## check conditions for the continuation
            if new_pars[self.par_index] < self.par_min or new_pars[self.par_index]>self.par_max:
                break
            ## adapt step size if needed
            if n_it<=4 and self.h < self.h_max:
                self.h *= 1.1
            if n_it > 15 and self.h > self.h_min:
                self.h /= 1.1
        
    def compute_stability(self, f_matrices, keep_n = 4, ignore_larger = 1e3, thr_zero = 1e-8):
        for st in self.points:
            st.compute_stability(f_matrices, keep_n, ignore_larger, thr_zero)
    def compute_measure(self, measure_function):
        """ applies the measure function to each point on the branch. """
        for st in self.points:
            st.compute_measure(measure_function)
    def get_point_measures(self):
        """ returns the point measures as an array. """
        muv = np.array([st.mu for st in self.points])
        return muv        
    def get_point_pars(self):
        """ returns the values of the parameter that is varied on this branch as an array """
        pv = np.array([st.pars[self.par_index] for st in self.points])
        return pv
    def get_point_stability(self):
        """returns a boolean array showing whether the points in the point list are stable or not."""
        stab_array = np.array([st.stable for st in self.points])
        return stab_array
    def get_tangents(self):
        """returns an array of arrays with the tangent vectors"""
        return np.array(self.tangent_vectors)
        
class Branch2Par():
    """ For 2-parameter bifurcation diagrams.
    Contains the methods to perform the pseudo-arclength continuation. """
    def __init__(self, name, par_index, starting_point, starting_tv):
        """ name is just to identify the branch.
        par_index is which parameters will be varied along the branch. 
        Assume that the starting point is already corrected. 
        The starting point should be of class BifurcationPoint. 
        Use method described by Seydel (2010)"""
        self.name = name
        self.rhs = starting_point.rhs
        self.par_index = par_index # this is a list with two numbers
        
        self.points = [starting_point]
        self.tangent_vectors = [starting_tv]
        
        # numerical parameters, can be adapted
        self.h = 0.1 # step size to start with. 
        self.h_min = 1e-3
        self.h_max = 1
        
        self.par_min = np.array([0,0])
        self.par_max = np.array([100,100])
         
        # define adapted rhs, where the two free parameters are appended to the X vector
        def rhs1(Y):
            pp = starting_point.pars.copy()
            pp[self.par_index[0]]=Y[-2]
            pp[self.par_index[1]]=Y[-1]
            return self.rhs(Y[:-2], pp)
        self.rhs_1 = rhs1
        
    def correct(self, Y0, v0):
        """do the correction of the X0 point with also the condition of orthogonality to vector v0.
        Also returns the number of iterations the solver needed. """
        def tosolve(y):
            return np.hstack((self.rhs_1(y), np.dot(v0, y-Y0)))
        Y1, infodict, _, _ =fsolve(tosolve, Y0, full_output=True, xtol=1e-14)
        return Y1, infodict['nfev']

    def continue_branch(self, Nsteps, force_tangent=None):
        """ use stepsize etc. from the numerical parameters 
        We use the notation X for the vector containing u (the values of the steady state) + the parameter that is varied"""
        for i in range(Nsteps):
            v_last = self.tangent_vectors[-1]
            P_last = self.points[-1] # object
            Y_last = np.hstack((P_last.X, P_last.pars[self.par_index])) # vector with parameters appended
            Y_predicted = Y_last + v_last*self.h
            # correction
            Y_corrected, n_it = self.correct(Y_predicted, v_last)
            # turn this into a new point
            new_pars = P_last.pars.copy()
            new_pars[self.par_index] = Y_corrected[-2:]
            new_P = BifurcationPoint(self.rhs, Y_corrected[:-2], new_pars, x=P_last.x)
            
            ## compute current tangent vector with normalization based on previous tangent vector
            J = approx_fprime(Y_corrected, self.rhs_1)
            M = np.zeros((J.shape[0]+1, J.shape[1]))
            M[:-1,:] = J
            M[-1,:] = v_last
            b = np.zeros(M.shape[0])
            b[-1] = 1
            v_new = np.linalg.solve(M, b)

            if force_tangent is not None: #### if we want to enforce the tangent vector going in a certain direction
                v_new[:-2]=0
            v_new /= np.linalg.norm(v_new)
            
            ## append
            self.points.append(new_P)
            self.tangent_vectors.append(v_new)
            ## check conditions for the continuation
            if np.any(new_pars[self.par_index] < self.par_min) or np.any(new_pars[self.par_index]>self.par_max):
                break
            ## adapt step size if needed
            if n_it<=4 and self.h < self.h_max:
                self.h *= 1.1
            if n_it > 15 and self.h > self.h_min:
                self.h /= 1.1
        
    def compute_measure(self, measure_function):
        """ applies the measure function to each point on the branch. """
        for st in self.points:
            st.compute_measure(measure_function)
    def get_point_measures(self):
        """ returns the point measures as an array. """
        muv = np.array([st.mu for st in self.points])
        return muv        
    def get_point_pars(self):
        """ returns the values of the parameters that are varied on this branch as two arrays """
        pv1 = np.array([st.pars[self.par_index[0]] for st in self.points])
        pv2 = np.array([st.pars[self.par_index[1]] for st in self.points])
        return pv1, pv2
    def get_tangents(self):
        """returns an array of arrays with the tangent vectors"""
        return np.array(self.tangent_vectors)
        
