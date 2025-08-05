# This code is associated to the quantum optimization benchmarking effort
#
# (C) Copyright IBM 2025.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

import numpy as np
from docplex.mp.model import Model


def birkhoff_constraints(n):
    "Function that returns a matrix with the Birkhoff constraints, i.e., sum columns and rows must be equal to 1. " 
    M = np.zeros((n*n, 2*n), dtype=int)

    # Fill the first part of matrix with ones
    for i in range(n):
        M[:,i][(i*n) : (i*n) + n] = np.ones((1,n))

    # Fill the second part of matrix with ones
    for k in range(0,n):
        for i in range(0,n):
            for j in range(0,n):
                if i == j:
                    M[k*n + i, n + j] = 1
                    
    return M

########## DOCPLEX ##########

class BirkCplex():

    def __init__(self, x_star) -> None:
        self.x_star = x_star
        self.n = int(np.sqrt(self.x_star.size))
        self.s = np.sum(x_star) / self.n
        self.time_limit = 3600

    def set_time_limit(self, time_limit):
        self.time_limit = time_limit
    
    def solve(self):       
        for k in range(1,(self.n-1)**2 + 2):
            status, p, w = self.solve_permutations_and_weights_fixed_size(k)
            if status == 'success':
                return status, p, w

    def solve_permutations_and_weights_fixed_size(self, k):
        
        """
        Function that computes a decomopsition of a doubly stochastic 
        matrix 'x_star' as the convex combinatino of 'k' permutation matrices. 
        """

        model = Model()
        p = model.binary_var_list(self.n*self.n*k, name='permutations')
        w = model.continuous_var_list(k, name='weights')
        w_aux = model.continuous_var_list(self.n*self.n*k, name='aux_weights')
        A = birkhoff_constraints(self.n)
        M = 100000

        # Constraints
        # sum rows and columns must be equal to 1
        for index in range(0,k):
            for j in range(0,2*self.n):
                model.add_constraint(model.sum(A[i,j]*p[self.n*self.n*index + i] for i in range(0,self.n*self.n)) == 1)
        # weights must sum to 1
        model.add_constraint(model.sum(w[index] for index in range(0,k)) == self.s)
        # weights need to be larger than 0
        for i in range(0,k):
            model.add_constraint(w[i] >= 0)
        # decomposition must be exact
        for i in range(0,self.n*self.n*k):
            model.add_constraint(w_aux[i] <= M*p[i])
        for i in range(0,k):
            for j in range(self.n*self.n):
                model.add_constraint(w_aux[self.n*self.n*i + j] <= w[i])
        for i in range(0,k):
            for j in range(self.n*self.n):
                model.add_constraint(w_aux[i*k + j] >= w[i] - M*(1-p[i*k + j]))
        # The last two lines change the equality to two inequalities and add slackness
        for i in range(0,self.n*self.n):
            model.add_constraint(model.sum(w_aux[i + j*self.n*self.n] for j in range(0,k)) == self.x_star[i])

        # Solve
        model.set_time_limit(self.time_limit)
        solution = model.solve()
        if (model.solve_details.status_code not in [101]):
                status = 'failed'
                return status, None, None
        
        ## Retrieve solutions
        status = 'success'
        perm = solution.get_value_list(p)
        c = solution.get_value_list(w)

        return status, perm, c