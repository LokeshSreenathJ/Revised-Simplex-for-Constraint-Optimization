# Revised-Simplex-for-Constraint-Optimization
Built a linear Constraints Optimization Algorithm from Scratch that works on any type of optimization problem

Built a Revised-Simplex with Two-Phase method.

PHASE-ONE
Let us follow the outline of the revised simplex method given on p.97 of Bertsimas and Tsitsiklis, "Introduction to Linear Optimization." However, we deviate from this approach in the following ways: 
Step1) The Phase I problem can never be unbounded, so we do not check for unboundedness after finding an improving direction.
Step2) Once an improving direction is found (i.e. once a nonbasic variable with a negative reduced cost is found (and decided upon) to enter the basis we force an artificial variable to leave the basis and immediately delete this artificial variable from the problem so that it cannot enter the basis in a subsequent iteration.

Inputs for Phase-One:
c = cost vector (n x I) 
A = constraint matrix (m x n). includes artificial variables for Phase I 
b = right-hand-side vector (m x I)
B = basis matrix (m x m)
B_inv = inverse of basis matrix B (m x m) 
bvar = row vector of basic variable indices (1 x m) 
nbvar = row vector of nonbasic variable indices (1 x (n-m)) 
enteringVar = the index in the set of the nonbasic variables of the nonbasic variables that will enter the basis. For example, if bvar = [1,2,3,4] (i.e ., x1,..,x4 are in the basis), and nbvar= [5,6,7,8] (i.e., x5, ... ,x8 are nonbasic variables), then enteringVar= 2 means that x5, the second index in the set of nonbasic variables will enter the basis. 

Output for Phase-One:
soln_status = -1 if infeasible, 0 otherwise 
A= possibly modified (m x n) constraint matrix; redundant rows may have been eliminated 
b = possibly modified (m x 1) right-hand-side vector; 
B_inv = inverse of basis matrix B (m x m) 
bvar = row vector of basic variable indices (1 x m) 
nbvar = row vector of nonbasic variable indices (1 x (n-m)) 

If any artificial varaiable present in the basis and its value is 0, then Remove all artificial variables from the set of nonbasic variables& Eliminate the corresponding columns in the A matrix as well.

If any artificial variable is in the basis but its value is not equals to 0, then try to pivot this artificial variable with any non-artificial variable that is in non-basis. If non-zero artificial variable doesnot goes out from basis using above two steps then the given solution is infeasible.

"We never get a unbounded condition in the Phase-one step because the cost function doesnot contain the non-artificial variables, so having exteme direction does not affect the main cost function".

If there are no artificial variables in the final Phase 1 basis, all artificial variables have been eliminated from the problem, so we may return the current basis to begin Phase 2 of the simplex method.

PHASE-TWO Method:
c = cost vector (n x 1) 
A= constraint matrix (m x n), includes artificial variables for Phase 1 
b = right-hand-side vector (m x 1)
c_B = cost vector of basic variables (m x 1) 
B = basis matrix (m x m) 
B_inv = inverse of basis matrix B (m x m) 
bvar = row vector of basic variable indices (1 x m) 
nbvar = row vector of nonbasic variable indices (1 x (n-m)) 

Output for the PHASE-TWO Method:
soln_status = 1 if unbounded, 2 if finite optimal solution 
z = objective function value: -inf if unbounded, finite otherwise 
x = (n x 1) optimal solution if soln_status = 2, a direction "u" of unboundedness otherwise 
p = an (m x 1) dual vector 
B_inv = inverse of basis matrix B (m x m) 
bvar = row vector of basic variable indices ( I x m) 
nbvar = row vector of nonbasic variable indices ( 1 x (n-m)) 

We return the inverse basis matrix, the index set of basic variables, and the index set of nonbasic variables so that we can perform a warm-start in the future. 

If the reduced costs of all nonbasic variables are nonnegative, the current basic feasible solution x_B is %optimal, and the algorithm terminates.

