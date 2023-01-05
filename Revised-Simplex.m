
% REVISED_PRIMAL_SIMPLEX_METHOD 
A = [1 1 -1 0 0 1 0 0; -1 1 0 -1 0 0 1 0; 0 1 0 0 1 0 0 1];
b = transpose([2 1 3]);
c = [1 -2 0 0 0];
revisedSimplexMethod(c,A,b)
function[soln_status,z,x,p] = revisedSimplexMethod(c,A,b)
%{ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%INPUT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

c = an (n x 1) cost vector 
A= an (m x n) constraint matrix 
b = an (m x 1) right-hand-side vector 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%OUTPUT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
soln_status = -1 if infeasible
				0 if feasible 
				1 if unbounded 
				2 if finite optimal solution

If soln_status = 2, the variables returned are 
	z = the optimal objective function value,
	x = a primal optimal solution,
	p = a dual optimal solution

If soln_status = 1, 
	z = inf (since we are solving a minimization problem) 
	x = a direction of unboundedness (an extreme ray)
	p = an infeasible dual vector
%}

[soln_status, A, b, B_inv, bvar, nbvar] = revisedSimplexPhaseOne(A,b);
if soln_status == -1
	% Original problem is infeasible 
	z = inf;
	x = []; 
	p=[];
	return;
else
    disp("Go to Phase2!")% soln_status == 0
    

	% Go to Phase II of the Simplex Method 
	[soln_status, z, x, p, B_inv, bvar, nbvar] = revisedSimplexPhaseTwo(c,A,b,B_inv,bvar,nbvar);
	if soln_status == 1 
		% Original problem is unbounded, x = u is a direction of unboundedness 
	else
		% Original problem has a finite optimal solution 
	end 
end
end
function [soln_status,A,b,B_inv,bvar,nbvar] = revisedSimplexPhaseOne(A_init,b_init)

[c,A,b,B_inv ,bvar,nbvar] = initializeSimplex(A_init,b_init); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% REVISED SIMPLEX METHOD for the Phase I Problem% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%{
The revised simplex method applied here solves the following problem: 
min y1 + y2 + ... + ym 
s.t. Ax+ y =b
x>= 0
y>=0

Let us follow the outline of the revised simplex method given on p.97 of Bertsimas and Tsitsiklis, "Introduction to Linear Optimization." However, we deviate from this approach in the following ways: 
I) The Phase I problem can never be unbounded, so we do not check for unboundedness after finding an improving direction.
2) Once an improving direction is found (i.e. once a nonbasic variable with a negative reduced cost is found (and decided upon) to enter the basis we force an artificial variable to leave the basis and immediately delete this artificial variable from the problem so that it cannot enter the basis in a subsequent iteration.

NOTES: 
enteringVar = the index in the set of the nonbasic variables of the nonbasic variables that will enter the basis. For example, if bvar = [1,2,3,4] (i.e ., x1,..,x4 are in the basis), and nbvar= [5,6,7,8] (i.e., x5, ... ,x8 are nonbasic variables), then enteringVar= 2 means that x5, the second index in the set of nonbasic variables will enter the basis. 
%} 
%{ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
INPUT: 
c = cost vector (n x I) 
A = constraint matrix (m x n). includes artificial variables for Phase I 
b = right-hand-side vector (m x I)
B = basis matrix (m x m)
B_inv = inverse of basis matrix B (m x m) 
bvar = row vector of basic variable indices (1 x m) 
nbvar = row vector of nonbasic variable indices (1 x (n-m)) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'%%%%%%%%%%%%%%%%%%%%
OUTPUT: 
soln_status = -1 if infeasible, 0 otherwise 
A= possibly modified (m x n) constraint matrix; redundant rows may have been eliminated 
b = possibly modified (m x 1) right-hand-side vector; 
B_inv = inverse of basis matrix B (m x m) 
bvar = row vector of basic variable indices (1 x m) 
nbvar = row vector of nonbasic variable indices (1 x (n-m)) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'%%%%%%%%%%%%%%%%%%%%
%}
soln_status = -2; 
n = length(A); 
m = length(b); 
numArtVars = length(find(c));
ArtVaridx = find(c);
bvar = ArtVaridx;% It gives the index of artificial variables
% PRECONDITIONS 
% initialize 
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf('\nPHASE ONE REVISED SIMPLEX METHOD\n');
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n'); 

phaseOneIter = -1; 

while soln_status < 0
	phaseOneIter = phaseOneIter+ 1;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	% STEP I: In a typical iteration of the revised simplex, we start with a basis B consisting of the basic %columns of A_[B(1)],....,A_[B(m)], an associated basic feasible solution x, and the inverse B_inv of the %basis matrix
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    B = A(:,bvar);
    B_inv = inv(B);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%STEP 2a: Compute the dual vector p'B=c_B' 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = c(bvar)*(B_inv);
	x = zeros(n, 1); 
	x(bvar) = (B_inv)*b;
    
 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	% STEP 2b: Compute the reduced costs of the nonbasic variables 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Choose as the entering variable "enteringVar" the *first* variable with a negative reduced cost 
    enteringVar = zeros(1,length(nbvar)+1);
    reduced_cost = zeros(1,length(nbvar));
    for j = 1:length(nbvar)       %(n-m)
		reduced_cost(j) = c(nbvar(j)) - p*A(:,nbvar(j)); 
        enteringVar(j+1) = reduced_cost(j); 
    end
    if min(enteringVar)== 0
        enteringVar = 0;
    else 
        enteringVar = find(enteringVar<0,1)-1;
    end
   
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% STEP 2c: Check necessary and sufficient optimality condition 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% If the reduced costs of all nonbasic variables are nonnegative,
	% the current basic feasible solution x_B is optimal, and the algorithm terminates. 
	if enteringVar==0
		% Current solution is OPTIMAL! 
		fprintf('Optimal PHASE ONE solution found in %d iterations', phaseOneIter); 
		break;
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% STEP 3: Compute the improving direction u = B_inv*A(:,j)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Bertsimas & Tsitsiklis call d_B = -B_inv*A(:,j) the improving direction.
	u = (B_inv)*A(:,nbvar(enteringVar));
   
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% STEP 4: Compute the stepsize (theta) in the improving direction u
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	x_B = x(bvar);
    leavingVar2 = 0;
	stepsize = 1000000; 
	countPosUjComponents = 0; 
	for i = 1:m
		if u(i)>0
			countPosUjComponents = countPosUjComponents+ 1;
			temp= x_B(i)/u(i);
		    if (temp < stepsize)
                leavingVar2 = [];
				stepsize = temp; 
				leavingVar2 = bvar(i);
                
            elseif temp == stepsize
                leavingVar2 = [leavingVar2; bvar(i)];
            end
    
	% Else if there is a tie for the leaving variable, force an artificial variable to leave the basis first 
			%TODO%
			
		end % if u(i)>0
    end
  if intersect(leavingVar2,ArtVaridx)~= 0
      leavingVar = min(intersect(leavingVar2, ArtVaridx));
  else
      leavingVar = leavingVar2(1);
  end  
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% STEP 5: Compute the new solution x_B(i) = x_B(i) – stepsize*u(i)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Form a new basis by replacing A_{B(leavingVar)} with A_{B(enteringVar)} 
	% lf y is the new basic feasible solution, the values of the new basic variables are 
	% y_{enteringVar) = stepsize and y_{B(i)} = x_{B(i)} – stepsize*u_i, for i ≠ leavingVar 
	%enteringVar

	%nbvar(enteringVar) 
    
	y = zeros(1,length(bvar));
    for i = 1: length(bvar)
        if bvar(i)~= enteringVar
            y(i) = x_B(i) - (stepsize*u(i));
        else
            y(i)= stepsize;
        end
    end
    x(bvar) = y;
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% STEP 6: Compute the new inverse basis matrix B_inv 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i= 1:length(bvar)
        if bvar(i)== leavingVar
            count = i;
        end
    end    
    B(:,count) = A(:,enteringVar);
    B_inv = inv(B);

	% UPDATE THE SET OF BASIC AND NONBASIC VARIABLES 
	temp1= bvar(count);
	bvar(count) = nbvar(enteringVar); 
	nbvar(enteringVar) = temp1;
   
end %while

% POST OPTIMALITY CHECK FOR PHASE 1 
z = dot(c(bvar),x(bvar)) ;
if z > 0
	soln_status = -1; % Original problem is INFEASIBLE!
	return;
else % z == 0, so the original problem is feasible
	soln_status = 0; % feasible
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Remove all artificial variables from the set of nonbasic variables
	% Eliminate the corresponding columns in the A matrix as well
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% KEY ASSUMPTION: an artificial variable was introduced in each row prior to starting Phase 1 
    
	nbArtIndices = intersect(nbvar,ArtVaridx); % find all nonbasic artificial variables 
    
    A(:,(nbArtIndices)) = []; 
	x((nbArtIndices)) = [];
	c((nbArtIndices)) = [];
    nbvar = setdiff(nbvar,nbArtIndices);
	

	numArtificialVarsInBasis = length(intersect(bvar,ArtVaridx)); 
    
    
	% If there are no artificial variables in the final Phase 1 basis, all artificial variables have been eliminated  
	%from the problem, so we may return the current basis to begin Phase 2 of the simplex method. 
	if numArtificialVarsInBasis == 0 
		return; 
	else
		% There is at least one artificial variable in the basis that must he removed. 
		bArtIndices = intersect(bvar,ArtVaridx); % find all basic artificial variables 
		rowsToRemove = []; 
		colsToRemove = []; 
        nbavr_only = setdiff(nbvar,ArtVaridx);
		AN= B_inv*A(:,nbavr_only); 
        
		for L = 1:length(bArtIndices)
			leaving_Var= bArtIndices(L); %: The l-th basic variable is an artificial var(Changed a line)
            for i = 1: length(bvar)
                if bvar(i)== bArtIndices(L)
                    count_1 = i;
                end
            end

			if ~any(AN(count_1,:)) == 1
                % if all entries in the leavingVar row arc zeros
				       % Remove the leavingVar row and the leavingVar-th basic variable 
				A(count_1,:) = [];
                A(:,count_1) = [];
                x(count_1) = [];
                c(count_1) = [];
                bvar(count_1) = [];
                b(count_1) = [];
				    % The basis matrix reduces from m x m to m-1 x m-1 
				    %TODO
                B(count_1,:) = [];
                B(:, count_1)= [];
                B_inv = inv(B);
            
            else

				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				% STEP 6: Compute the new inverse basis matrix B_inv 
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				% Find the first non-zero component in the leavingVar-th row of the basis %matrix A(:,bvar) 
				% NOTE: This entry may be negative, but the steps are same as in normal implementation of Phase 1. 
				 enteringVar = find(AN(count_1,:),1);
                 
				% Compute the direction vector u 
                u_new = AN(:,enteringVar);
               
			 % column vector (m x 1) 

				for i = 1:m
					if bvar(i) ~= bArtIndices(L)  %leavingVar = var index in the basis, not its original index
						B(:,bvar(i)) = A(:,bvar(i));
                    else
                        B(:,count_1) = A(:,enteringVar);
					end 
				end 
				B_inv = inv(B);
				 
				% UPDATE THE SET OF BASIC AND NONBASIC VARIABLES 
				temp= bvar(bArtIndices(L));
				bvar(bArtIndices(L)) = nbvar(enteringVar); 
				nbvar(enteringVar) = temp;
			end % else there is at least one nonzero entry in the leavingVar row

		end % 
        for L= 1: length(bArtIndices)
            A(rowsToRemove,:) = []; 
		    A(:,colsToRemove) = []; 
		
		    x(colsToRemove) = [];
		    c(colsToRemove) = [];
           
		    bvar(rowsToRemove) = []; 
		    B_inv = inv(A(:,bvar));

	    end %ELSE artificial vars still in basis 
end
end

end

function[soln_status,z,x,p,B_inv,bvar,nbvar] = revisedSimplexPhaseTwo(c,A,b,B_inv,bvar,nbvar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% REVISED SIMPLEX METHOD PHASE TWO 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%{ 
The revised simplex method applied here solves the following problem: 
min c'x
s.t. Ax=b
x>=0
 
whose dual problem is given by 
max p'b 
s.t. p'A <= c' 

We follow the outline given on p.97 of Bertsimas and Tsitsiklis, "Introduction to Linear Optimization" 
%}
 
%{ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INPUT: 
c = cost vector (n x 1) 
A= constraint matrix (m x n), includes artificial variables for Phase 1 
b = right-hand-side vector (m x 1)
 
c_B = cost vector of basic variables (m x 1) 
B = basis matrix (m x m) 
B_inv = inverse of basis matrix B (m x m) 
bvar = row vector of basic variable indices (1 x m) 
nbvar = row vector of nonbasic variable indices (1 x (n-m)) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OUTPUT: 
soln_status = 1 if unbounded, 2 if finite optimal solution 
z = objective function value: -inf if unbounded, finite otherwise 
x = (n x 1) optimal solution if soln_status = 2, a direction "u" of unboundedness otherwise 
p = an (m x 1) dual vector 
B_inv = inverse of basis matrix B (m x m) 
bvar = row vector of basic variable indices ( I x m) 
nbvar = row vector of nonbasic variable indices ( 1 x (n-m)) 

We return the inverse basis matrix, the index set of basic variables, and the index set of nonbasic variables so that we can perform a warm-start in the future. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
%[soln_status,A,b,B_inv,bvar,nbvar] = revisedSimplexPhaseOne(A_init,b_init);
n = length(A);
m = length(b);

soln_status = 0; 

% PRECONDITIONS 
% initialize 
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf('\nPHASE TWO REVISED SIMPLEX METHOD\n'); 
fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n'); 

phaseTwoIter = -1; 
while soln_status < 1
	phaseTwoIter = phaseTwoIter + 1;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%STEP I: In a typical iteration of the revised simplex, we start with a basis B consisting of the basic %columns A_{B(1)), .... A_{B(m)), an associated basic feasible solution x, and the inverse B_inv of the basis %matrix
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%STEP 2a: Compute the dual vector p'B=c_B'
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	p = c(bvar)*(B_inv);
	x = zeros(n, 1); 
	x(bvar) = (B_inv)*(b);
    

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%STEP 2b: Compute the reduced costs of the nonbasic variables
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Choose as the entering variable the first variable with a negative reduced cost
    enteringVar4 = zeros(1,length(nbvar)+1);
    reduced_cost = zeros(1,length(nbvar));
    for j = 1:length(nbvar)       %(n-m)
		reduced_cost(j+1) = c(nbvar(j)) - p*A(:,nbvar(j)); 
        enteringVar4(j) = reduced_cost(j); 
    end
   
    if min(enteringVar4)==0
        enteringVar = 0;
    else
        enteringVar = find(enteringVar4<0,1)-1; %This is the local index value in nbvar
    end
    
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%STEP 2c: Check necessary and sufficient optimality condition
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% If the reduced costs of all nonbasic variables are nonnegative, the current basic feasible solution x_B is %optimal, and the algorithm terminates. 
	if enteringVar == 0
		%Current solution is OPTIMAL!
		fprintf('Optimal PHASE TWO solution found in %d iterations', phaseTwoIter');
		z = dot(c(bvar),x(bvar));
        disp(z)
        break;
			
	end
    

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% STEP 3: Compute the improving direction u = B_inv*A(:,j) 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Bertsimas & Tsitsiklis call d_B = -B_inv*A(:,j) the improving direction.
    u = zeros(m,1);
	u = (B_inv)*A(:,nbvar(enteringVar));   % column vector (m x 1)
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% STEP 4: Compute the stepsize theta in the improving direction u
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	x_B = x(bvar);
    leavingVar = 0;
	stepsize = 1000000; 
	countPosUjComponents = 0; 
	for i = 1:m
		if u(i)>0
			countPosUjComponents = countPosUjComponents+ 1;
			temp= x_B(i)/u(i);
		    if (temp < stepsize)
                leavingVar = [];
				stepsize = temp; 
				leavingVar = bvar(i);
            end
        end
    end 
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% STEP 4b: Check unboundedness condition
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if countPosUjComponents == 0
		% LP is UNBOUNDED! 
		fprintf('Unbounded PHASE TWO solution found in %d iterations', phaseTwolter');
		z = -inf ;
		x = u; %u is a direction of unboundedness
			%Alternatively, we could output the column nbvar(enteringVar)
		return
	end


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% STEP 5: Compute the new solution x_B(i) = x_B(i) - stepsize*u(i)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Form a new basis by replacing A_[B(leavingVar)] with A_[B(enteringVar)] 
	% If y is the new basic feasible solution, the values of the new basic variables are y_{enteringVar) = %stepsize and y_{B(i)} = x_{ B(i)} - stepsize*u_i, for i \neq leavingVar 
	y = zeros(1,length(bvar));
    for i = 1: length(bvar)
        if i~= enteringVar
            y(i) = x_B(i) - (stepsize*u(i));
        else
            y(i)= stepsize;
        end
    end


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% STEP 6: Compute the new inverse basis matrix B_inv, similar to that in Phase 1 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	B = zeros(m,m);
    for i = 1:m
		if bvar(i) ~= leavingVar %leavingVar = var index in the basis, not its original index
			B(:,bvar(i)) = A(:,bvar(i));
        else
            B(:,leavingVar)= A(:,enteringVar);
		end
	end
	B_inv = inv(B);

	% UPDATE THE SET OF BASIC AND NONBASIC VARIABLES 
	temp= bvar(leavingVar); 
	bvar(leavingVar) = nbvar(enteringVar); 
	nbvar(enteringVar) = temp;
    
    
end


disp(soln_status)
disp(z)
disp(x(bvar))
disp(nbvar)
disp( bvar)
disp(B_inv)

end
function[c,A,b,B_inv ,bvar,nbvar] = initializeSimplex(A_init,b_init)
c = [0 0 0 0 0 1 1 1];
A = A_init;
b = b_init;
bvar = find(c);
var = 1: length(A);
nbvar = setdiff(var,bvar);
B_inv = inv(A(:,bvar));
end
