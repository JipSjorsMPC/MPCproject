function [H,h,kstar] = computeMaxConstraintAdmissibleSet(Ak,F,g)
% computeOutputAdmissibleSet computes the set of initial states x given by
%   H x <= h,
% such that the constraint
%   F x <= g
% is satisfied for all future k for any x(k) starting in this set under the autonomous system
%   x(k+1)= Ak x(k)

% Input: 
%   g       :   s x 1 column vector, where s are the number of constraints
%   F       :   s x nx matrix where s in the number of constraints and nx is the
%               number of states.
% Output:
%   H,h     :   the polytope that is maximal constraint admissible given by
%               H*x<=h

nx = size(Ak,1);
s = size(F,1);

% Algorithm implementation
kstarIsFound = false;
k=0;
while ~kstarIsFound 
    % Set the constraints for the optimization problem
     H=[];
     h=[];
     for t = 0:k
         Ht = [];
         for j = 1:s
            Ht(j,:) = F(j,:)*Ak^t;
         end
         H = [H; Ht];         
     end
     h = repmat(g,k+1,1);
     
    optVal = zeros(1,s);
    for i=1:s
         % Set the optimization problem: objective function and constraints
         fi = F(i,:)*Ak^(k+1);
         
         % Solve the optimization problem with CVX
         cvx_precision best
         cvx_begin 
            variable x(nx)
            maximize(fi*x-g(i))
            subject to
            H*x <= h;
        cvx_end
        
        % Save optimal value
        optVal(i) = cvx_optval;
    end
    
    % Evaluating optimality condition
    if (all(optVal<=0-eps) && ~strcmp(cvx_status,'Unbounded'))
        kstarIsFound = true;
        kstar = k;       
    else        
        k = k+1;
    end
end


end