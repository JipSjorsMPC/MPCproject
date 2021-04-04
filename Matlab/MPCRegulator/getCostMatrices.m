function [H,h,const] = getCostMatrices(T,S,x0,Q,R,P,N,M)
%%getCostMatrices outputs the matrices (H,h,C) in:
%       V_N(U_N,x0) = 1/2*U_N'*H*U_N+h'*U_N+const;
%                   where V_N(U_N,x0) is generated from the problem
%       V_N(X_N,U_N) = 1/2*sum_{k=0}^{k=N}[x(k)'*Q*x(k)+u(k)'*R*u(k)+2*x(k)'*M*u(k)]+x(N)'*P*x(N);
%                   with X_N eliminated by substitution as a function of U_N:
%       X_N = T*x0+S*U_N;
%                   where X_N = (x(0), x(1), ..., x(N)) and U_N=(u(0), u(1),..., u(N-1));
% Input:
%   T,S     :   prediction matrices T,S in X_N = T*x0+S*U_N;
%   x0      :   initial state
%   Q,R,P   :   weight matrices Q,R and the quadratic terminal cost P
%   N       :   prediction horizon
%   M       :   (optional) the weight matrix M penalizing crossterms
%               x(k)'M*u(k)

%Check if there are any crossterms
if ~exist('M','var') || isempty(M)    
    M = zeros(size(x0,1),size(S,2)/N);
end

%Build the weighting matrices
Qt = blkdiag(kron(eye(N),Q),P);
Rt = kron(eye(N),R);
Mt = vertcat(kron(eye(N),M),zeros(size(x0,1),size(S,2)));

%Compute the optimization matrices
H = Rt + S'*Qt*S+2*S'*Mt;
h = (Qt*S+Mt)'*T*x0;
const = 0.5*x0'*T'*Qt*T*x0;
end