function [H,h,const] = getCostMatrices(T,S,x0,Q,R,P,N)
%%getCostMatrices outputs the matrices H,h and const in:
%       Vn(Un,x0) = 1/2*U_N'*H*U_N+h'*U_N+const;
%                   where V_N(U_N,x0) is generated from the problem
%       Vn(Xn,Un) = 1/2*sum_{k=0}^{k=N}[x(k)'*Q*x(k)+u(k)'*R*u(k)]+x(N)'*P*x(N);
%                   with Xn eliminated by substitution as a function of Un:
%       Xn = T*x0+S*Un;
%                   where Xn = (x(0), x(1), ..., x(N)) and Un =(u(0), u(1),..., u(N-1));
% Input:
%   T,S     :   prediction matrices T,S in Xn = T*x0+S*Un;
%   x0      :   initial state
%   Q,R,P   :   weight matrices Q,R and the quadratic terminal cost P
%   N       :   prediction horizon
%   M       :   (optional) the weight matrix M penalizing crossterms
%               x(k)'M*u(k)

%Build the weighting matrices
Qt = blkdiag(kron(eye(N),Q),P);
Rt = kron(eye(N),R);

%Compute the optimization matrices
H = Rt + S'*Qt*S;
h = (Qt*S)'*T*x0;
const = 0.5*x0'*T'*Qt*T*x0;
end