function [T,S] = getStatePredictionMatrices(A,B,N,C)
% getStatePrediction matrices outputs the predictionmatrices T,S in 
% X = T*x0+S*U , where X =[x(0), x(1) .... x(N)]^T and U = [u(0), u(1), ... u(N-1)]^T. 
% If C specified, instead the output matrices will be defined according to
% Y = T*x0+S*U, where Y = [y(0), y(1), ..., y(N)];
%
% Input:
%   A,B,N   :   discrete matrices A(ni x ni),B(nx x ni) and prediction horizon N
%   C       :   (optional) C matrix if the output is required instead
% Ouput:
%   T       :   the prediction matrix T (n*(N+1) x n)
%   S       :   the prediction matrix S (n*(N+1) x (ni*N))

nx = size(A,1);
ni = size(B,2);

% C is an optional argument
if ~exist('C','var') || isempty(C)
    C = eye(size(A));   
else
    no = size(C,1);
    nx = no;
end

%Generate P matrix
T = [];
for i = 0:N
    T = [T; C*A^(i)];
end

%Generate the S matrix
S = zeros(nx,ni*N);
Srow = zeros(nx,ni*N);
for i = 1:N
 Srow = circshift(Srow,ni,2);
 Srow(:,1:ni) = C*A^(i-1)*B; 
 S = [S ; Srow]; 
end
    
end

