function [xr,ur] = getOTSReference(model,constraints,rsp,dh)
%Input:
%   model           : model of system and disturbance
%   constraints     : constraints on x and u
%   rsp             : desired setpoint
%   dh              : estimate of disturbance
%Output:
%   xr,ur           : optimal steady state target

%Unpack
A = model.A;
B = model.B;
Bd = model.Bd;
C = model.C;
Cd = model.Cd;
H = model.H;
Ax = constraints.Ax;
bx = constraints.bx;
Au = constraints.Au;
bu = constraints.bu;
nx = size(A,1);
ni = size(B,2);

%Build the constraints
Aeq = [eye(nx)-A -B; H*C zeros(size(H,1),ni)];
beq = [Bd*dh; rsp-H*Cd*dh];
Aineq = blkdiag(Ax,Au);
bineq = [bx; bu];

% Objective function ||(.)||2
Hp = eye(ni+nx);

%Solve OTS
cvx_begin
        variable z(nx+ni);
        minimize( z'*Hp*z);
        %Define constraints
        subject to        
        Aeq*z == beq;
        Aineq*z <= bineq;
cvx_end 

%Extract solution
xr = z(1:nx);
ur = z(nx+1:end);

end