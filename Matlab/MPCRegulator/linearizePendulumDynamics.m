function [A,B] = linearizePendulumDynamics(xref,uref)
%linearizePendulumDynamics obtains the A,B matrices around the reference
%   point (xref,uref). Let ex =  x- xref and eu = u-uref. Then the
%   matrices describe the linearized dynamics around the reference point
%   dexdt = A*ex+Beu
%       y = C*ex+D*eu
% Input:
%   model   : The model of the nonlinear pendulum to be linearized
%   xref    : a reference state
%   uref    : a reference input

% Output:   : the A,B matrices

global Lr mp Lp Rm kt km g Br Bp Jr Jp

syms x1 x2 x3 x4 u

% Define the manipulator form H(q)qdd+C(qd,q)qd+G(q)=Bu
H = [Jr+mp*Lr^2+1/4*mp*Lp^2*sin(x2)^2 -1/2*mp*Lr*Lp*cos(x2);...
    -1/2*mp*Lr*Lp*cos(x2) Jp+1/4*mp*Lp^2];
C = [1/4*mp*Lp^2*sin(2*x2)*x4+Br 0.5*mp*Lr*Lp*x4*sin(x2);...
-1/8*mp*Lp^2*x3*sin(2*x2) Bp];
G = [0; -0.5*mp*Lp*g*sin(x2)];
Bt = [1; 0];

%Include actuator dynamics
tau = kt*(u-km*x3)/Rm;

%Obtain dxdt = f(x,u,t);
dxdt = [x3; x4;...
    H\(Bt*tau-C*[x3; x4]-G)];

%Now linearize dxdt = f(x,u,t) w.r.t x and u;
dfdx =  simplify(jacobian(dxdt,[x1; x2; x3; x4]));
dfdu = simplify(jacobian(dxdt,u));

%Get the A = dfdx(xref,uref) and B = dfdu(xref,uref);
A = subs(dfdx,[x1 x2 x3 x4 u],[xref' uref]); 
B = subs(dfdu,[x1 x2 x3 x4 u],[xref' uref]);

end

