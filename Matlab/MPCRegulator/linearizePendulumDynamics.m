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

syms x [4 1], syms u

% Define the manipulator form H(q)qdd+C(qd,q)qd+G(q)=Bu
H = [Jr+mp*Lr^2+1/4*mp*Lp^2*sin(x(2))^2 -1/2*mp*Lr*Lp*cos(x(2));...
    -1/2*mp*Lr*Lp*cos(x(2)) Jp+1/4*mp*Lp^2];
C = [1/4*mp*Lp^2*sin(2*x(2))*x4+Br 0.5*mp*Lr*Lp*x(4)*sin(x(2));...
-1/8*mp*Lp^2*x(3)*sin(2*x(2)) Bp];
G = [0; -0.5*mp*Lp*g*sin(x(2))];
Bt = [1; 0];

%Include actuator dynamics
tau = kt*(u-km*x(3))/Rm;

%Obtain dxdt = f(x,u,t);
dxdt = [x(3); x(4);...
    H\(Bt*tau-C*[x(3); x(4)]-G)];

%Now linearize dxdt = f(x,u,t) w.r.t x and u;
dfdx =  simplify(jacobian(dxdt,x));
dfdu = simplify(jacobian(dxdt,u));

%Get the A = dfdx(xref,uref) and B = dfdu(xref,uref);
A = subs(dfdx,[x.' u],[xref' uref]); 
B = subs(dfdu,[x.' u],[xref' uref]);

end
