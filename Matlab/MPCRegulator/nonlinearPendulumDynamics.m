function dxdt = nonlinearPendulumDynamics(t,x)
%nonlinearPendulumDynamics outputs the dynamics of the rotary pendulum
%   in the form xdot = f(x,t); The state x will be defined as
%   x = [theta alpha thetadot alphadot]';
%
% Input:
%   t           :   time vector
%   x           :   state vector
% Output:
%   [dxdt]      :   the dynamics dxdt=f(x,t)

global Lr mp Lp Rm kt km g Br Bp Jr Jp
global u

% Define the manipulator form H(q)qdd+C(qd,q)qd+G(q)=Bu
H = [   Jr+mp*Lr^2+1/4*mp*Lp^2*sin(x(2))^2      -1/2*mp*Lr*Lp*cos(x(2));...
        -1/2*mp*Lr*Lp*cos(x(2))                 Jp+1/4*mp*Lp^2 ];
    
C = [   1/4*mp*Lp^2*sin(2*x(2))*x(4)+Br         0.5*mp*Lr*Lp*x(4)*sin(x(2));...
        -1/8*mp*Lp^2*x(3)*sin(2*x(2))           Bp ];
    
G = [   0;...
        -0.5*mp*Lp*g*sin(x(2)) ];
    
B = [   1;... 
        0 ];

%Include actuator dynamics
tau = kt*(u-km*x(3))/Rm;

%Define dxdt = f(x,u,t)
dxdt = [x(3); x(4);...
    H\(B*tau-C*[x(3); x(4)]-G)];

end

