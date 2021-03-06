clc, clear

%% Inverted rotary pendulum parameters
global mr Lr mp Lp mh rh Jh Rm kt km Jm Lm g Br Bp Jr Jp

mr = 0.095;         % Rotary arm mass
Lr = 0.085;         % Rotary arm length
mp = 0.024;         % Pendulum link mass
Lp = 0.129;         % Pendulum link length
mh = 0.016;         % Module attachment hub mass
rh = 0.0111;        % Module attachment hub radius
Jh = 0.6*10^-6;     % Module attachment hub moment of inertia
Rm = 8.4;           % Terminal resistance
kt = 0.042;         % Torque constant
km = 0.042;         % Motor back-emf constant
Jm = 4.0*10^-6;     % Rotor inertia
Lm = 1.16*10^-3;    % Rotor inductance
g = 9.81;          	% Gravitational constant

Br = 8.0508*10^-4;  % Damping coefficient rotary arm link
Bp = 1.4*10^-5;     % Damping coefficient pendulum link
Jr = 2.2923*10^-4;	% Pendulum arm moment of inertia
Jp = 1.2551*10^-4;	% Rotary arm moment of inertia

%Constraints
Vub = 10;           % Upperbound voltage
Vlb = -10;          % Lowerbound voltage

%Initial conditions
x0 = [0 pi 0 0]; 

tfin=10;
tspan=0:1e-2:tfin;       % To obtain solution at specific times

%% Matrices motion equations 
syms x dx [4 1] , syms u

H = [   Jr+mp*Lr^2+.25*mp*Lp^2*sin(x(2))^2     -.5*mp*Lr*Lp*cos(x(2));
        .5*mp*Lr*Lp*cos(x(2))                   Jp+.25*mp*Lp^2  ];

C = [	.25*mp*Lp^2*sin(2*x(2))*x(4)+Br         .5*mp*Lr*Lp*x(4)*sin(x(2));
       -.125*mp*Lp^2*x(3)*sin(2*x(2))           Bp  ];

G = [   0;
       -.5*mp*Lp*g*sin(x(1))];

B = [   1;
        0];
    
dx(1) = x(3);
dx(2) = x(4);
dx(3:4) = (-C*x(3:4) - G + B*u)\H;
