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
%Simulation without controller
u=@(t,x) 0; % Define the input, as a function of time and/or state; in this
            % case, u is the zero function
f=@(t,x)manipulator(t,x,u); % ODE definition    
[t,x]=ode45(f,tspan,x0);       % Solving ODE
