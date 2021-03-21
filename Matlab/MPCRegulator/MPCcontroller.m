close all;
clear all;
clc;

%% Define parameters
global Lr mp Lp Rm kt km g Br Bp Jr Jp;

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


%% Simulate nonlinear system



%% Linearize continous system using function

xref = [0; 0; 0; 0];
uref = 0;
[A,B] = linearizePendulumDynamics(xref,uref);

%% Verify lineariztion manually

%Define A,B
detH = Jp*mp*Lr^2+Jr*Jp+1/4*Jr*mp*Lp^2;
Ac = 1/detH * [0 0 detH 0; 0 0 0 detH;...
    0 1/4*mp^2*Lp^2*Lr*g -(Jp+1/4*mp*Lp^2)*Br -0.5*mp*Lp*Lr*Bp;...
    0 1/2*mp*Lp*g*(Jr+mp*Lr^2) -0.5*mp*Lp*Lr*Br -(Jr+mp^Lr^2)*Bp];

Bc = 1/detH * [0; 0; Jp+1/4*mp*Lp^2; 0.5*mp*Lp*Lr];
nx = size(A,1);
nu = size(B,2);

%Include the actuator dynamics
Ac = A-[zeros(nx,2) B*-kt*km/Rm zeros(nx,1)];
Bc = kt/Rm*B;

%% Simulate the system
