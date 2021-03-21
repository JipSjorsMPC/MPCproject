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


%% Simulate autonomous nonlinear system

%Set input to zero and simulate around small deviation from upright
%position
global u
u = 0;
x0 = [0; 0.25; 0; 0];
Tf = 5;
[t,x] = ode45(@nonlinearPendulumDynamics,[0 Tf],x0);

%Plot the results
subplot(2,2,1);
plot(t,x(:,1),'c');
title('$\theta(t)$','interpreter','latex');
xlabel('$t [s]$','interpreter','latex');
ylabel('$x_1$[rad]','interpreter','latex');
legend('$x_1(t)$','interpreter','latex');
subplot(2,2,2);
plot(t,x(:,2),'g');
title('$\alpha(t)$','interpreter','latex');
xlabel('$t [s]$','interpreter','latex');
ylabel('$x_2$[rad]','interpreter','latex');
legend('$x_2(t)$','interpreter','latex');
subplot(2,2,3);
plot(t,x(:,3),'b');
title('$\dot{\theta(t)}$','interpreter','latex');
xlabel('$t [s]$','interpreter','latex');
ylabel('$x_3$[rad/s]','interpreter','latex');
legend('$x_3(t)$','interpreter','latex');
subplot(2,2,4);
plot(t,x(:,4),'r');
title('$\dot{\alpha(t)}$','interpreter','latex');
xlabel('$t[s]$','interpreter','latex');
ylabel('$x_4$[rad/s]','interpreter','latex');
legend('$x_4(t)$','interpreter','latex');


%% Linearize continous system

%Linearize continuous dynamics around reference point
xref = [0; 0; 0; 0];
uref = 0;
[A,B] = linearizePendulumDynamics(xref,uref);

%Define the C and D matrix
nx = size(A,1);
ni = size(B,2);
C = eye(nx);
no = size(C,1);
D = zeros(no,ni);

%% Verify linearization manually
% %Define A,B
% detH = Jp*mp*Lr^2+Jr*Jp+1/4*Jr*mp*Lp^2;
% Am = 1/detH * [0 0 detH 0;...
%     0 0 0 detH;...
%     0 1/4*mp^2*Lp^2*Lr*g -(Jp+1/4*mp*Lp^2)*Br -0.5*mp*Lp*Lr*Bp;...
%     0 1/2*mp*Lp*g*(Jr+mp*Lr^2) -0.5*mp*Lp*Lr*Br -(Jr+mp*Lr^2)*Bp];
% 
% Bm = 1/detH * [0; 0; Jp+1/4*mp*Lp^2; 0.5*mp*Lp*Lr];
% nx = size(Am,1);
% nu = size(Bm,2);
% 
% %Include the actuator dynamics
% Am = Am-[zeros(nx,2) Bm*kt*km/Rm zeros(nx,1)];
% Bm = kt/Rm*Bm;

%Indeed checked and verified this gives same result
%% Discretizatize linearized system

%Discretize system using sample time Ts
h = 0.1;
sys = ss(A,B,C,D);
sysd = c2d(sys,h,'zoh');
[Ad,Bd,Cd,Dd] = ssdata(sysd);

%% Solve MPC problem

%Parameters
Q





%Get infinite control Riccati cost +control
[P,~,K] = dare(A,B,Q,R);
K = -K;




