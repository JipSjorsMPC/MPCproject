close all;
clear all;
clc;

%% Define parameters
global Lr mp Lp Rm kt km g Br Bp Jr Jp u;

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

%Lump motor, hub and rotary arm
Jr = Jr + Jh + Jm;

%% Linearize continous system

%Linearize continuous dynamics around reference point
xRef = [0; 0; 0; 0];
uRef = 0;
[Ac,Bc] = linearizePendulumDynamics(xRef,uRef);

%Define the C and D matrix
nx = size(Ac,1);
ni = size(Bc,2);
C = eye(nx);
no = size(C,1);
D = zeros(no,ni);

%Verify linearization manually
%Define A,B
detH = Jp*mp*Lr^2+Jr*Jp+1/4*Jr*mp*Lp^2;
Am = 1/detH * [0 0 detH 0;...
    0 0 0 detH;...
    0 1/4*mp^2*Lp^2*Lr*g -(Jp+1/4*mp*Lp^2)*Br -0.5*mp*Lp*Lr*Bp;...
    0 1/2*mp*Lp*g*(Jr+mp*Lr^2) -0.5*mp*Lp*Lr*Br -(Jr+mp*Lr^2)*Bp];

Bm = 1/detH * [0; 0; Jp+1/4*mp*Lp^2; 0.5*mp*Lp*Lr];
nx = size(Am,1);
nu = size(Bm,2);

%Include the actuator dynamics
Am = Am-[zeros(nx,2) Bm*kt*km/Rm zeros(nx,1)];
Bm = kt/Rm*Bm;

%Indeed checked and verified this gives same result
%% Discretizatize linearized system

%Discretize system using sample time Ts
Ts = 0.02;
sysd = c2d(ss(Ac,Bc,C,D),Ts,'zoh');
[A,B,C,D] = ssdata(sysd);

%% Solve MPC problem with MPT3 toolbox

%Weights
Q = blkdiag(10,100,1,5);
R = 0.01;
N = 20;

%Define model
model = LTISystem('A', A, 'B', B, 'C', C, 'D', D, 'Ts', Ts);

%State and input constraints
model.u.min = -1;
model.u.max = 1;
model.x.min = [-pi -pi/10 -20*2*pi -20*2*pi];
model.x.max = [pi pi/10 20*2*pi 20*2*pi];

%Set cost function
model.x.penalty = QuadFunction(Q);
model.u.penalty = QuadFunction(R);

%Terminal set and terminal penalty
model.x.with('terminalPenalty');
model.x.with('terminalSet');
model.x.terminalPenalty = model.LQRPenalty;
model.x.terminalSet = model.LQRSet;

%Get controller
mpc_controller = MPCController(model, N);
loop = ClosedLoop(mpc_controller, model);

%% Now simulate on the nonlinear system

%Prepare simulation
Tsim = 2;
Nsim = round(Tsim/Ts);
x0 = [0; 4*pi/180; 0; 0];
x = zeros(nx,Nsim+1);
x(:,1) = x0;
time = (0:Nsim+1)*Ts;

data = loop.simulate(x0,Nsim);

%Simulate the system using zero order hold
for k = 1:Nsim 
    uOpt = mpc_controller.evaluate(x(:,k)); 
    unonlin(:,k) = uOpt; 
    u = unonlin(:,k);
    [tout,x_interval] = ode45(@nonlinearPendulumDynamics,[time(k) time(k+1)],x(:,k));
    
    %Evaluate the next state
    x(:,k+1) = x_interval(end,:)';

    %Normalize angles
    x(1:2,k+1) = wrapToPi( x(1:2,k+1) );
end
%%

%Define colors
orange = [255 128 0]/256;
purple = [102 0 204]/256;

figure('Name','Linear vs Nonlinear','Position',[850 500 400 400])
subplot(3,1,1)
    stairs(Ts*(0:Nsim),data.X(1,:)'*180/pi,'Color',purple,'LineWidth',1); hold on
    stairs(Ts*(0:Nsim),x(1,:)'*180/pi,'m','LineWidth',1); hold on
    stairs(Ts*(0:Nsim),data.X(2,:)'*180/pi,'b','LineWidth',1); hold on
    stairs(Ts*(0:Nsim),x(2,:)'*180/pi,'c','LineWidth',1);
%     xlabel('t[s]');
    ylabel('$\theta, \alpha [deg]$','interpreter','latex','fontsize',12);   
    legend('$\theta(t)(Lin.)$','$\theta(t)(Nonlin.)$','$\alpha(t) (Lin.)$','$\alpha(t) (Nonlin.)$','interpreter','latex','fontsize',8);
%     title('Linear vs Nonlinear: state \theta and \alpha');
subplot(3,1,2);
    stairs(Ts*(0:Nsim),data.X(3,:)'*180/pi,'Color',purple,'LineWidth',1); hold on
    stairs(Ts*(0:Nsim),x(3,:)'*180/pi,'m','LineWidth',1); hold on
    stairs(Ts*(0:Nsim),data.X(4,:)'*180/pi,'b','LineWidth',1); hold on
    stairs(Ts*(0:Nsim),x(4,:)'*180/pi,'c','LineWidth',1); hold on
%     xlabel('t[s]');
    ylabel('$\dot{\theta}$,$\dot{\alpha}$ [deg/s]','interpreter','latex');
    legend('$\dot{\theta}(t)(Lin.)$','$\dot{\theta}(t)(Nonlin.)$','$\dot{\alpha}(t)(Lin.)$','$\dot{\alpha}(t)(Nonlin.)$','interpreter','latex','fontsize',8);
%     title('Linear vs Nonlinear: state \theta_d and \alpha_d');
subplot(3,1,3);
    stairs(Ts*(0:Nsim-1),data.U,'r','LineWidth',1); hold on
    stairs(Ts*(0:Nsim-1),unonlin','Color',orange,'LineWidth',1); hold on
    plot((0:Nsim-1)*Ts,[model.u.max;model.u.min]*ones(1,Nsim),'g--');
    xlabel('t[s]');
    ylabel('$u[V]$','interpreter','latex');
    ylim([-1.2 1.2]);
    legend('$V_m(t)(Lin.)$','$V_m(t) (Nonlin.)$','interpreter','latex','fontsize',8);
%     title('Linear vs Nonlinear: control input u(t)');
