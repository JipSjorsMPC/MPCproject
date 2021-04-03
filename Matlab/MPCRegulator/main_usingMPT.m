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
Ts = 0.05;
sysd = c2d(ss(Ac,Bc,C,D),Ts,'zoh');
[A,B,C,D] = ssdata(sysd);

%% Solve MPC problem with MPT3 toolbox

%Weights
Q = blkdiag(1,100,1,10);
R = 0.01;
N = 5;

%Define model
model = LTISystem('A', A, 'B', B, 'C', C, 'D', D, 'Ts', Ts);

%State and input constraints
model.u.min = -10;
model.u.max = 10;

model.x.min = [-pi/2 -pi/4 -5*2*pi -5*2*pi];
model.x.max = [pi/2 pi/4 5*2*pi 5*2*pi];

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

%Simulate
loop = ClosedLoop(mpc_controller, model);
x0 = [-5*pi/180; 2*pi/180; 0; 0];
e0 = x0-xRef;
Tsim = 10;
Nsim = Tsim/Ts;
data = loop.simulate(e0,Nsim);

%Plot the results
figure();
subplot(3,2,1)
stairs((0:Nsim)*Ts,data.X(1,:)'*180/pi);
xlabel('t');
ylabel('\theta[deg]');
title('x_1(t)');
subplot(3,2,2)
stairs((0:Nsim)*Ts,data.X(2,:)'*180/pi);
xlabel('t');
ylabel('\alpha[deg]');
title('x_2(t)');
subplot(3,2,3)
stairs((0:Nsim)*Ts,data.X(3,:)'*180/pi);
xlabel('t');
ylabel('\theta_d[deg/s]');
title('x_3(t)');
subplot(3,2,4)
stairs((0:Nsim)*Ts,data.X(4,:)'*180/pi);
xlabel('t');
ylabel('\alpha_d[deg/s]');
title('x_4(t)');
subplot(3,2,[5 6]);
stairs((0:Nsim-1)*Ts,data.U,'r');
title('Control input du(t)')
xlabel('t');
ylabel('u [V]');


%% Plot Xn and Xf

%Compute DARE gain and cost
[K, P] = dlqr(A, B, Q, R);
K = -K; 

% Calculate Xn and Xf (maximum LQR-invariant set) using normal penalty.
[Xn, V, Z] = findXn(A, B, K, N, model.x.min, model.x.max, model.u.min, model.u.max, 'lqr');
XN = Polyhedron(Xn{end}.A,Xn{end}.b);
Xf = Polyhedron(Xn{1}.A,Xn{1}.b);
XN.minHRep();
Xf.minHRep();

%Plot the ROA in 3D
XN_2D = XN.projection(1:3);
Xf_2D = Xf.projection(1:3);
plot(XN_2D,'color','g',Xf_2D,'color','r');
legend('$$X_N$$','$$X_f$$','Interpreter','latex');


%%











%% Now simulate on the nonlinear system

x = zeros(nx,Nsim+1);
e = zeros(nx,Nsim);
x(:,1) = x0;
uApl = zeros(ni,Nsim);
time = (0:Nsim-1)*Ts;
for k = 1:Nsim    
%Compute the optimal control input
e(:,k) = x(:,k)-xRef;
uOpt = mpc_controller.evaluate(e(:,k));

%Simulate the system using zero order hold
uApl(:,k) = uOpt + uRef;
u = uApl(:,k);
[tout,x_interval] = ode45(@nonlinearPendulumDynamics,[t(k) t(k+1)],x(:,k));

%Evaluate the next state
x(:,k+1) = x_interval(end,:)';

end

figure();
subplot(3,2,1)
stairs(Ts*(0:Nsim),x(1,:)'*180/pi);
xlabel('t[s]');
ylabel('$\theta(t)$ [deg]','interpreter','latex');
title('State x_1(t)');
subplot(3,2,2);
stairs(Ts*(0:Nsim),x(2,:)'*180/pi);
xlabel('t[s]');
ylabel('$\alpha(t)$ [deg]','interpreter','latex');
title('State x_2(t)');
subplot(3,2,3);
stairs(Ts*(0:Nsim),x(3,:)'*180/pi);
xlabel('t[s]');
ylabel('$\dot{\theta}(t)$ [deg/s]','interpreter','latex');
title('State x_3(t)');
subplot(3,2,4);
stairs(Ts*(0:Nsim),x(4,:)'*180/pi);
xlabel('t[s]');
ylabel('$\dot{\alpha}(t)$ [deg/s]','interpreter','latex');
title('State x_4(t)');
subplot(3,2,[5 6]);
stairs(Ts*(0:Nsim-1),uApl','r');
title(['Control input using MPC with N = ', num2str(N)]);
xlabel('t[s]');
ylabel('u[V]');





