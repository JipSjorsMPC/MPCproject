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
Ts = 0.05;
sysd = c2d(ss(Ac,Bc,C,D),Ts,'zoh');
[A,B,C,D] = ssdata(sysd);

%% Solve MPC problem with MPT3 toolbox

%Weights
Q = blkdiag(1,10,0.1,0.1);
R = 0.1;
N = 50;

%Define model
model = LTISystem('A', A, 'B', B, 'C', C, 'D', D, 'Ts', Ts);

%State and input constraints
model.u.min = -20;
model.u.max = 20;
model.x.min = [-pi -pi -inf -inf];
model.x.max = [pi pi inf inf];

%Set cost function
model.x.penalty = QuadFunction(Q);
model.u.penalty = QuadFunction(R);

%Terminal set and terminal penalty
model.x.with('terminalPenalty');
model.x.with('terminalSet');
P = model.LQRPenalty;
Xf = model.LQRSet;
model.x.terminalPenalty = P;
model.x.terminalSet = Xf;

%Get controller
mpc_controller = MPCController(model, N);

%Simulate
loop = ClosedLoop(mpc_controller, model);
x0 = [pi/10; pi/10; 0; 0];
e0 = x0-xRef;
Tsim = 2;
Nsim = Tsim/Ts;
data = loop.simulate(e0,Nsim);
figure();
subplot(2,1,1)
stairs(0:Nsim,data.X');
title('State error $\tilde{x}(k)$','interpreter','latex')
xlabel('k');
ylabel('$\tilde{x}$(k)','interpreter','latex');
legend('$\tilde{\theta}(k)$','$\tilde{\alpha}(k)$','$\tilde{\dot{\theta}}(k)$','$\tilde{\dot{\alpha}}(k)$','interpreter','latex');
subplot(2,1,2);
stairs(0:Nsim-1,data.U);
title('Input $\tilde{u}(k)$','interpreter','latex')
xlabel('k');
ylabel('u [V]','interpreter','latex');
legend('$\tilde{u}(k)$','interpreter','latex');

% %Plot Xn and Xf
% figure();
% Xn = model.invariantSet();
% %plot(Xf,'color','g',Xn,'color','r');
% plot(Xf,'color','g');
% legend('$$X_N$$','$$X_f$$','Interpreter','latex');


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
stairs(Ts*(0:Nsim),x(1,:)'');
xlabel('t[s]');
ylabel('$\theta(t)$ [rad]','interpreter','latex');
title('State x_1(t)');
subplot(3,2,2);
stairs(Ts*(0:Nsim),x(2,:)');
xlabel('t[s]');
ylabel('$\alpha(t)$ [rad]','interpreter','latex');
title('State x_2(t)');
subplot(3,2,3);
stairs(Ts*(0:Nsim),x(3,:)'');
xlabel('t[s]');
ylabel('$\dot{\theta}(t)$ [rad/s]','interpreter','latex');
title('State x_3(t)');
subplot(3,2,4);
stairs(Ts*(0:Nsim),x(4,:)'');
xlabel('t[s]');
ylabel('$\dot{\alpha}(t)$ [rad/s]','interpreter','latex');
title('State x_4(t)');
subplot(3,2,[5 6]);
stairs(Ts*(0:Nsim-1),uApl','r');
title(['Control input using MPC with N = ', num2str(N)]);
xlabel('t[s]');
ylabel('u[V]');



% %% Solve with methods from Class
% %Parameters
% Q = blkdiag(100,100,10,10);
% R = 1;
% N = 8;
% 
% %Get unconstrained infinite control Riccati cost and control
% [P,~,G] = dare(A,B,Q,R);
% K = -G;
% 
% %Constraints (see further below)
% % |x1(k)|<= pi/4    for all k
% % |u(k)|<= 5        for all k
% % |x2(k)| <= pi/8  for all k
% 
% %Compute the terminal set Xf
% %Build constraints 
% Ak = A+B*K;
% Ff = vertcat([1 0 0 0; -1 0 0 0],[0 1 0 0; 0 -1 0 0],[K ;-K]);
% gf = vertcat([pi/4; pi/4],[pi/8; pi/8],[5;5]);
% 
% %Compute Xf under u =Kx such that x satisfies input and state constraints
% [Hf,hf,kstar] = computeMaxConstraintAdmissibleSet(Ak,Ff,gf);
% 
% %Plot it
% figure();
% Xf = Polyhedron(Hf,hf);
% 
% %Remove redundancies
% Xf.minHRep();
% Hf = Xf.A;
% hf = Xf.b;
% 
% %Get Prediction matrices
% [T,S] = getStatePredictionMatrices(A,B,N);
% 
% %First build G,H and phi matrix in Gx+Hu+phi<=0
% F1h = kron(eye(N+1),[1 0 0 0; -1 0 0 0; 0 1 0 0; 0 -1 0 0]);
% g1h = repmat([pi/4;pi/4;pi/8;pi/8],N+1,1);
% F2h = kron(eye(N),[1;-1]);
% g2h = repmat(5*[1;1],N,1);
% G = [F1h*T;...
%  zeros(size(F2h,1),size(A,1));...
%  Hf*T(end-nx+1:end,:)];
% H = [F1h*S; F2h; Hf*S(end-nx+1:end,:)];
% phi = [-g1h; -g2h; -hf];
% 
% % Compute Xn, i.e. the region of attraction
% [PXn,gammaXn] = computeRegionOfAttraction(G,H,phi);
%  
% % get of Attraction (ROA)
% figure();
% Xn = Polyhedron(PXn,-gammaXn);
% plot(Xn,'color','g',Xf,'color','r');
% legend('$$X_N$$','$$X_f$$','Interpreter','latex');
% 
% 



