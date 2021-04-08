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

%% Obtain model

%Linearize continuous dynamics around reference point
xRef = [0; 0; 0; 0];
uRef = 0;
[Ac,Bc] = linearizePendulumDynamics(xRef,uRef);

%Define the C and D matrix
nx = size(Ac,1);
ni = size(Bc,2);
Cc = eye(nx);
no = size(Cc,1);
Dc = zeros(no,ni);

%Discretize system using sample time Ts
Ts = 0.05;
sysd = c2d(ss(Ac,Bc,Cc,Dc),Ts,'zoh');
[A,B,C,D] = ssdata(sysd);

%% Parameters
%Parameters
Q = blkdiag(10,100,1,5);
R = 0.01;
N = 20;
model.u.min = -2;
model.u.max = 2;
model.x.min = [-pi/2; -pi/10; -20*2*pi; -20*2*pi];
model.x.max = [pi/2; pi/10; 20*2*pi; 20*2*pi];
x0 = [0; 2*pi/180; 0; 0];

%% Build the MPC problem

% Get unconstrained infinite control Riccati cost and control
[P,~,G] = dare(A,B,Q,R);
K = -G;

% Define state and input polyhedrons
Au = [1; -1];
bu = [model.u.max; -model.u.min];
Ax = [eye(nx); -eye(nx)];
bx = [model.x.max; -model.x.min];

% Calculate Xn and Xf (maximum LQR-invariant set) 
[Xn, V, Z] = findXn(A, B, K, N, model.x.min, model.x.max, model.u.min, model.u.max, 'lqr');
XN = Polyhedron(Xn{end}.A,Xn{end}.b);
XN.minHRep();
Xf = Polyhedron(Xn{1}.A,Xn{1}.b);
Xf.minHRep();
plot(XN.projection(1:3),'color','g',Xf.projection(1:3),'color','r');
legend('$$X_N$$','$$X_f$$','Interpreter','latex');
xlabel('\theta');
ylabel('\alpha');
zlabel('$\dot{\theta}$','interpreter','latex');

%Get Prediction matrices
[T,S] = getStatePredictionMatrices(A,B,N);

%Build optimization constaint matrices Fu<=g(x(0))
F1h = blkdiag(kron(eye(N),Ax),Xf.A);
g1h = [repmat(bx,N,1); Xf.b];
F2h = kron(eye(N),Au);
g2h = repmat(bu,N,1);
F_con = [F1h*S; F2h];
g_con = [g1h-F1h*T*x0; g2h];

%% Simulate the online optimization problem for the linear system
%Prepare problem
Tsim = 1.8; %[s]
ksim = round(Tsim/Ts);
u_apl = zeros(ni,ksim);
x = zeros(nx,ksim);  
x(:,1) = x0;
for k = 1:ksim     
    
    %Recompute cost matrices in Vn(Un,x0), parametric in x0    
    [Ho,ho,~] = getCostMatrices(T,S,x(:,k),Q,R,P,N);
    
    %Define constraints Fu<=g(x(0))
    F_con = [F1h*S; F2h];
    g_con = [g1h-F1h*T*x(:,k); g2h];
        
    %Solve the problem
    cvx_begin
        cvx_precision high
        variable uvar(ni*N);
        minimize( 0.5*uvar'*Ho*uvar+ho'*uvar);        
        subject to
        F_con*uvar<=g_con;     
    cvx_end    
    uvar = uvar(:);
    u_apl(:,k) = uvar(1:ni);
    
    %Propagate system forward using the first input of the MPC solution
    x(:,k+1) = A*x(:,k)+B*u_apl(:,k);      
    
end

%% Plot the results
figure();
subplot(3,2,1:2)
stairs(Ts*(0:ksim),x(1,:)'*180/pi);
xlabel('t[s]');
ylabel('$\theta$,$\alpha$ [deg]','interpreter','latex');
hold on;
stairs(Ts*(0:ksim),x(2,:)'*180/pi);
legend('$\theta(t)$','$\alpha(t)$','interpreter','latex');
title('State x_1(t) and x_2(t)');
subplot(3,2,3:4);
stairs(Ts*(0:ksim),x(3,:)'*180/pi);
xlabel('t[s]');
ylabel('$\dot{\theta}$,$\dot{\alpha}$ [deg/s]','interpreter','latex');
hold on
stairs(Ts*(0:ksim),x(4,:)'*180/pi);
legend('$\dot{\theta}(t)$','$\dot{\alpha}(t)$','interpreter','latex');
title('State x_3(t) and x_4(t)');
subplot(3,2,5:6);
stairs(Ts*(0:ksim-1),u_apl','r');
title('Control input u(t)');
xlabel('t[s]');
ylabel('$u[V]$','interpreter','latex');
legend('$V_m(t)$','interpreter','latex');

%% Output MPC and disturbance rejection
C = [1 0 0 0; 0 1 0 0];
no = size(C,1);
Bdtrue = B;
Cdtrue = [0; 0];

%Assumption model
ndass = 2;
Cd = [0 1; 0 1];
Bd = [B zeros(nx,1)];

%Desired setpoint on theta
H = [1 0];
rsp = pi/4;

%Pick error dynamics
eigenplaces = [0.5; 0.4; 0.45; 0.6; 0.65; 0.55]; %working   


%Check theoretical conditions
rkObs = rank(ctrb(A',C'));
rkCond = rank([eye(nx)-A -Bd; C Cd]);   
    
%Build observer
Aobs = [A Bd; zeros(ndass,nx) eye(ndass)];
Bobs = [B; zeros(ndass,ni)]; 
Cobs = [C Cd];
Lobs = place(Aobs',Cobs',eigenplaces')'; 

%Prepare problem
model = struct('A',A,'B',B,'C',C,'Bd',Bd,'Cd',Cd,'H',H);
constraints = struct('Ax',Ax,'bx',bx,'Au',Au,'bu',bu);
u_apl = zeros(ni,ksim);
xr = zeros(nx,ksim+1);
ur = zeros(ni,ksim+1);
x = zeros(nx,ksim+1);
y = zeros(no,ksim);
zh = zeros(nx+ndass,ksim+1);
dh = zeros(ndass,ksim+1);
xh = zeros(nx,ksim+1);
v = zeros(no,ksim);

%Initial guess
d = 0.01;
xh(:,1) = zeros(4,1);
dh(:,1) = zeros(2,1);
zh(:,1) = [xh(:,1); dh(:,1)];

%Solve MPC in closed loop
for k = 1:ksim     
    
    % Solve OTS problem
    [xr(:,k),ur(:,k)] = getOTSReference(model,constraints,rsp,dh(:,k));
    
    %Obtain difference
    xt(:,k) = xh(:,k)-xr(:,k);
    
    %Recompute cost matrices in Vn(Un,x0), parametric in x0  
    [Ho,ho,~] = getCostMatrices(T,S,xt(:,k),Q,R,P,N);
    
    %Solve the problem
    cvx_begin
        variable du(ni*N);
        minimize( 0.5*du'*Ho*du+ho'*du);         
    cvx_end    
    du = du(:);
    
    %Propagate the true system forward using the first input of the MPC solution    
    u_apl(:,k) = du(1:ni)+ ur(:,k);                 
    x(:,k+1) = A*x(:,k) + B*u_apl(:,k) + Bdtrue*d;
    
    %Obtain measurement
    y(:,k) = C*x(:,k) + Cdtrue*d + v(:,k); 
    
    %Update our own belief of the state (i.e the observer)
    zh(:,k+1) = Aobs*zh(:,k)+Bobs*u_apl(:,k)+ Lobs*(y(:,k)-Cobs*zh(:,k)); 
    xh(:,k+1) = zh(1:nx,k+1);
    dh(:,k+1) = zh(nx+1:end,k+1);    
    
end

%% Plot the system dynamics
figure();
subplot(4,1,1)
stairs((0:ksim)*Ts,x(1:2,:)');
title('State x1 and x2');
xlabel('t[s]');
ylabel('x(k)');
legend('x_1(k)','x_2(k)');
subplot(4,1,2);
stairs((0:ksim)*Ts,x(3:4,:)');
title('State x3 and x4');
xlabel('t[s]');
ylabel('x(k)');
legend('x_3(k)','x_4(k)');
subplot(4,1,3);
stairs((0:ksim-1)*Ts,H*y);
hold on;
stairs((0:ksim-1)*Ts, rsp*ones(1,ksim),'--r');
title('Controlled variable r=\theta');
xlabel('t[s]');
yl = ylim;
ylim(yl+[-inf,1.05*rsp]);
ylabel('r(k)');
legend('r(k)','rsp');
subplot(4,1,4)
stairs((0:ksim-1)*Ts,u_apl(1,:),'r');
ylabel('u(k)');
title('MPC control input u(k)');
xlabel('t[s]');
legend('u_1(k)');

%Plot error dynamics
figure();
subplot(3,1,1);
stairs((0:ksim-1)*Ts,H*y-rsp);
title('Tracking Error');
xlabel('t[s]');
ylabel('\theta[rad]');
legend('r(k)');
subplot(3,1,2);
stairs((0:ksim)*Ts,dh');
title('Disturbance estimate');
xlabel('t[s]');
ylabel('d(k)');
legend('d_1(k)','d_2(k)');
subplot(3,1,3);
stairs((0:ksim)*Ts,(x-xh)');
title('State estimation error');
xlabel('t[s]');
ylabel('e_x(k)');
legend('e_1(k)','e_2(k)','e_3(k)','e_4(k)');
