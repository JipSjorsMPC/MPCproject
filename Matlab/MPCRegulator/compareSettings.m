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
% Q = blkdiag(10,100,1,5); %standard working
% R = 0.01;
% N = 20;
Q = blkdiag(1,5,1,1);
R = 1;
Nwaar = [2 10 20 30];
simulatiedata=cell(1,4);
for z=1:4
    N = Nwaar(z);
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
%     plot(XN.projection(1:3),'color','g',Xf.projection(1:3),'color','r');
%     legend('$$X_N$$','$$X_f$$','Interpreter','latex');
%     xlabel('\theta');
%     ylabel('\alpha');
%     zlabel('$\dot{\theta}$','interpreter','latex');

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
    simulatiedata{z}.x= x;
end

%% Plot and compare with nonlinear model
figure();
stairs((0:ksim)*Ts,simulatiedata{1}.x(1,:)'*180/pi,'c');
hold on;
stairs((0:ksim)*Ts,simulatiedata{2}.x(1,:)'*180/pi,'m');
hold on;
stairs((0:ksim)*Ts,simulatiedata{3}.x(1,:)'*180/pi,'b');
hold on;
stairs((0:ksim)*Ts,simulatiedata{3}.x(1,:)'*180/pi,'b');
xlabel('t[s]');
ylabel('\theta');
title('Different weights on Q');
legend('N = 5','N=10','N=20','N=30');

% figure();
% stairs((0:Nsim-1)*Ts,simulatiedata{1}.U','c');
% hold on;
% stairs((0:Nsim-1)*Ts,simulatiedata{2}.U','m');
% hold on;
% stairs((0:Nsim-1)*Ts,simulatiedata{3}.U','b');
% xlabel('t[s]');
% ylabel('U[V]');
% title('Different weights on R');
% legend('R= 0.01','R=1','R=10','R=100');




