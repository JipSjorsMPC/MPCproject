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

global u
u = 0;

%% Linearize continous system

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

%Indeed checked and verified this gives same result
%% Discretizatize linearized system

%Discretize system using sample time Ts
Ts = 0.05;
sysd = c2d(ss(Ac,Bc,Cc,Dc),Ts,'zoh');
[A,B,C,D] = ssdata(sysd);

%% Solve with methods from Class
%Parameters
Q = blkdiag(1,100,1,100);
R = 1;
N = 3;

%Get unconstrained infinite control Riccati cost and control
[P,~,G] = dare(A,B,Q,R);
K = -G;

Ak = A+B*K;
Qk = Q+K'*R*K;

%% Check if Linearized system is local CLF
%(A,B) stable?
if rank(ctrb(A,B))==length(A) disp('(A,B) is stabilizable'), else disp('(A,B) is <strong>not</strong> stabilizable'), end

%P positive definite?
if eig(P)>0 disp('P is pos. definite'), else disp('P is <strong>not</strong> pos. definite'), end 

%Ak stable?
if (trace(Ak)< 0) && (det(Ak)>0) disp('Ak is Hurwitz stable'), else disp('Ak is <strong>not</strong> Hurwitz stable'), end

%%

%Constraints (see further below)
% |x1(k)|<= pi/4    for all k
% |u(k)|<= 5        for all k
% |x2(k)| <= pi/8  for all k
x1_b = pi;
x2_b = pi;
u_b = 20;

%Compute the terminal set Xf
%Build constraints
Ff = vertcat([1 0 0 0; -1 0 0 0],[0 1 0 0; 0 -1 0 0],[K ;-K]);
gf = vertcat([x1_b; x1_b],[x2_b; x2_b],[u_b; u_b]);

%Compute Xf under u =Kx such that x satisfies input and state constraints
[Hf,hf,kstar] = computeMaxConstraintAdmissibleSet(Ak,Ff,gf);
Xf = Polyhedron(Hf,hf);
Xf.minHRep();

%%
%Project the 4 dimensional Xf onto 2 dimensional plane by removing the last
%2 columns
Hfmod = Hf;
hfmod = hf;
for i =1:2
    [Hfmod,hfmod] = fouriermotzkin(Hfmod,hfmod);
end
Xf2D = Polyhedron(Hfmod,hfmod);
Xf2D.minHRep();

%Get Prediction matrices
[T,S] = getStatePredictionMatrices(A,B,N);

%First build G,H and phi matrix in Gx+Hu+phi<=0
F1h = kron(eye(N+1),[eye(2) zeros(2); -eye(2) zeros(2)]);
g1h = repmat([x1_b; x2_b; x1_b; x2_b],N+1,1);
F2h = kron(eye(N),[1;-1]);
g2h = repmat(u_b*[1;1],N,1);
G = [F1h*T;...
 zeros(size(F2h,1),nx);...
 Hf*T(end-nx+1:end,:)];
H = [F1h*S; F2h; Hf*S(end-nx+1:end,:)];
phi = [-g1h; -g2h; -hf];

% Compute Xn, i.e. the region of attraction
[PXn,gammaXn] = computeRegionOfAttraction(G,H,phi);
AXn = PXn;
bXn = -gammaXn;
Xn = Polyhedron(AXn,bXn);
Xn.minHRep();

%Project Xn onto 2D by removing the last two columns
AXn_mod = Xn.A;
bXn_mod = Xn.b;
for i =1:2
    [AXn_mod,bXn_mod] = fouriermotzkin(AXn_mod,bXn_mod);
end
Xn2D = Polyhedron(AXn_mod,bXn_mod);
Xn2D.minHRep();

%Plot Xn and Xf
plot(Xn2D,'color','g',Xf2D,'color','r');
legend('$$X_N$$','$$X_f$$','Interpreter','latex');