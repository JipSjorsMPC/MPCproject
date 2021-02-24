%Inverted rotary pendulum parameters

mr = 0.095;   
Lr = 0.085; 
mp = 0.024;  
Lp = 0.129; 
mh = 0.016;  
rh = 0.0111;
Jh = 0.6*10^-6;
Rm = 8.4;
kt = 0.042;
km = 0.042;
Jm = 4.0*10^-6;
Lm = 1.16*10^-3; 
g = 9.81;

Br = 8.0508*10^-4;
Bp = 1.4*10^-5;
Jr = 2.2923*10^-4;
Jp = 1.2551*10^-4; 

%Constraints
Vub = 10;
Vlb = -10;

%Initial conditions
x0 = [0 pi 0 0]; 
