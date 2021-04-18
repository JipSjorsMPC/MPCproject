u001 = load('u001.mat','u_apl');
u1 = load('u1.mat','u_apl');
u10 = load('u10.mat','u_apl');
xa01 = load('xa01.mat','x');
xa10 = load('xa10.mat','x');
xa100 = load('xa100.mat','x');
xth01 = load('xth01.mat','x');
xth10 = load('xth10.mat','x');
xth100 = load('xth100.mat','x');
x0 = load('x0.mat','x');

Tsim = 2; %[s]
Ts = 0.05;
ksim = round(Tsim/Ts);

figure('Name','Different weights on Q and R','Position',[100 500 400 225]);
subplot(2,1,1)
    stairs(Ts*(0:ksim),xth01.x(1,:)'*180/pi,'b','LineWidth',1), hold on
    stairs(Ts*(0:ksim),xth10.x(1,:)'*180/pi,'m','LineWidth',1), hold on
    stairs(Ts*(0:ksim),xth100.x(1,:)'*180/pi,'c','LineWidth',1), hold on
%     xlabel('t[s]','interpreter','latex','fontsize',12);
    ylabel('$\theta$ [deg]','interpreter','latex','fontsize',12);
    legend('$Q_{\theta}$ = 0.1','$Q_{\theta}$ = 10','$Q_{\theta}$ = 100','interpreter','latex');
%     title('Different weights on Q_{\theta}'); 
% subplot(3,1,2)
%     stairs(Ts*(0:ksim),xa01.x(1,:)'*180/pi,'b'), hold on
%     stairs(Ts*(0:ksim),xa10.x(1,:)'*180/pi,'m'), hold on
%     stairs(Ts*(0:ksim),xa100.x(1,:)'*180/pi,'c'), hold on
%     stairs(Ts*(0:ksim),x0.x(1,:)'*180/pi,'r'), hold on
%     xlabel('t[s]');
%     ylabel('$\theta$ [deg]','interpreter','latex');
%     legend('$Q_{\alpha}$ = 0.1','$Q_{\alpha}$ = 10','$Q_{\alpha}$ = 100','interpreter','latex');
%     title('Different weights on Q_{\alpha}'); 
subplot(2,1,2);
    stairs(Ts*(0:ksim-1),u001.u_apl','b','LineWidth',1), hold on
    stairs(Ts*(0:ksim-1),u1.u_apl','m','LineWidth',1), hold on
    stairs(Ts*(0:ksim-1),u10.u_apl','c','LineWidth',1), hold on
%     title('Different weights on R');
    xlabel('$t [s]$','interpreter','latex','fontsize',12);
    ylabel('$u [V]$','interpreter','latex','fontsize',12);
    legend('$R = 0.01$','$R = 1$','$R = 10$','interpreter','latex','fontsize',8);
hold on