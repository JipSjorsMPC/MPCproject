xN5 = load('xN5.mat','x');
xN10 = load('xN10.mat','x');
xN15 = load('xN15.mat','x');
xN20 = load('xN20.mat','x');
xN25 = load('xN25.mat','x');

Tsim = 2; %[s]
Ts = 0.05;
Nsim = round(Tsim/Ts);

%Define colors
orange = [255 128 0]/256;
purple = [102 0 204]/256;

figure('Name','Different Prediction Horizon','Position',[100 500 400 225]);
subplot(2,1,1)
%     stairs(Ts*(0:Nsim),xN5.x(1,:),'b','LineWidth',1), hold on
    stairs(Ts*(0:Nsim),xN10.x(1,:),'m','LineWidth',1), hold on
%     stairs(Ts*(0:Nsim),xN15.x(1,:),'c','LineWidth',1), hold on
%     stairs(Ts*(0:Nsim),xN20.x(1,:),'g','LineWidth',1), hold on
%     stairs(Ts*(0:Nsim),xN25.x(1,:),'Color',orange,'LineWidth',1), hold on
%     stairs(Ts*(0:Nsim),data(6).X(1,:),'r','LineWidth',1), hold on
%     xlabel('t[s]');
    ylabel('$\theta [deg]$','interpreter','latex','fontsize',12);
    legend('$N = 5$','$N = 10$','$N = 15$','$N = 20$','$N = 25$','interpreter','latex','fontsize',8);
subplot(2,1,2)
%     stairs(Ts*(0:Nsim),xN5.x(2,:),'b','LineWidth',1), hold on
    stairs(Ts*(0:Nsim),xN10.x(2,:),'m','LineWidth',1), hold on
%     stairs(Ts*(0:Nsim),xN15.x(2,:),'c','LineWidth',1), hold on
%     stairs(Ts*(0:Nsim),xN20.x(2,:),'g','LineWidth',1), hold on
%     stairs(Ts*(0:Nsim),xN25.x(2,:),'Color',orange,'LineWidth',1), hold on
%     stairs(Ts*(0:Nsim),data(6).X(2,:),'r','LineWidth',1), hold on
    xlabel('$t [s]$','interpreter','latex','fontsize',12);
    ylabel('$\alpha [deg]$','interpreter','latex','fontsize',12);
