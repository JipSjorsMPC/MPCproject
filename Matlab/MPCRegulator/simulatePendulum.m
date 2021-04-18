function rp = simulatePendulum(x,Ts)
    
    global Lr Lp
    n = length(x);

    theta = x(1,:);%.*180/pi;
    alpha = x(2,:);%.*180/pi;

    link = [Lr*cos(theta);
            Lr*sin(theta);
            zeros(size(theta))]; %[x,y,z]
        
    tip = [Lr*cos(theta)+Lp*sin(alpha).*sin(theta);
           Lr*sin(theta)-Lp*sin(alpha).*cos(theta);
           Lp*cos(alpha)];  %[x,y,z]       
          
    figure('Name','Simulation MPC with swing-up','Position',[250 500 400 350]);
    plotcube([.1 .1 (Lp-.005)],[-.05 -.05 -(Lp+.005)],1,[32 32 32]/256); hold on %EDGES,ORIGIN,ALPHA,COLOR
    plotcylinder(0.012,0.01,[170 169 173]/256);
    pin = plot3(0,0,0,'Color',[32 32 32]/256,'linewidth',10); hold on
    rotary_arm = plot3(0,0,0,'Color',[170 169 173]/256,'linewidth',4); hold on
    pendulum_link = plot3(0,0,0,'Color',[135 22 20]/256,'linewidth',8); hold on
    traj = animatedline('Color','c');    
    set(pin,'XData',[0 0],'Ydata',[0 0],'Zdata',[-Lp 0]);
    axis([-Lp*1.2 Lp*1.2 -Lp*1.2 Lp*1.2 -Lp*1.2 Lp*1.2]);
    xlabel('x'), ylabel('y'), zlabel('z');
    hold on, grid on; 
    view(75,3)
    for i= 1:n
        set(rotary_arm,'XData',[0 link(1,i)],'Ydata',[0 link(2,i)],'Zdata',[0 link(3,i)]);
        set(pendulum_link,'XData',[link(1,i) tip(1,i)],'Ydata',[link(2,i) tip(2,i)],'Zdata',[link(3,i) tip(3,i)]);
        addpoints(traj,tip(1,i),tip(2,i),tip(3,i));
        drawnow;
        pause(Ts); %pause execution for to make animation visible    
    end
    hold off

return