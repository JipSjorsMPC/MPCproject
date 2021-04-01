function rp = simulate(x,Ts)
    
    global Lr Lp
    n = length(x);

    alpha = x(1,:);%.*180/pi;
    theta = x(2,:);%.*180/pi;
        
    link = [Lr*cos(theta);
            -Lr*sin(theta);
            zeros(size(theta))];
        
    tip = [Lr*cos(theta)+Lp*sin(alpha).*sin(theta);
           -Lr*sin(theta)-Lp*sin(alpha).*cos(theta);
           Lp*cos(alpha)];
    
    figure(5);
    rotary_arm = plot3(0,0,0,'b','linewidth',3); hold on 
    pendulum_link = plot3(0,0,0,'b','linewidth',6); hold on
    traj = animatedline('Color','r');    
    axis([-Lp Lp -Lp Lp -Lp Lp]);
    xlabel('x'), ylabel('y'), zlabel('z');
    hold on, grid on; 
    view(-75,30)
    for i= 1:n
        set(rotary_arm,'XData',[0 link(1,i)],'Ydata',[0 link(2,i)],'Zdata',[0 link(3,i)]);
        set(pendulum_link,'XData',[link(1,i) tip(1,i)],'Ydata',[link(2,i) tip(2,i)],'Zdata',[link(3,i) tip(3,i)]);
        addpoints(traj,tip(1,i),tip(2,i),tip(3,i));
        drawnow;
        pause(Ts); %pause execution for to make animation visible    
    end
    hold off

return