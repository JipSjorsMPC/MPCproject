function rp = simulate(x,Ts)
    
    global Lr Lp
    n = length(x);

%     x(1:2,:) = x(1:2,:).*180/pi;
    x(1:2,:) = [linspace(0,2*pi,n);linspace(0,2*pi,n)]
    rp = [Lr*cos(x(1,:))+0.5*Lp*sin(x(2,:)).*sin(x(1,:));
          Lr*cos(x(1,:))-0.5*Lp*sin(x(2,:)).*cos(x(1,:));
          0.5*Lp*cos(x(2,:))];
    
    link = [Lr*cos(x(1,:));
            Lr*sin(x(1,:));
            zeros(1,length(x))];
        
    tip = [Lr*cos(x(1,:))+Lp*sin(x(2,:)).*sin(x(1,:));
           Lr*cos(x(1,:))-Lp*sin(x(2,:)).*cos(x(1,:));
           Lp*cos(x(2,:))];
      
    figure(5);
    for i = 1:n
        plot3([0 link(1,i)],[0 link(2,i)],[0 link(3,i)],'r','linewidth',3),hold on
%         plot3([link(1,i) tip(1,i)],[link(2,i) tip(2,i)],[link(3,i) tip(3,i)],'r','linewidth',3);
%         plot3([link(1,i) rp(1,i)],[link(2,i) rp(2,i)],[link(3,i) rp(3,i)],'r','linewidth',3);
        pause(Ts), hold off
        plot3(rp(1,i),rp(2,i),rp(3,i),'or','MarkerSize',3,'MarkerFaceColor','c'), hold on;
        axis([-Lr Lr -Lr Lr -Lp Lp]);
        grid on
    end
    
%     plot3(rp(1,:),rp(2,:),rp(3,:))
%     
%     n = length(x);
%     XY = 10 * rand(2,n) - 5;
%     figure()
%     for i=1:n
%         plot3(link(1,:),link(2,:),link(3,:),'or','MarkerSize',5,'MarkerFaceColor','r')
%         axis(2*[-Lp Lp -Lp Lp -Lp Lp])
%         pause(.3)
%     end

return