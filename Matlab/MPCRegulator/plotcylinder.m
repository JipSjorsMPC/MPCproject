function plotcylinder(r,h,clr)
    clr
    [X,Y,Z] = cylinder(r,50);
    Z = Z*-h+ones(size(Z))*h/2;
    surf(X,Y,Z,'facecolor',clr,'EdgeColor','none','LineStyle','none');
    hold on
    fill3(X(1,:),Y(1,:),Z(1,:),clr,'EdgeColor',clr*1.01)
    fill3(X(2,:),Y(2,:),Z(2,:),clr,'EdgeColor',clr*1.01)
end