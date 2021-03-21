function dx = manipulator(dx,x,u)
global Lr mp Lp g Br Bp Jr Jp

    H = [   Jr+mp*Lr^2+.25*mp*Lp^2*sin(x(2))^2     -.5*mp*Lr*Lp*cos(x(2));
            .5*mp*Lr*Lp*cos(x(2))                   Jp+.25*mp*Lp^2  ];

    C = [	.25*mp*Lp^2*sin(2*x(2))*x(4)+Br         .5*mp*Lr*Lp*x(4)*sin(x(2));
           -.125*mp*Lp^2*x(3)*sin(2*x(2))           Bp  ];

    G = [   0;
           -.5*mp*Lp*g*sin(x(1))];

    B = [   1;
            0];
    
    dx(1) = x(3);
    dx(2) = x(4);
    dx(3:4) = H\(-C*x(3:4) - G + B*u);


end

