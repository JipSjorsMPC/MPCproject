function [P,gamma] = computeRegionOfAttraction(G,H,phi)
% computeRegionOfAttraction projects the constraints Z(x,u) onto x by
% means of eliminating u.
% Input:    Gx+Hu+phi<=0;
% Output:   Px+gamma<=0;

%Initialize algorithm
m = size(H,2);
j = m;
solutionIsNotFound = true;
Gj = [G H(:,1:m-1)];
Hj = H(:,m);
phij = phi;


while solutionIsNotFound == true

%Categorize the elements of Hj
I0 = find(Hj == 0);
Ip = find(Hj > 0);
In = find(Hj < 0);

%Compute Pj and gammaj
C = [Gj phij];
D0 = [];
D0 = C(I0,:);
% Dp = kron(Hj(Ip),C(In,:))-kron(C(Ip,:),Hj(In)); unstable for large
Dp = [];
for i = 1:length(Ip)
    Dpi = [];
    for z = 1:length(In)
        Dpi(z,:) = Hj(Ip(i))*C(In(z),:)-Hj(In(z))*C(Ip(i),:);        
    end
    Dp = [Dp ; Dpi];    
end

D = [D0; Dp];    
Pj = D(:,1:end-1);
gammaj = D(:,end);

%Remove redundant constraints
Poly = Polyhedron(Pj,-gammaj);
Poly.minHRep();
Pj = Poly.A;
gammaj = -Poly.b;


j = j-1;
 
    if j == 0
        P = Pj;
        gamma = gammaj;
        solutionIsNotFound = false;
    else        
        Gj = Pj(:,1:end-1);
        Hj = Pj(:,end);
        phij = gammaj;
    end
    
end

end
