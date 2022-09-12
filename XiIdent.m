function [xi,yi] = XiIdent(GradH,dx,dy,xGrid,yGrid,UT,dt)

%%% Hamiltonian and alpha0 s.t. H'(alpha0)=0
%H = @(p1,p2) (p1.^2 + p2.^2)/2;

%GradH = @(p1,p2) [p1,p2];

%%% Legendre transform of H
%L = @(p1,p2) (p1.^2 + p2.^2)/2;

[px,py] = gradient(UT,dx,dy);
Lapla = divergence(xGrid,yGrid,px,py);

xi=[];
yi=[];

for i = 1:length(xGrid(1,:))
    for j = 1:length(yGrid(:,1))
        if Lapla(j,i)>-1/dt
            X = [xGrid(1,i),yGrid(j,1)]-dt*GradH(px(j,i),py(j,i));
            xi = [xi,X(1)];
            yi = [yi,X(2)];
        end
    end
end