%%%%%%%%%%%%% EXAMPLE 4 %%%%%%%%%%%%%%%%
clear all

u0 = @(x,y) (1-sqrt((x+1).^2+y.^2)).*(sqrt((x+1).^2+y.^2)<1) - (1-sqrt((x-1).^2+y.^2)).*(sqrt((x-1).^2+y.^2)<1);
T = 4;
z0 = -1;
z1 = 1;

% u0 = @(x,y) (1-sqrt((x+1).^2+y.^2)).*(sqrt((x+1).^2+y.^2)<1) + 0.5*(1-sqrt((x-1).^2+y.^2)).*(sqrt((x-1).^2+y.^2)<1);
% T = 4;
% z0 = -1;
% z1 = 1;

%%% Hamiltonian
A = [1,0;0,1];
H = @(p1,p2) (A(1,1)*p1.^2 + 2*A(1,2)*p1.*p2 + A(2,2)*p2.^2)/2;
GradH = @(p1,p2) [A(1,1)*p1 + A(1,2)*p2,A(1,2)*p1 + A(2,2)*p2];

%%% Legendre transform of H
B = inv(A);
L = @(p1,p2) (B(1,1)*p1.^2 + 2*B(1,2)*p1.*p2 + B(2,2)*p2.^2)/2;

x0 = -4;
x1 =  4;
y0 = -3;
y1 =  3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Initial data %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 15_Vid7 and 16_Vid8 %%%%%%%%%%%%%%%%%%%

nx = 300;
ny = 200;

dx = (x1-x0)/(nx-1);
X = x0:dx:x1;
dy = (y1-y0)/(ny-1);
Y = y0:dy:y1;
[xGrid,yGrid] = meshgrid(X,Y);

%%% U0 %%%

U0 = [];

for i = X
    vi = [];
    for j = Y
        vij = u0(i,j);
        vi = [vi;vij];
    end
    U0 = [U0,vi];
end

subplot(1,3,1)
A = surf(xGrid,yGrid,U0,'LineStyle','none')
view(-12,10)
lightangle(20,20)
axis off
colormap winter
text(2.8,1,1,['$u_0$'],'interpreter','latex','FontSize', 20)
zlim([z0,z1])
saveas(A,['14_Fig41.png'])
hold off

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Evolution of Reconstruction scheme %%%%%%%%%
%%%%%%%% 15_Vid7 and 16_Vid8 %%%%%%%%%%%%%%%%%%%

% nt = 100;
% 
% nx = 200;
% ny = 150;
% 
% dx = (x1-x0)/(nx-1);
% X = x0:dx:x1;
% dy = (y1-y0)/(ny-1);
% Y = y0:dy:y1;
% [xGrid,yGrid] = meshgrid(X,Y);
% 
% %%% U0 %%%
% 
% U0 = [];
% 
% for i = X
%     vi = [];
%     for j = Y
%         vij = u0(i,j);
%         vi = [vi;vij];
%     end
%     U0 = [U0,vi];
% end
% 
% for n = 1:nt
%     dt = n*T/nt;
%     tic
%     UT = SL2d(L,X,Y,dt,U0);
%     U0tilde = SL2dBackward(L,X,Y,dt,UT);
%     toc
% 
%     tic
%     [xi,yi] = XiIdent(GradH,dx,dy,xGrid,yGrid,UT,dt);
%     toc
%     xi = xi';
%     yi = yi'; 
%     
%     subplot(1,3,1)
%     A = surf(xGrid,yGrid,UT,'LineStyle','none')
%     view(-12,10)
%     lightangle(20,20)
%     axis off
%     colormap winter
%     text(2.8,1,1,['$u_T$'],'interpreter','latex','FontSize', 20)
%     text(2.8,1,0.6,['$T=$' num2str(dt)],'interpreter','latex','FontSize', 15)
%     zlim([z0,z1])
%     hold off
%     
%     subplot(1,3,2)
%     surf(xGrid,yGrid,U0tilde,'LineStyle','none')
%     view(-12,10)
%     lightangle(20,20)
%     axis off
%     colormap winter
%     text(2.8,1,1,['$\tilde{u}_0$'],'interpreter','latex','FontSize', 20)
%     zlim([z0,z1])
%     hold off
%     
%     subplot(1,3,3)
%     %scatter3(xi,yi,0*xi,'filled','MarkerFaceColor',[0 0.2 0.35])
%     plot(xi,yi,'.','MarkerFace',[0 0.1 0.3])
%     %view(-12,10)
%     %lightangle(20,20)
%     axis off
%     %colormap winter
%     text(-3.8,-2.7,['$X_T(u_T)$'],'interpreter','latex','FontSize', 18)
%     zlim([z0,z1])
%     hold off
%     %saveas(A,['17_Vid9_' num2str(n) '.png'])
% end