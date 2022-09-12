%%%%%%%%%%%%% EXAMPLE 2 %%%%%%%%%%%%%%%%
clear all

% uT = @(x,y) (1-sqrt((x+2).^2+y.^2)).*(sqrt((x+2).^2+y.^2)<1) - (1-sqrt((x-2).^2+y.^2)).*(sqrt((x-2).^2+y.^2)<1);
% T = 1;
% z0 = -1;
% z1 = 1;

uT = @(x,y) (1-sqrt((x+1).^2+y.^2)).*(sqrt((x+1).^2+y.^2)<1) + 0.5*(1-sqrt((x-1).^2+y.^2)).*(sqrt((x-1).^2+y.^2)<1);
T = 1;
z0 = -1;
z1 = 1;

%%% Hamiltonian
A = [1,0;0,1];
H = @(p1,p2) (A(1,1)*p1.^2 + 2*A(1,2)*p1.*p2 + A(2,2)*p2.^2)/2;
GradH = @(p1,p2) [A(1,1)*p1 + A(1,2)*p2,A(1,2)*p1 + A(2,2)*p2];

%%% Legendre transform of H
B = inv(A);
L = @(p1,p2) (B(1,1)*p1.^2 + 2*B(1,2)*p1.*p2 + B(2,2)*p2.^2)/2;

x0 = -4;
x1 =  4;
y0 = -2;
y1 =  2;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Unreachable target and time horizon %%%
%%%%%%%% 07_Fig21 %%%%%%%%%%%%%%%%%%%%%%%%%

nx = 300;
ny = 200;

dx = (x1-x0)/(nx-1);
X = x0:dx:x1;
dy = (y1-y0)/(ny-1);
Y = y0:dy:y1;
[xGrid,yGrid] = meshgrid(X,Y);

%%% UT %%%

UT = [];

for i = X
    vi = [];
    for j = Y
        vij = uT(i,j);
        vi = [vi;vij];
    end
    UT = [UT,vi];
end

subplot(1,2,2)
A = surf(xGrid,yGrid,UT,'LineStyle','none')
view(-12,10)
lightangle(20,20)
axis off
colormap winter
text(3,1,0.8,['$u_T$'],'interpreter','latex','FontSize', 20)
zlim([z0,z1])
saveas(A,['07_Fig21.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Projection on reachable targets %%%%%%%
%%%%%%%% 08_Vid4 %%%%%%%%%%%%%%%%%%%%%%%%%%

% nx = 150;
% ny = 100;
% 
% dx = (x1-x0)/(nx-1);
% X = x0:dx:x1;
% dy = (y1-y0)/(ny-1);
% Y = y0:dy:y1;
% [xGrid,yGrid] = meshgrid(X,Y);
% 
% %%% UT %%%
% 
% UT = [];
% 
% for i = X
%     vi = [];
%     for j = Y
%         vij = uT(i,j);
%         vi = [vi;vij];
%     end
%     UT = [UT,vi];
% end
% 
% 
% nt = 1;
% 
% for n = 1:nt
%     dt = n*T/nt;
%     
%     w = SL2dBackward(L,X,Y,dt,UT);
%     A = surf(xGrid,yGrid,w,'LineStyle','none')
%     view(-12,10)
%     lightangle(20,20)
%     axis off
%     colormap winter
%     text(2.8,1,1,['$S_{T-t}^- u_T$'],'interpreter','latex','FontSize', 20)
%     text(2.8,1,0.75,['$t=$' num2str(T-dt)],'interpreter','latex','FontSize', 20)
%     zlim([z0,z1])
%     hold off
%     saveas(A,['08_Vid42_b' num2str(n) '.png'])
% end
% U0 = w;
% for n = 1:nt
%     dt = n*T/nt;
%     
%     u = SL2d(L,X,Y,dt,U0);
%     A = surf(xGrid,yGrid,u,'LineStyle','none')
%     view(-12,10)
%     lightangle(20,20)
%     axis off
%     colormap winter
%     text(2.8,1,1,['$S_t^+(S_{T}^- u_T)$'],'interpreter','latex','FontSize', 20)
%     text(2.8,1,0.75,['$t=$' num2str(dt)],'interpreter','latex','FontSize', 20)
%     zlim([z0,z1])
%     hold off
%     saveas(A,['08_Vid42_f' num2str(n) '.png'])
% end
% 
% A = surf(xGrid,yGrid,UT,'LineStyle','none')
% view(-12,10)
% lightangle(20,20)
% axis off
% colormap winter
% text(2.8,1,1,['$S_t^+(S_{T}^- u_T)$'],'interpreter','latex','FontSize', 20)
% text(2.8,1,0.75,['$t=$' num2str(dt)],'interpreter','latex','FontSize', 20)
% text(-2.8,1,1,['$u_T$'],'interpreter','latex','FontSize', 20)
% zlim([z0,z1])
% hold off
% saveas(A,['08_Vid42_b000.png'])
% 
% A = surf(xGrid,yGrid,u,'LineStyle','none')
% view(-12,10)
% lightangle(20,20)
% axis off
% colormap winter
% text(2.8,1,1,['$S_t^+(S_{T}^- u_T)$'],'interpreter','latex','FontSize', 20)
% text(2.8,1,0.75,['$t=$' num2str(dt)],'interpreter','latex','FontSize', 20)
% text(-2.8,1,1,['$u^\ast_T$'],'interpreter','latex','FontSize', 20)
% zlim([z0,z1])
% hold off
% saveas(A,['08_Vid41_f81.png'])
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reconstruction scheme %%%%%%%%%%%%
%%%%%%%% 09_Fig22 and 09_Fig23 %%%%%%%

% nx = 300;
% ny = 200;
% 
% dx = (x1-x0)/(nx-1);
% X = x0:dx:x1;
% dy = (y1-y0)/(ny-1);
% Y = y0:dy:y1;
% [xGrid,yGrid] = meshgrid(X,Y);
% 
% %%% UT %%%
% 
% UT = [];
% 
% for i = X
%     vi = [];
%     for j = Y
%         vij = uT(i,j);
%         vi = [vi;vij];
%     end
%     UT = [UT,vi];
% end
% 
% 
% U0 = SL2dBackward(L,X,Y,T,UT);
% UTReach = SL2d(L,X,Y,T,U0);
% 
% [xi,yi] = XiIdent(GradH,dx,dy,xGrid,yGrid,UTReach,T);
% xi = xi';
% yi = yi';
% 
% subplot(1,3,1)
% A = surf(xGrid,yGrid,UTReach,'LineStyle','none')
% view(-12,10)
% lightangle(20,20)
% axis off
% colormap winter
% text(2.8,1,1,['$u^\ast_T$'],'interpreter','latex','FontSize', 20)
% zlim([z0,z1])
% hold off
% 
% subplot(1,3,2)
% surf(xGrid,yGrid,U0,'LineStyle','none')
% view(-12,10)
% lightangle(20,20)
% axis off
% colormap winter
% text(2.8,1,1,['$\tilde{u}_0$'],'interpreter','latex','FontSize', 20)
% zlim([z0,z1])
% hold off
% 
% subplot(1,3,3)
% scatter3(xi,yi,0*xi,'filled','MarkerFaceColor',[0 0.2 0.35])
% view(-12,10)
% lightangle(20,20)
% axis off
% colormap winter
% text(2.8,1,0.6,['$X_T(u_T^\ast)$'],'interpreter','latex','FontSize', 20)
% zlim([z0,z1])
% hold off
% saveas(A,['09_Fig22.png'])