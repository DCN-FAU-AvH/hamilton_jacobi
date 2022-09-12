%%%%%%%%%%%%% EXAMPLE 1 %%%%%%%%%%%%%%%%

clear all
x0=-4;
x1=4;


% uT = @(x) (1-2*abs(x+1)).*(x<=0).*(x>-1.5) + (1-2*abs(x-1)).*(x>0).*(x<1.5);
% T = 0.3;
% y0 = -1;
% y1 = 1;

uT = @(x) (abs(x+1)-1).*(x<=0).*(x>-2) + (abs(x-1)-1).*(x>0).*(x<2);
T = 0.8;
y0 = -1;
y1 = 0.25;

%%% Hamiltonian
H = @(p) p.^2/2;
GradH = @(p) p;

%%% Legendre transform of H
L = @(p) p.^2/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Unreachable target and time horizon %%%
%%%%%%%% 01_Fig11 %%%%%%%%%%%%%%%%%%%%%%%%%



% subplot(1,2,2)
% A = plot(X,uT(X),'k','LineWidth',1.2)
% hold on
% ylim([y0,y1])
% legend('$u_T$','interpreter','latex','location','southeast','FontSize', 15)
% axis off
% saveas(A,['01_Fig11.png'])
% hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Projection on reachable targets %%
%%%%%%%% 02_Vid11 %%%%%%%%%%%%%%%%%%%%

% nx=500;
% 
% dx = (x1-x0)/(nx-1);
% 
% X = x0:dx:x1;
% 
% nt = 100;
% UT = uT(X);

%%% Backward
% for n = 1:nt
%     dt = n*T/nt;
%     w = Backward_1d(UT,dt,L,X);
%     
%     A = plot(X,w,'k','LineWidth',1.2)
%     hold on
%     plot(X,UT,':k','LineWidth',1.2)
%     ylim([y0,y1])
%     legend('$S_{T-t}^- u_T$','$u_T$','interpreter','latex','location','southeast','FontSize', 15)
%     text(3,1,['$t=$' num2str(T-dt)],'interpreter','latex','FontSize', 20)
%     axis off
%     saveas(A,['02_Vid11_b' num2str(n) '.png'])
%     hold off
% end

%%% Forward

% U0 = w;
% 
% for n = 1:nt
%     dt = n*T/nt;
%     u = Forward_1d(U0,dt,L,X);
%     
%     A = plot(X,u,'k','LineWidth',1.2)
%     hold on
%     plot(X,UT,':k','LineWidth',1.2)
%     ylim([y0,y1])
%     legend('$S_t^+(S_T^- u_T)$','$u_T$','interpreter','latex','location','southeast','FontSize', 15)
%     text(3,1,['$t=$' num2str(dt)],'interpreter','latex','FontSize', 20)
%     axis off
%     saveas(A,['02_Vid11_f' num2str(n) '.png'])
%     hold off
% end
% 
% A = plot(X,UT,'k','LineWidth',1.2)
% hold on
% plot(X,UT,':k','LineWidth',1.2)
% ylim([y0,y1])
% legend('$S_{T-t}^- u_T$','$u_T$','interpreter','latex','location','southeast','FontSize', 15)
% text(3,1,['$t=$' num2str(T)],'interpreter','latex','FontSize', 20)
% text(-4,1,['$u_T$'],'interpreter','latex','FontSize', 20)
% axis off
% saveas(A,['02_Vid11_b000.png'])
% hold off
% 
% 
% A = plot(X,u,'k','LineWidth',1.2)
% hold on
% plot(X,UT,':k','LineWidth',1.2)
% ylim([y0,y1])
% legend('$S_T^+(S_T^- u_T)$','$u_T$','interpreter','latex','location','southeast','FontSize', 15)
% text(3,1,['$t=$' num2str(T)],'interpreter','latex','FontSize', 20)
% text(-4,1,['$u^\ast_T = S_T^+(S_T^- u_T)$'],'interpreter','latex','FontSize', 20)
% axis off
% saveas(A,['02_Vid11_f101.png'])
% hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reconstruction scheme %%%%%%%%%%%%
%%%%%%%% 03_Fig12 %%%%%%%%%%%%%%%%%%%%

nx=500;

dx = (x1-x0)/(nx-1);

X = x0:dx:x1;

UT = uT(X);
U0 = Backward_1d(UT,T,L,X);
UTReach = Forward_1d(U0,T,L,X);

%%% Limits of interval in X_T

l_lims = x0;
u_lims = [];

tol = -10;

%%% Avoid discretization errors
parity = 0;

for i = 1:length(X)-2
    p1 = (UTReach(i+1)-UTReach(i))/dx;
    p2 = (UTReach(i+2)-UTReach(i+1))/dx;
    if (p2-p1)/dx<tol && parity == 0
        u_lims = [u_lims, X(i)-T*GradH(p1)];
        p3 = (UTReach(i+3)-UTReach(i+2))/dx;
        l_lims = [l_lims, X(i+1)-T*GradH(p3)];
        parity = 1;
    else
        parity = 0;
    end
end
u_lims = [u_lims, x1];

%%% Reconstruction Scheme Plots

% subplot(1,2,1)
% A = plot(X,UTReach,'k','LineWidth',1.2);
% hold on
% ylim([y0,y1])
% legend('$u^\ast_T$','interpreter','latex','location','southeast','FontSize', 15)
% axis off
% hold off
% 
% subplot(1,2,2)
% for n = 1:length(l_lims)
%     Xn = X(X>=l_lims(n) & X<=u_lims(n)+dx);
%     U0n = U0(X>=l_lims(n) & X<=u_lims(n)+dx);
%     plot(Xn,U0n,'r','LineWidth',1.2)
%     hold on
%     if n < length(l_lims)
%         Xn = X(X>=u_lims(n) & X<=l_lims(n+1)+dx);
%         U0n = U0(X>=u_lims(n) & X<=l_lims(n+1)+dx);
%         plot(Xn,U0n,'k','LineWidth',1.2)
%     end
% end
% ylim([y0,y1])
% legend('$\tilde{u}_0|_{X_T}$','$\tilde{u}_0$','interpreter','latex','location','southeast','FontSize', 15)
% axis off
% %saveas(A,['04_Fig13.png'])
% hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Examples of Initial Data %%%%%%%%%
%%%%%%%% 04_Vid2 %%%%%%%%%%%%%%%%%%%%

pt = x0;
fpt = 0;
if length(l_lims)>1
    for n = 1:length(l_lims)-1
        k = datasample(8:10,1);
        Xn = X(X>=u_lims(n)+dx & X<=l_lims(n+1));
        pt = [pt,u_lims(n), sort(datasample(Xn,k,'Replace',false)), l_lims(n+1)+dx];
        fpt = [fpt,0,0.1+0.4*rand(1,k),0];
    end
end
pt = [pt, x1];
fpt = [fpt, 0];

U0Ex = U0 + interp1q(pt',fpt',X')';

A = plot(X,U0Ex,'k','LineWidth',1.2);
hold on
plot(X,U0,':k','LineWidth',1.2)
ylim([y0,y1+0.5])
legend('$u_0$','$\tilde{u}_0$','interpreter','latex','location','southeast','FontSize', 15)
text(3,0.2,['$t=0$'],'interpreter','latex','FontSize', 20)
axis off
saveas(A,['06_Vid33_000.png'])
hold off
%%
nt = 100;

for n = 1:nt
    dt = n*T/nt;
    u = Forward_1d(U0Ex,dt,L,X);
    
    A = plot(X,u,'k','LineWidth',1.2);
    hold on
    uast = Forward_1d(U0,dt,L,X);
    plot(X,uast,':k','LineWidth',1.2)
    ylim([y0,y1+0.5])
    legend('$S_t^+u_0$','$S_t^+\tilde{u}_0$','interpreter','latex','location','southeast','FontSize', 15)
    text(3,0.2,['$t=$' num2str(dt)],'interpreter','latex','FontSize', 20)
    axis off
    saveas(A,['06_Vid33_' num2str(n) '.png'])
    hold off
end
