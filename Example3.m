%%%%%%%%%%%%% EXAMPLE 3 %%%%%%%%%%%%%%%%

clear all
x0=-4;
x1=4;


% u0 = @(x) (1-2*abs(x+1)).*(x<=0).*(x>-1.5) + (1-2*abs(x-1)).*(x>0).*(x<1.5);
% T = 4;
% y0 = -1;
% y1 = 1;

u0 = @(x) (abs(x+1)-1).*(x<=0).*(x>-2) + (abs(x-0.5)-0.5).*(x>0).*(x<1);
T = 4;
y0 = -1;
y1 = 0;

%%% Hamiltonian
H = @(p) p.^2/2;
GradH = @(p) p;

%%% Legendre transform of H
L = @(p) p.^2/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Initial data  %%%%%%%%%%%%%%%%%
%%%%%%%% 11_Fig31 %%%%%%%%%%%%%%%%%%%%%%%%%

% nx=500;
% 
% dx = (x1-x0)/(nx-1);
% 
% X = x0:dx:x1;
% 
% subplot(1,2,2)
% A = plot(X,u0(X),'k','LineWidth',1.2);
% hold on
% ylim([y0,y1])
% legend('$u_0$','interpreter','latex','location','southeast','FontSize', 15)
% axis off
% saveas(A,['11_Fig31.png'])
% hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Evolution of the target %%%%%%%%%%%%%%%
%%%%%%%% 12_Vid5 and 13_Vid6 %%%%%%%%%%%%%%

nx=500;

dx = (x1-x0)/(nx-1);

X = x0:dx:x1;
U0 = u0(X);

nt = 50;

for n = 1:nt
    dt = n*T/nt;
    UT = Forward_1d(U0,dt,L,X);
    
    l_lims = x0;
    u_lims = [];
    tol = -1/dt;
    %%% Avoid discretization errors
    parity = 0;
    for i = 1:length(X)-2
        p1 = (UT(i+1)-UT(i))/dx;
        p2 = (UT(i+2)-UT(i+1))/dx;
        if (p2-p1)/dx<tol && parity == 0
            u_lims = [u_lims, X(i)-dt*GradH(p1)];
            p3 = (UT(i+3)-UT(i+2))/dx;
            l_lims = [l_lims, X(i+1)-dt*GradH(p3)];
            parity = 1;
        else
            parity = 0;
        end
    end
    u_lims = [u_lims, x1];
    
    subplot(1,2,1)
    A = plot(X,UT,'k','LineWidth',1.2);
    hold on
    ylim([y0,y1])
    legend('$u_T$','interpreter','latex','location','southeast','FontSize', 15)
    text(3,y1-0.2,['$T=$' num2str(dt)],'interpreter','latex','FontSize', 20)
    axis off
    hold off
    
    U0tilde = Backward_1d(UT,dt,L,X);
    subplot(1,2,2)
    plot(X,U0,':k','LineWidth',1.2);
    hold on
    for m = 1:length(l_lims)
        Xn = X(X>=l_lims(m) & X<=u_lims(m)+dx);
        U0n = U0tilde(X>=l_lims(m) & X<=u_lims(m)+dx);
        plot(Xn,U0n,'r','LineWidth',1.2)
        hold on
        if m < length(l_lims)
            Xn = X(X>=u_lims(m) & X<=l_lims(m+1)+dx);
            U0n = U0tilde(X>=u_lims(m) & X<=l_lims(m+1)+dx);
            plot(Xn,U0n,'k','LineWidth',1.2)
        end
    end
    ylim([y0,y1])
    legend('$u_0$','$\tilde{u}_0|_{X_T}$','$\tilde{u}_0$','interpreter','latex','location','southeast','FontSize', 15)
    axis off
    saveas(A,['13_Vid6_'  num2str(n) '.png'])
    hold off
end