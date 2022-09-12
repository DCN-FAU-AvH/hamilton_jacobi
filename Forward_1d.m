function u= Forward_1d(U0,dt,L,X)
%%% U0 must be a vector of the same length as X

u = [];

for i = 1:length(X)
   [ui,j] = min(U0 + dt*L((X(i)-X)/dt));
   u = [u, ui];
end
