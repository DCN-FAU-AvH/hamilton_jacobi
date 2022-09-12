function w = Backward_1d(UT,dt,L,X)
%%% UT must be a vector of the same length as X

w = [];

for i = 1:length(X)
   wi = max(UT - dt*L((X-X(i))/dt));
   w = [w, wi];
end