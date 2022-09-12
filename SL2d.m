function V = SL2d(L,X,Y,dt,U0)

[xGrid,yGrid] = meshgrid(X,Y);

V = [];

for i = X
    vi = [];
    for j = Y
        vij = min(min(U0 + dt*L((i-xGrid)/dt,(j-yGrid)/dt)));
        vi = [vi;vij];
    end
    V = [V,vi];
end
