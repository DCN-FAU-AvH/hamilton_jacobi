function V = SL2dBackward(L,X,Y,dt,UT)

[xGrid,yGrid] = meshgrid(X,Y);

V = [];

for i = X
    vi = [];
    for j = Y
        vij = max(max(UT - dt*L((xGrid-i)/dt,(yGrid-j)/dt)));
        vi = [vi;vij];
    end
    V = [V,vi];
end