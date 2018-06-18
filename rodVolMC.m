function [ vol,mass,eqvDiam,AR ] = rodVolMC( nsphere,overlap,sphDiam,rho, nSim )
%DUMBELLVOL Calculate the volume and mass of a dumbell shape

d=sphDiam/2;
centres = zeros(3,nsphere);
centres(:,1) = [d sphDiam/2 sphDiam/2];
for i=2:nsphere
    d = d + (1-overlap)*sphDiam;
    centres(:,i) = [d sphDiam/2 sphDiam/2];
end

AR = (d+sphDiam/2)/sphDiam;

%define a bounding box, assume the rod's major axis in +x

xlow = 0;
xhigh = (d+sphDiam/2);
Lx = xhigh - xlow;

ylow = 0;
yhigh = sphDiam;
Ly = yhigh - ylow;

zlow = 0;
zhigh = sphDiam;
Lz = zhigh - zlow;

sampl = zeros(3,nSim);

sampl(1,:) = xlow + Lx*rand(nSim, 1);
sampl(2,:) = ylow + Ly*rand(nSim, 1);
sampl(3,:) = zlow + Lz*rand(nSim, 1);

acceptInds = zeros(nSim,1,'logical');

for i = 1:nsphere
    acceptInds = acceptInds | (sqrt( sum(bsxfun(@minus, sampl, centres(:,i)).^2) ) < sphDiam/2)';
end

vol = sum( acceptInds )/nSim*Lx*Ly*Lz;

mass=vol*rho;

eqvDiam = (6*vol/pi)^(1/3);

end

