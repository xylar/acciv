fileName = 'days_1_2/pass2/outGridVelocity.h5';
vAtr0 = 1.0;
r0 = 0.125;
imageScale = 1.0;
x0 = 0.5;
y0 = 0.5;

imageBounds = imageScale*[0, 1, 0, 1];
omega0 = vAtr0/r0*(sqrt(exp(1)))/(sqrt(exp(1))-1);

bounds = hdf5read(fileName, '/bounds');
vx = hdf5read(fileName, '/dataX');
vy = hdf5read(fileName, '/dataY');
xSize = size(vx,1);
ySize = size(vx,2);

[x,y] = ndgrid(linspace(bounds(1),bounds(2),xSize), linspace(bounds(3),bounds(4),ySize));

rSquared = (x-x0).^2 + (y-y0).^2;
indices = find(rSquared ~= 0);
invRSquared = zeros(size(rSquared));
invRSquared(indices) = 1./rSquared(indices);
vPhiOverR = omega0*r0^2.*invRSquared.*(1 - exp(-rSquared/(2*r0^2)));
vxExact = vPhiOverR.*(y-y0);
vyExact = -vPhiOverR.*(x-x0);

indices = find(rSquared  < (0.4*imageScale)^2);
vxDiff = abs(vxExact - vx);
vyDiff = abs(vyExact - vy);
maxVxDiff = max(vxDiff(indices))
maxVyDiff = max(vyDiff(indices))
errorVx = maxVxDiff/max(vxExact(indices))
errorVy = maxVyDiff/max(vyExact(indices))

rmsError = sqrt(mean(vxDiff(indices).^2 + vyDiff(indices).^2))

figure(1);
pcolor(x,y,vxDiff);
shading flat;
axis equal;
caxis([0, maxVxDiff]);

figure(2);
pcolor(x,y,vyDiff);
shading flat;
axis equal;
caxis([0, maxVyDiff]);

r = linspace(0, 0.4, 1024);
invR = zeros(size(r));
invR(2:end) = 1./r(2:end);
vPhi = omega0*r0^2.*invR.*(1 - exp(-r.^2/(2*r0^2)));

R = sqrt(rSquared);
invR = sqrt(invRSquared);
VPhi = (vx.*(y - y0) - vy.*(x - x0)).*invR;
VR =  (vx.*(x - x0) + vy.*(y - y0)).*invR;

figure(3);
plot(r, vPhi, 'r', R(:), VPhi(:), '.');

figure(4);
plot(R(:), VR(:), '.');


figure(5);
plot(R(:), VPhi(:) - vPhiOverR(:).*R(:), '.');