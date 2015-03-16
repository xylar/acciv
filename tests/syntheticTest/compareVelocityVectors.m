fileName = 'days_1_2/pass2/outScatteredVelocity.h5';

vAtr0 = 1.0;
r0 = 0.125;
imageScale = 1.0;
x0 = 0.5;
y0 = 0.5;

imageBounds = imageScale*[0, 1, 0, 1];
omega0 = vAtr0/r0*(sqrt(exp(1)))/(sqrt(exp(1))-1);

tiePointX = hdf5read(fileName, '/x');
tiePointY = hdf5read(fileName, '/y');
tiePointVx = hdf5read(fileName, '/dataX');
tiePointVy = hdf5read(fileName, '/dataY');

rSquared = (tiePointX-x0).^2 + (tiePointY-y0).^2;
indices = find(rSquared ~= 0);
invRSquared = zeros(size(rSquared));
invRSquared(indices) = 1./rSquared(indices);
vPhiOverR = omega0*r0^2.*invRSquared.*(1 - exp(-rSquared/(2*r0^2)));
vxExact = vPhiOverR.*(tiePointY-y0);
vyExact = -vPhiOverR.*(tiePointX-x0);


indices = find(rSquared  < (0.4*imageScale)^2);
vxDiff = abs(vxExact - tiePointVx);
vyDiff = abs(vyExact - tiePointVy);
maxVxDiff = max(vxDiff(indices))
maxVyDiff = max(vyDiff(indices))
vMax = max(max(vxExact(indices)), max(vyExact(indices)));
errorVx = maxVxDiff/vMax
errorVy = maxVyDiff/vMax

rmsError = sqrt(mean(vxDiff(indices).^2 + vyDiff(indices).^2))

rmsPercentError = rmsError/vMax

r = linspace(0, 0.4, 1024);
invR = zeros(size(r));
invR(2:end) = 1./r(2:end);
vPhi = omega0*r0^2.*invR.*(1 - exp(-r.^2/(2*r0^2)));

R = sqrt(rSquared);
invR = sqrt(invRSquared);
VPhi = (tiePointVx.*(tiePointY - y0) - tiePointVy.*(tiePointX - x0)).*invR;
VR =  (tiePointVx.*(tiePointX - x0) + tiePointVy.*(tiePointY - y0)).*invR;

scale = 1.5*max(maxVxDiff, maxVyDiff);

figure(1);
plot(r, vPhi, 'r', R(:), VPhi(:), '.');
axis([0 0.4 -scale 1.2]);

figure(2);
plot(R(:), VR(:), '.');
axis([0 0.4 -scale scale]);

figure(3);
plot(R(:), VPhi(:) - vPhiOverR(:).*R(:), '.');
axis([0 0.4 -scale scale]);
