function [vx, vy] = latLonPerSToV(lats, lonPerS, latPerS)

Re = 7.1492e7; %71,492 km
Rp = 6.6854e7; %66,854 km
epsilon = Rp/Re;

% convert to planetographic latitude
lats = 180/pi*atan(Re^2/Rp^2*tan(lats*pi/180));

t = atan(epsilon*tan(lats*pi/180));
dtdLatitude = pi/180*epsilon*(sec(lats*pi/180)).^2./(1 + epsilon^2*(tan(lats*pi/180)).^2);
dydLatitude = sqrt(Re^2*sin(t).^2 + Rp^2*cos(t).^2).*dtdLatitude;

dxdLongitude = -Re*cos(t)*pi/180;

vx = lonPerS.*dxdLongitude;
vy = latPerS.*dydLatitude;