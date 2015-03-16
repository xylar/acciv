clear variables;

%fileName = 'later/pass1/work/noCurvedPathScatteredVelocity.h5';
fileName = 'full/pass3/outScatteredVelocity.h5';
%fileName = 'full/pass3/work/noCurvedPathScatteredVelocity.h5';

lons = hdf5read(fileName, '/x');
lats = hdf5read(fileName, '/y');
vx = hdf5read(fileName, '/vx');
vy = hdf5read(fileName, '/vy');

figure(1);

maxPoints = 100000; %5000;
if(length(lons) > maxPoints)
    indices = ceil(rand(1,maxPoints)*length(lons));
else
    indices = 1:length(lons);
end
quiver(lons(indices), lats(indices), -vx(indices), vy(indices));
set(gca, 'xdir', 'reverse');

%fileName = 'later/pass1/work/noCurvedPathGriddedVelocity.h5';
fileName = 'full/pass3/outGridVelocity.h5';
%fileName = 'full/pass3/work/noCurvedPathGriddedVelocity.h5';

bounds = hdf5read(fileName, '/bounds');
gridVx = hdf5read(fileName, '/vx');
gridVy = hdf5read(fileName, '/vy');
xSize = size(gridVx,1);
ySize = size(gridVx,2);

[gridLons,gridLats] = ndgrid(linspace(bounds(1),bounds(2),xSize), linspace(bounds(3),bounds(4),ySize));

figure(2);
pcolor(gridLons',gridLats',gridVx'); shading flat; colorbar;
set(gca, 'xdir', 'reverse');

figure(3);
pcolor(gridLons',gridLats',gridVy'); shading flat; colorbar;
set(gca, 'xdir', 'reverse');

fileName = 'image001.h5';
imageData = hdf5read(fileName, '/data');
figure(4);
pcolor(gridLons',gridLats',imageData'); shading flat; colormap('gray');
set(gca, 'xdir', 'reverse');
% 
% fileName = 'image002.h5';
% imageData = hdf5read(fileName, '/data');
% figure(5);
% pcolor(gridLons',gridLats',imageData'); shading flat; colormap('gray');
% set(gca, 'xdir', 'reverse');
% 
% fileName = 'image003.h5';
% imageData = hdf5read(fileName, '/data');
% figure(6);
% pcolor(gridLons',gridLats',imageData'); shading flat; colormap('gray');
% set(gca, 'xdir', 'reverse');

