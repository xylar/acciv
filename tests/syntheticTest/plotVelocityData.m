clear variables;

fileName = 'days_1_2/pass2/outScatteredVelocity.h5';
%fileName = 'day1/pass1/work/noCurvedPathScatteredVelocity.h5';

maxPoints = 5000;

lons = hdf5read(fileName, '/x');
lats = hdf5read(fileName, '/y');
vx = hdf5read(fileName, '/dataX');
vy = hdf5read(fileName, '/dataY');
indexFractions = rand(1,maxPoints);
if(length(lons) > maxPoints)
    indices = ceil(indexFractions*length(lons));
else
    indices = 1:length(x);
end
lons = lons(indices);
lats = lats(indices);
vx = vx(indices);
vy = vy(indices);


figure(1);

quiver(lons, lats, vx, vy);
axis equal; axis tight;

fileName = 'days_1_2/pass2/outGridVelocity.h5';
%fileName = 'day1/pass1/work/noCurvedPathGriddedVelocity.h5';

bounds = hdf5read(fileName, '/bounds');
gridVx = hdf5read(fileName, '/dataX');
gridVy = hdf5read(fileName, '/dataY');
xSize = size(gridVx,1);
ySize = size(gridVx,2);

[gridLons,gridLats] = ndgrid(linspace(bounds(1),bounds(2),xSize), linspace(bounds(3),bounds(4),ySize));

skip = 1;
figure(4);
pcolor(gridLons(1:skip:end,1:skip:end)',gridLats(1:skip:end,1:skip:end)',gridVx(1:skip:end,1:skip:end)'); shading flat; colorbar;
axis equal; axis tight;

figure(5);
pcolor(gridLons(1:skip:end,1:skip:end)',gridLats(1:skip:end,1:skip:end)',gridVy(1:skip:end,1:skip:end)'); shading flat; colorbar;
axis equal; axis tight;

