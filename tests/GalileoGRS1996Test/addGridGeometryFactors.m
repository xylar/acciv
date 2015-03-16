fileName = 'gridGeometryFactors.h5';
bounds = hdf5read(fileName, '/imageBounds');
imageSize = hdf5read(fileName, '/imageSize');
nx = double(imageSize(1));
ny = double(imageSize(2));

[gridLons,gridLats] = ndgrid(linspace(bounds(1),bounds(2),nx), linspace(bounds(3),bounds(4),ny));

vxFactor = ones(nx,ny);
vyFactor = ones(nx,ny);

[vxFactor, vyFactor] = latLonPerSToV(gridLats, vxFactor, vyFactor);

hdf5write(fileName, '/bounds', bounds, 'WriteMode', 'overwrite');
%hdf5write(fileName, '/imageSize', imageSize, 'WriteMode', 'append');
hdf5write(fileName, '/dataX', vxFactor, 'WriteMode', 'append');
hdf5write(fileName, '/dataY', vyFactor, 'WriteMode', 'append');
