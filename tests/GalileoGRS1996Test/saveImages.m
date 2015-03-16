fileName = 'grid.h5';
imageData = imread(sprintf('grs%iman.tif', 1));

nx = size(imageData,2);
ny = size(imageData,1);

%The map scale is 20 km/pix, defined at the equator (R=71492 km).  In 
%degrees per pixel it is:

degpix = 20./(2. * 71492. * pi) * 360.;  %= 0.0160286

%The conversion is:

lats = (450 - [ny,1]) * degpix - 20;
lons = (600 - [1,nx]) * degpix + 318.5;

% If we did it correctly, (1,1) should be at (12.80S, 328.10W).
% The timing is:
% grs1man.tif = 1996, day of year 178,  4:19:38.9
% grs2man.tif = 1996, day of year 178, 13:19:34.9
% grs3man.tif = 1996, day of year 178, 14:30:21.5
times = [4,19,38.9;
    13,19,34.9;
    14,30,21.5];

times = 3600*times(:,1)+60*times(:,2)+times(:,3);

bounds = [lons, lats];
hdf5write(fileName, '/imageSize', int32([nx, ny]), 'WriteMode', 'overwrite');
hdf5write(fileName, '/imageBounds', bounds, 'WriteMode', 'append');

for(imageIndex = 1:3)
    imageData = double(imread(sprintf('grs%iman.tif', imageIndex)))/255;
    imageData = imageData(end:-1:1,:);
    mask = imageData ~= 0;
    fileName = sprintf('image%03i.h5', imageIndex);
    hdf5write(fileName, '/data', imageData', 'WriteMode', 'overwrite');
    hdf5write(fileName, '/bounds', bounds, 'WriteMode', 'append');
    hdf5write(fileName, '/mask', uint8(mask)', 'WriteMode', 'append');
    hdf5write(fileName, '/time', times(imageIndex), 'WriteMode', 'append');
end