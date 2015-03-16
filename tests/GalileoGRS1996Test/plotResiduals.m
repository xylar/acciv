clear variables;

pointCount = 0;

maxPoints = 5000;

indexFractions = rand(1,maxPoints);

for index = 0:9
    
    %fileName = 'full/pass1/work/noCurvedPathScatteredVelocity.h5';
    fileName = sprintf('full/pass3/work/curvedPathScatteredVelocity_%i.h5',index);

    x = hdf5read(fileName, '/x');
    y = hdf5read(fileName, '/y');
    residuals = hdf5read(fileName, '/residuals');

    medianRes = median(residuals)

    figure(1);
    hist(residuals/medianRes,100);

    maxPoints = 5000;
    if(length(x) > maxPoints)
        indices = ceil(indexFractions*length(x));
    else
        indices = 1:length(x);
    end
    figure(2);
    plot3(x(indices),y(indices),residuals(indices)/medianRes, '.');

    figure(3);
    plot(x,y, '.', 'MarkerSize', 1);

end