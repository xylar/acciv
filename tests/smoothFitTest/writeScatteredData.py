#!/usr/bin/python
import numpy
import numpy.random
import h5py
import matplotlib.pyplot as plt

def velX(x,y,bounds):
  x0 = 0.5*(bounds[0]+bounds[1])
  y0 = 0.5*(bounds[2]+bounds[3])
  denom = (x0-bounds[0])*(x0-bounds[1])*(y0-bounds[2])*(y0-bounds[3])
  result = (x-bounds[0])*(x-bounds[1])*(y-bounds[2])*(y-bounds[3])/denom
  return result

outFileName = 'parabolic2/scatteredField.h5'

pointCount = 10000

vx0 = 0.2
vy0 = 0.0

bounds = [3., -1., 0., 2.]

x = 4*numpy.random.rand((pointCount))-1;
y = 2*numpy.random.rand((pointCount));

#vx = vx0*2*(x-0.5)+0.01*(2*numpy.random.rand((pointCount))-1)
#vy = vy0*numpy.ones(y.shape)+0.01*(2*numpy.random.rand((pointCount))-1)

vx = velX(x,y,bounds)
vy = vy0*numpy.ones(y.shape)

vx += 0.01*(2*numpy.random.rand((pointCount))-1);

  
h5File = h5py.File(outFileName, 'w')
dataset = h5File.create_dataset("dataX", data=vx)
dataset = h5File.create_dataset("dataY", data=vy)
dataset = h5File.create_dataset("x", data=x)
dataset = h5File.create_dataset("y", data=y)
h5File.close()

nx = 401;
ny = 201;
(xGrid, yGrid) = numpy.meshgrid(numpy.linspace(bounds[0],bounds[1],nx),
                                numpy.linspace(bounds[2],bounds[3],ny))
vxGrid = velX(xGrid,yGrid,bounds)
vyGrid = vy0*numpy.ones(yGrid.shape)

outFileName = 'parabolic2/grid.h5'
h5File = h5py.File(outFileName, 'w')
dataset = h5File.create_dataset("dataX", data=vxGrid)
dataset = h5File.create_dataset("dataY", data=vyGrid)
dataset = h5File.create_dataset("bounds", data=bounds)
h5File.close() 

plt.figure(1)
plt.plot(x,y,'.')
plt.show()

