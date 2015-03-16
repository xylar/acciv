#!/usr/bin/python
import numpy
import h5py
import matplotlib.pyplot as plt

from basemapInterp2D import basemapInterp2D


inImageFileName = '../syntheticCuspVortex/image001.h5'
outFolder= '.'

h5File = h5py.File(inImageFileName, 'r')
bounds = h5File["bounds"][...]
inImageData = h5File["data"][...]
inMask = numpy.array(h5File["mask"][...],float)
h5File.close()


inx = numpy.linspace(bounds[0],bounds[1],inImageData.shape[1])
iny = numpy.linspace(bounds[2],bounds[3],inImageData.shape[0])

(inX,inY) = numpy.meshgrid(inx,iny)

vx = 0.4*(inY-0.5)
vy = numpy.zeros(inImageData.shape)


vmax = numpy.amax(inImageData)
vmin = numpy.amin(inImageData)
  

ts = numpy.linspace(0.,1.,4)

for tIndex in range(len(ts)):
  t = ts[tIndex]
  
  # find the position backwards in time
  outX = -vx*t+inX
  outY = -vy*t+inY
  
  maskValue = -1e8
  outImageData = basemapInterp2D(inImageData,inx,iny,outX,outY,masked=maskValue, order=1)
  outMask = basemapInterp2D(inMask,inx,iny,outX,outY,masked=maskValue,order=1)
  outImageData[outImageData == maskValue] = 0.0
  outMask = numpy.array(outMask > 0.99,numpy.uint8)
  outImageData *= outMask
  
  h5File = h5py.File('%s/image%03i.h5'%(outFolder,tIndex+1), 'w')
  dataset = h5File.create_dataset("bounds", data=bounds)
  dataset = h5File.create_dataset("data", data=outImageData)
  dataset = h5File.create_dataset("mask", data=outMask)
  dataset = h5File.create_dataset("time", data=numpy.array(t))
  h5File.close() 

  plt.imsave('%s/image%03i.png'%(outFolder,tIndex+1), outImageData[::-1,:],
             vmin=vmin, vmax=vmax, cmap='gray')

#  fig = plt.figure(1)
#  ax = fig.add_subplot(111)
#  plt.imshow(inMask, extent=[bounds[0],bounds[1],bounds[3],bounds[2]], cmap='gray')
#  ax.set_ylim(ax.get_ylim()[::-1])
#  plt.axis('tight')
#
#  fig = plt.figure(2)
#  ax = fig.add_subplot(111)
#  plt.imshow(outMask, extent=[bounds[0],bounds[1],bounds[3],bounds[2]], cmap='gray')
#  ax.set_ylim(ax.get_ylim()[::-1])
#  plt.axis('tight')
#  
#  plt.show()

