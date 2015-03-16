#!/usr/bin/python
import numpy
import h5py
from optparse import OptionParser
import matplotlib
import os.path

parser = OptionParser()
parser.add_option("--folder", type="string", default='full/pass1/', dest="folder", help="folder of the output data to be plotted")
parser.add_option("--imageFileName", type="string", default='image001.h5', dest="imageFileName", help="folder of the output data to be plotted")
parser.add_option("--savePlots", action="store_true", dest="savePlots", help="include this flag save plots to files instead of displaying them")
parser.add_option("--gridFileName", type="string", default='outGridVelocity.h5', dest="gridFileName")
parser.add_option("--scatterFileName", type="string", default='outScatteredVelocity.h5', dest="scatterFileName")
parser.add_option("--figurePrefix", type="string", default='fig', dest="figurePrefix")
parser.add_option("--tiePointsFolder", type="string", default='_work', dest="tiePointsFolder")

options, args = parser.parse_args()


if options.savePlots:
  matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.colors as colors

folder = options.folder
scatterFileName = '%s/%s'%(folder,options.scatterFileName)
if(not os.path.exists(scatterFileName)):
  #print "not found:", scatterFileName
  exit()
gridFileName = '%s/%s'%(folder,options.gridFileName)
if(not os.path.exists(gridFileName)):
  print "not found:", gridFileName
  exit()
tiePointsFileName = '%s/%s/combinedCorrelationTiePoints.h5'%(folder,options.tiePointsFolder)
imageFileName = options.imageFileName
if(not os.path.exists(imageFileName)):
  print "not found:", imageFileName
  exit()
# plot a velocity vector every "skip" pixels from the gridded velocity data
skip = 8

# width and height of each figure (inches)
width = 12
height = 6.75

# number of points to be sampled from the scattered data
maxPoints = 10000

# the locations of the major and minor axes to plot
x0 = 319
y0 = -20.5
# the width around each axis to take points from when plotting axes
dx = 0.5
dy = 0.5

maxImageStd = 3.0 # Image data is clamped to be within 3 std. dev. from the mean.
                  # Decrease this number to increase image contrast.

# the data is cropped to be within this box [leftLon,rightLon,botLat,topLat]
cropBounds = [328.1011223432933,
308.8828490718264,
-27.212863196130208,
-12.803165388750081]
#cropBounds = [320,295,25,39]
#ellipseBounds = [310,302,33,40]
ellipseBounds = cropBounds

# decrease these numbers to increase the length of vectors, and visa versa
scatterVectorScale = 400.0
gridVectorScale = 800.0

h5File = h5py.File(imageFileName, 'r')
bounds = h5File["bounds"][...]
imageData = h5File["data"][...]
imageMask = numpy.array(h5File["mask"][...],bool)
h5File.close()


h5File = h5py.File(scatterFileName, 'r')
x = h5File["x"][...]
y = h5File["y"][...]
vx = h5File["vx"][...]
vy = h5File["vy"][...]
pixelVx = h5File["dataX"][...]
pixelVy = h5File["dataY"][...]
h5File.close()


h5File = h5py.File(tiePointsFileName, 'r')
deltaTs = h5File["deltaTs"][...]
residualsFound = "correlationVelocityResiduals" in h5File
if residualsFound:
  correlationVelocityResiduals = h5File["correlationVelocityResiduals"][...]
  correlationLocationResiduals = h5File["correlationLocationResiduals"][...]
h5File.close()

maxDeltaT = numpy.amax(deltaTs)

#print numpy.amax(numpy.abs(vx))
#print numpy.amax(numpy.abs(vy))

h5File = h5py.File(gridFileName, 'r')
gridVx = h5File["vx"][...]
gridVy = h5File["vy"][...]
h5File.close()
gx = numpy.linspace(bounds[0],bounds[1],gridVx.shape[1])
gy = numpy.linspace(bounds[2],bounds[3],gridVx.shape[0])
#print numpy.amax(numpy.abs(gridVx))
#print numpy.amax(numpy.abs(gridVy))

dLon = gx[1]-gx[0]
dLat = gy[1]-gy[0]

pixelVx = pixelVx/dLon*maxDeltaT
pixelVy = pixelVy/dLat*maxDeltaT

xMin = numpy.argmin(numpy.abs(gx-cropBounds[0]))
xMax = numpy.argmin(numpy.abs(gx-cropBounds[1]))+1
yMin = numpy.argmin(numpy.abs(gy-cropBounds[2]))
yMax = numpy.argmin(numpy.abs(gy-cropBounds[3]))+1

imageCropped = imageData[yMin:yMax,xMin:xMax]
maskCropped = imageMask[yMin:yMax,xMin:xMax]
boundsCropped = [gx[xMin],gx[xMax-1],gy[yMin],gy[yMax-1]]

# crop the gridded velocity
gridVx = gridVx[yMin:yMax,xMin:xMax]
gridVy = gridVy[yMin:yMax,xMin:xMax]
[gridX, gridY] = numpy.meshgrid(gx[xMin:xMax],gy[yMin:yMax])

# crop the scattered velocity
mask = numpy.logical_and(
  numpy.logical_and(x >= boundsCropped[1],x <= boundsCropped[0]),
  numpy.logical_and(y >= boundsCropped[2],y <= boundsCropped[3]))
x = x[mask]
y = y[mask]
vx = vx[mask]
vy = vy[mask]
pixelVx = pixelVx[mask]
pixelVy = pixelVy[mask]

vMag = numpy.sqrt(vx**2+vy**2)


xc = 0.5*(ellipseBounds[0]+ellipseBounds[1])
yc = 0.5*(ellipseBounds[2]+ellipseBounds[3])
xr = 0.5*(ellipseBounds[0]-ellipseBounds[1])
yr = 0.5*(ellipseBounds[3]-ellipseBounds[2])

ellipseMask = ((x-xc)/xr)**2 + ((y-yc)/yr)**2 <= 1.0

mask = numpy.logical_and(ellipseMask,numpy.abs(y-y0) < dy)
vxMean = numpy.mean(vx[mask])
print "mean vx along the major axis:", vxMean
vxMean = 0.0

imageMean = numpy.mean(imageCropped[maskCropped])
imageStd = numpy.std(imageCropped[maskCropped])

imageCropped *= maskCropped
imageCropped = numpy.maximum(imageMean-maxImageStd*imageStd,
   numpy.minimum(imageMean+maxImageStd*imageStd,imageCropped))

if(x.size > maxPoints):
  indices = numpy.array(numpy.random.rand(maxPoints)*x.size,int)
else:
  indices = numpy.array(numpy.linspace(0,x.size-1,x.size),int)

colorList1 = numpy.array(((0.0,0.0,0.0),
	(1.0,0.0,0.0),
	(1.0,1.0,0.0),
	(1.0,1.0,1.0)))

colorList2 = numpy.array(((0.5,0.5,0.5),
	(0.67,0.67,0.67),
	(0.83,0.83,0.83),
	(1.0,1.0,1.0)))
alpha = 0.25

colorList = alpha*colorList1 + (1-alpha)*colorList2
colorList = tuple(tuple(x) for x in colorList)

colormap = colors.LinearSegmentedColormap.from_list('my_map',colorList,N=256)


maskXAxis = numpy.abs(y-y0) < dy
maskYAxis = numpy.abs(x-x0) < dx
xAxisIndices = indices[maskXAxis[indices]]
yAxisIndices = indices[maskYAxis[indices]]

xAxisGridIndex = numpy.argmin(numpy.abs(gy-y0))
yAxisGridIndex = numpy.argmin(numpy.abs(gx-x0))

fig = plt.figure(1, figsize=[width,height])
#fig.subplots_adjust(left=0.075, right=0.975, bottom=0.05, top=0.95, wspace=0.2, hspace=0.25)
ax = fig.add_subplot(111)
plt.imshow(imageCropped, extent=(boundsCropped[0],boundsCropped[1],boundsCropped[3],boundsCropped[2]), cmap=colormap)
ax.set_ylim(ax.get_ylim()[::-1])
plt.quiver(x[indices], y[indices], vx[indices]-vxMean, vy[indices], color='k', pivot='mid',  scale_units='xy', scale=scatterVectorScale)
plt.quiver(x[xAxisIndices], y[xAxisIndices], vx[xAxisIndices]-vxMean, vy[xAxisIndices], color='r', pivot='mid',  scale_units='xy', scale=scatterVectorScale)
plt.quiver(x[yAxisIndices], y[yAxisIndices], vx[yAxisIndices]-vxMean, vy[yAxisIndices], color='b', pivot='mid',  scale_units='xy', scale=scatterVectorScale)
plt.title('a sample of %i scattered velocity vectors'%(indices.size))
ax.set_aspect('equal')
ax.autoscale(tight=True)


#fig = plt.figure(12, figsize=[width,height])
##fig.subplots_adjust(left=0.075, right=0.975, bottom=0.05, top=0.95, wspace=0.2, hspace=0.25)
#ax = fig.add_subplot(111)
#plt.imshow(imageCropped, extent=(boundsCropped[0],boundsCropped[1],boundsCropped[3],boundsCropped[2]), cmap=colormap)
#ax.set_ylim(ax.get_ylim()[::-1])
#plt.axis('tight')


fig = plt.figure(2, figsize=[width,height])
#fig.subplots_adjust(left=0.075, right=0.975, bottom=0.05, top=0.95, wspace=0.2, hspace=0.25)
ax = fig.add_subplot(111)
plt.imshow(imageCropped, extent=(boundsCropped[0],boundsCropped[1],boundsCropped[3],boundsCropped[2]), cmap=colormap)
ax.set_ylim(ax.get_ylim()[::-1])
plt.quiver(gridX[::skip,::skip], gridY[::skip,::skip], gridVx[::skip,::skip]-vxMean, gridVy[::skip,::skip], color='k', pivot='mid',  scale_units='xy', scale=gridVectorScale)
plt.title('gridded velocity vector (skip = %i)'%skip)
ax.set_aspect('equal')
ax.autoscale(tight=True)

fig = plt.figure(3, figsize=[width,height])
#fig.subplots_adjust(left=0.075, right=0.975, bottom=0.05, top=0.95, wspace=0.2, hspace=0.25)
ax = fig.add_subplot(111)
plt.imshow(gridVx-vxMean, extent=(boundsCropped[0],boundsCropped[1],boundsCropped[3],boundsCropped[2]), cmap=plt.get_cmap('jet'))
ax.set_ylim(ax.get_ylim()[::-1])
plt.colorbar()
plt.title('vx')
ax.set_aspect('equal')
ax.autoscale(tight=True)

fig = plt.figure(4, figsize=[width,height])
#fig.subplots_adjust(left=0.075, right=0.975, bottom=0.05, top=0.95, wspace=0.2, hspace=0.25)
ax = fig.add_subplot(111)
plt.imshow(gridVy, extent=(boundsCropped[0],boundsCropped[1],boundsCropped[3],boundsCropped[2]), cmap=plt.get_cmap('jet'))
ax.set_ylim(ax.get_ylim()[::-1])
plt.colorbar()
plt.title('vy')
ax.set_aspect('equal')
ax.autoscale(tight=True)


fig = plt.figure(5, figsize=[width,height])
#fig.subplots_adjust(left=0.075, right=0.975, bottom=0.05, top=0.95, wspace=0.2, hspace=0.25)
ax = fig.add_subplot(111)
plt.imshow(numpy.sqrt((gridVx-vxMean)**2 + gridVy**2), extent=(boundsCropped[0],boundsCropped[1],boundsCropped[3],boundsCropped[2]), cmap=plt.get_cmap('jet'))
ax.set_ylim(ax.get_ylim()[::-1])
plt.colorbar()
plt.title('|v|')
ax.set_aspect('equal')
ax.autoscale(tight=True)

weights = numpy.ones(vx.shape)/vx.size

fig = plt.figure(6, figsize=[width,height])
#fig.subplots_adjust(left=0.075, right=0.975, bottom=0.05, top=0.95, wspace=0.2, hspace=0.25)
ax = fig.add_subplot(111)
plt.hist(vMag,100,weights=weights,histtype='step')
plt.hist(vx,100,weights=weights,histtype='step')
plt.hist(vy,100,weights=weights,histtype='step')
plt.xlabel('velocity')
plt.ylabel('tie point fraction')
plt.title('velocity histograms')
plt.axis('tight')
plt.legend(['|v|','vx','vy'])

fig = plt.figure(7, figsize=[width,height])
#fig.subplots_adjust(left=0.075, right=0.975, bottom=0.05, top=0.95, wspace=0.2, hspace=0.25)
ax = fig.add_subplot(111)
plt.hist(numpy.sqrt(pixelVx**2+pixelVy**2),100,weights=weights,histtype='step')
plt.hist(pixelVx,100,weights=weights,histtype='step')
plt.hist(pixelVy,100,weights=weights,histtype='step')
plt.xlabel('velocity*maxDeltaT (pixels)')
plt.ylabel('tie point fraction')
plt.title('pixel offset histograms (search range)')
plt.axis('tight')
plt.legend(['|v|','vx','vy'])

#x0 = 0.5*(boundsCropped[0]+boundsCropped[1])
#y0 = 0.5*(boundsCropped[2]+boundsCropped[3])
#dx = 0.02*(boundsCropped[1]-boundsCropped[0])
#dy = 0.02*(boundsCropped[3]-boundsCropped[2])

fig = plt.figure(8, figsize=[width,height])
#fig.subplots_adjust(left=0.075, right=0.975, bottom=0.05, top=0.95, wspace=0.2, hspace=0.25)
ax = fig.add_subplot(111)
plt.plot(x[maskXAxis], vy[maskXAxis], '.k',gx,gridVy[xAxisGridIndex,:],'r')
plt.title('vy along x axis within dy = %.1f of y = %.1f'%(dy,y0))
plt.xlabel('x')
plt.ylabel('vy')
plt.axis('tight')
ax.set_xlim(ax.get_xlim()[::-1])

fig = plt.figure(9, figsize=[width,height])
#fig.subplots_adjust(left=0.075, right=0.975, bottom=0.05, top=0.95, wspace=0.2, hspace=0.25)
ax = fig.add_subplot(111)
plt.plot(vx[maskYAxis]-vxMean, y[maskYAxis], '.k',gridVx[:,yAxisGridIndex],gy,'b')
plt.title('vx along y axis within dx = %.1f of x = %.1f'%(dx,x0))
plt.xlabel('vx')
plt.ylabel('y')
plt.axis('tight')


if residualsFound:
  maxVal = 6.0*numpy.median(correlationVelocityResiduals)
  fig = plt.figure(10, figsize=[width,height])
  #fig.subplots_adjust(left=0.075, right=0.975, bottom=0.05, top=0.95, wspace=0.2, hspace=0.25)
  ax = fig.add_subplot(111)
  weights = numpy.ones(correlationVelocityResiduals.shape)/correlationVelocityResiduals.size
  plt.hist(correlationVelocityResiduals,100,range=[0.0,maxVal], weights=weights,histtype='step')
  plt.xlabel('correlation velocity uncertainty')
  plt.ylabel('tie point fraction')
  
  maxVal = 6.0*numpy.median(correlationLocationResiduals)/1000
  fig = plt.figure(11, figsize=[width,height])
  #fig.subplots_adjust(left=0.075, right=0.975, bottom=0.05, top=0.95, wspace=0.2, hspace=0.25)
  ax = fig.add_subplot(111)
  weights = numpy.ones(correlationLocationResiduals.shape)/correlationLocationResiduals.size
  plt.hist(correlationLocationResiduals/1000,100,range=[0.0,maxVal], weights=weights,histtype='step')
  plt.xlabel('correlation location uncertainty (km)')
  plt.ylabel('tie point fraction')


plt.draw()
if options.savePlots:
  lastFig = 9
  if residualsFound:
    lastFig = 11
  for index in range(1,lastFig+1):
    outFileName = '%s/%s%03i.png'%(folder,options.figurePrefix, index)
    plt.figure(index)
    plt.savefig(outFileName)
else:
  plt.show()

