#!/usr/bin/python
import numpy
import h5py
from optparse import OptionParser
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os.path
parser = OptionParser()
parser.add_option("--inFolder", type="string", default='.', dest="inFolder", help="folder where the images are")
parser.add_option("--outFolder", type="string", default='.', dest="outFolder", help="folder to write to")
parser.add_option("--imageCount", type="int", default=4, dest="imageCount", help="the number of images")


options, args = parser.parse_args()

maxImageStd = 3.0 # Image data is clamped to be within 3 std. dev. from the mean.
                  # Decrease this number to increase image contrast.

# the data is cropped to be within this box [leftLon,rightLon,botLat,topLat]
cropBounds = [328.1011223432933,
308.8828490718264,
-27.212863196130208,
-12.803165388750081]

for imageIndex in range(options.imageCount):
  imageFileName = '%s/earlierAdvectedImage_%i.h5'%(options.inFolder,imageIndex+1)
  if not(os.path.exists(imageFileName)):
    continue
  print imageFileName
  h5File = h5py.File(imageFileName, 'r')
  bounds = h5File["bounds"][...]
  imageData = h5File["data"][...]
  imageMask = numpy.array(h5File["mask"][...],bool)
  h5File.close()

  gx = numpy.linspace(bounds[0],bounds[1],imageData.shape[1])
  gy = numpy.linspace(bounds[2],bounds[3],imageData.shape[0])
    
  xMin = numpy.argmin(numpy.abs(gx-cropBounds[0]))
  xMax = numpy.argmin(numpy.abs(gx-cropBounds[1]))+1
  yMin = numpy.argmin(numpy.abs(gy-cropBounds[2]))
  yMax = numpy.argmin(numpy.abs(gy-cropBounds[3]))+1
  
  imageCropped = imageData[yMin:yMax,xMin:xMax]
  maskCropped = imageMask[yMin:yMax,xMin:xMax]

  imageMean = numpy.mean(imageCropped[maskCropped])
  imageStd = numpy.std(imageCropped[maskCropped])
  
  imageCropped *= maskCropped
  vmin = imageMean-maxImageStd*imageStd
  vmax = imageMean+maxImageStd*imageStd
  
  plt.imsave('%s/earlierAdvectedImage_%i.png'%(options.outFolder,imageIndex+1), imageCropped[::-1,:],
             vmin=vmin, vmax=vmax, cmap=cm.get_cmap('gray'))


