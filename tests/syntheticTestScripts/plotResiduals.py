#!/usr/bin/python
import numpy
import h5py
from optparse import OptionParser
import matplotlib

parser = OptionParser()
parser.add_option("--folder", type="string", default='full/pass1/', dest="folder", help="folder of the output data to be plotted")
parser.add_option("--savePlots", action="store_true", dest="savePlots", help="include this flag save plots to files instead of displaying them")
parser.add_option("--scatterFileName", type="string", default='outScatteredVelocity.h5', dest="scatterFileName")
parser.add_option("--figurePrefix", type="string", default='fig', dest="figurePrefix")


options, args = parser.parse_args()


if options.savePlots:
  matplotlib.use('Agg')

import matplotlib.pyplot as plt

folder = options.folder
scatterFileName = '%s/%s'%(folder,options.scatterFileName)

# width and height of each figure (inches)
width = 12
height = 10


h5File = h5py.File(scatterFileName, 'r')
x = h5File["x"][...]
y = h5File["y"][...]
resX = h5File["residualX"][...]
resY = h5File["residualY"][...]
h5File.close()


fig = plt.figure(1, figsize=[width,height])
fig.subplots_adjust(left=0.075, right=0.975, bottom=0.05, top=0.95, wspace=0.2, hspace=0.25)
ax = fig.add_subplot(111, aspect='equal')
plt.plot(x, resX, '.k')
plt.xlabel('x')
plt.ylabel('resX')
plt.axis('tight')

fig = plt.figure(2, figsize=[width,height])
fig.subplots_adjust(left=0.075, right=0.975, bottom=0.05, top=0.95, wspace=0.2, hspace=0.25)
ax = fig.add_subplot(111, aspect='equal')
plt.plot(y, resX, '.k')
plt.xlabel('y')
plt.ylabel('resX')
plt.axis('tight')

fig = plt.figure(3, figsize=[width,height])
fig.subplots_adjust(left=0.075, right=0.975, bottom=0.05, top=0.95, wspace=0.2, hspace=0.25)
ax = fig.add_subplot(111, aspect='equal')
plt.plot(x, resY, '.k')
plt.xlabel('x')
plt.ylabel('resY')
plt.axis('tight')

fig = plt.figure(4, figsize=[width,height])
fig.subplots_adjust(left=0.075, right=0.975, bottom=0.05, top=0.95, wspace=0.2, hspace=0.25)
ax = fig.add_subplot(111, aspect='equal')
plt.plot(y, resY, '.k')
plt.xlabel('y')
plt.ylabel('resY')
plt.axis('tight')


plt.draw()
if options.savePlots:
  outFileName = '%s/%s_x_resX.png'%(folder,options.figurePrefix)
  plt.figure(1)
  plt.savefig(outFileName)
  outFileName = '%s/%s_y_resX.png'%(folder,options.figurePrefix)
  plt.figure(2)
  plt.savefig(outFileName)
  outFileName = '%s/%s_x_resY.png'%(folder,options.figurePrefix)
  plt.figure(3)
  plt.savefig(outFileName)
  outFileName = '%s/%s_y_resY.png'%(folder,options.figurePrefix)
  plt.figure(4)
  plt.savefig(outFileName)
else:
  plt.show()

