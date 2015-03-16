import numpy
import h5py
image0 = numpy.loadtxt('image0.csv', delimiter=',')
image1 = numpy.loadtxt('image1.csv', delimiter=',')
mask = numpy.ones(image0.shape,numpy.uint8)
bounds = [0.0,image0.shape[1],0.0,image0.shape[0]]
time = 0.0
h5File = h5py.File('image0.h5', 'w')
dataset = h5File.create_dataset("bounds", data=numpy.array(bounds))
dataset = h5File.create_dataset("data", data=image0)
dataset = h5File.create_dataset("mask", data=mask)
dataset = h5File.create_dataset("time", data=numpy.array(time))
h5File.close()
time = 1.0
h5File = h5py.File('image1.h5', 'w')
dataset = h5File.create_dataset("bounds", data=numpy.array(bounds))
dataset = h5File.create_dataset("data", data=image1)
dataset = h5File.create_dataset("mask", data=mask)
dataset = h5File.create_dataset("time", data=numpy.array(time))
h5File.close()
