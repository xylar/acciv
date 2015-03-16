import numpy
from Vicar import readVicar
import h5py

import dateutil.parser
import datetime

import matplotlib.pyplot as plt

outFolder = '../../Omid'
#inFolder = '../../Omid/1673pho2cb2cl2gmaps'
#fileNames = ['W1673461218_1.map', 'W1673465178_1.map', 'W1673466365_1.map', 'W1673470325_1.map']

inFolder = '../../Omid'
fileNames = ['W1673470325_1_for_tracking.map']
#fileNames = ['W1585476794_1.map']

for index in range(len(fileNames)):
  fileName = '%s/%s'%(inFolder,fileNames[index])  
  (image,metadata) = readVicar(fileName)

  image = numpy.array(image,float)
  # flip the image over in y so north is now up
  image = image[::-1,:]

  fig = plt.figure(1)
  ax = fig.add_subplot(111, aspect='equal')
  plt.imshow(image*(image > 0), cmap='gray',interpolation='nearest')
  plt.colorbar()
  plt.title('image with bounds in pixels')
  ax.set_ylim(ax.get_ylim()[::-1])

  print 'image size:', image.shape
  (ny,nx) = image.shape
  
  mask = numpy.logical_not(image <= 0.0)
  print 'image intensity range:', numpy.amin(image[mask]), numpy.amax(image[mask])
  mask = numpy.array(mask,numpy.uint8)
  
  date = metadata.IMAGE_MID_TIME
  day = int(date[5:8])
  date = date[0:5]+'01-01'+date[8:]
  imageDate = dateutil.parser.parse(date)+datetime.timedelta(days=day)
  if(index == 0):
    firstDate = imageDate
  time = (imageDate-firstDate).total_seconds()
  
  Re = float(metadata.A_AXIS_RADIUS)*1000. # in m
  Rp = float(metadata.C_AXIS_RADIUS)*1000. # in m
  epsilon = Rp/Re
  kmPerPixel = float(metadata.MAP_SCALE)*1000. # in m
  
  lonPerPixel = 180./numpy.pi*(kmPerPixel/Re)
  latPerPixel = lonPerPixel


  if(hasattr(metadata,'CENTLAT')): # cropped
    xc = int(metadata.CENTSAMP)
    yc = int(metadata.CENTLINE) 
    lat0 = float(metadata.CENTLAT)
    lon0 = float(metadata.CENTLON)
  elif((nx == 3601) and (ny == 1801)): # full 360 by 180 map

    lat0 = 0.0 #float(metadata.CENTER_LATITUDE)
    lon0 = 0.0 #float(metadata.CENTER_LONGITUDE)
    xc = (nx-1)/2
    yc = (ny-1)/2
    lonPerPixel = 0.1 # fix rounding error?
    latPerPixel = 0.1
  else:
    print "Unknown map format -- couldn't find CENTLAT and does not appear to be a full map (3601 by 1801)"
    exit()

  print 'degrees per pixel:', lonPerPixel

  #kmPerPixel_y = epsilon**2*kmPerPixel_x
  
  planetocentric = (metadata.COORDINATE_SYSTEM_NAME == 'PLANETOCENTRIC')
  print 'Is the projection planetocentric?', planetocentric
  
  lon = numpy.arange(nx-1,-1,-1)*lonPerPixel
  lon += lon0 - lon[xc]
  
  lat = numpy.arange(ny)*latPerPixel
  lat += lat0 - lat[yc]
    
  print lon[xc], lon0
  print lat[yc], lat0
  bounds = [lon[0],lon[-1],lat[0],lat[-1]]
  print 'bounds in lon/lat:', bounds
  
  extent = [bounds[0],bounds[1],bounds[3],bounds[2]]
  fig = plt.figure(2)
  ax = fig.add_subplot(111, aspect='equal')
  plt.imshow(image*mask, extent = extent, cmap='gray',interpolation='nearest')
  plt.colorbar()
  plt.title('image with bounds in degrees')
  ax.set_ylim(ax.get_ylim()[::-1])
  
  h5File = h5py.File('%s/image%03i.h5'%(outFolder,index+1), 'w')
  dataset = h5File.create_dataset("bounds", data=numpy.array(bounds))
  dataset = h5File.create_dataset("data", data=image)
  dataset = h5File.create_dataset("mask", data=mask)
  dataset = h5File.create_dataset("time", data=numpy.array(time))
  h5File.close() 
  
  if(index == 0): 
    #write out the geometry factors, too
    (lon, lat) = numpy.meshgrid(lon*numpy.pi/180.,lat*numpy.pi/180.)
    if(planetocentric):
      lat = numpy.arctan(numpy.tan(lat)/epsilon**2)
    secLat = 1.0/numpy.cos(lat)
    tanLat = numpy.tan(lat)
    t = numpy.arctan(epsilon*tanLat)
    dtdLatitude = numpy.pi/180.*epsilon*(secLat**2/(1 + epsilon**2*tanLat**2))
    dydLatitude = numpy.sqrt((Re*numpy.sin(t))**2 + (Rp*numpy.cos(t))**2)*dtdLatitude
    dxdLongitude = -Re*numpy.cos(t)*numpy.pi/180.
    h5File = h5py.File('%s/gridGeometryFactors.h5'%(outFolder), 'w')
    dataset = h5File.create_dataset("bounds", data=numpy.array(bounds))
    dataset = h5File.create_dataset("dataX", data=dxdLongitude)
    dataset = h5File.create_dataset("dataY", data=dydLatitude)
    h5File.close() 

    fig = plt.figure(3)
    ax = fig.add_subplot(111, aspect='equal')
    plt.imshow(dxdLongitude, extent = extent, interpolation='nearest')
    plt.colorbar()
    plt.title('geometry factor dx/dlon')
    ax.set_ylim(ax.get_ylim()[::-1])

    fig = plt.figure(4)
    ax = fig.add_subplot(111, aspect='equal')
    plt.imshow(dydLatitude, extent = extent, interpolation='nearest')
    plt.colorbar()
    plt.title('geometry factor dy/dlat')
    ax.set_ylim(ax.get_ylim()[::-1])
  
  plt.show()

