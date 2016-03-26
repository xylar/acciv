#!/usr/bin/env/python

import os,re
import numpy as np
import h5py

from attrdict import AttrMap

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors


import scipy as scp

import skimage.exposure as exposure
from skimage.io import imread
from skimage.filters import *
from skimage.morphology import *
from skimage import data
from skimage import img_as_float, img_as_bool
from skimage.morphology import disk
from skimage import measure


images = ['images/image000.png', 'images/image001.png']


for idx,path in enumerate(images):
  
  image = img_as_float(imread(path))
  
  # bounds represent orthonormal standart canonical coordinate system 
  #  | y
  #  |
  # -----------> x
  # bounds = [ x_min, x_max, y_min, y_max ]
  
  bounds = np.array([ 0 , image.shape[1] , 0 , image.shape[0] ] )
  
  mask = np.ones(image.shape,np.uint8)
  time = idx * 1 # 1 second between images -> no scaling
  
  h5File = h5py.File('image%03i.h5'%idx, 'w')
  dataset = h5File.create_dataset("bounds", data=bounds)
  dataset = h5File.create_dataset("data", data=np.flipud(images))
  dataset = h5File.create_dataset("mask", data=np.flipud(mask))
  dataset = h5File.create_dataset("time", data=numpy.array(time))
  h5File.close() 



