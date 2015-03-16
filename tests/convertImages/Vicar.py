#!/usr/bin/python

import numpy

# based on the vicar2png module by Jessica McKellar (jesstess at mit.edu)
# substantial modifications have been made to the code.  However for
# thoroughness, I am including her Copyright under the MIT License below:

'''
The MIT License (MIT)

Copyright (c) 2012-2013 Jessica McKellar

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.
'''

class VICARMetadata(object):
    """
    Contains VICAR metadata accessible as uppercase class attributes,
    e.g.:

    vicar.RECSIZE
    vicar.FORMAT
    """
    def __init__(self, metadata):
        """
        metadata: A dictionary of VICAR label/value pairs.
        """
        for key, value in metadata.iteritems():
            if value.isdigit():
                value = int(value)
            setattr(self, key.upper(), value)


def addMetadataToDict(metadata, metadata_dict):
    gettingTag = True
    has_lparen = False
    has_lquote = False
    tag_buf = ''
    for char in metadata:
      if gettingTag:
        if char == '=':
          tag = tag_buf
          tag_buf = ''
          
          gettingTag = False
          has_lparen = False
          has_lquote = False
        elif char != ' ':
          tag_buf += char
          
      else: # getting value
        if char == "'":
          has_lquote = not has_lquote
          if has_lparen:
            tag_buf += char
        elif char == "(" and not has_lquote:
          has_lparen = True
          tag_buf += char
        elif char == ")" and not has_lquote:
          has_lparen = False
          tag_buf += char
        elif char == " " and tag_buf and not (has_lquote or has_lparen):
          # We have a full value, save it.
          value = tag_buf
          metadata_dict[tag] = value
          
          gettingTag = True
          has_lparen = False
          has_lquote = False
          tag_buf = ""
        elif char == " " and not has_lquote:
          continue
        else:
            tag_buf += char
    return metadata_dict
            
def process_metadata(metadata_fd):
    # A VICAR file must start with 'LBLSIZE=<integer label size>'.
    lblsize_field = metadata_fd.read(len("LBLSIZE="))
    if lblsize_field.upper() != "LBLSIZE=":
        raise ValueError("Malformed VICAR file: doesn't start with LBLSIZE.")

    lblsize = ""
    while True:
        char = metadata_fd.read(1)
        if char == " ":
            break
        else:
            lblsize += char
    try:
        lblsize = int(lblsize)
    except ValueError:
        raise ValueError("Malformed VICAR file: contains non-integer LBLSIZE.")

    # Read in the rest of the VICAR metadata.
    metadata_fd.seek(0)
    metadata = metadata_fd.read(lblsize)

    metadata_dict = {}
    
    metadata_dict = addMetadataToDict(metadata, metadata_dict)
    
    vicar = VICARMetadata(metadata_dict)
    
    if(hasattr(vicar, 'EOL')):
      if vicar.EOL == 1:
        if vicar.FORMAT == 'BYTE':
          byteCount = 1
        elif vicar.FORMAT == 'HALF':
          byteCount = 2
        elif vicar.FORMAT == 'FULL':
          byteCount = 4
        elif vicar.FORMAT == 'REAL':
          byteCount = 4
        elif vicar.FORMAT == 'DOUB':
          byteCount = 8
        else:
          raise ValueError('Unrecognized Vicar FORMAT: %s in file: %s'%(vicar.FORMAT,metadata_fd.name)) 

        # Read in the VICAR metadata from the end of the file
        metadata_fd.seek(vicar.LBLSIZE + vicar.NLB * vicar.RECSIZE 
          + byteCount*vicar.N1*vicar.N2*vicar.N3)
        # A VICAR file must start with 'LBLSIZE=<integer label size>'.
        lblsize_field = metadata_fd.read(len("LBLSIZE="))
        if lblsize_field.upper() != "LBLSIZE=":
            raise ValueError("Malformed VICAR file: EOL doesn't start with LBLSIZE.")
    
        lblsize = ""
        while True:
            char = metadata_fd.read(1)
            if char == " ":
                break
            else:
                lblsize += char
        try:
            lblsize = int(lblsize)
        except ValueError:
            raise ValueError("Malformed VICAR file: contains non-integer LBLSIZE.")

        metadata_fd.seek(vicar.LBLSIZE + vicar.NLB * vicar.RECSIZE 
          + byteCount*vicar.N1*vicar.N2*vicar.N3)
        metadata = metadata_fd.read(lblsize)
  
        metadata_dict = addMetadataToDict(metadata, metadata_dict)

    metadata_fd.close()
    
    return VICARMetadata(metadata_dict)

def extract_image(vicar, image_fd):
    image_fd.seek(vicar.LBLSIZE + vicar.NLB * vicar.RECSIZE)


    if vicar.FORMAT == 'BYTE':
      outType = numpy.int8
    elif vicar.FORMAT == 'HALF':
      outType = numpy.int16
    elif vicar.FORMAT == 'FULL':
      outType = numpy.int32
    elif vicar.FORMAT == 'REAL':
      outType = numpy.float32
    elif vicar.FORMAT == 'DOUB':
      outType = numpy.float64
    else:
      raise ValueError('Unrecognized Vicar FORMAT: %s in file: %s'%(vicar.FORMAT,image_fd.name)) 
 
    if vicar.ORG != 'BSQ':
      raise ValueError('Vicar ORG: %i is not supported.'%vicar.ORG)
      
    if vicar.NB > 1:
      print 'Reading only the first image of %i images in the file'%vicar.NB
    
    nx = vicar.NS
    ny = vicar.NL
    image = numpy.fromfile(image_fd,dtype=outType,count=nx*ny).reshape(ny,nx)
    return image

def readVicar(infile):
  
    metadata_fd = open(infile, "r")

    vicar_metadata = process_metadata(metadata_fd)

    image_fd = open(infile, "rb")

    image = extract_image(vicar_metadata, image_fd)
    
    return (image,vicar_metadata)
