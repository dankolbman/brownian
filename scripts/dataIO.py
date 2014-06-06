"""
  Read and write operations for python scripts.
  
  readPos - Read position data
  readGr - Read g(r) data
  readConf - Read system configuration data
  
  Dan Kolbman 2014
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import sys
import re

def readPos(filen):
  """ readPos : String -> float[] float[]
  Read position data from a file and return x and y lists
  """

  xpos=[]
  ypos=[]
  try:
    f = open(filen)
    for line in f:
      l = line.split()
      if len(l) < 2: break
      xpos.append(float(l[0]))
      ypos.append(float(l[1]))
  except IOError as e:
    print('IO Error!', e.strerror)
  f.close()
  return xpos,ypos

def readGr(filen):
  """ readGr : String -> float[] float[]
  Reads radial distribution function data from file
  Expected file format:
  run iteration time r g(r)
  """
  r=[]
  gr=[]
  try:
    f = open(filen)
    for line in f:
      l = line.split()
      r.append(float(l[3]))
      gr.append(float(l[4]))
  except IOError as e:
    print('IO Error!', e.strerror)
  f.close()
  return r,gr

def readConf(filen):
  """ readConf : String -> Dict
  Reads a system parameter file and returns a dictionary with param val keys
  """
  conf=dict() 
  try:
    f = open(filen)
    for line in f:
      # Search for a line with a variable assignment
      if re.search('[A-Za-z0-9]+([A-Za-z0-9 ]+)?(?= =)',line):
        # Match variable name
        mat = re.search('[A-Za-z0-9]+([A-Za-z0-9 ]+)?(?= =)', line)
        key = mat.group(0)
        # Remove whitespace
        key = "".join(key.split())
        key = key.lower()
        # Match the variable's value
        mat = re.search('(?<![A-Za-z])[0-9\.]{1,12}([E+-]{2}[0-9]{1,4})?', line)
        val = mat.group(0)
        conf[key] = float(val)
  except IOError as e:
    print('Could not process file', e.strerror)
  
  return conf

if __name__ == '__main__':
  print('This is a file containing only helper functions. See source')
