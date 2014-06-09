"""
  Statistical package.
  
  Dan Kolbman 2014
"""
import numpy as np
import sys
import os
import fnmatch

import dataIO

def avgGr(prefix, fout):
  """ avgGr
  Average radial distribution functions.
  """
  r = []
  grsum = []
  nfiles = 0
  
  for filen in os.listdir('.'):
    if fnmatch.fnmatch(filen, prefix+'*.dat'):
      r,gr = dataIO.readGr(filen)
      if(grsum != []): np.add(gr,grsum)
      else: grsum = gr
      nfiles += 1


  # Write
  f = open(fout, 'w')
  for i in range(len(grsum)):
    grsum[i] = grsum[i]/nfiles
    # Need to follow g(r) format
    f.write('*\t*\t*\t'+str(r[i]) + '\t' + str(grsum[i]) + '\n')
  f.close()
  

if __name__ == "__main__":
  print('This file only contains statistics functions')

