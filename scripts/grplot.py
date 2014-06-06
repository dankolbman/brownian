"""
  Utility for plotting radial distribution functions
  
  Dan Kolbman 2014
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import sys
import re

import dataIO

def plotGr(conf, arg):
  """ plotGr : Dict -> None
  Plots position data
  """
  fig = plt.figure()
  colors = ['r', 'g', 'b']
  ax = fig.add_subplot(1,1,1)

  for i in range(1,len(arg)):
    r, gr = dataIO.readGr(arg[i])
    # Line
    ax.plot(r, gr, color=colors[(i-2)%3], label=str(i))
    # Dots
    #ax.plot(r, gr, 'o', color=colors[(i-2)%3], label=str(i))

  plt.show()


if __name__ == '__main__':
  if(len(sys.argv) < 3):
    print('Correct usage: grplot sysparam.dat gr1.dat gr2.dat ...')
  else:
    # Get the configuration variables
    conf = dataIO.readConf(sys.argv[1])
    plotGr(conf, sys.argv)
  
