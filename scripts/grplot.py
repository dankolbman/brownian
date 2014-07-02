"""
  Utility for plotting radial distribution functions.
  
  Dan Kolbman 2014
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import sys
import re

import dataIO

def plotGr(conf, arg):
  """ plotGr : Dict String[] -> None
  Plots the specified radial distribution files on a single plot
  """
  colors = ['#E82C2C', '#245BFF', 'c', 'm']
  for i in range(0,len(arg)):
    r, gr = dataIO.readGr(arg[i])
    # Line
    plt.plot(r, gr, color=colors[(i)%3], label=str(i))
    # Dots
    #plot(r, gr, 'o', color=colors[(i)%3], label=str(i))

  plt.xlabel('r')
  plt.ylabel('g(r)')
  plt.savefig('gr.png')


"""
  If called directly, only plot g(r)
"""
if __name__ == '__main__':
  if(len(sys.argv) < 3):
    print('Correct usage: python grplot.py sysparam.dat gr1.dat gr2.dat ...')
  else:
    # Get the configuration variables
    conf = dataIO.readConf(sys.argv[1])
    plotGr(conf, sys.argv[2:])
    plt.show()
  
