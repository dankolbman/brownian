"""
  Utility for plotting mean square displacements
  
  Dan Kolbman 2014
"""
import matplotlib.pyplot as plt
import sys

import dataIO

def plotMsd(conf, arg):
  """ plotMsd : Dict String[] -> None
  Plots mean square displacement data
  """
  colors = ['r', 'g', 'b']
  for i in range(0,len(arg)):
    r, gr = dataIO.readMsdave(arg[i])
    # Line
    plt.plot(r, gr, color=colors[(i)%3], label=str(i))
    # Dots
    #plt.plot(r, gr, 'o', color=colors[(i)%3], label=str(i))

"""
  If called directly, only plot msd
"""
if __name__ == '__main__':
  if(len(sys.argv) < 3):
    print('Correct usage: python msdplot.py sysparam.dat msd1.dat msd2.dat ...')
  else:
    # Get the configuration variables
    conf = dataIO.readConf(sys.argv[1])
    plotMsd(conf, sys.argv[2:])
    plt.show()
  
