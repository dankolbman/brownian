"""
  Utility for plotting mean square displacements
  
  Dan Kolbman 2014
"""
import matplotlib.pyplot as plt
import numpy as np
import sys

import dataIO

def plotMsd(conf, arg):
  """ plotMsd : Dict String[] -> None
  Plots mean square displacement data
  """
  colors = ['r', 'b', 'm', 'c']
  for i in range(0,len(arg)):
    t, msd = dataIO.readMsdave(arg[i])
    # Line
    plt.loglog(t, msd, color=colors[(i)%3], label=str(i))
    # Dots
    #plt.plot(t, msd, 'o', color=colors[(i)%3], label=str(i))
    # Current axes
    ax = plt.gcf().gca()
    # Linear fit
    slope,intercept=np.polyfit(t,msd,1)
    # Put fit on graph
    plt.text(0.1, 0.9-i*0.06,\
      'Slope: '+str(int(slope)),\
      transform = ax.transAxes)
  # Titles
  plt.gcf().gca().set_title('Mean Square Displacement')
  plt.xlabel('Time')
  plt.ylabel('MSD')
  plt.savefig('msd.png')


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
  
