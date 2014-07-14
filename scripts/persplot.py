"""
  Utility for plotting persistance
  
  Dan Kolbman 2014
"""
import matplotlib.pyplot as plt
import numpy as np
import sys

import dataIO

def plotPers(conf, arg):
  """ plotPers : Dict String[] -> None
  Plots persistance data
  """
  colors = ['#E82C2C', '#245BFF', 'c', 'm']
  for i in range(0,len(arg)):
    t, msd = dataIO.readPers(arg[i])
    # Line
    #plt.loglog(t, msd, color=colors[(i)%3], label=str(i))
    plt.plot(t, msd, color=colors[(i)%3], label=str(i))
    # Dots
    #plt.plot(t, msd, 'o', color=colors[(i)%3], label=str(i))
    # Current axes
    ax = plt.gcf().gca()
    # Linear fit
    #slope,intercept=np.polyfit(t,msd,1)
    # Put fit on graph
    #plt.text(0.1, 0.9-i*0.06,\
    #  'Slope: '+str(int(slope)),\
    #  transform = ax.transAxes)
  # Titles
  plt.gcf().gca().set_title('Persistance')
  plt.xlabel('Lag Time')
  plt.ylabel('Correlation')
  plt.ylim((-1,1))
  plt.axhline(color='black')
  plt.savefig('pers.png')


"""
  If called directly, only plot msd
"""
if __name__ == '__main__':
  if(len(sys.argv) < 3):
    print('Correct usage: python msdplot.py sysparam.dat msd1.dat msd2.dat ...')
  else:
    # Get the configuration variables
    conf = dataIO.readConf(sys.argv[1])
    plotPers(conf, sys.argv[2:])
    plt.show()
  
