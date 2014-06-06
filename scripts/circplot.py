"""
  Plot a multispecies system with circular boundaries.

  Dan Kolbman 2014
"""
import matplotlib.pyplot as plt
import sys

import dataIO

def plotSys(conf, arg):
  """ plotPos : Dict -> None
  Plots position data
  """
  colors = ['r', 'g', 'b']

  for i in range(0,len(arg)):
    xpos, ypos = dataIO.readPos(arg[i])
    s = [conf['diameter']**2/4*3.1415 for i in range(len(xpos))]
    plt.scatter(xpos, ypos, s=s,color=colors[(i)%3])
  
  # Radius of boundary
  rad = conf['circlerad']

  # Draw everything
  #ax = plt.add_subplot(1, 1, 1)
  plt.ylim( (-rad*1.1, rad*1.1) )
  plt.xlim( (-rad*1.1, rad*1.1) )
  circ = plt.Circle((0, 0), radius=conf['circlerad'])
  circ.fill = False
  # Add to current figure
  fig = plt.gcf()
  fig.gca().add_artist(circ)

"""
  If called directly, only show the position plot
"""
if __name__ == '__main__':
  if(len(sys.argv) < 3):
    print('Correct usage: python plotpos.py sysparam.dat pos1.dat pos2.dat ...')
  else:
    # Get the configuration variables
    conf = dataIO.readConf(sys.argv[1])
    plotSys(conf, sys.argv[2:])
    plt.show()
  
