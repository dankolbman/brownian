"""
  Plot positions of a multispecies system.

  Dan Kolbman 2014
"""
import matplotlib.pyplot as plt
import sys

import dataIO

def plotSys(conf, arg):
  """ plotPos : Dict -> None
  Plots position data
  """
  colors = ['#E82C2C', '#245BFF', 'c', 'm']

  for i in range(0,len(arg)):
    fig = plt.gcf()
    ax = fig.gca()
    xpos, ypos = dataIO.readPos(arg[i])
    circleScatter(xpos, ypos, ax,\
            radius=conf['diameter']/2,\
            color=colors[i%len(colors)])
    #s = [conf['diameter']**2/4*3.1415 for i in range(len(xpos))]
    #plt.scatter(xpos, ypos, s=s,color=colors[(i)%3])

  plotBounds(conf, plt.gcf().gca())
  plt.savefig('fpos.png', figsize=(1,1), dpi=100)


def plotBounds(conf, axes):
  """ plotBounds : Dict Axes -> True
  Makes a rectangle or circle depending on the geometry of the boundary
  """
  # Rectangle
  if(conf['circ'] == 0):
    width = conf['boxwidth']
    height = conf['boxheight']
    plt.ylim( (0-height*0.1, height*1.1) )
    plt.xlim( (0-width*0.1, width*1.1) )
    shape = plt.Rectangle((0,0), width, height)
  else:
    rad = conf['circlerad']
    plt.ylim( (-rad*1.1, rad*1.1) )
    plt.xlim( (-rad*1.1, rad*1.1) )
    shape = plt.Circle((0, 0), radius=conf['circlerad'])

  shape.fill = False
  axes.add_artist(shape)
  return True
  

def circleScatter(xpos, ypos, axes, **kwargs):
  """ circleScatter : float[] float[] -> True
  Creates a scatter plot of circles
  """
  for x,y in zip(xpos, ypos):
    circle = plt.Circle((x,y), **kwargs)
    axes.add_patch(circle)

  return True


"""
  If called directly, only show the position plot
"""
if __name__ == '__main__':
  if(len(sys.argv) < 3):
    print('Correct usage: python posplot.py sysparam.dat pos1.dat pos2.dat ...')
  else:
    # Get the configuration variables
    conf = dataIO.readConf(sys.argv[1])
    plt.gcf().add_subplot(111, aspect='equal')
    plotSys(conf, sys.argv[2:])
    plt.show()
  
