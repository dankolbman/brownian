"""
  Plots a 2x2 grid plot with initial positions in top left, final positions
  in bottom left, g(r) in top right, and avg msd in bottom right.
  
  Dan Kolbman 2014
"""
import matplotlib.pyplot as plt
import sys

import dataIO
import grplot
import msdplot
import circplot
import stats

if __name__ == '__main__':
  if(len(sys.argv) < 2):
    print('Correct Usage:\n\tpython fourplot.py sysparam.dat\n\t\
python fourplot.py sysparam.dat initpos.dat gr.dat finalpos.dat msd.dat')
    exit()
  elif(len(sys.argv) == 2):
    ipos = [ 'config11.dat' ]
    fpos = [ 'fpos11.dat', 'fpos21.dat' ]
    gr = [ 'fgr11' ]
    msd = [ 'msdave1.dat' ]
  else:
    ipos = [ sys.argv[2] ]
    fpos = [ sys.argv[3] ]
    gr = [ sys.argv[4] ]
    msd = [ sys.argv[5] ]

  # Compute average g(r)
  stats.avgGr(gr[0], 'avggr.dat')
  gr = [ 'avggr.dat' ]

  conf = dataIO.readConf(sys.argv[1])
  
  fig = plt.figure()

  plt.suptitle(str(int(conf['npart1']))+' of Type 1, '\
              +str(int(conf['npart2']))+' of Type 2',\
              fontsize=18)

  # Initial positions
  fig.add_subplot(2,2,1)
  circplot.plotSys(conf, ipos)
  fig.gca().set_title('Initial Configuration')

  # Radial Distribution
  fig.add_subplot(2,2,2)
  grplot.plotGr(conf, gr)
  fig.gca().set_title('Radial Distribution')

  # Final Positions
  fig.add_subplot(2,2,3)
  circplot.plotSys(conf, fpos)
  fig.gca().set_title('Final Configuration')
  
  # Mean square displacement
  fig.add_subplot(2,2,4)
  msdplot.plotMsd(conf, msd)
  fig.gca().set_title('Mean Square Displacement')
  
  plt.show()
