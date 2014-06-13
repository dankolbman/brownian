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
import posplot
import stats

if __name__ == '__main__':
  if(len(sys.argv) < 2):
    print('Correct Usage:\n\tpython fourplot.py sysparam.dat\n\t\
python fourplot.py sysparam.dat initpos.dat gr.dat finalpos.dat msd.dat')
    exit()
  elif(len(sys.argv) == 2):
    ipos = [ 'config11.dat','config21.dat' ]
    fpos = [ 'fpos11.dat', 'fpos21.dat' ]
    gr = [ 'fgr11', 'fgr22' ]
    msd = [ 'msdave1.dat', 'msdave2.dat' ]
  else:
    ipos = [ sys.argv[2] ]
    fpos = [ sys.argv[3] ]
    gr = [ sys.argv[4] ]
    msd = [ sys.argv[5] ]

  # Compute average g(r)
  stats.avgGr(gr[0], 'avggr1.dat')
  stats.avgGr(gr[1], 'avggr2.dat')
  gr = [ 'avggr1.dat', 'avggr2.dat' ]

  conf = dataIO.readConf(sys.argv[1])
  
  fig = plt.figure()
  #fig.subplots_adjust(wspace=0.4,hspace=0.2)

  plt.suptitle(str(int(conf['npart1']))+' of Type 1, '\
              +str(int(conf['npart2']))+' of Type 2, '\
              +str(int(conf['nrun']*conf['ncor']))+' iterations, '\
              +str(int(conf['nitn']))+' trials',\
              fontsize=18)

  # Initial positions
  fig.add_subplot(2,2,1,aspect='equal')
  posplot.plotSys(conf, ipos)
  fig.gca().set_title('Initial Configuration')

  # Radial Distribution
  fig.add_subplot(2,2,2)
  grplot.plotGr(conf, gr)
  fig.gca().set_title('Final Average Radial Distribution')

  # Final Positions
  fig.add_subplot(2,2,3,aspect='equal')
  posplot.plotSys(conf, fpos)
  fig.gca().set_title('Final Configuration')
  
  # Mean square displacement
  fig.add_subplot(2,2,4)
  msdplot.plotMsd(conf, msd)
  fig.gca().set_title('Mean Square Displacement')
  
  fig.tight_layout()
  fig.subplots_adjust(top=0.88)
  
  plt.savefig('circ.png')
  plt.show()
