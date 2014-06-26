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
    ipos = [ 'config11.dat', 'config21.dat' ]
    fpos = [ 'fpos11.dat', 'fpos21.dat']
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

  """
  plt.suptitle(str(int(conf['npart1']))+' of Type 1, '\
              #+str(int(conf['npart2']))+' of Type 2, '\
              +str(int(conf['nrun']*conf['ncor']))+' iterations, \n'\
              +'Adhesion: '+str(conf['adh1'])+', '\
              +'Contact: '+str(conf['contact'])+', \n'\
              +'Propulsion: '+str(conf['vprop1'])+', '\
              +'Repulsion: '+str(conf['repul1'])+', \n'\
              +str(int(conf['nitn']))+' trials',\
              fontsize=14)
  """

  # Initial positions
  fig.add_subplot(2,3,1,aspect='equal')
  posplot.plotSys(conf, ipos)
  fig.gca().set_title('Initial Configuration')

  # Radial Distribution
  fig.add_subplot(2,3,2)
  grplot.plotGr(conf, gr)
  fig.gca().set_title('Final Average Radial Distribution')

  # Final Positions
  fig.add_subplot(2,3,4,aspect='equal')
  posplot.plotSys(conf, fpos)
  fig.gca().set_title('Final Configuration')
  
  # Mean square displacement
  fig.add_subplot(2,3,5)
  msdplot.plotMsd(conf, msd)
  fig.gca().set_title('Mean Square Displacement')

  #fig.gca().axis([-0.05, 1.05, -0.1, 2.05])
  #plt.legend(loc='lower left', bbox_to_anchor=(1.02, 0), borderaxespad=0)

  # Parameter text
  textstr = 'Parameters:\nIterations: '+str(int(conf['nrun']*conf['ncor']))+\
    '\nTrials: '+str(int(conf['nitn']))+\
    '\nType 1: '+str(int(conf['npart1']))+\
    '\n  Prop: '+str(conf['vprop1'])+\
    '\n  Rep: '+str(conf['repul1'])+\
    '\n  Adh: '+str(conf['adh1'])+\
    '\n\nType 2: '+str(int(conf['npart2']))+\
    '\n  Prop: '+str(conf['vprop2'])+\
    '\n  Rep: '+str(conf['repul2'])+\
    '\n  Adh: '+str(conf['adh2'])

  fig.gca().text(1.3, 2.3, textstr, transform=fig.gca().transAxes, fontsize=14,
        verticalalignment='top')
  
  fig.tight_layout()
  #fig.subplots_adjust(top=0.85)
  
  fig.set_size_inches(11,7)
  plt.savefig('fourplot.png', figsize=(10,1), dpi=100)
  plt.show()
