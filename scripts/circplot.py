import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import sys
import re

def readPos(filen):
  """ readPos : String -> float[] float[]
  Read position data from a file and return x and y lists
  """

  xpos=[]
  ypos=[]
  try:
    f = open(filen)
    for line in f:
      l = line.split()
      if len(l) < 2: break
      xpos.append(float(l[0]))
      ypos.append(float(l[1]))
  except IOError as e:
    print('IO Error!', e.strerror)
  f.close()
  return xpos,ypos

def readConf(filen):
  """ readConf : String -> Dict
  Reads a system parameter file and returns a dictionary with param val keys
  """
  conf=dict() 
  try:
    f = open(filen)
    for line in f:
      if re.search('[A-Za-z0-9]+([A-Za-z0-9 ]+)?(?= =)',line):
        #print(re.search('[A-Za-z]+[A-Za-z ]*[=]{1}[ ]*[0-9\.E+-]*',line)).group(0)
        #abc = re.search('[A-Za-z]+[A-Za-z ]', line)
        # Variable name
        abc = re.search('[A-Za-z0-9]+([A-Za-z0-9 ]+)?(?= =)', line)
        key = abc.group(0)
        # Remove whitespace
        key = "".join(key.split())
        key = key.lower()
        print(key)
        #abc = re.search('[0-9\.]+[E+-]{2}[0-9]+', line)
        # Variable value
        abc = re.search('(?<![A-Za-z])[0-9\.]{1,12}([E+-]{2}[0-9]{1,4})?', line)
        val = abc.group(0)
        #print(val)
        conf[key] = float(val)
  except IOError as e:
    print('Could not process file', e.strerror)
  
  return conf

def plotSys(conf, arg):
  """ plotPos : Dict -> None
  Plots position data
  """
  fig = plt.figure()

  colors = ['r', 'g', 'b']

  for i in range(2,len(arg)):
    xpos, ypos = readPos(arg[i])
    s = [conf['diameter']**2/4*3.1415 for i in range(len(xpos))]
    plt.scatter(xpos, ypos, s=s,color=colors[(i-2)%3])


  # Radius of boundary
  rad = conf['circlerad']

  # Draw everything
  ax = fig.add_subplot(1, 1, 1)
  plt.ylim( (-rad*1.1, rad*1.1) )
  plt.xlim( (-rad*1.1, rad*1.1) )
  circ = plt.Circle((0, 0), radius=conf['circlerad'])
  circ.fill = False
  ax.add_patch(circ)

  plt.show()


if __name__ == '__main__':
  if(len(sys.argv) < 3):
    print('Correct usage: plotpos sysparam.dat pos1.dat pos2.da ...')
  else:
    # Get the configuration variables
    conf = readConf(sys.argv[1])
    plotSys(conf, sys.argv)
  
