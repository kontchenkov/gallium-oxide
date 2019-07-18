import sys
import re
import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib import cm

results_dir = Path(sys.argv[1])
for d in results_dir.glob('*'):
  if not d.is_dir():
    continue
  title = d.name + ' K'
  field = []
  velocity = []
  rates = []
  for f in sorted(d.glob('*.info'), key=lambda x: float(x.name.split('.')[0])):
    field.append(float(f.name.split('.')[0]))
    velocity.append(-float(re.search(r'Average velocity: { (\S*), \S*, \S* }', f.read_text()).group(1)))
    s = re.search(r'Scattering rates: \[(.*)\]', f.read_text()).group(1)
    rates.append([float(x) for x in s.split(',')])
  plt.cla()
  plt.title(title)
  plt.plot(field, velocity, 'k+')
  plt.savefig(str(d) + '_vah.png')
  plt.cla()
  plt.title(title)
  colors = cm.tab20(np.linspace(0, 1, 20))
  plt.stackplot(field, np.array(rates).T, colors=colors)
  plt.savefig(str(d) + '_rates.png')
  
