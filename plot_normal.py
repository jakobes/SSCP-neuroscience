import pickle
import pylab
import numpy
import sys
from itertools import chain

data = pickle.load(open('control.sav', 'r'))

# Extract data vectors
data_time = data[14][0]
data_soma = data[15][0]
data_dend = data[16][0]

# Initiate figure
fig = pylab.figure()
ax = fig.add_subplot(111)

ax.plot(data_time, data_soma)
ax.plot(data_time, data_dend)
pylab.show()
