import pickle
import pylab
import numpy
import sys
from itertools import chain

data = pickle.load(open('control.sav', 'r'))

def myPlot(time, *yvals, figname="foo.png"):
    """ plot multiple potentials vs time.

    Arguments:
        time: arraylike
        yvals: list of (data, label)

    plot(time, data, label=label)
    """
    fig = plt.figue()
    ax = fig.add_subplot(111)

    ax.xlabel("Time (ms)")
    ax.ylabel("Transmembrane Potential (mV)")

    for v in yvals:
        plt.plot(time, yvals)
