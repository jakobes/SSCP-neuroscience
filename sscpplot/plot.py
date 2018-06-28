import numpy as np

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import seaborn as sns
sns.set()


def myPlot(time, *yvals, figname="foo.png"):
    """ plot multiple potentials vs time.

    Arguments:
        time: arraylike
        yvals: list of (data, label)

    plot(time, data, label=label)
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Transmembrane Potential (mV)")
    ax.set_title("Cell potentials vs time")

    for v, l in yvals:
        plt.plot(time, v, label=l)

    ax.legend()
    fig.savefig(figname)


if __name__ == "__main__":
    time = np.linspace(0, 1, 11)
    y = 2*time

    myPlot(time, (y, "AAAAAy"))
