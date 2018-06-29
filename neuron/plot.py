import numpy as np

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import seaborn as sns
sns.set()


def myPlot(*vals, indices=None, figname="foo.png"):
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

    for t, v, style, l in vals:
        if indices is not None:
            t = t[indices]
            v = v[indices]
        plt.plot(t, v, style, label=l)

    ax.legend()
    fig.savefig(figname)


if __name__ == "__main__":
    time = np.linspace(0, 1, 11)
    y = 2*time

    myPlot(time, (y, "AAAAAy"))
