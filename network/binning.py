"""A module for fitting parameters to a fI curve."""

import time

from math import sqrt
from itertools import takewhile

import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import seaborn as sns

from collections import namedtuple


def grouper(mylist, interval):
    bucket_list = []
    mylist = list(mylist)       # Copy
    T = max(mylist)

    threshold = interval
    while threshold <= T + interval:
        bucket = []
        while len(mylist) > 0:
            if mylist[0] < threshold:
                bucket.append(mylist[0])
                mylist.pop(0)
            else:
                break
        bucket_list.append(bucket)
        threshold += interval
    return bucket_list


def compute_frequency(spike_map, dt):
    """Compute the frequency for each row in the spike_map. Assumes dt in ms.

    Args:
        spike_map (np.ndarray): A NxM matrix. Compte the frequency for each of the N channels.
        dt (float): timestep.

    Returns:
        array of shape (N,) with freqiencies in Hz
    """
    return spike_map.sum(1)/spike_map.shape[1]/dt*1000


def rolling_window(a, window):
    a = list(a)
    sumlist = []
    while len(a) > 0:
        sumlist.append(len(list(takewhile(lambda x: x - a[0] <= window, a))))
        a.pop(0)
    return sumlist


def tsodyks_solver(
        stimulus,
        tau,
        threshold,
        refractory_period,
        icv = 13.5,
        T = 11000,
        tau1 = 3,
        tau_rec = 800,
        tau_f = 1000,
        U = 0.5,
        dt = 0.05,
        num_n = 512,
        R = 1,
        syn_frac = 0.1,
        syn_weight = 1.8,
        noise_scale = 1,
        ie_frac = 0,
        seed = None
):
    """
    Solve the Tsodyks model for either a singlkke neuron or a newtwork.

    Arguments:
        stimulus (float): External stimulus.
        tau (float): Membrane potential time parameter.
        threshold (float): Firing threshold.
        refractory_period (float): Period of quiescence after having fired.
        icv (float): Initial membrane potential.
        T (float): End time of simulation (ms).
        tau1 (float): Synaptic current time parameter.
        tau_rec (float): Synaptic current recovery time parameter.
        U (float): Spike increse in synaptic resource usage.
        dt (float): Time step.
        num_n (int): Number of neurons.
        R (float): synaptic resistance.
        syn_frac (float): Fraction of neurons that are connected.
        syn_weight (float): Synaptic weight.
    """
    np.random.seed(seed)

    # Creeate connectivity matrix
    random_matrix = np.random.random((num_n, num_n))
    A = np.ones(shape=(num_n, num_n))
    num_ex = int(num_n*(1 - ie_frac))
    A[:num_ex, :num_ex] = 1.8*syn_weight    # ee
    A[:num_ex, num_ex:] = 5.4*syn_weight    # ei
    A[num_ex:, :num_ex] = -7.2*syn_weight   # ie
    A[num_ex:, num_ex:] = 0                 # ii

    # Define different parameters for excitatory and inhibitory
    _U = np.zeros(num_n)
    _U[:num_ex] = U
    _U[num_ex:] = 0.4

    _tau_rec = np.zeros(num_n)
    _tau_rec[:num_ex] = tau_rec
    _tau_rec[num_ex:] = 100

    _refractory_period = np.zeros(num_n)
    _refractory_period[:num_ex] = refractory_period
    _refractory_period[num_ex:] = 2

    _tau1 = np.zeros(num_n)
    _tau1[:num_ex] = tau1
    _tau1[num_ex:] = 10

    # Sparse connectivity matrix
    A[random_matrix >= syn_frac] = 0
    np.fill_diagonal(A, 0)

    # Initialise solution arrays
    V_sol = np.ones(num_n)*icv
    x_sol = np.ones_like(V_sol)
    y_sol = np.zeros_like(V_sol)
    z_sol = np.zeros_like(V_sol)
    u_sol = np.ones_like(V_sol)*_U

    # Spike map
    N = int(T/dt + 1)
    # spike_map = np.zeros((N, num_n))
    SpikeTuple = namedtuple("SpikeTuple", ["T", "id"])
    fire_list = []

    # Solution loop
    spike_times = -np.ones(num_n)*2*_refractory_period
    random_component = np.random.random(size=num_n)
    Ib = dt*np.ones(shape=num_n)*R*stimulus

    x_list = np.empty(shape=(N, num_n))
    for i, t in enumerate(range(N)):
        t *= dt     # Scale time to real world (ms)

        # Compute indices of active neurons (That is, not resting)
        refractory_idx = (t - spike_times >= _refractory_period)   # 3 ms refractory time

        # Synaptic currents
        I_syn = A@y_sol

        # Transmembrane potential
        dv = (R*I_syn - V_sol)/tau

        # Synaptic kinetics
        dx = z_sol/_tau_rec
        dy = -y_sol/_tau1
        dz = y_sol/_tau1 - z_sol/_tau_rec

        # Synaptic resources
        du = -u_sol/tau_f

        # Background noise
        Ib_noise = noise_scale*R*np.random.normal(
            0,
            scale=sqrt(dt),
            size=num_n
        )

        # Update solutions
        V_sol += (dt*dv + (Ib + Ib_noise)/tau)*refractory_idx
        x_sol += dt*dx
        y_sol += dt*dy
        z_sol += dt*dz
        # u_sol += dt*du

        # Spiking and not resting
        update_idx = (V_sol > threshold) & refractory_idx
        ex_update = np.where(update_idx[:num_ex])[0]       # excitatory
        in_update = num_ex + np.where(update_idx[num_ex:])[0]       # inhibitory

        y_sol[ex_update] += u_sol[ex_update]*x_sol[ex_update]
        y_sol[in_update] += u_sol[in_update]*(1 - y_sol[in_update])

        x_sol[ex_update] -= u_sol[ex_update]*x_sol[ex_update]
        x_sol[num_ex:] = 1 - y_sol[num_ex:] - z_sol[num_ex:]
        # u_sol[update_idx] += (_U*(1 - u_sol))[update_idx]

        V_sol[update_idx] = icv
        spike_times[update_idx] = t

        # Find indices of 
        for j in np.where(update_idx)[0]:
            fire_list.append(SpikeTuple(T=i, id=j))

        x_list[i] = x_sol
    return fire_list, np.asarray(x_list)



import numpy as np

def nonzero_runs(a):
    # Create an array that is 1 where a is 0, and pad each end with an extra 0.
    iszero = np.concatenate(([0], np.greater(a, 0).view(np.int8), [0]))
    absdiff = np.abs(np.diff(iszero))
    # Runs start and end where absdiff is 1.
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    return ranges


def binning():
    DT = 0.5
    T = 11000
    NUM_N = 500
    tick = time.clock()
    fire_list, x_map = tsodyks_solver(
        tau=41.4,                 # 30
        threshold=12.6,           # 15
        icv=-87.1,                # 13.5
        refractory_period=5,      # 3
        R=42,                     # 1
        dt=DT,
        num_n=NUM_N,
        T=T,
        stimulus=0.0,
        syn_frac=0.1,
        syn_weight=5e-1,
        U=0.5,                   # 0.5
        tau_f=1000,              # 1000
        tau_rec=80,              # 800
        tau1=3,                  # 3
        seed=42,
        noise_scale=1,
        ie_frac=0.2
    )
    tock = time.clock()
    print("Time: ", tock - tick)

    import operator
    id_vector = np.fromiter(map(operator.attrgetter("id"), fire_list), dtype="f4")
    time_array = np.fromiter(map(operator.attrgetter("T"), fire_list), dtype="f4")

    # Reduce the id_array into bins of 5 ms
    bins = list(map(len, grouper(time_array, int(5/DT))))

    res = 1
    # sns.set()
    fig = plt.figure()

    ax = fig.add_subplot(311)
    ax.set_title("Average frequency: {:.2f}".format(len(fire_list)/(T/1000)/NUM_N))
    ax.plot(np.linspace(0, T, len(bins)), list(map(lambda x: x/NUM_N, bins)))

    ax = fig.add_subplot(312)
    ax.plot(np.linspace(0, T, len(x_map)), x_map.sum(1)/NUM_N)

    ax = fig.add_subplot(313)
    ax.plot([0] + list(time_array*DT), [0] + list(id_vector), "x", markersize=0.2)
    fig.savefig("foo.png")

    cumsum = np.asarray(rolling_window(time_array, 5/DT))
    cumsum[cumsum < 20] = 0
    bursts = nonzero_runs(cumsum)
    print("Length of bursts: ", (bursts[:, -1] - bursts[:, 0])*DT)
    print("Burst frequency: ", bursts.shape[0]/T*1000)


if __name__ == "__main__":
    binning()
