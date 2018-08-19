"""A module for fitting parameters to a fI curve."""

import time
import operator

from math import sqrt
from itertools import takewhile

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from collections import namedtuple


def num_spikes(fire_list, dt):
    burst_dict = {}
    for fire in map(operator.attrgetter("T"), fire_list):
        if fire in burst_dict:
            burst_dict[fire] += 1
        else:
            burst_dict[fire] = 1

    xvals = np.fromiter((i*dt for i in sorted(burst_dict.keys())), dtype="f8")
    yvals = np.fromiter((burst_dict[i] for i in burst_dict.keys()), dtype="f8")
    return xvals, yvals


def nonzero_runs(a):
    # Create an array that is 1 where a is 0, and pad each end with an extra 0.
    iszero = np.concatenate(([0], np.greater(a, 0).view(np.int8), [0]))
    absdiff = np.abs(np.diff(iszero))
    # Runs start and end where absdiff is 1.
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    return ranges


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


    a = list(a)
    sumlist = []
    while len(a) > 0:
        sumlist.append(len(list(takewhile(lambda x: abs(x - a[0]) <= window, a))))
        a.pop(0)
    return sumlist


def find_bursts(fire_list, interval, cutoff):
    burst_dict = {0: []}

    prev_spike = 0
    burst_id = 0
    for fire in map(operator.attrgetter("T"), fire_list):
        if fire - prev_spike < interval:
            burst_dict[burst_id].append(fire)
        else:
            burst_id += 1
            burst_dict[burst_id] = []
        prev_spike = fire

    new_dict = {}
    for k, v in burst_dict.items():
        if len(v) > cutoff:
            new_dict[k] = (min(v), max(v))
    return new_dict


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
        ie_frac = 0.2,
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
    A[:num_ex, num_ex:] = -1.8*syn_weight    # ie -- 7.2
    A[num_ex:, :num_ex] = 2.4*syn_weight   # ei -- 5.4
    A[num_ex:, num_ex:] = -0*syn_weight                 # ii -- -2.2

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


def binning(tau, threshold, icv, refractory_period, R, name=None, syn_weight=5e-1, seed=42):
    DT = 0.05
    T = 11000
    NUM_N = 500
    tick = time.clock()
    fire_list, x_map = tsodyks_solver(
        tau=tau,                 # 30
        threshold=12.6,           # 15
        icv=-87.1,                # 13.5
        refractory_period=5,      # 3
        R=42,                     # 1
        dt=DT,
        num_n=NUM_N,
        T=T,
        stimulus=0.0,
        syn_frac=0.1,
        syn_weight=syn_weight,
        U=0.5,                   # 0.5
        tau_f=1000,              # 1000
        tau_rec=80,              # 800
        tau1=3,                  # 3
        seed=seed,
        noise_scale=1,
        ie_frac=0.2,     # 0.4
    )
    tock = time.clock()
    print("Time: ", tock - tick)

    id_vector = np.fromiter(map(operator.attrgetter("id"), fire_list), dtype="f4")
    time_array = np.fromiter(map(operator.attrgetter("T"), fire_list), dtype="f4")

    # Reduce the id_array into bins of 5 ms
    bins = list(map(len, grouper(time_array, int(5/DT))))

    sns.set()
    # fig = plt.figure(figsize=(20, 20))
    # fig = plt.figure()
    fig = plt.figure(figsize=(11,7))

    # ax = fig.add_subplot(111)
    # ax.set_title("Average frequency: {:.2f}".format(len(id_vector)/(T/1000)/NUM_N), fontsize=48)
    # ax.plot(np.linspace(0, T, len(bins)), list(map(lambda x: x/(NUM_N), bins)), linewidth=2)
    # ax.tick_params(axis='both', which='major', labelsize=26)
    # ax.tick_params(axis='both', which='minor', labelsize=24)

    x, y = [0] + list(time_array*DT), [0] + list(id_vector)
    # np.savetxt("foo.txt", (x, y))
    # ax = fig.add_subplot(111)
    # ax.plot([0] + list(time_array*DT), [0] + list(id_vector), "x", markersize=0.6, mew=0.8)
    # ax.tick_params(axis='both', which='major', labelsize=26)
    # ax.tick_params(axis='both', which='minor', labelsize=24)

    cumulative_bursts = find_bursts(fire_list, 5/DT, 10)
    z = np.fromiter(map(lambda x: sum(x)/2*DT, cumulative_bursts.values()), "f8")
    # import IPython
    # IPython.embed()
    # ax.plot(z, -5*np.ones_like(z), "ro", markersize=10, mew=1)
    # fig.tight_layout()
    # print(z)

    ax = fig.add_subplot(111)
    ax.plot(x, y, "x", markersize=0.6, mew=0.8)
    ax.tick_params(axis='both', which='major', labelsize=28)
    ax.tick_params(axis='both', which='minor', labelsize=26)

    ax.set_xlabel("time [ms]", fontsize=28)
    ax.set_ylabel("neuron id", fontsize=28)
    ax.set_title("Nework activity", fontsize=28)

    ax.plot(z, -5*np.ones_like(z), "ro", markersize=5, mew=1)
    fig.tight_layout()
    fig.savefig("network.png")


    if name is None:
        name = "foo"
    # fig.savefig("images/{}.pdf".format(name))
    fig.savefig("images/{}.png".format(name))
    plt.close(fig)

    BurstTuple = namedtuple("BurstTuple", ["durations", "frequency"])
    return BurstTuple(
        list(map(lambda x: (x[1] - x[0])*DT, cumulative_bursts.values())),
        len(cumulative_bursts)
    )


if __name__ == "__main__":
    ParameterTuple = namedtuple(
        "ParameterTuple",
        ["tau", "threshold", "icv", "refractory_period", "R"]
    )
    spec_dict = {
        "mut1_100": ParameterTuple(tau=43.6, threshold=20.1, icv=-98.2, refractory_period=5, R=66.8),
        "mut1_95": ParameterTuple(tau=37.4, threshold=11.6, icv=-80.1, refractory_period=5, R=38.7),
        "mut1_90": ParameterTuple(tau=48.9, threshold=14.6, icv=-81.7, refractory_period=5, R=48.8),
        "mut1_85": ParameterTuple(tau=43.8, threshold=10.9, icv=-83.0, refractory_period=5, R=36.4),
        # "mut1_80": ParameterTuple(tau=41.9, threshold=9.0, icv=-85.6, refractory_period=5, R=32.3),
        # "mut1_75": ParameterTuple(tau=42.6, threshold=8.0, icv=-81.7, refractory_period=5, R=29.7),

        "mut2_100": ParameterTuple(tau=43.6, threshold=20.1, icv=-98.2, refractory_period=5, R=66.8),
        "mut2_95": ParameterTuple(tau=46.9, threshold=21.2, icv=-99.8, refractory_period=5, R=70.5),
        "mut2_90": ParameterTuple(tau=52.9, threshold=15.0, icv=-78.7, refractory_period=5, R=49.8),
        "mut2_85": ParameterTuple(tau=60.3, threshold=20.8, icv=-84.6, refractory_period=5, R=69.2),
        "mut2_80": ParameterTuple(tau=61.4, threshold=20.6, icv=-88.3, refractory_period=5, R=68.7),
        # "mut2_80": ParameterTuple(tau=61.4, threshold=20.6, icv=-88.3, refractory_period=5, R=68.7),
    }

    for seed in np.random.randint(0, 100, 3):
        mutant_dict = {}
        for fraq, spec in spec_dict.items():
            print("Computing fraq: ", fraq)
            mutant_dict[fraq] = binning(
                tau=spec.tau,
                threshold=spec.threshold,
                icv=spec.icv,
                refractory_period=spec.refractory_period,
                R=spec.R,
                name=fraq,
                syn_weight=10,
                seed=seed
            )
            break

        print("seed: ", seed)
        for k, v in mutant_dict.items():
            print("{}: ({}, {})".format(k, v.frequency, v.durations))
        break   # Seems seed does nothing
