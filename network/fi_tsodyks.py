"""A module for fitting parameters to a fI curve."""

import time

from math import sqrt

import numpy as np
import matplotlib.pyplot as plt


try:
    from numba import jit
except ModuleNotFoundError:
    print("Could not import numba. Using unity decorator")


    def jit(func, **kwargs):
        """Unity decorator."""
        return func



@jit(cache=True, nopython=True, nogil=True)
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
        dt = 0.1,
        num_n = 512,
        R = 1,
        syn_frac = 0.1,
        syn_weight = 1.8,
        white_noise = True
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
        white_noise (bool): If True, use Euler-Maruyama else Euler.
    """
    # Creeate connectivity matrix
    random_matrix = np.random.random((num_n, num_n))
    A = syn_weight*(0.5 + np.random.uniform(0, 1, size=(num_n, num_n)))

    # Because numba is stupid
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if random_matrix[i, j] >= syn_frac or i == j:
                A[i, j] = 0

    # Use this if not numba
    # np.fill_diagonal(A, 0)      # No self-connections

    # Initialise solution arrays
    V_sol = icv + np.random.random(num_n)*1.2       # V13.5?
    x_sol = np.ones_like(V_sol)
    y_sol = np.zeros_like(V_sol)
    z_sol = np.zeros_like(V_sol)
    u_sol = np.ones_like(V_sol)

    # Spike map
    N = int(T/dt + 1)
    V_array = np.zeros((N, 4))
    spike_map = np.zeros((N, num_n))

    # Solution loop
    spike_times = np.zeros(num_n)
    for i, t in enumerate(range(N)):
        t *= dt     # Scale time to real world (ms)

        # Compute indices of active neurons (That is, not resting)
        refractory_idx = (t - spike_times > refractory_period)   # 3 ms refractory time

        I_syn = A@y_sol     # Synaptic currents

        # Transmembrane potential
        dv = (R*I_syn - V_sol)/tau

        # Synaptic kinetics
        dx = z_sol/tau_rec
        dy = -y_sol/tau1
        dz = y_sol/tau1 - z_sol/tau_rec

        # Synaptic resources
        du = -u_sol/tau_f

        # Background noise
        Ib = np.ones(shape=num_n)*dt*R*stimulus
        if white_noise:
            Ib[:] = R*stimulus*np.random.normal(
                0,
                scale=sqrt(dt),
                size=num_n
            )

        # Update solutions
        V_sol = V_sol + (dt*dv + Ib)*refractory_idx
        x_sol = x_sol + dt*dx
        y_sol = y_sol + dt*dy
        z_sol = z_sol + dt*dz
        u_sol = u_sol + dt*du

        # Spiking and not resting
        update_idx = (V_sol > 15) & refractory_idx
        x_sol[update_idx] -= u_sol[update_idx]*x_sol[update_idx]
        y_sol[update_idx] -= u_sol[update_idx]*x_sol[update_idx]
        u_sol[update_idx] += (U*(1 - u_sol))[update_idx]

        # Again, because numba is stupid
        for j, t in enumerate(update_idx):
            if t:
                spike_map[i, j] = 1

        # "plot" solutions
        V_array[i] = V_sol[:4]
        V_sol[update_idx] = icv     # 13.5
        spike_times[update_idx] = t
    return spike_map, V_array


if __name__ == "__main__":
    DT = 0.05
    synaptic_potential = np.arange(16)      # milli Volts = nano Ampere times Mega Ohm

    threshold = 15
    tau = 30
    refractory_period = 3
    T = 11000

    frequency_list = []
    for stim in synaptic_potential:
        tick = time.clock()
        spike_map, V_array = tsodyks_solver(
            stimulus=stim,
            tau=30,
            threshold=15,
            refractory_period=refractory_period,
            dt = DT,
            num_n = 1
        )
        tock = time.clock()
        frequency_list.append(spike_map.sum()/(T/1000))
        print("Time: ", tock - tick)
    print(frequency_list)

    avg_freq = (spike_map.sum(0)/spike_map.shape[0]).sum()/spike_map.shape[1]
    fig = plt.figure(figsize=(10, 10), dpi=93)

    ax = fig.add_subplot(111)
    ax.set_title(r"Avg spike frequency = %.3f Hz" % (avg_freq*1000))
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Neuron number")
    ax.imshow(spike_map[::1].T, cmap="binary")

    labels = map(lambda x: int(DT*float(x.get_text()[1:-1])), ax.get_xticklabels())
    ax.set_xticklabels(labels)

    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    ax.set_aspect(0.5*abs(x1 - x0)/abs(y1 - y0))

    fig.savefig("network_tsodyks.png")
    import seaborn as sns

    fig = plt.figure(figsize=(10, 10), dpi=93)

    ax = fig.add_subplot(111)
    ax.set_title("Four excitatory neurons")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Membrane potential (mV)")
    ax.plot(V_array)

    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    ax.set_aspect(0.5*abs(x1 - x0)/abs(y1 - y0))

    fig.savefig("V_tsodyks.png")
