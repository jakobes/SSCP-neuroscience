"""Recreating the model from tsodyks 2000."""

from SDE import SDEsolver

import numpy as np


class TsodyksNetwork:
    """ """
    def __init__(
            self,
            tau,
            tau1,
            tau_rec,
            tau_facil,
            R,
            U,
            threshold,
            num_n,
            syn_freq,
            syn_weight
    ):
        """Store parameters and compute connectivity matrix.

        Arguments:
            tau (float): Membrane potential time parameter.
            tau1 (float): synaptic weight time parameter.
            tau_rec (float): Recovery time parameter.
            tau_facil (float): Facilitating time parameter.
            U (float): Spike synaptic resourtce consumption.
            threshold (float): Threshold for classifying spike.
            num_n (int): Number of neurons.
            syn_freq (float): Probability of synnaptic connection.
            syn_weight (float): Weight of synaptic connection.
        """

        self.R = R
        self.tau = tau
        self.tau1 = tau1
        self.tau_rec = tau_rec
        self.tau_facil = tau_facil
        self.U = U
        self.threshold = threshold

        self.A = np.random.random((num_n, num_n))
        self.A[self.A >= syn_freq] = 0
        self.A[self.A < syn_freq] = syn_weight
        np.fill_diagonal(self.A, 0)      # No self-connections

    def a(self, Y):
        """a(Xt)dt.

        Arguments:
            Y (array-like): The solution at the previous time step.
        """
        V, x, y, z, u = Y

        I_syn = self.A@y
        dv = (self.R*I_syn - V)/self.tau

        dx = z/self.tau_rec
        dy = -y/self.tau1
        dz = y/self.tau1 - z/self.tau_rec
        du = -u/self.tau_facil
        return np.asarray((dv, dx, dy, dz, du))

    def b(self, y):
        return self.R*self.threshold


T = 400
dt = 0.1
N = int(T/dt + 1)
num_n = 512

rhs = TsodyksNetwork(30, 3, 800, 1000, 0.1, 0.5, 15, num_n, 0.1, 1.8)


class SpikeMap:
    def __init__(self, U, num_n, N, threshold, reset=-70.6):
        self.U = U
        self.reset = reset
        self.threshold = threshold
        self.spike_map = np.zeros((N, num_n))
        self.V_array = np.zeros((N, 4))
        self.counter = 0

    def __call__(self, Y):
        V, x, y, z, u = Y
        spike_idx = V > self.threshold

        x[spike_idx] -= u[spike_idx]*x[spike_idx]
        y[spike_idx] -= u[spike_idx]*x[spike_idx]
        u[spike_idx] += self.U*(1 - u[spike_idx])

        self.spike_map[self.counter, spike_idx] = 1
        self.counter += 1

        V[spike_idx] = self.reset    # 13.5?
        return V, x, y, z, u


V = -70 + np.random.random(num_n)*1.2       # V13.5?
x = np.ones_like(V)
y = np.zeros_like(V)
z = np.zeros_like(V)
u = np.ones_like(V)
sde_solver = SDEsolver(rhs.a, rhs.b, (V, x, y, z, u))
spike_mapper = SpikeMap(rhs.U, num_n, N, 13.5)

solution = sde_solver.solve((0, T), dt, spike_mapper, True)
spike_map = spike_mapper.spike_map


import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt


spike_avg = spike_map.sum()/spike_map.size
# fig = plt.figure(figsize=(10, 10), dpi=93)
fig = plt.figure()

ax = fig.add_subplot(111)
ax.set_title(r"Avg spike frequency = %2.f (\(ms^{-1}\))" % spike_avg)
ax.set_xlabel("Time (ms)")
ax.set_ylabel("Neuron number")
ax.imshow(spike_map[::1].T, cmap="binary")

x0, x1 = ax.get_xlim()
y0, y1 = ax.get_ylim()
ax.set_aspect(0.5*abs(x1 - x0)/abs(y1 - y0))

fig.savefig("network_tsodyks.png")


import seaborn as sns

# fig = plt.figure(figsize=(10, 10), dpi=93)
fig = plt.figure()

ax = fig.add_subplot(111)
ax.set_title("Four excitatory neurons")
ax.set_xlabel("Time (ms)")
ax.set_ylabel("Membrane potential (mV)")
ax.plot(solution[0, :4])

x0, x1 = ax.get_xlim()
y0, y1 = ax.get_ylim()
ax.set_aspect(0.5*abs(x1 - x0)/abs(y1 - y0))

fig.savefig("V_tsodyks.png")
