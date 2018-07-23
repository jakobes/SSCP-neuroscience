import numpy as np

import random
from math import sqrt

import time

try:
    from tqdm import tqdm
except ModuleNotFoundError:
    tqdm = lambda x: x


import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import math

np.random.seed(42)


class RHS:
    def __init__(
            self,
            tau,
            tau1,
            tau_rec,
            tau_f,
            U,
            dt,
            num_n,
            threshold,
            R,
            syn_frac,
            syn_weight,
            random_values = False
    ):
        self.threshold = threshold
        self.R = R

        self.tau = tau
        self.tau1 = tau1
        self.tau_rec = tau_rec
        self.U = U
        self.tau_f = tau_f
        if random_values:
            self.tau = tau/2 + tau*np.random.uniform(0, 1, size=num_n)
            self.tau1 = tau1/2 + tau1*np.random.uniform(0, 1, size=num_n)
            self.tau_rec = tau_rec/2 + tau_rec*np.random.uniform(0, 1, size=num_n)
            self.U = U/2 + U*np.random.uniform(0, 1, size=num_n)
            self.tau_f = tau_f/2 + tau_f*np.random.uniform(0, 1, size=num_n)

        random_matrix = np.random.random((num_n, num_n))
        self.A = syn_weight*(0.5 + np.random.uniform(0, 1, size=(num_n, num_n)))
        self.A[random_matrix >= syn_frac] = 0
        np.fill_diagonal(self.A, 0)      # No self-connections
        self.num_n = num_n
        self.dt = dt

    def a(self, y):
        V, x, y, z, u = y
        I_syn = self.A@y

        dv = (self.R*I_syn - V)/self.tau

        dx = z/self.tau_rec
        dy = -y/self.tau1
        dz = y/self.tau1 - z/self.tau_rec
        du = -u/self.tau_f
        return dv, dx, dy, dz, du

    def b(self):
        return self.R*self.threshold*np.random.normal(0, scale=sqrt(self.dt), size=self.num_n)


T = 4300
dt = 0.1
icv = -70
num_n = 512
threshold = 15

rhs = RHS(30, 3, 800, 1000, 0.5, dt, num_n, threshold, 1, 0.5, 1.8)

spike_times = np.zeros(num_n)

def solve(icv, T, dt, num_n, a, b, threshold, refractory_period=3):
    N = int(T/dt + 1)

    V_sol = icv + np.random.random(num_n)*1.2       # V13.5?
    x_sol = np.ones_like(V_sol)
    y_sol = np.zeros_like(V_sol)
    z_sol = np.zeros_like(V_sol)
    u_sol = np.ones_like(V_sol)

    V_array = np.zeros((N, 4))
    spike_map = np.zeros((N, num_n))

    for i, t in tqdm(enumerate(range(N))):
        t *= dt
        refractory_idx = (t - spike_times > refractory_period)   # 3 ms refractory time
        dV, dx, dy, dz, du = rhs.a((V_sol, x_sol, y_sol, z_sol, u_sol))
        V_sol = V_sol + (dt*dV + rhs.b())*refractory_idx
        x_sol = x_sol + dt*dx
        y_sol = y_sol + dt*dy
        z_sol = z_sol + dt*dz
        u_sol = u_sol + dt*du

        update_idx = (V_sol > 15) & refractory_idx
        x_sol[update_idx] -= u_sol[update_idx]*x_sol[update_idx]
        y_sol[update_idx] -= u_sol[update_idx]*x_sol[update_idx]
        u_sol[update_idx] += (rhs.U*(1 - u_sol))[update_idx]
        # u_sol[update_idx] += rhs.U*(1 - u_sol[update_idx])
        spike_map[i, update_idx] = 1

        V_array[i] = V_sol[:4]
        V_sol[update_idx] = 13.5
        spike_times[update_idx] = t
    return spike_map, V_array


spike_map, V_array = solve(icv, T, dt, num_n, rhs.a, rhs.b, threshold)

avg_freq = (spike_map.sum(0)/spike_map.shape[0]).sum()/spike_map.shape[1]

fig = plt.figure(figsize=(10, 10), dpi=93)

ax = fig.add_subplot(111)
ax.set_title(r"Avg spike frequency = %.3f Hz" % (avg_freq*1000))
ax.set_xlabel("Time (ms)")
ax.set_ylabel("Neuron number")
ax.imshow(spike_map[::1].T, cmap="binary")

labels = map(lambda x: int(dt*float(x.get_text()[1:-1])), ax.get_xticklabels())
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

