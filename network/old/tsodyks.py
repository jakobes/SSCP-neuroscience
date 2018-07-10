import numpy as np

from scipy.integrate import odeint
import random
from math import sqrt
from tqdm import tqdm

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import math

np.random.seed(42)


class RHS:
    def __init__(self, tau, tau1, tau_rec, tau_f, R, U, dt):
        self.R = R
        self.tau = tau
        self.tau1 = tau1
        self.tau_rec = tau_rec
        self.U = U
        self.dt = dt
        self.tau_f = tau_f

    def __call__(self, y, t, I):
        V, x, y, z, u = y

        dv = (self.R*I - V)/self.tau

        dx = z/self.tau_rec
        dy = -y/self.tau1
        dz = y/self.tau1 - z/self.tau_rec
        du = -u/self.tau_f
        return dv, dx, dy, dz, du

    def delta(self, x, a=3):
        return 1/(a*math.sqrt(math.pi))*np.exp(-(x/a)**2)


T = 400
dt = 1
N = int(T/dt + 1)

num_n = 512
V_sol = -70 + np.random.random(num_n)*1.2       # V13.5?
x_sol = np.ones_like(V_sol)
y_sol = np.zeros_like(V_sol)
z_sol = np.zeros_like(V_sol)
u_sol = np.ones_like(V_sol)

syn_fraq = 0.1
A = np.random.random((num_n, num_n))
A[A >= syn_fraq] = 0
A[A < syn_fraq] = 1.8
np.fill_diagonal(A, 0)      # No self-connections

rhs = RHS(30, 3, 800, 1000, 0.1, 0.5, dt)
V_array = np.zeros((N, 4))
spike_map = np.zeros((N, num_n))


for i, t in tqdm(enumerate(range(N))):
    # Scale = standard deviation
    I_b = 15*np.random.normal(0, scale=sqrt(dt), size=num_n)
    I_syn = A@y_sol

    dV, dx, dy, dz, du = rhs(
        (V_sol, x_sol, y_sol, z_sol, u_sol),
        t,
        I_syn
    )
    V_sol = V_sol + dt*dV + rhs.R*I_b
    x_sol = x_sol + dt*dx
    y_sol = y_sol + dt*dy
    z_sol = z_sol + dt*dz
    u_sol = u_sol + dt*du

    spike_idx = V_sol > 15

    x_sol[spike_idx] -= u_sol[spike_idx]*x_sol[spike_idx]
    y_sol[spike_idx] -= u_sol[spike_idx]*x_sol[spike_idx]
    u_sol[spike_idx] += rhs.U*(1 - u_sol[spike_idx])
    spike_map[i, spike_idx] = 1

    V_array[i] = V_sol[:4]
    # V_sol[spike_idx] = 13.5
    V_sol[spike_idx] = -70.6


spike_avg = spike_map.sum()/spike_map.size
fig = plt.figure(figsize=(10, 10), dpi=93)

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

fig = plt.figure(figsize=(10, 10), dpi=93)

ax = fig.add_subplot(111)
ax.set_title("Four excitatory neurons")
ax.set_xlabel("Time (ms)")
ax.set_ylabel("Membrane potential (mV)")
ax.plot(V_array)

x0, x1 = ax.get_xlim()
y0, y1 = ax.get_ylim()
ax.set_aspect(0.5*abs(x1 - x0)/abs(y1 - y0))

fig.savefig(f"V_tsodyks{dt}.png")

