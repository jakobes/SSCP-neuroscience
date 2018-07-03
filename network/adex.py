import numpy as np
import random

random.seed(42)

R = 1
U = 1
tau = 30
tau1 = 3
tau_rec = 800
tau_f = 1000

I_b = 13.475 + random.random()*0.05

NT = 4000
time, dt = np.linspace(0, 300, NT, retstep=True)

N = 512
Vsol = -70 + np.random.random(N)*1.2
wsol = np.zeros_like(Vsol)

# A = np.zeros(shape=(N, N))
# A_idx = np.random.random((N, N))
syn_fraq = 0.1
# A[A_idx < syn_fraq] = 1
A = np.random.random((N, N)) < syn_fraq
np.fill_diagonal(A, 0)      # No self-connections

SW = 1
syn_weights = np.random.random(N)*SW/(N*syn_fraq)


def rhs(V, w, EL=-70.6, C=281, gl=30, DT=2, VT=-50.4, a=4, tau_w=144):
    dv = 1./C*(gl*DT*np.exp((V - VT)/DT) - w - gl*(V - EL))
    dw = 1./tau_w*(a*(V - EL) - w)
    return dv, dw


spike_map = np.zeros((NT, N))


for i, t in enumerate(time[:-1]):
    Ib = 1*np.random.random(N)
    I_syn = A@(syn_weights*wsol)

    dv, dw = rhs(Vsol, wsol)
    Vsol = Vsol + dt*dv + Ib*dt + dt*I_syn
    wsol = wsol + dt*dw

    for j, V in enumerate(Vsol):
        if V > 20:
            Vsol[j] = -70.6
            wsol[j] += 0.0805
            spike_map[i, j] = 1


import matplotlib as mpl
mpl.use("Agg")
mpl.rc('text', usetex=True)

import matplotlib.pyplot as plt

spike_avg = spike_map.sum()/spike_map.size

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title(r"Avg spike frequency = %2.f (\(ms^{-1}\))" % spike_avg)
ax.set_xlabel("Time (ms)")
ax.set_ylabel("Neuron number")
ax.imshow(spike_map[::4].T, cmap="binary")

fig.savefig("adex_network.png")
