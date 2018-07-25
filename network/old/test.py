import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt


def solve(v0, tau, I, R, delta, dt, theta, T):
    N = int(T/dt) + 1
    V = np.zeros(N)
    V[0] = v0

    spike_counter = 0
    last_spike = -2*delta
    for i in range(N - 1):
        if V[i] > theta:
            spike_counter += 1
            last_spike = i*dt
            V[i + 1] = v0
        elif i*dt - last_spike > delta:
            V[i + 1] = V[i] + dt*(R*I - V[i])/tau
    return V, spike_counter


if __name__ == "__main__":
    T  = 1000
    DT = 0.01
    V, spike_counter = solve(0, 10, 1.5, 1, 0, DT, 1, T)
    print(spike_counter/T)
    plt.plot(V)
    plt.savefig("foo.png")
