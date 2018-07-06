"""Module for solving ito stochastic differential eqwuations (SDEs)."""

from math import sqrt

try:
    from tqdm import tqdm
except ModuleNotFoundError:
    tqdm = lambda x: x


import numpy as np


class SDEsolver:
    """Class for solving ito stochastic differential equations (SDEs)."""

    def __init__(self, a, b, x0):
        """Solver an ito SDE on the form 'dXt = a(Xt)dt + b(Xt)dWt' with Euler Maruyama.

        Arguments:
            a (callable): a(Xt)dt
            b (callable): b(Xt)dWt
            x0 (array-like): Initial condition
        """
        self.a = a
        self.b = b
        self.x0 = np.asarray(x0)

    def _step(self, Y, dt):
        """Euler Matuyama."""
        dWt = np.random.normal(0, scale=sqrt(dt), size=self.x0.shape[-1])      # mean 0, std sqrt(dt)
        _Y = Y + self.a(Y)*dt
        _Y[0] += self.b(Y)*dWt
        return _Y

    def solve(self, interval, dt, transform=None, save=False):
        """Solve the equation in `interval` with the time step `dt`.

        Arguments:
            interval (tuple(float, float)): (start_time, stop_time)
            dt (float): time step
            transform (callable): `transform` is called on the solution every timestep
                before it is stored. This enables to solution of e.g. equations with a
                Dirac delta function on the right hand side.
        """
        t_start, t_end = interval

        if transform is None:
            transform = lambda x: x

        N = int((t_end - t_start)/dt + 1)

        Y = self.x0
        if save:
            solution = np.zeros(([N] + list(Y.shape)))
            solution[0] = Y

        for i in tqdm(range(1, N)):
            Y = transform(self._step(Y, dt))
            if save:
                solution[i] = Y

        if save:
            return solution
        return Y
