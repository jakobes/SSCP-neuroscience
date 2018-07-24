"""A function for creating a spike map from a system o tsodyks neurons."""

import time
import operator

from multiprocessing import Pool
from functools import partial
from fi_tsodyks import tsodyks_solver
from scipy.optimize import minimize

import numpy as np


try:
    from numba import jit
except ModuleNotFoundError:
    print("Could not import numba. Using unity decorator")


    def jit(func, **kwargs):
        """Unity decorator"""
        return func


def objective_function(x, reference):
    """
    Compute the MSE between tsodyks fI curve and the reference one.

    NB! Multiprocess apparently does not like the function signature. Look into repplying numba

    Args:
        x tuple(float, float, float): tau, threshold, refractory_period
        reference (array like): The target fI curve

    Returns float
    """
    tau, threshold, refractory_period = x
    DT = 0.05       # Seems to change based on DT
    T = 11000       # End time

    # Input current
    synaptic_potential = np.arange(16)      # milli Volts = nano Ampere times Mega Ohm
    frequency_list = []
    for stim in synaptic_potential:     # Compute frequency for each stimulus amplitude
        spike_map, _ = tsodyks_solver(
            stimulus=stim,
            tau=tau,
            threshold=threshold,
            refractory_period=refractory_period,
            dt=DT,
            num_n=1
        )
        frequency_list.append(spike_map.sum()/T/1000)   # compute frequency (Hz ~ [1/s])

    # Compute error
    fl = np.array(frequency_list)
    return x, np.power(refractory_period - fl, 2).sum()


def nelderMead(objective, x0=(30.9762, 15.3912, 3.0544), maxiter=150):
    """
    Uses scipy to optimise the objective

    Args:
        objective (callable): The objective function with an unknown gradient.
            Takes x0 as an argument
        maxiter (int): The maximum number of iterations (default = 150).

    Returns scipy result object # TODO: Fix this
    """
    tick = time.time()
    foo = minimize(
        objective,
        x0,
        method="Nelder-Mead",
        tol=1e-6,
        options={
            "disp": True,
            "maxiter": maxiter
        }
    )
    tock = time.time()
    print("Time spent: {}".format(tock - tick))
    return foo


def grid_search(objective, N=9):
    """Serial grid search.

    # TODO: add more parameters for parameter space

    Args:
        objective (callable): The objective function, takes (float, float, float) as only argument.
        N (int): Number of points in the grid.
    returns (np.ndarray, float)
    """
    parameter_space = np.zeros(shape=(N, 3))
    parameter_space[:, 0] = np.random.randint(15, 45, N)
    parameter_space[:, 1] = np.random.randint(5, 20, N)
    parameter_space[:, 2] = np.random.randint(0, 10, N)

    results = {}
    for i, p in enumerate(map(tuple, parameter_space)):
        print("Iteration: ", i)
        if p not in results:
            results[p] = objective(p, reference)

        if i == 0:
            best_key = p
        elif results[p] < results[best_key]:
            best_key = p
            print(best_key, results[best_key])

    k = min(results, key=results.get)
    return k, results[k]


def grid_search_mp(objective, N=300, num_p=4):
    """
    Parallel random grid search.

    # TODO: add more parameters for parameter space

    Args:
        objective (callable): The objective function, takes (float, float, float) as only argument.
        N (int): Number of points in the grid.
        num_p (int): Number of processes.

    returns (np.ndarray, float)
    """
    # Create parameter space
    parameter_space = np.zeros(shape=(N, 3))

    # Use random points, no permutations
    parameter_space[:, 0] = np.random.randint(15, 45, N)
    parameter_space[:, 1] = np.random.randint(5, 20, N)
    parameter_space[:, 2] = np.random.randint(0, 10, N)

    # Run process pool
    processPool = Pool(num_p)
    result = processPool.map(objective, parameter_space)
    return min(result, key=operator.itemgetter(1))

if __name__ == "__main__":
    np.random.seed(42)
    reference = np.load("fi_data.npy")

    result = grid_search_mp(partial(objective_function, reference=reference), N=10)
    # result = nelderMead(
    #     lambda x: partial(jit(cache=True, nopython=True, nogil=True)(objective_function), reference=reference)(x)[-1],  # TODO: Really ugly
    #     x0=(33, 13, 0)
    # )
    print(result)       # 33, 13, 0 is best so far
