# integrate-and-fire_fit.py
import numpy as np
import random
import matplotlib as mpl
from scipy.optimize import curve_fit
mpl.use("Agg")
mpl.rc('text', usetex=True)

import matplotlib.pyplot as plt
#import pickle

##### Functions #####
def unpack_infile(infilename):
    infile = open(infilename, 'r')
    lines = infile.readlines()
    
    #print(infilename)
    
    Is = []
    fs = []
    for line in lines:
        words = line.split()
        Is.append(float(words[0]))
        fs.append(float(words[1]))
    infile.close()
    
    return Is, fs

def if_freq(I0, tau_m, factor, refr):
    # The frequency of spikes in the integrate- and fire-method
    NI = len(I0)
    answer = np.zeros(NI)
    for i in range(NI):
        #print("i =", i)
        #print("I0[i] =", I0[i])
        logarithm = np.log(I0[i]/(I0[i]-factor))
        answer[i] = 1./(refr+tau_m*logarithm)
        #print(answer[i])
    #logarithm = np.log(R*I0/(R*I0-threshold))
    #divisor = tau_m*logarithm
    #return 1./divisor
    #return 1./(tau_m*np.log(R*I0/(R*I0-threshold)))
    return answer

def get_rms(data, testfunction):
    # Only works if the two arrays are evaluated at the same points.
    Nd = len(data)
    rms = 0
    for i in range(Nd):
        rms += (data[i]-testfunction[i])**2
    rms /= Nd
    return rms

# Want to write the best fit parameters to a file
# Fraction of unmutated channels:
fraction = 0.75
fr = str(int(round(fraction*100)))
# Input and output files
infilename     = 'fI_Q1481K_fraction'+str(int(round(fraction*100)))+'_percent_long_plaintext.txt'
sim_Is, sim_fs = unpack_infile(infilename)
NI             = len(sim_Is)

startat = 2
I0 = sim_Is[startat:]
sf = sim_fs[startat:]
NI2 = len(I0)

# Fitting factors
tau_m = 0.413#0.83
refr = 0
factors = np.linspace(0.19,0.2,5)

plt.figure()
plt.plot(sim_Is, sim_fs, 'o', label='Single neuron model')
for factor in factors:
    test = if_freq(sim_Is[2:], tau_m, factor, refr)
    print("Root mean square deviation, factor =", factor, ": ", get_rms(sf,test))
    plt.plot(sim_Is[2:], test, label='factor = %f' % factor)

plt.title(r"Fits of spike frequency vs input current, $\tau_m=$%f" % tau_m)# % spike_avg)
plt.xlabel("Current [nA]")
plt.ylabel("Spike frequency [Hz]")
plt.legend(loc='best')
plt.tight_layout()
plt.show()
plt.savefig("AtestAtest_factoronly_Integrate_and_fire_gridsearch_single_neuron_f-I.png")


# Fiting tau_m
tau_ms = np.linspace(0.41, 0.415, 5) #0.413#0.83
print("tau_ms:", tau_ms)
refr = 0
factor = 0.1925 #np.linspace(0.17,0.28,5)

plt.figure()
plt.plot(sim_Is, sim_fs, 'o', label='Single neuron model')
for tau_m in tau_ms:
    test = if_freq(sim_Is[2:], tau_m, factor, refr)
    print("Root mean square deviation, tau_m =", tau_m, ": ", get_rms(sf,test))
    plt.plot(sim_Is[2:], test, label=r'$tau_m$ = %f' % tau_m)

plt.title(r"Fits of spike frequency vs input current, factor = %f" % factor)# % spike_avg)
plt.xlabel("Current [nA]")
plt.ylabel("Spike frequency [Hz]")
plt.legend(loc='best')
plt.tight_layout()
plt.show()
plt.savefig("AtestAtest_tau_m_only_Integrate_and_fire_gridsearch_single_neuron_f-I.png")

