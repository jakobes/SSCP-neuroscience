import pickle
import pylab
import matplotlib as mpl
import numpy as np
mpl.use('Agg')
import matplotlib.pyplot as plt

# Frequency has units Hz (from dendrite_stimulation_fI_generator)
# The input current has units nA.

fractions = [0.75, 0.8, 0.85, 0.9, 0.95, 1.00]

for fraction in fractions:
    fr = int(round(fraction*100))
    loadfilename = 'fI_Q1481K_fraction'+str(int(round(fraction*100)))+'_percent_long.sav'
    plotfilename = 'Apical_stim_several_f-I_Q1481K_farend_fraction'+str(int(round(fraction*100)))+'_percent'
    test = pickle.load(open(loadfilename, 'r'))
    input_curr = test[0]
    spikefreq = test[1]
    
    plt.figure()
    plt.plot(input_curr, spikefreq)
    plt.title(r'Spike freq $f$ vs $I$, %i' % fr + '%', fontsize=18)
    plt.xlabel(r'Input Current $I$ [nA]', fontsize=16)
    plt.ylabel(r'Frequency $f$ [Hz]', fontsize=16)
    plt.tight_layout()
    plt.savefig(plotfilename)

# Plotting all f-I curves in the same plot

plotfilename_all = 'Apical_stim_several_f-I_Q1481K_farend_all'

plt.figure()
plt.plot(input_curr, spikefreq, label='100%') # This is 100 %
plt.hold('on')
for fraction in fractions[:-1]:
    fr = int(round(fraction*100))
    loadfilename = 'fI_Q1481K_fraction'+str(int(round(fraction*100)))+'_percent_long.sav'
    
    test = pickle.load(open(loadfilename, 'r'))
    input_curr = test[0]
    spikefreq = test[1]

    plt.plot(input_curr, spikefreq,'--', label='%i' % fr + '%')

plt.title(r'Spike freq $f$ vs $I$', fontsize=18)
plt.xlabel(r'Input Current $I$', fontsize=16)
plt.ylabel(r'Frequency $f$', fontsize=16)
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(plotfilename_all)

