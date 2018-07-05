import pickle
import pylab
import matplotlib as mpl
import numpy as np
mpl.use('Agg')
import matplotlib.pyplot as plt

# Frequency has units Hz (from dendrite_stimulation_fI_generator)
# The input current has units nA.

# Unmutated neuron data
fraction = 1.0
test2 = pickle.load(open('Q1481K_fraction'+str(int(round(fraction*100)))+'_percent_f-I_curve.sav', 'r'))

spikefreq_unmut  = test2[0][0]
input_curr_unmut = test2[19]

# Mutated neuron data
fractions = [0.75, 0.8, 0.85, 0.9, 0.95]

for fraction in fractions:
    fr = int(round(fraction*100))
    loadfilename = 'Q1481K_fraction'+str(int(round(fraction*100)))+'_percent_f-I_curve.sav'
    plotfilename = 'Apical_stim_f-I_Q1481K_farend_fraction'+str(int(round(fraction*100)))+'_percent'
    test = pickle.load(open(loadfilename, 'r'))
    spikefreq_mut = test[0][0]
    input_curr_mut = test2[19]
    
    # Plotting 100 % and one mutated in each plot
    plt.figure()
    plt.plot(input_curr_unmut, spikefreq_unmut, label='100%')
    plt.hold('on')
    plt.plot(input_curr_mut, spikefreq_mut,'--', label='95%')
    plt.title(r'Spike freq $f$ vs $I$', fontsize=18)
    plt.xlabel(r'Input Current $I$', fontsize=16)
    plt.ylabel(r'Frequency $f$', fontsize=16)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(plotfilename)

# Plotting all f-I curves in the same plot
plotfilename_all = 'Apical_stim_f-I_Q1481K_farend_all'

plt.figure()
plt.plot(input_curr_unmut, spikefreq_unmut, label='100%')
plt.hold('on')
for fraction in fractions:
    fr = int(round(fraction*100))
    loadfilename = 'Q1481K_fraction'+str(int(round(fraction*100)))+'_percent_f-I_curve.sav'
    
    test = pickle.load(open(loadfilename, 'r'))
    spikefreq_mut = test[0][0]
    input_curr_mut = test2[19]
        
    plt.plot(input_curr_mut, spikefreq_mut,'--', label='%i' % fr + '%')

plt.title(r'Spike freq $f$ vs $I$', fontsize=18)
plt.xlabel(r'Input Current $I$', fontsize=16)
plt.ylabel(r'Frequency $f$', fontsize=16)
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(plotfilename_all)

