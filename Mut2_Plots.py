
import pickle
import pylab
import numpy
import sys
from itertools import chain
import matplotlib.pyplot as plt

data_list = ['control_mutant_control.sav', 'control_mutant_95.sav', 'control_mutant_90.sav','control_mutant_85.sav', 'control_mutant_80.sav', 'control_mutant_75.sav']
title_list = ['Control', '95% Mutated', '90% Mutated', '85% Mutated', '80% Mutated', '75% Mutated']

fig, axs = plt.subplots(6,1, figsize=(6,15), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = .5, wspace = .05)
axs = axs.ravel()

for i in range(0,6):
    print(i)
    file = data_list[i]
    data_mutant = pickle.load(open(file, 'r'))
    data_time_mutant = data_mutant[14][0]
    data_soma_mutant = data_mutant[15][0]
    data_dend_mutant = data_mutant[16][0]
    axs[i].plot(data_time_mutant, data_soma_mutant, label = 'soma_mutant')
    axs[i].plot(data_time_mutant, data_dend_mutant, label = 'dendrite_mutant')
    axs[i].set_title(title_list[i])
    axs[i].set_ylabel('Voltage (mV)')
    if i == 5:
        axs[i].set_xlabel('Time (mS)')

fig.suptitle('Spiking Behavior of Mutant 2', fontsize=20)
pylab.show()