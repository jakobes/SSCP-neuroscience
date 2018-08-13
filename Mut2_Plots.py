
import pickle
import pylab
import numpy
import sys
from itertools import chain
import matplotlib.pyplot as plt

data_list = ['control_mutant_control.sav', 'control_mutant_95.sav', 'control_mutant_90.sav','control_mutant_85.sav', 'control_mutant_80.sav', 'control_mutant_75.sav']
title_list = ['Control', '95%', '90%', '85%', '80%', '75%']

fig, axs = plt.subplots(6,1, figsize=(7,15), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = .5, wspace = .05)
axs = axs.ravel()

for i in range(0,6):
    print(i)
    file = data_list[i]
    data_mutant = pickle.load(open(file, 'r'))
    data_time_mutant = data_mutant[14][0]
    data_soma_mutant = data_mutant[15][0]
    data_dend_mutant = data_mutant[16][0]
    axs[i].plot(data_time_mutant, data_soma_mutant, label = 'soma_mutant', lw=4, color="blue", alpha=0.9)
    axs[i].plot(data_time_mutant, data_dend_mutant, label = 'dendrite_mutant', lw=4, color="red", alpha=0.9)
    axs[i].set_title(title_list[i])
    #axs[i].set_ylabel('Voltage (mV)')
	# Remove the plot frame lines. They are unnecessary chartjunk.  
    axs[i].spines["top"].set_visible(False)    
    axs[i].spines["bottom"].set_visible(False)    
    axs[i].spines["right"].set_visible(False)    
    axs[i].spines["left"].set_visible(False)
	# Ensure that the axis ticks only show up on the bottom and left of the plot.    
	# Ticks on the right and top of the plot are generally unnecessary chartjunk.    
    axs[i].get_xaxis().tick_bottom()    
    axs[i].get_yaxis().tick_left()  
	# Limit the range of the plot to only where the data is.    
	# Avoid unnecessary whitespace.    
    plt.ylim(-81, 41)    
    plt.xlim(0, 401)
	# Remove the tick marks;
    #plt.ytick_params(axis="both", which="both", bottom="off", top="off", labelbottom="on", left="off", right="off", labelleft="on")  
	# Provide tick lines across the plot to help your viewers trace along    
# the axis ticks. Make sure that the lines are light and small so they    
# don't obscure the primary data lines. 
    #axs[i].yticks(range(-80, 41, 40))
    #plt.xticks(range(0, 100, 400))
    for y in range(-80, 41, 40):    
        axs[i].plot(range(0, 401), [y] * len(range(0, 401)), "--", lw=0.5, color="black", alpha=0.4) 
    #if i == 5:
    #    axs[i].set_xlabel('Time (mS)')
    

#fig.text(0.5, 0.04, 'Time (mS)', ha='center', fontsize=20)
#fig.text(0.04, 0.5, 'Voltage (mV)', va='center', rotation='vertical', fontsize=20)

fig.suptitle('Spiking Behavior at Different Proportions of Mutated Channels', fontsize=20)

#fig.tight_layout()
pylab.show()