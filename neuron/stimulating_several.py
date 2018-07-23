from neuron import h
import matplotlib
matplotlib.use('Agg')
from pylab import *
import mytools
import pickle
import time
import sys
import random

random.seed(1) # Give a seed for the random number generator to produce a fixed distribution of synapses

# Give path to morphologies
morphology_file = "morphologies/cell1.asc"
biophys_file = "models/L5PCbiophys3.hoc"
template_file = "models/L5PCtemplate.hoc"
# Set quantities
v0 = -80
ca0 = 0.0001
proximalpoint = 400 # This probably has something to do with the distance from the soma or something
distalpoint = 620
#fs = 8 # We don't use this...
# Simulation parameters
ITERS = 20
tstop = 11000.0
tstop_s = tstop/1000.      # tstop in seconds
squareDur = tstop          # Duration of the (square) amplitude

# Loading file with information about the synapse locations
unpicklefile = open('synlocs.sav', 'rb')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
Nsyns =  unpickledlist[0]
maxSynsPerSeg = unpickledlist[1]
maxLens = unpickledlist[2]
synlocsAll = unpickledlist[3]

fractions = [0.75, 0.8, 0.85, 0.9, 0.95, 1.00]

icell = 0
morphology_file = "morphologies/cell"+str(icell+1)+".asc"
biophys_file = "models/L5PCbiophys3.hoc"
template_file = "models/L5PCtemplate.hoc"

for fraction in fractions: # Only want one version of the data
    synlocs = synlocsAll[icell]
    frequency = []
    vsoma_pickle = []
    times_pickle = [] 
    #theseCoeffsAll = theseCoeffsAllAll[icell] # Extracting coeffs. of some kind
    # Initial call to h:
    h("""
load_file("stdlib.hoc")
load_file("stdrun.hoc")
objref cvode
cvode = new CVode()
cvode.active(1)
load_file("import3d.hoc")
objref L5PC
load_file(\""""+biophys_file+"""\")
load_file(\""""+template_file+"""\")
L5PC = new L5PCtemplate(\""""+morphology_file+"""\")
access L5PC.soma
objref vsoma,sl,syn1,tvec, syns["""+str(Nsyns)+"""]
vsoma = new Vector()
tvec = new Vector()
sl = new List()
double siteVec[2]
L5PC.soma cvode.record(&v(0.5),vsoma,tvec)
""")
    for istim in range(0,Nsyns):
        h("""
siteVec[0] = """+str(synlocs[istim][0])+"""
siteVec[1] = """+str(synlocs[istim][1])+"""
L5PC.apic[siteVec[0]] {
  syns["""+str(istim)+"""] = new AlphaSynapse(siteVec[1])
  syns["""+str(istim)+"""].e = 0
  syns["""+str(istim)+"""].tau = 5
  syns["""+str(istim)+"""].onset = 10000
}
""")
    # We change the fraction of unmutated channels here, both for the soma and the apical dendrites:
    h('L5PC.soma gNaTa_tbar_NaTa_t = '+str(fraction*2.04))
    h('L5PC.soma gNaTa_tbar_NaTa_modt = '+str((1-fraction)*2.04))
    h('forsec L5PC.apical gNaTa_tbar_NaTa_t = '+str(fraction*0.0213))
    h('forsec L5PC.apical gNaTa_tbar_NaTa_modt = '+str((1-fraction)*0.0213))
    Is = [0.1/Nsyns*x for x in range(0,16)] # Input in one synapse
    Is_tot = [0.1*x for x in range(0,16)]   # Total input current
    spikfreqs = len(Is)*[0]
    for iI in range(0,len(Is)): # Loop over 0.1*x
        squareAmp = Is[iI]        # Increasing the (square) amplitude
        # Loop over g's (or equivalent) somewhere...
        thisCa = h.L5PC.soma[0].minCai_CaDynamics_E2
        # Testing for the kind of cell model we have
        # to determine the conductances
        hasErred = 0
        # Performing the simulation
        # Looping over synapses?
        for istim in range(0,Nsyns): # Stimulating the neuron at several locations at once
            # Feeding more information into h
            #h("syns["+str(istim)+"].gmax = "+str(thisg))
            #access syns["""+str(istim)+"""]
            h("""objref st"""+str(istim)+"""
            st"""+str(istim)+""" = new IClamp(0.0)
            syns["""+str(istim)+"""] st"""+str(istim)+"""
            st"""+str(istim)+""".amp = """+str(squareAmp)+"""
            st"""+str(istim)+""".dur = """+str(squareDur)+"""
            st"""+str(istim)+""".del = 0
            """)
        h("""
        tstop = """+str(tstop)+"""
        cai0_ca_ion = """+str(thisCa)+"""
        v_init = """+str(v0)+"""
        """)
        h.init()             # Initializing the system
        try:                 # Running if we can
            h.run()
        except RuntimeError:
            hasErred = 1
            print("Too large g!") # If we get an error, our conductance is apparently too large
            continue
        
        # Making arrays
        times=np.array(h.tvec)
        Vsoma=np.array(h.vsoma)
        # Testing for spikes
        spikes = mytools.spike_times(times,Vsoma,-35,-45)
        nSpikes1 = len(spikes)
        freq_this = nSpikes1/tstop_s # The frequency
        frequency.append(freq_this)
        # storing one vsoma:
        if iI==10:
            vsoma_pickle = Vsoma
            times_pickle = times
  
    outfilename = 'Plots/fI_Q1481K_fraction'+str(int(round(fraction*100)))+'_percent_long.sav'
    pickle_file = open(outfilename, "wb")
    picklelist = [Is_tot, frequency, times_pickle, vsoma_pickle]
    pickle.dump(picklelist, pickle_file)
    pickle_file.close()
  
  
