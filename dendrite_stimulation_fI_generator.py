# runcontrols
# A script for determining the control neuron F-I curve and limit cycle.
#
# The input code for the hoc-interface is based on BAC_firing.hoc by Etay Hay (2011)
#
# Tuomo Maki-Marttunen, Jan 2015
# (CC BY)
# And some modifications by Kine :

### About the code
# This program gives the frequency of spikes for several different amplitudes of a constant input current, together with the rest of the output from runcontrols.py.

# We apply a constant input current close to the end of the apical dendrite. The input current lasts for the entirety of the simulation. The input current has units of nA, and the frequency is given in Hz.

# The frequency arrays can be accessed from the output files of this script by
# frequencies = myFile[0][0]
# While the current amplitude arrays can be accessed by
# currents = myFile[19]

from neuron import h
import mytools
import pickle
import numpy as np

# Giving different fractions of unmutated channels to run for:
fractions = [0.75, 0.8, 0.85, 0.9, 0.95, 1.0]

# Filenames and parameters that we don't need to define inside the loop
icell = 0 # Using only one cell model
morphology_file = "morphologies/cell"+str(icell+1)+".asc"
biophys_file = "models/L5PCbiophys3.hoc"
template_file = "models/L5PCtemplate.hoc"
v0 = -80
ca0 = 0.0001

proximalpoint = 400
distalpoint = 620
BACdt = 5.0

tstop = 4000   # End simulation at this time. This is also the duration of the input current. In ms
tstop_s = float(tstop/1000.) # Time in seconds. A float in order to avoid integer division

for fraction in fractions:
    # Resetting lists
    spikfreqsAll = []
    timescAll = []
    VsomacAll = []
    VDerivcAll = []
    VDcoeffAll = []
    VdendcAll = []
    VdDcoeffAll = []
    VdDerivcAll = []
    CasomacAll = []
    CaDerivcAll = []
    CaDcoeffAll = []
    CadendcAll = []
    CadDerivcAll = []
    CadDcoeffAll = []
    times_controlAll = []
    Vsoma_controlAll = []
    Vdend_controlAll = []
    Casoma_controlAll = []
    Cadend_controlAll = []
    # Giving name of output file
    outfilename = 'Plots/Q1481K_fraction'+str(int(round(fraction*100)))+'_percent_f-I_curve.sav'
    # First input to h. st1 is a single pulse current injection, given at the apical dendrite (but we probably have more than one apical dendrite...). These are the lines:
    #access L5PC.apic
    #objref st1
    #st1 = new IClamp(0.5) # The argument is the location of the stimulus. 0.5 is the middle of the selected compartment
    #L5PC.apic st1
    # From online manual, https://www.neuron.yale.edu/neuron/static/docs/help/neuron/stdrun/electrod.html:
    # Electrode: A current injection electrode inserted in the middle of the current section which can be switched between current and voltage clamp modes and can do simple voltage clamp families.
    # IClamp: Switches the Electrode to single pulse current injection. Uses IClamp point process.
    h("""
load_file("stdlib.hoc")
load_file("stdrun.hoc")
objref cvode
cvode = new CVode()
cvode.active(1)
cvode.atol(0.0002)
load_file("import3d.hoc")
objref L5PC
load_file(\""""+biophys_file+"""\")
load_file(\""""+template_file+"""\")
L5PC = new L5PCtemplate(\""""+morphology_file+"""\")
access L5PC.apic
objref st1
st1 = new IClamp(0.9)
L5PC.apic st1
objref vsoma, vdend, recSite, vdend2, isoma, cadend, cadend2, casoma
vsoma = new Vector()
casoma = new Vector()
vdend = new Vector()
cadend = new Vector()
vdend2 = new Vector()
cadend2 = new Vector()
objref sl,ns,tvec
tvec = new Vector()
sl = new List()
double siteVec[2]
sl = L5PC.locateSites("apic","""+str(distalpoint)+""")
maxdiam = 0
for(i=0;i<sl.count();i+=1){
  dd1 = sl.o[i].x[1]
  dd = L5PC.apic[sl.o[i].x[0]].diam(dd1)
  if (dd > maxdiam) {
    j = i
    maxdiam = dd
  }
}
siteVec[0] = sl.o[j].x[0]
siteVec[1] = sl.o[j].x[1]
print "distalpoint gCa_HVA: ", L5PC.apic[siteVec[0]].gCa_HVAbar_Ca_HVA
print "distalpoint gCa_LVA: ", L5PC.apic[siteVec[0]].gCa_LVAstbar_Ca_LVAst
L5PC.apic[siteVec[0]] cvode.record(&v(siteVec[1]),vdend,tvec)
L5PC.apic[siteVec[0]] cvode.record(&cai(siteVec[1]),cadend,tvec)
L5PC.soma cvode.record(&v(0.5),vsoma,tvec)
L5PC.soma cvode.record(&cai(0.5),casoma,tvec)
""")

    # We change the fraction of unmutated channels here, both for the soma and the apical dendrites:
    h('L5PC.soma gNaTa_tbar_NaTa_t = '+str(fraction*2.04))
    h('L5PC.soma gNaTa_tbar_NaTa_modt = '+str((1-fraction)*2.04))
    h('forsec L5PC.apical gNaTa_tbar_NaTa_t = '+str(fraction*0.0213))
    h('forsec L5PC.apical gNaTa_tbar_NaTa_modt = '+str((1-fraction)*0.0213))
  
    Is = [0.1*x for x in range(0,16)]
    spikfreqs = len(Is)*[0]
    for iI in range(0,len(Is)): # Loop over 0.1*iI
        squareAmp = Is[iI]      # Increasing the (square) amplitude
        # The duration of the (square) amplitude is tstop
        h("""
tstop = """+str(tstop)+"""
v_init = """+str(v0)+"""
cai0_ca_ion = """+str(ca0)+"""
st1.amp = """+str(squareAmp)+"""
st1.dur = """+str(tstop)+"""
st1.del = 0
""")                          # Sending in more parameters
        h.init()              # Initiate
        h.run()               # Run

        # Making h-objects into numpy arrays
        times=np.array(h.tvec)
        Vsoma=np.array(h.vsoma)
        Vdend=np.array(h.vdend)
        Casoma=np.array(h.casoma)
        Cadend=np.array(h.cadend)
        spikes = mytools.spike_times(times,Vsoma,-35,100)          # Find the time of the spikes using Tuomo's module
        #print "i=", iI, ": ",len(spikes), " spikes"#: ", spikes
        spikfreqs[iI] = len(spikes)/tstop_s # Frequency in Hz
        
        if abs(Is[iI]-1.0) < 0.0001: # A tolerance test # Will always kick in at i=10, it seems
            #print "if-test satisfied, i =", iI
            Vsoma_control = Vsoma
            Casoma_control = Casoma
            Vdend_control = Vdend
            Cadend_control = Cadend
            times_control = times
            spikes_control = spikes

    #print "length, spikes_control:", len(spikes_control)
    spts = spikes_control[len(spikes_control)-3:len(spikes_control)]
    istart = next((i for i,x in enumerate(times_control) if x > spts[0]))

    iend = next((i for i,x in enumerate(times_control) if x > spts[1]))+4
    nsteps = iend-istart-1
    tdiff = [y-x for x,y in zip(times_control[istart:iend-1],times_control[istart+1:iend])]
    cadiff = [y-x for x,y in zip(Casoma_control[istart:iend-1],Casoma_control[istart+1:iend])]
    caddiff = [y-x for x,y in zip(Cadend_control[istart:iend-1],Cadend_control[istart+1:iend])]
    caderiv1 = [y/x for x,y in zip(tdiff[0:nsteps-1],cadiff[0:nsteps-1])]
    caderiv2 = [y/x for x,y in zip(tdiff[1:nsteps],cadiff[1:nsteps])]
    caderiv = [(x+y)/2.0 for x,y in zip(caderiv1,caderiv2)]
    cadderiv1 = [y/x for x,y in zip(tdiff[0:nsteps-1],caddiff[0:nsteps-1])]
    cadderiv2 = [y/x for x,y in zip(tdiff[1:nsteps],caddiff[1:nsteps])]
    cadderiv = [(x+y)/2.0 for x,y in zip(cadderiv1,cadderiv2)]
    vdiff = [y-x for x,y in zip(Vsoma_control[istart:iend-1],Vsoma_control[istart+1:iend])]
    vddiff = [y-x for x,y in zip(Vdend_control[istart:iend-1],Vdend_control[istart+1:iend])]
    vderiv1 = [y/x for x,y in zip(tdiff[0:nsteps-1],vdiff[0:nsteps-1])]
    vderiv2 = [y/x for x,y in zip(tdiff[1:nsteps],vdiff[1:nsteps])]
    vderiv = [(x+y)/2.0 for x,y in zip(vderiv1,vderiv2)]
    vdderiv1 = [y/x for x,y in zip(tdiff[0:nsteps-1],vddiff[0:nsteps-1])]
    vdderiv2 = [y/x for x,y in zip(tdiff[1:nsteps],vddiff[1:nsteps])]
    vdderiv = [(x+y)/2.0 for x,y in zip(vdderiv1,vdderiv2)]

    Vsomac = Vsoma_control[istart+1:iend-1]
    Vdendc = Vdend_control[istart+1:iend-1]
    Casomac = Casoma_control[istart+1:iend-1]
    Cadendc = Cadend_control[istart+1:iend-1]
    timesc = times_control[istart+1:iend-1]
    VDerivc = vderiv[:]
    VDcoeff =  mytools.limitcyclescaledv(Vsomac,VDerivc,Vsomac,VDerivc)
    VdDerivc = vdderiv[:]
    VdDcoeff =  mytools.limitcyclescaledv(Vdendc,VdDerivc,Vdendc,VdDerivc)
    CaDerivc = caderiv[:]
    CaDcoeff =  mytools.limitcyclescaledv(Casomac,CaDerivc,Casomac,CaDerivc)
    CadDerivc = cadderiv[:]
    CadDcoeff =  mytools.limitcyclescaledv(Cadendc,CadDerivc,Cadendc,CadDerivc)

    spikfreqsAll.append(spikfreqs[:])
    timescAll.append(timesc[:])
    VsomacAll.append(Vsomac[:])
    VDerivcAll.append(VDerivc[:])
    VDcoeffAll.append(VDcoeff)
    VdendcAll.append(Vdendc[:])
    VdDerivcAll.append(VdDerivc[:])
    VdDcoeffAll.append(VdDcoeff)
    CasomacAll.append(Casomac[:])
    CaDerivcAll.append(CaDerivc[:])
    CaDcoeffAll.append(CaDcoeff)
    CadendcAll.append(Cadendc[:])
    CadDerivcAll.append(CadDerivc[:])
    CadDcoeffAll.append(CadDcoeff)
    times_controlAll.append(times_control[:])
    Vsoma_controlAll.append(Vsoma_control[:])
    Vdend_controlAll.append(Vdend_control[:])
    Casoma_controlAll.append(Casoma_control[:])
    Cadend_controlAll.append(Cadend_control[:])

    picklelist = [spikfreqsAll,
    	          timescAll,
    	          VsomacAll,
    	          VDerivcAll,
    	          VDcoeffAll,
    	          VdendcAll,
    	          VdDerivcAll,
    	          VdDcoeffAll,
    	          CasomacAll,
    	          CaDerivcAll,
    	          CaDcoeffAll,
                  CadendcAll,
                  CadDerivcAll,
                  CadDcoeffAll,
                  times_controlAll,
                  Vsoma_controlAll,
                  Vdend_controlAll,
                  Casoma_controlAll,
                  Cadend_controlAll,
                  Is] # Spike frequency: Index 0; Input current: Index 19
    file = open(outfilename, 'w')
    pickle.dump(picklelist,file)
    file.close()


