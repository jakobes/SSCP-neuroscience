# scalemutations.py
# A script for scaling down the effects of genetic variants displayed in the list MT returned by mutation_stuff.getMT()
# For each variant MT[i][j], a coefficient c is looked for which fulfils the following condition: If the variant is 
# scaled with c-eps for small eps > 0, then the scaled variant neuron obeys the conditions I-V (it behaves similarly to
# the control neuron). If the variant is scaled with c+eps for small eps, then the scaled variant neuron does not obey
# some of the conditions I-V (it behaves in a different manner than the control neuron). If c is larger than 1.99, then
# such a factor c was not found on the range [0,2], meaning that the variant need not be downscaled. The resulting scaling
# factors are saved to files 'scalings_cs'+str(icell)+'_'+str(counter)+'.sav', where icell refers to the cell morphology
# and counter refers to the variant number. Note that if variant MT[i][j] reports ranges of different parameter changes,
# then the scaling is performed for each combination of the end points of these ranges.
#
# The input code for the hoc-interface is based on BAC_firing.hoc by Etay Hay (2011)
#
# Tuomo Maki-Marttunen, Jan 2015
# (CC BY)
from neuron import h
import mytools
from pylab import *
import pickle
import sys
import mutation_stuff


# Get the table of variants and the default values for the model parameters
MT = mutation_stuff.getMT()
defVals = mutation_stuff.getdefvals()
keyList = defVals.keys()
for idefval in range(0,len(keyList)):
  if type(defVals[keyList[idefval]]) is not list:
    defVals[keyList[idefval]] = [defVals[keyList[idefval]], defVals[keyList[idefval]]] #make the dictionary values [somatic, apical]
updatedVars = ['somatic','apical','basal'] # the possible classes of segments that defVals may apply to
whichDefVal = [0,1,0]                      # use the defVal[0] for somatic and basal segments and defVal[1] for apical segments


# Get the properties of a control neuron
unpicklefile = open('control.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
spikfreqs_control_All = unpickledlist[0] #List of spiking frequencies
timesc_control_All = unpickledlist[1]    #Times for limit cycle
Vsomac_control_All = unpickledlist[2]    #V for limit cycle
VDerivc_control_All = unpickledlist[3]   #Derivative for limit cycle
VDcoeff_control_All = unpickledlist[4]   #Coefficient for the derivatives
Is_control = unpickledlist[19]           #List of current amplitudes used for calculating spikfreqs_control_All
Is = [0.2,0.4,0.6,0.8,1.0,1.2,1.4]       

theseCoeffsAllAll = [] # Collect here the scaling coefficients

#The somatic current amplitudes [A1,A2,A3] for both cell morphologies in three conditions: Somatic step alone (A1), Synaptic alone (A2), and Somatic+Synaptic (A3):
squareAmps = [[0.696,0,1.137],[0.872,0,0.993]]  

#The synaptic current amplitudes [S1,S2,S3] for both cell morphologies in three conditions: Somatic step alone (S1), Synaptic alone (S2), and Somatic+Synaptic (S3):
epsp_gmaxs = [[0,0.0612,0.100],[0,0.455,0.518]]

for icell in range(0,2): # icell goes through the two cell morphologies (cell #3 wouldn't produce spikes)
  spikfreqs_control = mytools.interpolate(Is_control,spikfreqs_control_All[icell],Is)
  Vsomac_control = Vsomac_control_All[icell]
  VDerivc_control = VDerivc_control_All[icell]
  VDcoeff_control = VDcoeff_control_All[icell]
  timesc_control = timesc_control_All[icell]
  print spikfreqs_control

  morphology_file = "morphologies/cell"+str(icell+1)+".asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"
  v0 = -80
  ca0 = 0.0001
  proximalpoint = 400
  distalpoint = 620
  BACdt = 5.0

  #Initialize the model
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
objref st1
st1 = new IClamp(0.5)
L5PC.soma st1
L5PC.distribute_channels("apic","gIhbar_Ih",2,-0.8696,3.6161,0.0,1.0*2.0870,0.0002)
L5PC.distribute_channels("apic","gCa_HVAbar_Ca_HVA",3,1.0,0.1,685.0,885.0,1.0*0.000555)
L5PC.distribute_channels("apic","gCa_LVAstbar_Ca_LVAst",3,1.0,0.01,685.0,885.0,1.0*0.0187)
objref sl,st2,ns,syn1,con1,isyn, tvec
isyn = new Vector()
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
L5PC.apic[siteVec[0]] st2 = new IClamp(siteVec[1])
st2.amp = 0
L5PC.apic[siteVec[0]] {
  st2
  syn1 = new AlphaSynapse(siteVec[1])
  syn1.e = 0
  syn1.tau = 5
  syn1.onset = 200 + """+str(BACdt)+""" 
  cvode.record(&syn1.i,isyn,tvec)
}
objref vsoma, vdend, recSite, vdend2, isoma, cadend, cadend2, casoma
vsoma = new Vector()
casoma = new Vector()
vdend = new Vector()
cadend = new Vector()
vdend2 = new Vector()
cadend2 = new Vector()
L5PC.soma cvode.record(&v(0.5),vsoma,tvec)
L5PC.soma cvode.record(&cai(0.5),casoma,tvec)
L5PC.apic[siteVec[0]] cvode.record(&v(siteVec[1]),vdend,tvec)
L5PC.apic[siteVec[0]] cvode.record(&cai(siteVec[1]),cadend,tvec)
sl = new List()
sl = L5PC.locateSites("apic","""+str(proximalpoint)+""")
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
L5PC.apic[siteVec[0]] recSite = new IClamp(siteVec[1])
recSite.amp = 0
L5PC.apic[siteVec[0]] {
        recSite
}
L5PC.apic[siteVec[0]] cvode.record(&v(siteVec[1]),vdend2,tvec)
L5PC.apic[siteVec[0]] cvode.record(&cai(siteVec[1]),cadend2,tvec)
L5PC.soma isoma = new Vector()
L5PC.soma cvode.record(&st1.i,isoma,tvec)
""")

  ITERS = 20
  theseCoeffsAll = []  #Save here the downscaling coefficients
  theseMutValsAll = [] #Save here the changed parameter values
  theseMutVarsAll = [] #Save here the names of changed parameters

  counter = -1
  for igene in range(0,len(MT)): # Go through the genes
   theseCoeffsGene = []
   for imut in range(0,len(MT[igene])): # Go through the variants of this gene
    theseCoeffsMut = []
    nVals = len(MT[igene][imut])*[0]
    thesemutvars = []
    for imutvar in range(0,len(MT[igene][imut])): # Go through the parameters changed by the variants
      thesemutvars.append(MT[igene][imut][imutvar][0])
      if type(MT[igene][imut][imutvar][1]) is int or type(MT[igene][imut][imutvar][1]) is float:
        MT[igene][imut][imutvar][1] = [MT[igene][imut][imutvar][1]] # If the type of the entry is not list, make it a list of one entry
      nVals[imutvar] = len(MT[igene][imut][imutvar][1]) # nVals shows how many end points there are for each parameter variation
    cumprodnVals = cumprod(nVals) # take the cumulative product to help determining the combinations of end points
    allmutvars = cumprodnVals[len(MT[igene][imut])-1]*[thesemutvars[:]]
    allmutvals = []
    for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]): # Reserve memory for allmutvals (the values for parameter value change of each combination of end points)
      allmutvals.append([0]*len(thesemutvars))
    for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]): # Set the values for parameter value change of each combination of end points
      for imutvar in range(0,len(MT[igene][imut])): 
        if imutvar==0:
          allmutvals[iallmutval][imutvar] = MT[igene][imut][imutvar][1][iallmutval%nVals[imutvar]]
        else:
          allmutvals[iallmutval][imutvar] = MT[igene][imut][imutvar][1][(iallmutval/cumprodnVals[imutvar-1])%nVals[imutvar]]
    theseMutValsAll.append(allmutvals[:])
    theseMutVarsAll.append(allmutvars[:])
    for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]): # Go through all combinations of end points of the parameter ranges corresponding to this variant
      # Allow parallelization if possible. The variants go from 0 to 51; the index to be considered can be given as the first argument, otherwise, all variants are considered one after another
      counter = counter + 1
      if len(sys.argv) > 1 and int(float(sys.argv[1])) != counter: 
        continue
      
      # Use bisection method to find a threshold c
      nextCoeffs = [0.0,2.0,1.0]    # nextCoeffs[0] shows the lower bound for c, nextCoeffs[1] the upper bound, and nextCoeffs[2] the next one to iterate
      for iter in range(0,ITERS+2): # Iterate ITERS times to find the threshold factor c when the downscaled variant violates some of the conditions I-V
        thisCoeff = nextCoeffs[min(iter,2)] #First try 0.0 and 2.0, and then start iterating
   
        mutText = "" # Save here a compact string presenting the effects of the scaled variant
        # Apply the scaling c for all parameter changes presented by this variant
        for imutvar in range(0,len(MT[igene][imut])): # Go through all model parameters
          if imutvar > 0 and imutvar%2==0:
            mutText = mutText+"\n"
          mutvars = allmutvars[iallmutval][imutvar]
          if type(mutvars) is str:
            mutvars = [mutvars]
          mutText = mutText + str(mutvars) + ": "
          mutvals = allmutvals[iallmutval][imutvar]
          for kmutvar in range(0,len(mutvars)): # Go through the parameters that should be changed at once
            if mutvars[kmutvar].find('offm') > -1 or mutvars[kmutvar].find('offh') > -1 or mutvars[kmutvar].find('ehcn') > -1: # Apply linear change to offset and reversal potentials
              newVal = [x+thisCoeff*mutvals for x in defVals[mutvars[kmutvar]]] # Add the scaled increment to the default value
              if mutvals >= 0 and kmutvar==0:
                mutText = mutText + "+" + str(mutvals*thisCoeff) +" mV"
              elif kmutvar==0:
                mutText = mutText  + str(mutvals*thisCoeff) +" mV"
            else: # Apply logarithmic change to all other types of model parameters
              newVal = [x*(mutvals**thisCoeff) for x in defVals[mutvars[kmutvar]]] # Multiply the default value with a scaled factor
              if kmutvar==0:
                mutText = mutText + "*" + str(mutvals**thisCoeff)
            if kmutvar < len(mutvars)-1:
              mutText = mutText + ", "
            if mutvars[kmutvar].find('_Ih') > -1: # For Ih, apply the changes to both somatic, apical and basal segments
              updateThese = [1,1,1]
            elif mutvars[kmutvar].find('_Ca_HVA') > -1 or mutvars[kmutvar].find('_Ca_LVAst') > -1 or mutvars[kmutvar].find('_SKv3.1') > -1 or mutvars[kmutvar].find('_Ca_HVA') > -1 or mutvars[kmutvar].find('_SK_E2') > -1 or mutvars[kmutvar].find('_NaTa_t') > -1 or mutvars[kmutvar].find('_CaDynamics_E2') > -1: # For these channels, apply the change to somatic and apical segments
              updateThese = [1,1,0]
            elif mutvars[kmutvar].find('_K_Pst') > -1 or mutvars[kmutvar].find('_K_Tst') > -1 or mutvars[kmutvar].find('_Nap_Et2') > -1: # For these channels, apply the change to somatic segments
              updateThese = [1,0,0]
            elif mutvars[kmutvar].find('_Im') > -1: # For these channels, apply the change to apical segments
              updateThese = [0,1,0]
            else:
              print "Error: str=" + str(mutvars[kmutvar])
              updatedVars = [0,0,0]
            for iupdated in range(0,3):
              if updateThese[iupdated]:
                print """forsec L5PC."""+str(updatedVars[iupdated])+""" {
"""+mutvars[kmutvar]+""" = """+str(newVal[whichDefVal[iupdated]])+"""
}"""
                h("""forsec L5PC."""+str(updatedVars[iupdated])+""" {
"""+mutvars[kmutvar]+""" = """+str(newVal[whichDefVal[iupdated]])+"""
}""")
        print mutText
    
        ############################################# Condition 1: Short burst #############################################
        tstop = 500.0
        squareAmp = squareAmps[icell][0]
        squareDur = 150.0
        epsp_gmax = 0.0
        h("""
tstop = """+str(tstop)+"""
v_init = """+str(v0)+"""
cai0_ca_ion = """+str(ca0)+"""
st1.amp = """+str(squareAmp)+"""
st1.del = 200
st1.dur = """+str(squareDur)+"""
syn1.gmax = """+str(epsp_gmax)+"""
syn1.onset = 200 + """+str(BACdt)+""" 
""")
        h.init()
        h.run()

        times=np.array(h.tvec)
        Vsoma=np.array(h.vsoma)
        spikes = mytools.spike_times(times,Vsoma,-35,-45)
        nSpikes1 = len(spikes)

        ############################################# Condition 2: Distal EPSC #############################################
        squareAmp = 0.0
        epsp_gmax = epsp_gmaxs[icell][1]
        h("""
tstop = """+str(tstop)+"""
v_init = """+str(v0)+"""
cai0_ca_ion = """+str(ca0)+"""
st1.amp = """+str(squareAmp)+"""
st1.dur = """+str(squareDur)+"""
syn1.gmax = """+str(epsp_gmax)+"""
syn1.onset = 200 + """+str(BACdt)+""" 
""")
        h.init()
        h.run()
        times=np.array(h.tvec)
        Vsoma=np.array(h.vsoma)
        spikes = mytools.spike_times(times,Vsoma,-35,-45)
        nSpikes2 = len(spikes)

        ############################################# Condition 3: Somatic stim + EPSC #############################################
        squareAmp = squareAmps[icell][2]
        squareDur = 10
        epsp_gmax = epsp_gmaxs[icell][2]
        h("""
tstop = """+str(tstop)+"""
v_init = """+str(v0)+"""
cai0_ca_ion = """+str(ca0)+"""
st1.amp = """+str(squareAmp)+"""
st1.dur = """+str(squareDur)+"""
syn1.gmax = """+str(epsp_gmax)+"""
syn1.onset = 200 + """+str(BACdt)+""" 
""")
        h.init()
        h.run()
        times=np.array(h.tvec)
        Vsoma=np.array(h.vsoma)
        spikes = mytools.spike_times(times,Vsoma,-35,-45)
        nSpikes3 = len(spikes)

        ############################################# Condition 4: IF curve #############################################
        spikfreqs = len(Is)*[0]
        for iI in range(0,len(Is)):
          tstop = 4000.0
          squareAmp = Is[iI]
          squareDur = 3800.0
          epsp_gmax = 0.0
          h("""
tstop = """+str(tstop)+"""
v_init = """+str(v0)+"""
cai0_ca_ion = """+str(ca0)+"""
st1.amp = """+str(squareAmp)+"""
st1.del = 200
st1.dur = """+str(squareDur)+"""
syn1.gmax = """+str(epsp_gmax)+"""
syn1.onset = 200 + """+str(BACdt)+""" 
""")
          h.init()
          h.run()

          times=np.array(h.tvec)
          Vsoma=np.array(h.vsoma)
          spikes = mytools.spike_times(times,Vsoma,-35,100)
          spikfreqs[iI] = sum([1 for x in spikes if x >= 500.0])/3.5
          if iI==4: # use the memb. pot. time course of 1.0nA for the limit cycle
            times_lc = times[:]
            Vsoma_lc = Vsoma[:]
            spikes_lc = spikes[:]

        spikfreqdiffsum = sum([abs(x-y) for x,y in zip(spikfreqs,spikfreqs_control)])
        spikfreqdiffrel = spikfreqdiffsum/sum(spikfreqs_control)

        ############################################# Condition 5: Limit cycle #############################################
        if len(spikes_lc) < 3:
          lcdiff = 1e6
        else:
          spts = spikes_lc[len(spikes_lc)-3:len(spikes_lc)]
          istart = next((i for i,x in enumerate(times_lc) if x > spts[0]))
          iend = next((i for i,x in enumerate(times_lc) if x > spts[1]))+4
          nsteps = iend-istart-1
          Vsomac = Vsoma_lc[istart:iend]
          timesc = times_lc[istart:iend]
          VDerivc = mytools.membpotderivs(timesc,Vsomac)
          VDcoeff =  mytools.limitcyclescaledv(Vsomac,VDerivc,Vsomac,VDerivc)
          lcdiff1 = mytools.limitcyclediff(Vsomac[1:nsteps-1],VDerivc,Vsomac_control,VDerivc_control,VDcoeff_control)
          lcdiff2 = mytools.limitcyclediff(Vsomac_control,VDerivc_control,Vsomac[1:nsteps-1],VDerivc,VDcoeff_control)
          lcdiff = 0.5*(lcdiff1+lcdiff2)

        #Print the parameters and their default values:
        for idefval in range(0,len(defVals.keys())):
          thisdefval = defVals.keys()[idefval]
          if thisdefval.find('_Im') > -1:
            h('print "L5PC.apic[0].'+thisdefval+' = ", L5PC.apic[0].'+thisdefval+', "Default = ", '+str(defVals[thisdefval][1]))
          else:
            h('print "L5PC.soma[0].'+thisdefval+' = ", L5PC.soma[0].'+thisdefval+', "Default = ", '+str(defVals[thisdefval][0]))

        # Determine if the variant obeys (isChanged = 0) or violates (isChanged = 1) the conditions I-V:
        isChanged = nSpikes1 != 4 or nSpikes2 != 1 or nSpikes3 != 2 or spikfreqdiffrel > 0.1 or lcdiff > 600.0
        print isChanged
        if iter==0 and isChanged:
          print "Even null mutation causes different spiking!! igene="+str(igene)+", imut="+str(imut)+", iallmutval="+str(iallmutval)
          continue
        if iter==1 and not isChanged:
          print "This mutation effect does not alter spiking even when doubled!! igene="+str(igene)+", imut="+str(imut)+", iallmutval="+str(iallmutval)
          continue
        if iter>=2:
          if isChanged:
            nextCoeffs = [nextCoeffs[0],nextCoeffs[2],0.5*nextCoeffs[0]+0.5*nextCoeffs[2]]
          else:
            nextCoeffs = [nextCoeffs[2],nextCoeffs[1],0.5*nextCoeffs[1]+0.5*nextCoeffs[2]]

      #Restore default values:
      for imutvar in range(0,len(MT[igene][imut])):
        mutvars = allmutvars[iallmutval][imutvar]
        if type(mutvars) is str:
          mutvars = [mutvars]
        mutvals = allmutvals[iallmutval][imutvar]
        for kmutvar in range(0,len(mutvars)):
          newVal = defVals[mutvars[kmutvar]]
          if mutvars[kmutvar].find('_Ih') > -1:
            updateThese = [1,1,1]
          elif mutvars[kmutvar].find('_Ca_HVA') > -1 or mutvars[kmutvar].find('_Ca_LVAst') > -1 or mutvars[kmutvar].find('_SKv3.1') > -1 or mutvars[kmutvar].find('_Ca_HVA') > -1 or mutvars[kmutvar].find('_SK_E2') > -1 or mutvars[kmutvar].find('_NaTa_t') > -1 or mutvars[kmutvar].find('_CaDynamics_E2') > -1:
            updateThese = [1,1,0]
          elif mutvars[kmutvar].find('_K_Pst') > -1 or mutvars[kmutvar].find('_K_Tst') > -1 or mutvars[kmutvar].find('_Nap_Et2') > -1: 
            updateThese = [1,0,0]
          elif mutvars[kmutvar].find('_Im') > -1:
            updateThese = [0,1,0]
          else:
            print "Error: str=" + str(mutvars[kmutvar])
            updatedVars = [0,0,0]
          for iupdated in range(0,3):
            if updateThese[iupdated]:
              h("""forsec L5PC."""+str(updatedVars[iupdated])+""" {
"""+mutvars[kmutvar]+""" = """+str(newVal[whichDefVal[iupdated]])+"""
}""")
      theseCoeffsMut.append(nextCoeffs[0]+0.0)
      # Save the results to be used by the drawfig1.py and others
      picklelist = [nextCoeffs[0]+0.0,igene,imut,iallmutval,counter,MT]
      file = open('scalings_cs'+str(icell)+'_'+str(counter)+'.sav', 'w')
      pickle.dump(picklelist,file)
      file.close()

    theseCoeffsGene.append(theseCoeffsMut[:])
   theseCoeffsAll.append(theseCoeffsGene[:])
  theseCoeffsAllAll.append(theseCoeffsAll[:])

# After all scalings_csX_Y.sav have been saved, they should be collected to a single file by collectscalings.py

