# findthresholddistalamps.py
# A script for determining the threshold conductance for eliciting a spike with a stimulus that
# is widely distributed along the apical dendrite
#
# The input code for the hoc-interface is based on BAC_firing.hoc by Etay Hay (2011)
#
# Tuomo Maki-Marttunen, Jan 2015
# (CC BY)
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

morphology_file = "morphologies/cell1.asc"
biophys_file = "models/L5PCbiophys3.hoc"
template_file = "models/L5PCtemplate.hoc"
v0 = -80
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620
fs = 8
ITERS = 20
tstop = 11000.0

unpicklefile = open('synlocs.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
Nsyns =  unpickledlist[0]
maxSynsPerSeg = unpickledlist[1]
maxLens = unpickledlist[2]
synlocsAll = unpickledlist[3]

import mutation_stuff
MT = mutation_stuff.getMT()
defVals = mutation_stuff.getdefvals()
keyList = defVals.keys()
for idefval in range(0,len(keyList)):
  if type(defVals[keyList[idefval]]) is not list:
    defVals[keyList[idefval]] = [defVals[keyList[idefval]], defVals[keyList[idefval]]] #make the dictionary values [somatic, apical]
updatedVars = ['somatic','apical','basal'] # the possible classes of segments that defVals may apply to
whichDefVal = [0,1,0]                      # use the defVal[0] for somatic and basal segments and defVal[1] for apical segments
unpicklefile = open('scalings.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
theseCoeffsAllAll = unpickledlist[0]
theseMutValsAllAll = unpickledlist[2]

gsAllAll = []

for icell in range(0,2):
  synlocs = synlocsAll[icell]
  gsAll = []
  morphology_file = "morphologies/cell"+str(icell+1)+".asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"

  theseCoeffsAll = theseCoeffsAllAll[icell]
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
""") #"""

  coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]

  counter = -1
  
  for igene in range(0,len(MT)):
   gsThisGene = []
   for imut in range(0,len(MT[igene])):
    gsThisMut = []
    nVals = len(MT[igene][imut])*[0]
    thesemutvars = []
    theseCoeffs = theseCoeffsAll[igene][imut]
    for imutvar in range(0,len(MT[igene][imut])):
      thesemutvars.append(MT[igene][imut][imutvar][0])
      if type(MT[igene][imut][imutvar][1]) is int or type(MT[igene][imut][imutvar][1]) is float:
        MT[igene][imut][imutvar][1] = [MT[igene][imut][imutvar][1]]
      nVals[imutvar] = len(MT[igene][imut][imutvar][1])
    cumprodnVals = cumprod(nVals)
    allmutvars = cumprodnVals[len(MT[igene][imut])-1]*[thesemutvars]
    allmutvals = []
    for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
      allmutvals.append([0]*len(thesemutvars))
    for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
      for imutvar in range(0,len(MT[igene][imut])):
        if imutvar==0:
          allmutvals[iallmutval][imutvar] = MT[igene][imut][imutvar][1][iallmutval%nVals[imutvar]]
        else:
          allmutvals[iallmutval][imutvar] = MT[igene][imut][imutvar][1][(iallmutval/cumprodnVals[imutvar-1])%nVals[imutvar]]
      
    for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
      counter = counter + 1                                                                                                                                                               
      if len(sys.argv) > 1 and int(float(sys.argv[1])) != counter:
        continue
      gsThisMutVal = []
      close("all")
      f, axarr = plt.subplots(2, 2)
      maxCac = 0
      maxCadc = 0
      for iter in [0, 2, 5, 6, 8, -1]:
        gsThisIter = []
        if iter >= 0:
          thisCoeff = coeffCoeffs[iter][0]*theseCoeffs[iallmutval] + coeffCoeffs[iter][1]*(1.0 - 0.5*theseCoeffs[iallmutval])
        else:
          thisCoeff = 0
        if iter == -1 and (igene > 0):
          continue # do the control only once!
        print "iter="+str(iter)+", thisCoeff="+str(thisCoeff)
          
        mutText = ""
        for imutvar in range(0,len(MT[igene][imut])):
          if imutvar > 0 and imutvar%2==0:
            mutText = mutText+"\n"
          mutvars = allmutvars[iallmutval][imutvar]
          mutvals = allmutvals[iallmutval][imutvar]
          if type(mutvars) is str:
            mutvars = [mutvars]
          mutText = mutText + str(mutvars) + ": "
          for kmutvar in range(0,len(mutvars)):
            mutvar = mutvars[kmutvar]
            if mutvar.find('offm') > -1 or mutvar.find('offh') > -1 or mutvar.find('ehcn') > -1:
              newVal =  [x+mutvals*thisCoeff for x in defVals[mutvar]]
              if mutvals >= 0 and kmutvar==0:
                mutText = mutText + "+" + str(mutvals) +" mV"
              elif kmutvar==0:
                mutText = mutText  + str(mutvals) +" mV"
            else:
              newVal = [x*(mutvals**thisCoeff) for x in defVals[mutvar]]
              if kmutvar==0:
                mutText = mutText + "*" + str(mutvals)
            if kmutvar < len(mutvars)-1:
              mutText = mutText + ", "
            if mutvar.find('_Ih') > -1:
              updateThese = [1,1,1]
            elif mutvar.find('_Ca_HVA') > -1 or mutvar.find('_Ca_LVAst') > -1 or mutvar.find('_SKv3.1') > -1 or mutvar.find('_Ca_HVA') > -1 or mutvar.find('_SK_E2') > -1 or mutvar.find('_NaTa_t') > -1 or mutvar.find('_CaDynamics_E2') > -1:
              updateThese = [1,1,0]
            elif mutvar.find('_K_Pst') > -1 or mutvar.find('_K_Tst') > -1 or mutvar.find('_Nap_Et2') > -1:
              updateThese = [1,0,0]
            elif mutvar.find('_Im') > -1:
              updateThese = [0,1,0]
            else:
              print "Error: str=" + str(mutvar)
              updatedThese = [0,0,0]
            for iupdated in range(0,3):
              if updateThese[iupdated]:
                print """forsec L5PC."""+str(updatedVars[iupdated])+""" {
"""+mutvar+""" = """+str(newVal[whichDefVal[iupdated]])+"""
}"""
                h("""forsec L5PC."""+str(updatedVars[iupdated])+""" {
"""+mutvar+""" = """+str(newVal[whichDefVal[iupdated]])+"""
}""")
        print mutText
        thisCa = h.L5PC.soma[0].minCai_CaDynamics_E2
        if icell==0:
          nextgs = [0.00,0.003,0.0015]
        if icell==1:
          nextgs = [0.00,0.06,0.03]
        hasSpiked = 0
        hasErred = 0
        for iterg in range(0,ITERS+2):
            thisg = nextgs[min(iterg,2)]
            for istim in range(0,Nsyns):
              h("syns["+str(istim)+"].gmax = "+str(thisg))
            h("""
tstop = """+str(tstop)+"""
cai0_ca_ion = """+str(thisCa)+"""
v_init = """+str(v0)+"""
""")
            h.init()
            try:
              h.run()
            except RuntimeError:
              hasErred = 1
              print "Too large g!"
              if iterg == 1:
                nextgs = [0.0,4.0,3.0]
                continue
              else:
                nextgs = [nextgs[0],nextgs[2],nextgs[0]+nextgs[2]]
              continue
  
            times=np.array(h.tvec)
            Vsoma=np.array(h.vsoma)
            spikes = mytools.spike_times(times,Vsoma,-35,-45)
            nSpikes1 = len(spikes)
            hasSpiked = hasSpiked or (nSpikes1 > 0)
  
            print "iterg="+str(iterg)+" done, g="+str(thisg)+", "+str(nSpikes1)+" spikes"
            if iterg==0 and nSpikes1 > 0:
              print "Even zero g causes spiking!! igene="+str(igene)+", imut="+str(imut)+", iallmutval="+str(iallmutval)+", iter="+str(iter)
              nextgs = [0.0,0.0,0.0]
              break
            if iterg==1 and not hasSpiked:
              print "No spiking with iterg==1, adding 25% to the current! igene="+str(igene)+", imut="+str(imut)+", iallmutval="+str(iallmutval)
              nextgs = [nextgs[0],2.0*nextgs[1],1.25*nextgs[min(iterg,2)]]
              continue

            if iterg>=2 and iterg < ITERS+2:
              if nSpikes1 > 0:
                nextgs = [nextgs[0],nextgs[2],0.5*nextgs[0]+0.5*nextgs[2]]
              else:
                nextgs = [nextgs[2],nextgs[1],0.5*nextgs[1]+0.5*nextgs[2]]

        #Print the parameters and their default values:
        for idefval in range(0,len(defVals.keys())):
          thisdefval = defVals.keys()[idefval]
          if thisdefval.find('_Im') > -1:
            h('print "L5PC.apic[0].'+thisdefval+' = ", L5PC.apic[0].'+thisdefval+', "Default = ", '+str(defVals[thisdefval][1]))
          else:
            h('print "L5PC.soma[0].'+thisdefval+' = ", L5PC.soma[0].'+thisdefval+', "Default = ", '+str(defVals[thisdefval][0]))
  
        #Restore default values:
        for imutvar in range(0,len(MT[igene][imut])):
          mutvars = allmutvars[iallmutval][imutvar]
          mutvals = allmutvals[iallmutval][imutvar]
          if type(mutvars) is str:
            mutvars = [mutvars]
          for kmutvar in range(0,len(mutvars)):
            mutvar = mutvars[kmutvar]
            newVal = defVals[mutvar]
            if mutvar.find('_Ih') > -1:
              updateThese = [1,1,1]
            elif mutvar.find('_Ca_HVA') > -1 or mutvar.find('_Ca_LVAst') > -1 or mutvar.find('_SKv3.1') > -1 or mutvar.find('_Ca_HVA') > -1 or mutvar.find('_SK_E2') > -1 or mutvar.find('_NaTa_t') > -1 or mutvar.find('_CaDynamics_E2') > -1:
              updateThese = [1,1,0]
            elif mutvar.find('_K_Pst') > -1 or mutvar.find('_K_Tst') > -1 or mutvar.find('_Nap_Et2') > -1:
              updateThese = [1,0,0]
            elif mutvar.find('_Im') > -1:
              updateThese = [0,1,0]
            else:
              print "Error: str=" + str(mutvar)
              updatedThese = [0,0,0]
            for iupdated in range(0,3):
              if updateThese[iupdated]:
                print """forsec L5PC."""+str(updatedVars[iupdated])+""" {
"""+mutvar+""" = """+str(newVal[whichDefVal[iupdated]])+"""
}"""
                h("""forsec L5PC."""+str(updatedVars[iupdated])+""" {
"""+mutvar+""" = """+str(newVal[whichDefVal[iupdated]])+"""
}""")
        gsThisMutVal.append(nextgs[2])
      gsThisMut.append(gsThisMutVal[:])
      picklelist = [theseCoeffsAll,gsThisMutVal,MT]
      file = open('thresholddistalamp_cs'+str(icell)+'_'+str(counter)+'.sav', 'w')
      pickle.dump(picklelist,file)
      file.close()
    gsThisGene.append(gsThisMut[:])
   gsAll.append(gsThisGene[:])
  gsAllAll.append(gsAll[:])
  
  
