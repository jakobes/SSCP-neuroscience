# collectthresholddistalamps.py
# A script for collecting all variants' threshold conductances for spike
# generation as a response to distal stimulus
#
# Tuomo Maki-Marttunen, Jan 2015
# (CC BY)
from neuron import h
import matplotlib
matplotlib.use('Agg')
from pylab import *
import mytools
import pickle

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

gsAllAll = []

for icell in range(0,1):
  gsAll = []
  theseCoeffsAll = theseCoeffsAllAll[icell]

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
      try:
        unpicklefile = open('thresholddistalamp_cs'+str(icell)+'_'+str(counter)+'.sav', 'r')
        unpickledlist = pickle.load(unpicklefile)
        unpicklefile.close()
        gsThisMut.append(unpickledlist[1])
      except:
        gsThisMut.append([])        
    gsThisGene.append(gsThisMut[:])
   gsAll.append(gsThisGene[:])
  gsAllAll.append(gsAll[:])
  
picklelist = [theseCoeffsAllAll,gsAllAll,MT]
file = open('thresholddistalamp.sav', 'w')
pickle.dump(picklelist,file)
file.close()
  
