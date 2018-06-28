# savesynapselocations.py
# A script for randomly picking locations along the apical dendrite (from 300um on)
# A maximum of 50 synapses per segment are allowed, while the total number of synapses
# is 3000.
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

random.seed(1)

proximalpoint = 400
distalpoint = 620
BACdt = 5.0
fs = 8
maxLens = [1300,1185]

maxSynsPerSeg = 50
Nsyns = 3000

synlocsAll = []

for icell in range(0,2):
  morphology_file = "morphologies/cell"+str(icell+1)+".asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"
  synlocs = []

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
objref sl,syn1,con1,isyn, tvec, syns["""+str(Nsyns)+"""]
tvec = new Vector()
sl = new List()
double siteVec[2]
""")
  synsInSegs = [0]*len(h.L5PC.apic)
  for istim in range(0,Nsyns):
    myiseg = -1
    while myiseg == -1:
      x = 300.0+(maxLens[icell]-300)*random.random()
      h("""sl = L5PC.locateSites("apic","""+str(x)+""")
Nsegs_x = sl.count()
""")
      iseg = random.randint(0,h.Nsegs_x-1)
      h("iseg = sl.o["+str(iseg)+"].x[0]")
      if synsInSegs[int(h.iseg)] < maxSynsPerSeg:
        myiseg = int(h.iseg)
        break
      print "istim = "+str(istim)+", x = "+str(x)+", continue searching for iseg..."
    synsInSegs[myiseg] = synsInSegs[myiseg] + 1
    h("""
siteVec[0] = sl.o[j].x[0]
siteVec[1] = sl.o[j].x[1]
L5PC.apic[siteVec[0]] {
  syns["""+str(istim)+"""] = new AlphaSynapse(siteVec[1])
  syns["""+str(istim)+"""].e = 0
  syns["""+str(istim)+"""].tau = 5
  syns["""+str(istim)+"""].onset = 10000 + """+str(BACdt)+""" 
}
""")
    synlocs.append([h.siteVec[0],h.siteVec[1]])
  
  synlocsAll.append(synlocs[:])
picklelist = [Nsyns,maxSynsPerSeg,maxLens,synlocsAll]
file = open('synlocs.sav', 'w')
pickle.dump(picklelist,file)
file.close()
