# drawfig4
# A script for plotting the neuron responsiveness to combination of apical and
# somatic stimuli during "up" and "down" states
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

useLatex = False
if useLatex:
  Vmlabel = 'max $V_m$ (mV)'
else:
  Vmlabel = 'max V_m (mV)'

v0 = -80
ca0 = 0.0001
proximalpoint = 200
distalpoint = 760
tstop = 5000.0
longSquareAmps = [[0,0.42],[0,0.42]] # see Hay 2011. The values should be readjusted for cell #2
shortSquareAmps = [[1.8,0.5],[1.8,0.5]]  # see Hay 2011. The values should be readjusted for cell #2
epspAmps = [[0.5,0.5],[0.5,0.5]]         # see Hay 2011. The values should be readjusted for cell #2
epspdts = [0.25*x for x in range(-80,81)]

styles = ['b-','b-','b-','b-','b-','b-','b-','b-','b-','b-']
downstyles = ['b--','b--','b--','b--','b--','b--','b--','b--','b--','b--']
cols = ['#444444','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#009999','#772277','#00cc00']
col_control = '#2222ff'
coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]
lw = 0.6
fs = 10

import mutation_stuff
MT = mutation_stuff.getMT()
geneNames = mutation_stuff.getgenenames()
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

variants = [[0,1,0],[1,2,14],[3,0,1],[6,1,0],[8,0,0],[11,0,0]]
f, axarr23 = plt.subplots(2, 3)
for iy in range(0,2):
  for ix in range(0,3):
    axarr23[iy,ix].set_position([0.1+0.19*ix, 0.1+0.26*(1-iy), 0.19, 0.2])
axarr = [axarr23[0,0],axarr23[0,1],axarr23[0,2],axarr23[1,0],axarr23[1,1],axarr23[1,2]]

for ivar in range(0,len(variants)):
  icell = 0
  igene = variants[ivar][0]
  imut = variants[ivar][1]
  theseCoeffsAll = theseCoeffsAllAll[icell]
  longSquareAmp = longSquareAmps[icell]
  shortSquareAmp = shortSquareAmps[icell]
  epspAmp = epspAmps[icell]
  morphology_file = "morphologies/cell"+str(icell+1)+".asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"

  theseCoeffsAll = theseCoeffsAllAll[icell]
  times_control = [[[],[],[]],[[],[],[]]]
  Vsoma_control = [[[],[],[]],[[],[],[]]]
  Vdend_control = [[[],[],[]],[[],[],[]]]

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
access L5PC.soma
objref st1, st2
L5PC.soma st2 = new IClamp(0.5)
st2.amp = 0
st2.del = 1000
st2.dur = 5
objref vsoma, vdend, vdend2, cadend, cadend2, casoma
vsoma = new Vector()
casoma = new Vector()
vdend = new Vector()
cadend = new Vector()
vdend2 = new Vector()
cadend2 = new Vector()
objref sl,syn1,tvec
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
L5PC.apic[siteVec[0]] cvode.record(&v(siteVec[1]),vdend,tvec)
L5PC.apic[siteVec[0]] cvode.record(&cai(siteVec[1]),cadend,tvec)
L5PC.apic[siteVec[0]] {
  syn1 = new epsp(siteVec[1])
  syn1.tau0 = 0.5
  syn1.tau1 = 5
  syn1.onset = 1000
  syn1.imax = 0
}
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
L5PC.apic[siteVec[0]] st1 = new IClamp(siteVec[1])
st1.amp = 0
st1.del = 900
st1.dur = 200
L5PC.soma cvode.record(&v(0.5),vsoma,tvec)
L5PC.soma cvode.record(&cai(0.5),casoma,tvec)
""") #"""

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
    
  iallmutval = variants[ivar][2]

  VsomaupThisMutVal = []
  VsomadownThisMutVal = []
  VdendupThisMutVal = []
  VdenddownThisMutVal = []

  iters = [0,2,6,8,-1]
  for iiter in range(0,len(iters)):
    iter = iters[iiter]
    if iter >= 0:
      thiscol = cols[iter]
      thisstyle = styles[iter]
      thisdownstyle = downstyles[iter]
    else:
      thiscol = col_control
      thisstyle = 'b-'
      thisdownstyle = 'b--'
      
    VsomaupThisIter = []
    VsomadownThisIter = []
    VdendupThisIter = []
    VdenddownThisIter = []

    if iter >= 0:
      thisCoeff = coeffCoeffs[iter][0]*theseCoeffs[iallmutval] + coeffCoeffs[iter][1]*(1.0 - 0.5*theseCoeffs[iallmutval])
    else:
      thisCoeff = 0
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

    if iter==-1:
      filename = 'fig4_cs'+str(icell)+'_control.sav'
    else:
      filename = 'fig4_cs'+str(icell)+'_ivar'+str(ivar)+'_iter'+str(iter)+'.sav'
    try: # If the simulation has already been made, don't bother rerun it
      unpicklefile = open(filename, 'r')
      unpickledlist = pickle.load(unpicklefile)
      unpicklefile.close()
      VsomadownThisIter = unpickledlist[0]
      VdenddownThisIter = unpickledlist[1]
      VsomaupThisIter = unpickledlist[2]
      VdendupThisIter = unpickledlist[3]
    except:
      for idt in range(0,len(epspdts)):
        for iup in range(0,2):
          h("st1.amp = "+str(longSquareAmp[iup]))
          h("st2.amp = "+str(shortSquareAmp[iup]))
          h("syn1.imax = "+str(epspAmp[iup]))
          h("""
tstop = """+str(tstop)+"""
cai0_ca_ion = """+str(thisCa)+"""
v_init = """+str(v0)+"""
syn1.onset = """+str(1000+epspdts[idt])+"""
""")
          h.init()
          h.run()

          times=np.array(h.tvec)
          Casoma=np.array(h.casoma)
          Vsoma=np.array(h.vsoma)
          Cadend=np.array(h.cadend)
          Vdend=np.array(h.vdend)

          if iup == 0:
            VsomadownThisIter.append(max(Vsoma))
            VdenddownThisIter.append(max(Vdend))
          else:
            VsomaupThisIter.append(max(Vsoma))
            VdendupThisIter.append(max(Vdend))
      picklelist = [VsomadownThisIter,VdenddownThisIter,VsomaupThisIter,VdendupThisIter]
      file = open(filename, 'w')
      pickle.dump(picklelist,file)
      file.close()

    axarr[ivar].plot(epspdts,VdendupThisIter, thisstyle, color=thiscol, linewidth=lw)
    axarr[ivar].plot(epspdts,VdenddownThisIter, thisdownstyle, color=thiscol, linewidth=lw)[0].set_dashes([2,1])

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

  axarr[ivar].set_title(geneNames[igene],fontsize=fs+2)
  if ivar % 3 == 0:
    axarr[ivar].set_xticks([-20,-10,0,10,20])
  else:
    axarr[ivar].set_xticks([-10,0,10,20])
  axarr[ivar].set_xlim([-20,20])
  axarr[ivar].set_yticks([-40,-20,0,20])
  axarr[ivar].set_ylim([-50,23])
  if ivar < 3:
    axarr[ivar].set_xticklabels(['']*(4+(ivar%3==0)))
  else:
    axarr[ivar].set_xlabel('ISI (ms)',fontsize=fs+2)
  if ivar % 3 > 0:
    axarr[ivar].set_yticklabels(['','',''])
  else:
    axarr[ivar].set_ylabel(Vmlabel,fontsize=fs+2)
  for tick in axarr[ivar].yaxis.get_major_ticks()+axarr[ivar].xaxis.get_major_ticks():
    tick.label.set_fontsize(fs)
  axarr[ivar].xaxis.set_ticks_position('bottom')
  axarr[ivar].yaxis.set_ticks_position('left')

  if useLatex:
    params = {'text.latex.preamble': [r"\usepackage{upgreek}"],
              'text.usetex': True}
    plt.rcParams.update(params)
  f.savefig("figure4.eps")
