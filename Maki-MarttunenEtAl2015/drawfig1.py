# drawfig1
# A script for plotting the result figure with [Ca] responses to short somatic step-current input.
# For CACNA1I variants, also the membrane potential and [Ca] time courses are illustrated.
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
  tlabel = '$t$ (ms)'
  Vmlabel = '$V_m$ (mV)'
  xlabel_Ca = '[Ca$^{2+}$] ($\upmu$M)'  
  ylabel_Ca = '$d$[Ca$^{2+}$]/$dt$ ($\upmu$M/ms)'
else:
  tlabel = 't (ms)'
  Vmlabel = 'V_m (mV)'
  xlabel_Ca = '[Ca2+] (uM)'  
  ylabel_Ca = 'd[Ca2+]/dt (uM/ms)'

v0 = -80
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620
fs = 8
tstop = 13000.0
squareAmps = [1.626,1.6204] # induces stably one spike (halfway from threshold of 1 to threshold of 2 spikes) (not very stable for the second cell!)

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

CasomasAllAll = []
VsomasAllAll = []
timesAllAll = []

styles = ['g-','g-','g-','g-','g-','g-','g-','g-','g-']
cols = ['#444444','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#009999','#772277','#00cc00']
col_control = '#2222ff'
coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]
lw = 0.6

variants = [[0,1,0],[1,2,14],[3,0,1],[6,1,0],[8,0,0],[11,0,0]]
close("all")
f, axarr25 = plt.subplots(2, 5)
for iy in range(0,2):
  for ix in range(0,3):
    axarr25[iy,ix+2].set_position([0.42+0.19*ix, 0.1+0.4*(1-iy), 0.19, 0.34])
  axarr25[iy,0].set_position([0.1, 0.1+0.4*(1-iy), 0.19, 0.34])
  axarr25[iy,1].set_position([0.14, 0.29+0.4*(1-iy), 0.135, 0.12])
axarr = [axarr25[0,2],axarr25[0,3],axarr25[0,4],axarr25[1,2],axarr25[1,3],axarr25[1,4]]
axtimecourse = [axarr25[0,0],axarr25[1,0]]
axzoom = [axarr25[0,1],axarr25[1,1]]
axtimecourse[0].add_patch(Rectangle((9995,-100),30,200,edgecolor='None', facecolor='#CCCCCC'))
axtimecourse[1].add_patch(Rectangle((9995,-1),30,2,edgecolor='None', facecolor='#CCCCCC'))
axzoom[0].plot([10000,10005],[-85,-85], 'b-', color='#FF0000', linewidth=lw)
axzoom[1].plot([10000,10005],[0.0975,0.0975], 'b-', color='#FF0000', linewidth=lw)

for ivar in range(0,len(variants)):
  icell = 0
  squareAmp = squareAmps[icell]

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
cvode.atol(0.0002)
load_file("import3d.hoc")
objref L5PC
load_file(\""""+biophys_file+"""\")
load_file(\""""+template_file+"""\")
L5PC = new L5PCtemplate(\""""+morphology_file+"""\")
access L5PC.soma
objref st1,sl,tvec
st1 = new IClamp(0.5)
L5PC.soma st1
objref vsoma, vdend, recSite, vdend2, isoma, cadend, cadend2, casoma
vsoma = new Vector()
casoma = new Vector()
vdend = new Vector()
cadend = new Vector()
vdend2 = new Vector()
cadend2 = new Vector()
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

  counter = -1

  igene = variants[ivar][0]
  imut = variants[ivar][1]

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

  iters = [0, 2, 6, 8, -1]
  for iiter in range(0,len(iters)):
    iter = iters[iiter]
    if iter >= 0:
      thiscol = cols[iter]
    else:
      thiscol = col_control

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
            mutText = mutText + "+" + str(mutvals*thisCoeff) +" mV"
          elif kmutvar==0:
            mutText = mutText  + str(mutvals*thisCoeff) +" mV"
        else:
          newVal = [x*(mutvals**thisCoeff) for x in defVals[mutvar]]
          if kmutvar==0:
            mutText = mutText + "*" + str(mutvals**thisCoeff)
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

    if iter==-1:
      filename = 'fig1_cs'+str(icell)+'_control.sav'
    else:
      filename = 'fig1_cs'+str(icell)+'_ivar'+str(ivar)+'_iter'+str(iter)+'.sav'
    try: # If the simulation has already been made, don't bother rerun it
      unpicklefile = open(filename, 'r')
      unpickledlist = pickle.load(unpicklefile)
      unpicklefile.close()
      times = unpickledlist[0]
      Vsoma = unpickledlist[1]
      Casoma = unpickledlist[2]
    except: # Otherwise, run the simulation and save the results
      thisCa = h.L5PC.soma[0].minCai_CaDynamics_E2
      h("""
tstop = """+str(tstop)+"""
cai0_ca_ion = """+str(thisCa)+"""
v_init = """+str(v0)+"""
st1.amp = """+str(squareAmp)+"""
st1.del = 10000
st1.dur = 5
""")
      h.init()
      h.run()

      times=np.array(h.tvec)
      Vsoma=np.array(h.vsoma)
      Casoma=np.array(h.casoma)
      picklelist = [times,Vsoma,Casoma]
      file = open(filename, 'w')
      pickle.dump(picklelist,file)
      file.close()

    if ivar == 2: # Draw the time courses for CACNA1I variant
      axtimecourse[0].plot(times,Vsoma, 'b-', color=thiscol, linewidth=lw)
      axtimecourse[1].plot(times,[1000.0*x for x in Casoma], 'b-', color=thiscol, linewidth=lw)
      axzoom[0].plot(times,Vsoma, 'b-', color=thiscol, linewidth=lw)
      axzoom[1].plot(times,[1000.0*x for x in Casoma], 'b-', color=thiscol, linewidth=lw)

    times2 = [x for x in times if x > 9900]
    Vsoma2 = [Vsoma[i] for i,x in enumerate(times) if x > 9900]
    Casoma2 = [Casoma[i] for i,x in enumerate(times) if x > 9900]
    itoremove = min([i for i,x in enumerate(times2) if x > 10005])
    times2.pop(itoremove)
    Vsoma2.pop(itoremove)
    Casoma2.pop(itoremove)

    istart = 2
    iend = len(Vsoma2)-1
    nsteps = iend-istart-1
    tdiff = [y-x for x,y in zip(times2[istart:iend-1],times2[istart+1:iend])]
    cadiff = [y-x for x,y in zip(Casoma2[istart:iend-1],Casoma2[istart+1:iend])]
    caderiv1 = [y/x for x,y in zip(tdiff[0:nsteps-1],cadiff[0:nsteps-1])]
    caderiv2 = [y/x for x,y in zip(tdiff[1:nsteps],cadiff[1:nsteps])]
    caderiv = [(x+y)/2.0 for x,y in zip(caderiv1,caderiv2)]

    axarr[ivar].plot([1000.0*x for x in Casoma2[istart+1:iend-1]],[1000.0*x for x in caderiv], color=thiscol, linewidth=lw)

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

  axarr[ivar].set_title(geneNames[igene])
  axarr[ivar].set_xticks([0.1,0.12,0.14])
  axarr[ivar].set_xlim([0.095,0.155])
  axarr[ivar].set_yticks([0,0.01,0.02])
  axarr[ivar].set_ylim([-0.001,0.022])
  if ivar < 3:
    axarr[ivar].set_xticklabels(['','',''])
  else:
    axarr[ivar].set_xlabel(xlabel_Ca)
  if ivar % 3 > 0:
    axarr[ivar].set_yticklabels(['','',''])
  else:
    axarr[ivar].set_ylabel(ylabel_Ca)

for ix in range(0,5):
  for iy in range(0,2):
    axarr25[iy,ix].xaxis.set_ticks_position('bottom')
    axarr25[iy,ix].yaxis.set_ticks_position('left')

axtimecourse[0].set_xticks([10000,11000,12000])
axtimecourse[0].set_xticklabels(['','',''])
axtimecourse[0].set_xlim([9900,12100])
axtimecourse[0].set_yticks([-80,-40,0,40])
axtimecourse[0].set_ylim([-90,50])
axtimecourse[0].set_ylabel(Vmlabel)
axtimecourse[0].set_xlabel(tlabel)
axtimecourse[1].set_xticks([10000,11000,12000])
axtimecourse[1].set_xticklabels(['0','1000','2000'])
axtimecourse[1].set_xlim([9900,12100])
axtimecourse[1].set_yticks([0.1,0.12,0.14])
axtimecourse[1].set_ylim([0.095,0.155])
axtimecourse[1].set_ylabel(xlabel_Ca)
axtimecourse[1].set_xlabel(tlabel)
axzoom[0].set_xticks([10000,10020])
axzoom[0].set_xticklabels(['0','20'])
axzoom[0].set_xlim([9995,10025])
axzoom[0].set_yticks([])
axzoom[0].set_ylim([-90,50])
axzoom[1].set_xticks([])
axzoom[1].set_xticklabels(['0','20'])
axzoom[1].set_xlim([9995,10025])
axzoom[1].set_xticks([10000,10020])
axzoom[1].set_yticks([])
axzoom[1].set_ylim([0.095,0.155])

f.text(0.01, 0.81, 'A', fontsize=33)
f.text(0.01, 0.41, 'B', fontsize=33)
f.text(0.31, 0.81, 'C', fontsize=33)

if useLatex:
  params = {'text.latex.preamble': [r"\usepackage{upgreek}"],
            'text.usetex': True}
  plt.rcParams.update(params)
f.savefig("figure1.eps")



