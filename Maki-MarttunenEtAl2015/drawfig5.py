# drawfig5
# A script for plotting the threshold conductances for a second distal stimulus
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
  ylabel_g = '$g_{\mathrm{th}}$ ($\upmu$S)'  
  ylabel_cg = 'threshold $c_g$'
else:
  ylabel_g = 'g_th (uS)'  
  ylabel_cg = 'threshold c_g'

v0 = -80
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620
fs = 8
xs = range(700,1150,50);
tstop = 11000.0
currCoeff = 1.1 # use the threshold current g_th*1.1 for inducing the first spike, and g_th*1.1*c for the second spike, where c saved in thresholddistalamp.sav 
ITERS = 20
PPIdts = range(0,500,2) # use a fine temporal resolution
maxLens = [1300,1185]

barxs = [1,2,3,-1,-2,0]
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

unpicklefile = open('thresholddistalamp.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
gsAllAll = unpickledlist[1]
gs_control = gsAllAll[0][0][1][0][5]

unpicklefile = open('synlocs.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
Nsyns = unpickledlist[0]
synlocsAll = unpickledlist[3]

variants = [[0,1,0],[1,2,14],[3,0,1],[6,1,0],[8,0,0],[11,0,0]]
f, axarr23 = plt.subplots(2, 3)
axnew23 = [[0,0,0],[0,0,0]]
for iy in range(0,2):
  for ix in range(0,3):
    axarr23[iy,ix].set_position([0.1+0.19*ix, 0.1+0.26*(1-iy), 0.19, 0.2])
    axnew23[iy][ix] = f.add_axes([0.22+0.19*ix, 0.22+0.26*(1-iy), 0.06, 0.06],axisbg='w')
axarr = [axarr23[0,0],axarr23[0,1],axarr23[0,2],axarr23[1,0],axarr23[1,1],axarr23[1,2]]
axnew = [axnew23[0][0],axnew23[0][1],axnew23[0][2],axnew23[1][0],axnew23[1][1],axnew23[1][2]]

for ivar in range(0,len(variants)):
  icell = 0
  igene = variants[ivar][0]
  imut = variants[ivar][1]

  synlocs = synlocsAll[icell]
  theseCoeffs = theseCoeffsAllAll[icell][igene][imut]

  morphology_file = "morphologies/cell"+str(icell+1)+".asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"

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
objref vsoma, sl, tvec, syns["""+str(2*Nsyns)+"""]
vsoma = new Vector()
tvec = new Vector()
sl = new List()
double siteVec[2]
L5PC.soma cvode.record(&v(0.5),vsoma,tvec)
""") #"""
  for istim in range(0,Nsyns):
    h("""
siteVec[0] = """+str(synlocs[istim][0])+"""
siteVec[1] = """+str(synlocs[istim][1])+"""
L5PC.apic[siteVec[0]] {
  syns["""+str(istim)+"""] = new AlphaSynapse(siteVec[1])
  syns["""+str(istim)+"""].e = 0
  syns["""+str(istim)+"""].tau = 5
  syns["""+str(istim)+"""].onset = 10000
  syns["""+str(istim+Nsyns)+"""] = new AlphaSynapse(siteVec[1])
  syns["""+str(istim+Nsyns)+"""].e = 0
  syns["""+str(istim+Nsyns)+"""].tau = 5
  syns["""+str(istim+Nsyns)+"""].onset = 10000
}
""")

  nVals = len(MT[igene][imut])*[0]
  thesemutvars = []
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
  gs = gsAllAll[icell][igene][imut][iallmutval]

  iters = [0, 2, 5, 6, 8, -1]
  for iiter in range(0,len(iters)):
    iter = iters[iiter]
    if iter >= 0:
      thiscol = cols[iter]
      thisCoeff = coeffCoeffs[iter][0]*theseCoeffs[iallmutval] + coeffCoeffs[iter][1]*(1.0 - 0.5*theseCoeffs[iallmutval])
      gsThisIter = gs[iiter]        
    else:
      thiscol = col_control
      thisCoeff = 0
      gsThisIter = gs_control
    if iter==5: # Disregard this one (the one corresponding to unscaled variant), but keep it in the list of iters as it was there also in thresholddistalamp.sav
      continue
    gCoeffsThisIter = []
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
      filename = 'fig5_cs'+str(icell)+'_control.sav'
    else:
      filename = 'fig5_cs'+str(icell)+'_ivar'+str(ivar)+'_iter'+str(iter)+'.sav'
    try: # If the simulation has already been made, don't bother rerun it
      unpicklefile = open(filename, 'r')
      unpickledlist = pickle.load(unpicklefile)
      unpicklefile.close()
      PPIdts = unpickledlist[0]
      gCoeffsThisIter = unpickledlist[1]
    except:
      for iPPI in range(0,len(PPIdts)):
        PPIdt = PPIdts[iPPI]

        nextCoeffs = [0,10.0,5.0]
        hasSpiked = 0
        for iterI in range(0,ITERS+2):
          for istim in range(0,Nsyns):
            h("syns["+str(istim)+"].gmax = "+str(gsThisIter*currCoeff))
            h("syns["+str(istim+Nsyns)+"].gmax = "+str(gsThisIter*currCoeff*nextCoeffs[min(iterI,2)]))
            h("syns["+str(istim+Nsyns)+"].onset = "+str(10000+PPIdt))
          h("""
tstop = """+str(tstop)+"""
cai0_ca_ion = """+str(thisCa)+"""
v_init = """+str(v0)+"""
""")
          h.init()
          try:
            h.run()
          except RuntimeError:
            print "Too large I!"
            if iterI > 1:
              nextCoeffs = [nextCoeffs[0],nextCoeffs[2],0.5*(nextCoeffs[0]+nextCoeffs[2])]
            continue
          times=np.array(h.tvec)
          Vsoma=np.array(h.vsoma)
          nSpikes_total = len(mytools.spike_times(times,Vsoma,-35,-37.5))
          print "nextCoeffs="+str(nextCoeffs)+", "+str(nSpikes_total)+" spikes"
          if iterI==0:
            nSpikes_normal = nSpikes_total
          hasSpiked = hasSpiked or (nSpikes_total > nSpikes_normal)
          if iterI > 0 and not hasSpiked:
            nextCoeffs = [nextCoeffs[1],2*nextCoeffs[1],1.5*nextCoeffs[1]]
            continue
          if iterI > 1 and nSpikes_total > nSpikes_normal:
            nextCoeffs = [nextCoeffs[0],nextCoeffs[2],0.5*(nextCoeffs[0]+nextCoeffs[2])]
          if iterI > 1 and nSpikes_total <= nSpikes_normal:
            nextCoeffs = [nextCoeffs[2],nextCoeffs[1],0.5*(nextCoeffs[2]+nextCoeffs[1])]

        gCoeffsThisIter.append(nextCoeffs[2])
      picklelist = [PPIdts,gCoeffsThisIter]
      file = open(filename, 'w')
      pickle.dump(picklelist,file)
      file.close()

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

    axarr[ivar].plot(PPIdts,[x/1.1 for x in gCoeffsThisIter], color=thiscol, linewidth=lw)
    axnew[ivar].bar(barxs[iiter], 1000*gsThisIter, 0.8, color=thiscol)

  axarr[ivar].set_title(geneNames[igene],fontsize=fs+2)
  if ivar % 3 == 0:
    axarr[ivar].set_xticks([0,50,100,150,200])
  else:
    axarr[ivar].set_xticks([50,100,150,200])
  axarr[ivar].set_xlim([0,200])
  axarr[ivar].set_yticks([0,2,4])
  axarr[ivar].set_ylim([0,5.7])
  if ivar < 3:
    axarr[ivar].set_xticklabels(['']*(4+(ivar%3==0)))
  else:
    axarr[ivar].set_xlabel('ISI (ms)',fontsize=fs+2)
  if ivar % 3 > 0:
    axarr[ivar].set_yticklabels(['','',''])
  else:
    axarr[ivar].set_ylabel(ylabel_cg,fontsize=fs+2)
  for tick in axarr[ivar].yaxis.get_major_ticks()+axarr[ivar].xaxis.get_major_ticks():
    tick.label.set_fontsize(fs)
  axarr[ivar].xaxis.set_ticks_position('bottom')
  axarr[ivar].yaxis.set_ticks_position('left')

  axnew[ivar].set_xlim([-3,3.8])
  axnew[ivar].set_ylim([0,0.037])
  axnew[ivar].set_xticks([])
  axnew[ivar].set_yticks([0,0.01,0.02,0.03])
  for tick in axnew[ivar].yaxis.get_major_ticks()+axnew[ivar].xaxis.get_major_ticks():
    tick.label.set_fontsize(fs-2)
  axnew[ivar].xaxis.set_ticks_position('bottom')
  axnew[ivar].yaxis.set_ticks_position('left')
  axnew[ivar].set_xlabel(ylabel_g,fontsize=fs)

  if useLatex:
    params = {'text.latex.preamble': [r"\usepackage{upgreek}"],
              'text.usetex': True}
    plt.rcParams.update(params)
  f.savefig("figure5.eps")

