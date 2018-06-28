# drawfig2
# A script for plotting the F-I curves
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
  flabel = '$f$ (Hz)'
  ilabel = '$I$ (nA)'
else:
  flabel = 'f (Hz)'
  ilabel = 'I (nA)'

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

styles = ['b-','b-','b-','b-','b-','b-','b-','b-','b-','b-']
cols = ['#444444','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#009999','#772277','#00cc00']
col_control = '#2222ff'
coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]
lw = 0.6
fs = 10
ispDef = 1 # Consider every local maximum above -35mV a spike
Is = [0.35+0.05*x for x in range(0,22)]

v0 = -80
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620

variants = [[0,1,0],[1,2,14],[3,0,1],[6,1,0],[8,0,0],[11,0,0]]
f, axarr23 = plt.subplots(2, 3)
for iy in range(0,2):
  for ix in range(0,3):
    axarr23[iy,ix].set_position([0.1+0.19*ix, 0.1+0.26*(1-iy), 0.19, 0.2])
axarr = [axarr23[0,0],axarr23[0,1],axarr23[0,2],axarr23[1,0],axarr23[1,1],axarr23[1,2]]

for ivar in range(0,len(variants)):
  icell = 0
  theseCoeffsAll = theseCoeffsAllAll[icell]
  igene = variants[ivar][0]
  imut = variants[ivar][1]

  morphology_file = "morphologies/cell"+str(icell+1)+".asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"

  h("""
load_file("stdlib.hoc")                                                                                                                                                                                
load_file("stdrun.hoc")                                                                                                                                                                                
objref cvode                                                                                                                                                                                           
cvode = new CVode()                                                                                                                                                                                    
cvode.active(1)                                                                                                                                                                                        
cvode.atol(0.001)                                                                                                                                                                                      
load_file("import3d.hoc")                                                                                                                                                                              
objref L5PC                                                                                                                                                                                            
load_file(\""""+biophys_file+"""\")
load_file(\""""+template_file+"""\")
L5PC = new L5PCtemplate(\""""+morphology_file+"""\")
access L5PC.soma
objref st1, tvec, sl
L5PC.soma st1 = new IClamp(0.5)
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
objref vsoma, vdend, vdend2, cadend, cadend2, casoma
vsoma = new Vector()
casoma = new Vector()
vdend = new Vector()
cadend = new Vector()
vdend2 = new Vector()
cadend2 = new Vector()
cvode.record(&v(0.5),vsoma,tvec)
cvode.record(&cai(0.5),casoma,tvec)
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
L5PC.apic[siteVec[0]] cvode.record(&v(siteVec[1]),vdend2,tvec)
L5PC.apic[siteVec[0]] cvode.record(&cai(siteVec[1]),cadend2,tvec)
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

    if iter==-1:
      filename = 'fig2_cs'+str(icell)+'_control.sav'
    else:
      filename = 'fig2_cs'+str(icell)+'_ivar'+str(ivar)+'_iter'+str(iter)+'.sav'
    try: # If the simulation has already been made, don't bother rerun it
      unpicklefile = open(filename, 'r')
      unpickledlist = pickle.load(unpicklefile)
      unpicklefile.close()
      Is = unpickledlist[0]
      nSpikes = unpickledlist[1]
    except:
      nSpikes = []
      for iI in range(0,len(Is)):
        tstop = 8000.0
        squareAmp = Is[iI]
        squareDur = 7800.0
        h("""
tstop = """+str(tstop)+"""
v_init = """+str(v0)+"""  
cai0_ca_ion = """+str(ca0)+"""  
st1.amp = """+str(squareAmp)+"""
st1.del = 200                   
st1.dur = """+str(squareDur)+"""
""")
        h.init()
        h.run()

        times=np.array(h.tvec)
        Vsoma=np.array(h.vsoma)
        Vdend=np.array(h.vdend)
        spikes = mytools.spike_times(times,Vsoma,-35,-45)
        nSpikes.append(sum([1 for x in spikes if x >= 500]))

      picklelist = [Is,nSpikes]
      file = open(filename, 'w')
      pickle.dump(picklelist,file)
      file.close()

    axarr[ivar].plot(Is,[x/7.5 for x in nSpikes], color=thiscol, linewidth=lw)
    #Print the parameters and their default values:                                                                                                                                                
    for idefval in range(0,len(defVals.keys())):
      thisdefval = defVals.keys()[idefval]
      if thisdefval.find('_Im') > -1:
        h('print "L5PC.apic[0].'+thisdefval+' = ", L5PC.apic[0].'+thisdefval+', "Default = ", '+str(defVals[thisdefval][1]))
        #) #+" (def="+str(defVals[thisdefval])+")"                                                                                                                                                 
      else:
        h('print "L5PC.soma[0].'+thisdefval+' = ", L5PC.soma[0].'+thisdefval+', "Default = ", '+str(defVals[thisdefval][0]))
        #h('print L5PC.soma[0]."+thisdefval) #+" (def="+str(defVals[thisdefval])+")"                                                                                                               

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
  axarr[ivar].set_xticks([0.4,0.8,1.2])
  axarr[ivar].set_xlim([0.2,1.4])
  axarr[ivar].set_yticks([0,10,20])
  axarr[ivar].set_ylim([0,20])
  if ivar < 3:
    axarr[ivar].set_xticklabels(['','',''])
  else:
    axarr[ivar].set_xlabel(ilabel,fontsize=fs+2)
  if ivar % 3 > 0:
    axarr[ivar].set_yticklabels(['','',''])
  else:
    axarr[ivar].set_ylabel(flabel,fontsize=fs+2)
  for tick in axarr[ivar].yaxis.get_major_ticks()+axarr[ivar].xaxis.get_major_ticks():
    tick.label.set_fontsize(fs)
  axarr[ivar].xaxis.set_ticks_position('bottom')
  axarr[ivar].yaxis.set_ticks_position('left')

if useLatex:
  params = {'text.latex.preamble': [r"\usepackage{upgreek}"],
            'text.usetex': True}
  plt.rcParams.update(params)
f.savefig("figure2.eps")
