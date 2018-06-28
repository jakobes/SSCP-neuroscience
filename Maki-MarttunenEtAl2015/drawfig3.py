# drawfig3
# A script for plotting the steady-state firing properties
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
fs = 12

v0 = -80
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620
variants = [[0,1,0],[1,2,14],[3,0,1],[6,1,0],[8,0,0],[11,0,0]]

f, axarr25 = plt.subplots(2, 5)
for iy in range(0,2):
  for ix in range(0,3):
    axarr25[iy,ix+2].set_position([0.42+0.19*ix, 0.1+0.4*(1-iy), 0.19, 0.34])
  axarr25[iy,0].set_position([0.1, 0.1+0.4*(1-iy), 0.19, 0.34])
  axarr25[iy,1].set_position([0.17, 0.16+0.4*(1-iy), 0.105, 0.14])
axarr = [axarr25[0,2],axarr25[0,3],axarr25[0,4],axarr25[1,2],axarr25[1,3],axarr25[1,4]]
axtimecourse = [axarr25[0,0],axarr25[1,0]]
axzoom = [axarr25[0,1],axarr25[1,1]]

for ivar in range(0,len(variants)):
      icell = 0
      theseCoeffsAll = theseCoeffsAllAll[icell]

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
objref st1,sl,tvec
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
L5PC.apic[siteVec[0]] cvode.record(&v(siteVec[1]),vdend2,tvec)
L5PC.apic[siteVec[0]] cvode.record(&cai(siteVec[1]),cadend2,tvec)
""") #"""

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
                mutText = mutText + "+" + "{0:.3f}".format(mutvals*thisCoeff) +" mV"
              elif kmutvar==0:
                mutText = mutText  + "{0:.3f}".format(mutvals*thisCoeff) +" mV"
            else:
              newVal = [x*(mutvals**thisCoeff) for x in defVals[mutvar]]
              if kmutvar==0:
                mutText = mutText + "*" + "{0:.3f}".format(mutvals**thisCoeff)
            if kmutvar < len(mutvars)-1:
              mutText = mutText + ", "
            if mutvar.find('_Ih') > -1:
              updateThese = [1,1,1]
            elif mutvar.find('_Ca_HVA') > -1 or mutvar.find('_Ca_LVAst') > -1 or mutvar.find('_SKv3.1') > -1 or mutvar.find('_Ca_HVA') > -1 or mutvar.find('_SK_E2') > -1 or mutvar.find('_NaTa_t') > \
-1 or mutvar.find('_CaDynamics_E2') > -1:
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
          print "newVal="+str(newVal[0])+","+str(newVal[1])
        print geneNames[igene]+", iter="+str(iter)+", mutText: "+mutText
        print mutText                                                                                                                                                                                  

        if iter==-1:
          filename = 'fig3_cs'+str(icell)+'_control.sav'
        else:
          filename = 'fig3_cs'+str(icell)+'_ivar'+str(ivar)+'_iter'+str(iter)+'.sav'
        try: # If the simulation has already been made, don't bother rerun it          
          unpicklefile = open(filename, 'r')
          unpickledlist = pickle.load(unpicklefile)
          unpicklefile.close()
          times = unpickledlist[0]
          Vsoma = unpickledlist[1]
          Casoma = unpickledlist[2]
          spikes = unpickledlist[3]
        except:
          tstop = 4000.0                                                                                                                                                                                 
          squareAmp = 1.2                                                                                                                                                                                
          squareDur = 3800.0                                                                                                                                                                             
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
          Casoma=np.array(h.casoma)                        
          spikes = mytools.spike_times(times,Vsoma,-35,-45)

          picklelist = [times,Vsoma,Casoma,spikes]
          file = open(filename, 'w')
          pickle.dump(picklelist,file)
          file.close()

        spTimesThisCoeff = spikes[:]                     
        nSpikes1 = len(spikes)
        if nSpikes1 > 5:                                                                                                            
          spts = spikes[nSpikes1-3:nSpikes1]                                                                                  
          istart = next((i for i,x in enumerate(times) if x > spts[0]))                                                             
          iend = next((i for i,x in enumerate(times) if x > spts[1]))+4                                                             
          nsteps = iend-istart-1                                                                                                    
          tdiff = [y-x for x,y in zip(times[istart:iend-1],times[istart+1:iend])]                                                   
          cadiff = [y-x for x,y in zip(Casoma[istart:iend-1],Casoma[istart+1:iend])]                                                
          caderiv1 = [y/x for x,y in zip(tdiff[0:nsteps-1],cadiff[0:nsteps-1])]                                                     
          caderiv2 = [y/x for x,y in zip(tdiff[1:nsteps],cadiff[1:nsteps])]                                                         
          caderiv = [(x+y)/2.0 for x,y in zip(caderiv1,caderiv2)]                                                                   
          axarr[ivar].plot([1000.0*x for x in Casoma[istart+1:iend-1]], [1000.0*x for x in caderiv], color=thiscol, linewidth=lw)

        if ivar==0: # Draw the membrane potential and [Ca] time course for CACNA1C variants
          axtimecourse[0].plot(times,Vsoma, color=thiscol,linewidth=lw)
          axtimecourse[1].plot(times,[1000.0*x for x in Casoma], color=thiscol,linewidth=lw)
          axzoom[0].plot(times,Vsoma, color=thiscol,linewidth=lw)
          axzoom[1].plot(times,[1000.0*x for x in Casoma], color=thiscol,linewidth=lw)

      axarr[ivar].set_title(geneNames[igene],fontsize=fs+2.4)
      axarr[ivar].set_xticks([0.23,0.24,0.25])
      axarr[ivar].set_xlim([0.228,0.258])
      axarr[ivar].set_yticks([0,0.002,0.004,0.006])
      axarr[ivar].set_ylim([-0.0006,0.0073])
      if ivar < 3:
        axarr[ivar].set_xticklabels(['','',''])
      else:
        axarr[ivar].set_xlabel(xlabel_Ca)
      if ivar % 3 > 0:
        axarr[ivar].set_yticklabels(['','',''])
      else:
        axarr[ivar].set_ylabel(ylabel_Ca)
      for tick in axarr[ivar].yaxis.get_major_ticks()+axarr[ivar].xaxis.get_major_ticks():
        tick.label.set_fontsize(fs)

for ix in range(0,5):
  for iy in range(0,2):
    axarr25[iy,ix].xaxis.set_ticks_position('bottom')
    axarr25[iy,ix].yaxis.set_ticks_position('left')

axtimecourse[0].set_xticks([200,400,600])
axtimecourse[0].set_xticklabels(['0','200','400'])
axtimecourse[0].set_xlim([190,700])
axtimecourse[0].set_yticks([-100,-50,0])
axtimecourse[0].set_ylim([-300,40])
axtimecourse[0].set_ylabel(Vmlabel)

axtimecourse[1].set_xticks([200,400,600])
axtimecourse[1].set_xticklabels(['0','200','400'])
axtimecourse[1].set_xlim([190,700])
axtimecourse[1].set_yticks([0.1,0.15,0.2,0.25])
axtimecourse[1].set_ylim([0.09,0.27])
axtimecourse[1].set_ylabel(xlabel_Ca)
axtimecourse[1].set_xlabel(tlabel)

axzoom[0].set_xticks([3600,3800])
axzoom[0].set_xticklabels(['3400','3600'])
axzoom[0].set_xlim([3580,3862])
axzoom[0].set_yticks([-50,0])
axzoom[0].set_ylim([-70,25])

axzoom[1].set_xticks([3600,3800])
axzoom[1].set_xticklabels(['3400','3600'])
axzoom[1].set_xlim([3580,3862])
axzoom[1].set_yticks([0.24,0.25])
axzoom[1].set_ylim([0.234,0.255])

f.text(0.01, 0.81, 'A', fontsize=33)
f.text(0.01, 0.41, 'B', fontsize=33)
f.text(0.31, 0.81, 'C', fontsize=33)

if useLatex:
  params = {'text.latex.preamble': [r"\usepackage{upgreek}"],
            'text.usetex': True}
  plt.rcParams.update(params)
f.savefig("figure3.eps")
