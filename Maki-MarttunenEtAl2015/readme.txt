# Python-scripts for reproducing Figures 2-6 in the article 
# "Functional effects of schizophrenia-linked genetic variants on intrinsic single-neuron excitability: A modeling study."
#
# The scripts can be run using Python 2.7.8 with an interface to NEURON 7.3.
#
# If LaTeX installed, change the variable 'useLatex' to True for cleaner figure texts.
#
# All scripts are open for distribution and reuse through Creative Commons Attribution 2.0 Generic license
# (CC BY 2.0). The NEURON interface commands use the cell model and simulation of BAC firing by Etay Hay (Hay et al. 2011)
# The included .asc and .hoc files are copied as such from ModelDB entry 139653 with the author's permission, and the .mod
# files are copied with minor changes (hard-coded parameters are replaced with parameters whose values may be changed).
#
# Tuomo Maki-Marttunen, Jan 2015
#
# To run the python scripts, save them to the same directory with the .mod files, and make sure directories
# 'morphologies' and 'models' exist (in more detail, make sure that files 'morphologies/cell1.asc',
# 'morphologies/cell2.asc', 'models/L5PC_template.hoc', and 'models/L5PCbiophys3.hoc' are accessible). Then,
# follow the procedures below:
#
# 1) Run the .mod file complier:

nrnivmodl

# Before it is possible to draw any of figures, the mutations have to be scaled. To do this,
# do one of the following:
# 
# 2a)
#
# python runcontrols.py
# python scalemutations.py
# 
# Running runcontrols.py calculates the f-I curve and membrane potential limit cycle for the control neuron.
# Running scalemutations.py finds the scaling factors for each variant separately.
# This is a computationally heavy task, so if you can parallelize it to different computing nodes,
# you can alternatively run on each machine the following:
# 
# 2b)
#
# python scalemutations.py i
# 
# where i=0,...,57.
#
# After 2a) or 2b), run
# 
# python collectscalings.py
#
# to collect the scalings into a single file, scalings.sav
# 
# 2c) Alternatively, download the file scalings.sav, where the scalings are calculated in advance.
# 
# 3) You are now ready to draw Figures 2-5 from the article:

python drawfig1.py
python drawfig2.py
python drawfig3.py
python drawfig4.py

# On a tested standard computer, Figure 2 took around one minute to calculate, Figure 3
# was finished in 10 hours, Figure 4 in 15 minutes, and Figure 5 in five hours.
# 
# 4) To draw Figure 6, first calculate the spiking thresholds for uniformly distributed
# synaptic inputs on apical dendrite:

python savesynapselocations.py
python findthresholddistalamps.py
python collectthresholddistalamps.py

# (The second script takes some time to calculate, it can be distributed to different
# computers as was done in 2b with scalemutations.py)
#
# 5) Now, Figure 6 can be drawn:

python drawfig5.py

# This is the figure with heaviest computational load, due to the fine resolution of ISI.
# It took around 6 days to finish the Figure 6 on the tested machine.

