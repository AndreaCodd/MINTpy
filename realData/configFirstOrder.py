import numpy as np
"""
testing Real Data from EDI
"""
##############################################
#  DATA
edi_path= "~/MINTpy/realData/EDIfiles"
name = "TEST/"
# computed data directory
Dname = "data/"

# data file names
rsEastsName = Dname+"rsEasts.npy"
EastsName   = Dname+"Easts.npy"
ObsName     = Dname+"ObsName.txt" 
freqsName   = Dname+"frequencies.npy"
RMWSsqName  = Dname+"RMWSsq.npy"
SensorsName = Dname+"sensors.npy"
dataZxyName = Dname+"dataZxy.npy"
dataZyxName = Dname+"dataZyx.npy"
wxyName     = Dname+"wxy.npy"
wyxName     = Dname+"wyx.npy"
ExyName     = Dname+"Exy.npy"
EyxName     = Dname+"Eyx.npy"
phasexyName = Dname+"phasexy.npy"
appRxyName  = Dname+"appRxy.npy"
phaseyxName = Dname+"phaseyx.npy"
appRyxName  = Dname+"appRyx.npy"

# Signs for the TE and TM data, mu and appropriate data scale
signTE = -1.
signTM = -1.
mu = 4*np.pi*1e-7
dataScale = 1000.*mu

##############################################
#  MESH file names
geofile = "mesh/shallow.geo"
flyfile = "mesh/shallow.fly"

#  MESH specifications depend on the data

# CORE
# spans (1 + 2*smx)*sensor span
smx = 0.1
# make core multiples of sqx horizontally and sqy vertically
sqx = 1000.
sqy = 1000.
# depth and ground level
cy = 8000.
grLevel = 0.0

# Buffer
# scales and air height
sbx = 3 #8.
sby = 79 #149.
airLayer = 600000.

# element size 
# scales sensors span 
# core
cg = 40.
cb = 200.

# buffer
bb = 80000. 
bg = 20.
ba = 80000.

#sensors, scale for tolerance for location  
sg = cg 

##############################################
#  First order INVERSION
w0 = 0.
w1 = 1.e6
wa = 1000.
m_tol = 1.e-3

sigmaBG = 0.1

fixBC = False
withEps = False
eps=0.005
Reduced = False
firstS = 0
lastS = -1
stepS = 1
firstF = 0 
lastF = -1
stepF = 20

name = "shallow"

if fixBC:
    name = name+"BC"
else:
    name = name+"1D"
if not withEps:
    name = name+"Errors"

silofile = "silos/"+name+"_%f_%f_"%(np.log10(w1),int(np.log10(wa)))
FinalZxyName = "output/Zxy_"+name+"_%f_%f.npy"%(np.log10(w1),int(np.log10(wa)))
FinalZyxName = "output/Zyx_"+name+"_%f_%f.npy"%(np.log10(w1),int(np.log10(wa)))
InitZxyName = "output/Zxy0_"+name+"_%f_%f.npy"%(np.log10(w1),int(np.log10(wa)))
InitZyxName = "output/Zyx0_"+name+"_%f_%f.npy"%(np.log10(w1),int(np.log10(wa)))










