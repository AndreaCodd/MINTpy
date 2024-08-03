#!/usr/bin/python3
import numpy as np
import importlib, sys, os
import argparse

'''
Makes a mesh with the core region with horizontal extent greater than the observations.  
Name for the observation points is in the config file.
Must have observation points FIRST.  
'''

parser = argparse.ArgumentParser(description='creates a GMSH geo file', epilog="a.codd@uq.edu.au")
parser.add_argument(dest='config', metavar='configfile', type=str, help='py thon setting configuration')
args = parser.parse_args()
config = importlib.import_module(args.config)
print("configuration "+args.config+" imported.")

# get points and ground level
easts = np.load(config.rsEastsName)
snum = len(easts)
grLevel = config.grLevel

# core spans sensors + a bit more
Sxmin = min(easts)
Sxmax = max(easts)
print(Sxmin,Sxmax)
sspan = Sxmax-Sxmin

# Core dimensions
sidex = config.sqx
sidey = config.sqy

Cxmin = np.floor((Sxmin - sspan*config.smx)/sidex)*sidex 
Cxmax = np.ceil((Sxmax + sspan*config.smx)/sidex)*sidex
Cymax = config.grLevel
Cymin = np.floor((Cymax - config.cy)/sidey)*sidey

Cxspan = Cxmax-Cxmin
Cyspan = Cymax-Cymin

print("corex", Cxmin, Cxmax, Cxspan) 
print("corey", Cymin, Cymax, Cyspan)


# core element sizes
cg = config.cg  # ground level
cb = config.cb  # bottom of core

# buffer 
Bx   = sspan*config.sbx
By   = config.cy*config.sby
Bair = config.airLayer

# buffer element sizes
bb = config.bb
bg = config.bg
ba = config.ba

# sensor element size
sg = config.sg
cs = 5.#config.cs

# assume core x is span of sensors
out = "Cxmin =%e;\n"%Cxmin
out+= "Cxmax =%e;\n"%Cxmax
out+= "Cymin =%e;\n"%Cymin
out+= "Cymax =%e;\n"%Cymax  
out+= "cg =%e;\n"%cg
out+= "cb =%e;\n"%cb
#out+= "cs =%e;\n"%cs

out+= "Bx = %e;\n"%Bx
out+= "By = %e;\n"%By
out+= "Bair = %e;\n"%Bair
out+= "bb = %e;\n"%bb
out+= "bg = %e;\n"%bg
out+= "ba = %e;\n"%ba
out+= "snum = %e;\n"%snum
out+= "sspan = %e;\n"%sspan
out+= "sg = %e;\n"%sg



out+="""
sx= sspan/(snum-1);
s1= sspan/2.;


// 2D
// CORE
// core Points

Point(1) = {Cxmin, Cymin, 0, cb};
Point(2) = {Cxmax, Cymin, 0, cb};
Point(3) = {Cxmax, Cymax, 0, cg};
Point(4) = {Cxmin, Cymax, 0, cg};

Point(5) = {Cxmin - Bx , Cymin-By, 0.0, bb};
Point(6) = {Cxmax + Bx , Cymin-By, 0.0, bb};
Point(7) = {Cxmax + Bx , Cymax, 0.0, bg};
Point(8) = {Cxmax + Bx , Bair, 0.0, ba};
Point(9) = {Cxmin - Bx , Bair, 0.0, ba};
Point(10) = {Cxmin - Bx , Cymax, 0.0, bg};

// mesh nodes 
k=newp;
""" 

for s in range(snum):
    out+="Point(k+%d)={%e, Cymax, 0.0, sg};\n"%(s,easts[s])

out+="""
// core lines, ground, bottom, sides
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,k};
For i In {1:snum-1}
    Line(3+i) = {k+i-1,k+i};
EndFor
l=newl;
Line(l) = {k+snum-1,4};
Line(l+1) = {4,1};
Line Loop(1) = {1:l+1};
Plane Surface(1) = {1};
Physical Surface("core") = {1};
m=newl;
Line(m) = {5,6};
Line(m+1) = {6,7};
Line(m+2) = {7,8};
Line(m+3) = {8,9};
Line(m+4) = {9,10};
Line(m+5) = {10,5};
Line(m+6) = {10,4};
Line(m+7) = {3,7};
Line Loop(2) = {m,m+1,-m-7,-2,-1,-l-1,-m-6,m+5};
Plane Surface(2) = {2};
Physical Surface("buffer") = {2};

Line Loop(3) = {-m-7,3:l,-m-6,-m-4,-m-3,-m-2};
Plane Surface(3) = {3};
Physical Surface("air") = {3};
//MeshSize{k:k+snum} = cs;


"""

open(config.geofile,'w').write(out)
print("geofile created written to "+config.geofile)

print("made geo file")

