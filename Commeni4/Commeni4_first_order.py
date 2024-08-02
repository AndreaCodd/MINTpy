from esys.escript import *
from esys.escript import length
import numpy as np
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions
from esys.escript.pdetools import Locator, ArithmeticTuple, MaskFromTag
from esys.finley import ReadMesh, ReadGmsh
from esys.weipa import saveSilo
from MTmodel import MT2Dmodel 
import argparse
import importlib, sys, os
import os.path
import logging
#from time import time
from esys.downunder import MinimizerLBFGS, MeteredCostFunction
sys.path.append(os.getcwd())

class MTInv(MeteredCostFunction):
    provides_inverse_Hessian_approximation=True
    def __init__(self, domain, sigmaBG, sensors, periods, m0, Dxys=[], Dyxs=[],
                       w0=0., w1=1.,wa=1., mu=4*np.pi*1e-7, eps=0.05):

        super(MTInv, self).__init__()

        self.domain = domain
        XF=ReducedFunction(self.domain).getX()
        self.GrF = interpolate(whereNegative(XF[1]), Function(domain))
        self.GrF.expand()
        self.mu = mu 
        self.w0 = w0
        self.w1 = w1 
        self.wa = wa
        self.sigma = None
        self.rho = None
        self.fixB = fixB
        self.loc = Locator(ReducedFunction(self.domain),[s for s in sensors])
        if isinstance(sigmaBG, float) or isinstance(sigmaBG, int):
            self.sigmaBG = sigmaBG*self.GrF
        else:
            self.sigmaBG = interpolate(sigmaBG, Function(domain))
        self.setSigmaRho(m0)
              
        self.periods=periods
        self.numP = len(periods)

        self.TEpde = self.createPDE()
        self.adjTEpde =self.createPDE()
        self.TEpde.setValue(A=kronecker(self.domain.getDim()))
        self.adjTEpde.setValue(A=kronecker(self.domain.getDim()))

        self.TMpde = self.createPDE()
        self.adjTMpde =self.createPDE()        

        z = self.domain.getX()[self.domain.getDim()-1]
        t = sup(z)
        b = inf(z)
        air= wherePositive(z)

        self.TEpde.setValue(q=whereZero(z-t)+whereZero(z-b), r=(z-b)/(t-b))
        self.TMpde.setValue(q=air+whereZero(z-b), r=air) 
             
        self.adjTEpde.setValue(q=self.TEpde.getCoefficient('q'), r=Data())
        self.adjTMpde.setValue(q=self.TMpde.getCoefficient('q'), r=Data())      


        rr = 100
        adfn=[]
        self.AF = Scalar(0., Function(self.domain))
        xyz=ReducedFunction(self.domain).getX()
        for ind1 in range(len(sensors)):
            locPr = Locator(ReducedFunction(self.domain),sensors[ind1])
            RR = (length(xyz-locPr.getValue(xyz)))**2
            af=interpolate(whereNegative(RR-rr), Function(self.domain))
            adfn.append(af)
            self.AF += af
        intAF = integrate(self.AF)
        self.wDxys=[] 
        self.wDyxs=[] 
        self.scxy=[]
        self.scyx=[]
        for ind1 in range(self.numP):
            Dxy = Dxys[ind1]
            Dyx = Dyxs[ind1]
            wDxy = ComplexScalar(0.+0.j, Function(self.domain))
            wDyx = ComplexScalar(0.+0.j, Function(self.domain))
            for ind2 in range(len(sensors)):
                wDxy += Dxy[ind2]*adfn[ind2] 
                wDyx += Dyx[ind2]*adfn[ind2]
            self.wDxys.append(wDxy)
            self.wDyxs.append(wDyx)
            a1 = eps**2*integrate(length(wDxy)**2)
            b1 = eps**2*integrate(length(wDyx)**2)        
            self.scxy.append(1./a1/self.numP)#/intAF)
            self.scyx.append(1./b1/self.numP)#/intAF)
        # Hessian pde
        self.Hpde=LinearSinglePDE(self.domain, isComplex=False)
        self.Hpde.setSymmetryOn()
        A= Data(0.,(2,2),Function(self.domain))
        A[0,0] = Scalar(self.w1*self.wa, Function(self.domain))
        A[1,1] = Scalar(self.w1, Function(self.domain))
        self.Hpde.setValue(A=A, D = self.w0, q=whereZero(z-t))
        optionsH=self.Hpde.getSolverOptions()
        optionsH.setTolerance(1.e-2)  
        optionsH.setPackage(SolverOptions.TRILINOS)
        optionsH.setSolverMethod(SolverOptions.PCG)
        optionsH.setPreconditioner(SolverOptions.AMG)

        self.getValue(m0)

    def createPDE(self):
        pde=LinearSinglePDE(self.domain, isComplex=True)
        options = pde.getSolverOptions()
        options.setSolverMethod(SolverOptions.DIRECT)
        options.setTolerance(1.e-10)
        pde.setSymmetryOn()
        return pde

    def _getDualProduct(self, m, r):
        rmY = integrate(r[0]*interpolate(m,r[0].getFunctionSpace()))             
        rmX = integrate(inner(r[1],grad(m)))  
        rmy = integrate(inner(r[2],interpolate(m,r[2].getFunctionSpace())))
        #print("dual", rmY,rmX,rmy)     
        return rmY+rmX+rmy 

    def _getNorm(self, x):
        return Lsup(x)

    def setSigmaRho(self,m):
        mF = clip(interpolate(m, Function(self.domain)),minval=-20., maxval = 20. )
        mF = mF*self.GrF
        self.sigma = self.sigmaBG*exp(mF)
        air = (1.-self.GrF)
        self.rhoG = safeDiv(1., self.sigmaBG)*self.GrF*exp(-mF) 
        self.rho = self.rhoG + air*1.e10
        return self

    def _getArguments(self, m):
        self.setSigmaRho(m)
        Exs=[]
        Hxs=[]
        Zxys = []
        Zyxs = []

        z = FunctionOnBoundary(self.domain).getX()[self.domain.getDim()-1]
        b = inf(z)
        t = sup(z)

        for p in self.periods:
            #print("p",p)
            om=2.*np.pi/p
            
            self.TEpde.setValue(D=1j*om*self.mu*self.sigma)
            Ex = self.TEpde.getSolution() 
            gEx = grad(Ex, Function(self.domain))
            Hy = -1./(1j*om*self.mu)*gEx[self.domain.getDim()-1]
            Zxy = interpolate(Ex, Hy.getFunctionSpace())/Hy #safeDiv(interpolate(Ex, Hy.getFunctionSpace()), Hy)
            Exs.append(Ex)
            Zxys.append(Zxy)
            self.TMpde.setValue(A=self.rho*kronecker(self.domain.getDim()))
            self.TMpde.setValue(D=1j*om*self.mu)
            Hx = self.TMpde.getSolution()
            gHx = grad(Hx, Function(self.domain))
            Ey = self.rho*gHx[self.domain.getDim()-1] 
            Zyx = safeDiv(Ey,interpolate(Hx, Ey.getFunctionSpace()))
            Hxs.append(Hx)
            Zyxs.append(Zyx)
        return Exs, Hxs, Zxys, Zyxs

    def _getValue(self, m, *args):       
        if len(args)==0:
            args=self.getArguments(m)
        Exs, Hxs, Zxys, Zyxs = args
        normTE = 0.
        normTM = 0.
        
        for ip in range(self.numP):
            DEFxy = Zxys[ip] - self.wDxys[ip]
            normTE += integrate(self.AF*length(DEFxy)**2)*self.scxy[ip]            
            DEFyx = Zyxs[ip] - self.wDyxs[ip] 
            normTM += integrate(self.AF*length(DEFyx)**2)*self.scyx[ip]
        gm = grad(m)     
        A1 = self.wa*self.w1*integrate(gm[0]**2) +  self.w1*integrate(gm[0]**2)    
        A0=self.w0*integrate(interpolate(m, Function(self.domain))**2) 
        print ("reg, mf_xy, mf_yx, totmf =",A0/2.+A1/2.,normTE/2.,normTM/2.,normTE/2.+ normTM/2.)
        
        return (A0+A1+normTE+normTM)/2. 

    def _getGradient(self, m, *args):
        if len(args)==0: 
            args=self.getArguments(m)
        Exs, Hxs, Zxys, Zyxs = args
        gm = grad(m) 
        FS = gm.getFunctionSpace()
        X= Data(0.,(2,),FS)
        X[0] = gm[0]*self.w1*self.wa
        X[1] = gm[1]*self.w1 
        #X = self.w1*grad(m)
        Y = self.w0*interpolate(m,FS)
        y = Scalar(0., FunctionOnBoundary(self.domain))
        self.adjTMpde.setValue(A=self.rho*kronecker(self.domain.getDim()))

        #z = self.domain.getX()[self.domain.getDim()-1]
        z = FunctionOnBoundary(self.domain).getX()[self.domain.getDim()-1]
        b = inf(z)

        self.adjTMpde.setValue(A=self.rho*kronecker(self.domain.getDim()))

        for ip in range(self.numP):
            om = 2*np.pi/self.periods[ip] 

            self.adjTEpde.setValue(D=1j*om*self.mu*self.sigma)
            DEFxy = Zxys[ip] - self.wDxys[ip]
            gEx = grad(Exs[ip], Function(self.domain))                
            hatDxy = self.scxy[ip]*safeDiv(DEFxy.conjugate(), gEx[self.domain.getDim()-1])
            YTE = -1j*om*self.mu*self.AF*hatDxy
            XTE = -self.AF*hatDxy*Zxys[ip]*kronecker(2)[1]
            self.adjTEpde.setValue(Y=YTE, X=XTE) 
            ExStar = self.adjTEpde.getSolution()
            Y += (-1j*om*self.mu*self.sigma*interpolate(Exs[ip],FS)*interpolate(ExStar,FS)).real()

            self.adjTMpde.setValue(D=1j*om*self.mu)                 
            DEFyx = Zyxs[ip] - self.wDyxs[ip] 
            hatDyx =self.scyx[ip]*safeDiv(DEFyx.conjugate(), interpolate(Hxs[ip], DEFyx.getFunctionSpace()))
            YTM = -self.AF*hatDyx*Zyxs[ip]
            XTM = self.rho*self.AF*hatDyx*kronecker(2)[1]
            self.adjTMpde.setValue(Y=YTM, X=XTM) 
            HxStar = self.adjTMpde.getSolution()
            g1=grad(Hxs[ip], FS)
            g2=grad(HxStar, FS)
            gg = inner(g1,g2)  
            bob = self.AF*hatDyx*g1[1]  
            Y += (self.rhoG*(gg - bob)).real()
        return ArithmeticTuple(Y,X,y) 

    def _getInverseHessianApproximation(self, m, r, *args):
        self.Hpde.setValue(X=r[1], Y=r[0])
        return self.Hpde.getSolution()



## -------------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------------
print()
print("** 2D isotropic MT BFGS inversion Testing**")
print()

#############################################
#############################################
# domain, sensors, periods,
domain = ReadMesh("comm4Mesh.fly", numDim=2)
sensors = [(-20000., 0., 0.), (-10000., 0., 0.), (-7000., 0., 0.), 
           (-6000., 0., 0.), (-5000., 0., 0.), (0., 0., 0.), 
           (2000., 0., 0.), (5000., 0., 0.), (8000., 0., 0.), 
           (16000., 0., 0.)]
periods = np.logspace(-2, 2, num=40, endpoint=True, base=10.0, dtype=float)
mu = 4*np.pi*1e-7

######
# inversion variables
sigmaBG = 1./5.  
w0 = 0.
w1 = 1.e-5
wa = 10.
fixB = True   #config.fixBottom
m0 = Scalar(0., ContinuousFunction(domain)) 
eps = 0.05

#############################################
#############################################
# run Forward problem to get data

# ground 
XF=ReducedFunction(domain).getX() 
groundF = whereNegative(XF[1])
groundF = interpolate(groundF, Function(domain)) 

# set sigma and rho
rhoBase   = 5.
rhoMiddle = 1000.
rhoBlob1  = 10.
rhoBlob2  = 2.5
rhoTop    = 25.
rhoAir    = 1.e10

trueSigma = Scalar(0., Function(domain))
trueSigma.setTaggedValue("Base", 1./rhoBase) 
trueSigma.setTaggedValue("Middle", 1./rhoMiddle)
trueSigma.setTaggedValue("Blob1", 1./rhoBlob1)
trueSigma.setTaggedValue("Blob2", 1./rhoBlob2)
trueSigma.setTaggedValue("Top", 1./rhoTop)
trueSigma.expand()

trueRho = Scalar(0., Function(domain))
trueRho.setTaggedValue("Base", rhoBase) 
trueRho.setTaggedValue("Middle", rhoMiddle)
trueRho.setTaggedValue("Blob1", rhoBlob1)
trueRho.setTaggedValue("Blob2", rhoBlob2)
trueRho.setTaggedValue("Top", rhoTop)
trueRho.setTaggedValue("air", rhoAir)
trueRho.expand()

# use esys-escript models
model = MT2Dmodel(domain, sigmaBG, sensors, periods, mu=mu, fixBottom=True, airLayer = 0.0)
model.setSigmaRhoForward(trueSigma, trueRho, Ground = groundF)
sigma1, rho1 = model.getSigmaRho()
Exs, Hxs, Zxy, Zyx = model.getExHxZs()
print("finished forward\n\n")

#######################################################
#######################################################
# Inversion

if True:             
    print("w1", w1, "sigmaBG", sigmaBG, "wa", wa)
    print("inversion 2D MT fixed base")
    costFn = MTInv(domain, sigmaBG, sensors, periods, m0, 
                       Dxys = Zxy, Dyxs = Zyx, w0 = w0, w1 = w1, wa=wa, mu = mu,
                       eps = eps)
    solver = MinimizerLBFGS(J=costFn, imax=1000)
    solver.run(m0)
    m=solver.getResult()
    print("m",m)
    sigma=costFn.sigma
    rho = costFn.rho
    saveSilo("silos/Commeni4_sigmaBG_%f_w1_%f_wa_%f"%(sigmaBG,np.log10(w1),wa), m=m, sigma=sigma, rho=rho)

print("finished inversion")


