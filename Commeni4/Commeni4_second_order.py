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
from esys.downunder import MinimizerLBFGS, MeteredCostFunction
sys.path.append(os.getcwd())

class MTInv2(MeteredCostFunction):
    provides_inverse_Hessian_approximation=True
    def __init__(self, domain, sigmaBG, sensors, periods, m0, Dxys=[], Dyxs=[],
                       w0=1., a=1., b=1., mu=4*np.pi*1e-7, doTE=True, doTM = True,
                       eps=0.05):

        super(MTInv2, self).__init__()

        self.domain = domain
        XF=ReducedFunction(self.domain).getX()
        self.GrF = interpolate(whereNegative(XF[1]), Function(domain))
        self.GrF.expand()
        
        self.count=0
        self.mu = mu 
        self.w0 = w0
        self.a=a
        self.b=b
        self.sigma = None
        self.rho = None
        if isinstance(sigmaBG, float) or isinstance(sigmaBG, int):
            self.sigmaBG = sigmaBG*self.GrF
            self.sigmaBGCF = sigmaBG
        else:
            self.sigmaBG = interpolate(sigmaBG, Function(domain))
            self.sigmaBGCF = interpolate(sigmaBG, ContinuousFunction(domain))
        self.setSigmaRho(m0)

        self.periods = periods
        self.numP = len(periods) 
        self.doTE = doTE
        self.doTM = doTM 
        if self.doTE:
            self.TEpde = self.createPDE()
            self.adjTEpde =self.createPDE()
            self.TEpde.setValue(A=kronecker(self.domain.getDim()))
            self.adjTEpde.setValue(A=kronecker(self.domain.getDim()))

        if self.doTM:
            self.TMpde = self.createPDE()
            self.adjTMpde =self.createPDE()        

        z = self.domain.getX()[self.domain.getDim()-1]
        t = sup(z)
        b = inf(z)
        air=wherePositive(z)

        if self.doTE:
            self.TEpde.setValue(q=whereZero(z-t)+whereZero(z-b), r=(z-b)/(t-b))
        if self.doTM:
            self.TMpde.setValue(q=air+whereZero(z-b), r=air) 
                     
        if self.doTE:
            self.adjTEpde.setValue(q=self.TEpde.getCoefficient('q'), r=Data())
        if self.doTM:
            self.adjTMpde.setValue(q=self.TMpde.getCoefficient('q'), r=Data())      

        rr = 100
        adfn=[]
        self.AF = Scalar(0., Function(self.domain))
        xyz=ReducedFunction(self.domain).getX()
        count=1
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
            self.scxy.append(1./a1/self.numP))
            self.scyx.append(1./b1/self.numP))
      
        # Hessian pde
        self.Hpde=LinearSinglePDE(self.domain, isComplex=False)
        self.Hpde.setSymmetryOn()
        A= Data(0.,(2,2),Function(self.domain))
        A[0,0] = Scalar(self.w0*self.a*self.b, Function(self.domain))
        A[1,1] = Scalar(self.w0*self.a, Function(self.domain))                
        self.Hpde.setValue(A=A, D = self.w0, q=whereZero(z-t))
        optionsH=self.Hpde.getSolverOptions()
        optionsH.setTolerance(1.e-2)  
        optionsH.setPackage(SolverOptions.TRILINOS)
        optionsH.setSolverMethod(SolverOptions.PCG)
        optionsH.setPreconditioner(SolverOptions.AMG)


    def createPDE(self):
        pde=LinearSinglePDE(self.domain, isComplex=True)
        options = pde.getSolverOptions()
        options.setSolverMethod(SolverOptions.DIRECT)
        options.setTolerance(1.e-10)
        pde.setSymmetryOn()
        return pde

    def _getDualProduct(self, M, R):
        My = M[0]
        Mz = M[1]
        m  = M[2]
        gMy = grad(My)
        gMz = grad(Mz)      
        gm = grad(m)
        R0 = R[0]
        R1 = R[1]
        R2 = R[2]  
        FS = R0.getFunctionSpace()
        iMy = interpolate(My, FS)
        iMz = interpolate(Mz, FS)
        im  = interpolate(m , FS)
        dpY = integrate( R0[0]*iMy + R0[1]*iMz + R0[2]*im)
        FS = R1.getFunctionSpace()
        igMy = interpolate(gMy, FS)
        igMz = interpolate(gMz, FS)
        igm  = interpolate(gm , FS)
        dpX = integrate(R1[0,0]*igMy[0]+R1[0,1]*igMy[1]+R1[1,0]*igMz[0]+R1[1,1]*igMz[1]+R1[2,0]*igm[0]+R1[2,1]*igm[1])
        dpy = integrate(inner(R[2],interpolate(m,R[2].getFunctionSpace()))) 
        return dpY+dpX+dpy 

    def _getNorm(self, x):
        return Lsup(x)

    def setSigmaRho(self,M):
        mF = clip(interpolate(M[2], Function(self.domain)),minval=-20., maxval = 20. )
        self.sigma = self.GrF*self.sigmaBG*exp(mF)
        air = (1-self.GrF)
        self.rhoG = (1./(self.sigmaBG+1e-10*air))*self.GrF*exp(-mF) 
        self.rho = interpolate(self.rhoG + air*1.e10, Function(self.domain))
        return self

    def _getArguments(self, M):
        self.setSigmaRho(M)
        Exs=[]
        Hxs=[]
        Zxys = []
        Zyxs = []
        z = FunctionOnBoundary(self.domain).getX()[self.domain.getDim()-1]
        b = inf(z)
        t = sup(z)
        for p in self.periods:
            om=2.*np.pi/p            
            if self.doTE:
                self.TEpde.setValue(D=1j*om*self.mu*self.sigma)
                Ex = self.TEpde.getSolution() 
                gEx = grad(Ex, Function(self.domain))
                Hy = -1./(1j*om*self.mu)*gEx[self.domain.getDim()-1]
                Zxy = safeDiv(interpolate(Ex, Hy.getFunctionSpace()), Hy)
                Exs.append(Ex)
                Zxys.append(Zxy)
            if self.doTM:
                self.TMpde.setValue(A=self.rho*kronecker(self.domain.getDim()))
                self.TMpde.setValue(D=1j*om*self.mu)
                Hx = self.TMpde.getSolution()
                gHx = grad(Hx, Function(self.domain))
                Ey = self.rho*gHx[self.domain.getDim()-1] 
                Zyx = safeDiv(Ey,interpolate(Hx, Ey.getFunctionSpace()))
                Hxs.append(Hx)
                Zyxs.append(Zyx)
        return Exs, Hxs, Zxys, Zyxs

    def _getValue(self, M, *args):
        My = M[0]
        Mz = M[1]
        m  = M[2]
        gMy = grad(My)
        gMz = grad(Mz)      
        gm = grad(m)
        FS = gm.getFunctionSpace()
        iMy = interpolate(My, FS)
        iMz = interpolate(Mz, FS)
        im  = interpolate(m , FS)
        divM = self.b*gMy[0] + gMz[1]
        curlM = self.b*gMz[0] - gMy[1]

        if len(args)==0:
            args=self.getArguments(M)
        Exs, Hxs, Zxys, Zyxs = args
        normTE = 0.
        normTM = 0.
        for ip in range(self.numP):
            if self.doTE:
                DEFxy = Zxys[ip] - self.wDxys[ip]
                normTE += integrate(self.AF*length(DEFxy)**2)*self.scxy[ip]            
            if self.doTM:
                DEFyx = Zyxs[ip] - self.wDyxs[ip] 
                normTM += integrate(self.AF*length(DEFyx)**2)*self.scyx[ip]
        A0 = self.w0*integrate( (iMy - self.a*self.b*gm[0] )**2 + (iMz - self.a*gm[1] )**2 )
        A1 = self.w0*integrate((im - self.a*divM)**2) 
        A2 = self.w0*self.a**2*integrate((curlM)**2)

        print ("reg, mf_xy, mf_yx, totmf =", A0/2.+A1/2.+A2/2., normTE/2.,normTM/2.,normTE/2.+ normTM/2.)
        return (A0+A1+A2+normTE+normTM)/2. 

    def _getGradient(self, M, *args):
        if len(args)==0: 
            args=self.getArguments(M)
        Exs, Hxs, Zxys, Zyxs = args
        My = M[0]
        Mz = M[1]
        m  = M[2]
        gMy = grad(My)
        gMz = grad(Mz)      
        gm  = grad(m)
        FS = gm.getFunctionSpace()
        im  = interpolate(m , FS)
        iMy = interpolate(My, FS)
        iMz = interpolate(Mz, FS)
        divM = self.b*gMy[0] + gMz[1]
        curlM = self.b*gMz[0] - gMy[1]

        Y = Data(0.,(3,),FS)
        Y[0] = self.w0*( iMy - self.a*self.b*gm[0] )
        Y[1] = self.w0*( iMz - self.a*gm[1] )
        Y[2] = self.w0*(  im - self.a*divM  )

        X = Data(0.,(3,2),FS)
        X[0,0] = - self.a*self.w0 *(im - self.a*divM)
        X[0,1] = - self.a**2*self.w0*curlM
        X[1,0] =   self.a**2*self.w0*curlM
        X[1,1] = - self.a*self.w0 *(im - self.a*divM)
        X[2,0] = - self.a*self.w0*(iMy - self.a*self.b*gm[0])
        X[2,1] = - self.a*self.w0*(iMz - self.a*gm[1])        
        
        y = Data(0., FunctionOnBoundary(self.domain))
        if self.doTM:
            self.adjTMpde.setValue(A=self.rho*kronecker(self.domain.getDim()))
        z = FunctionOnBoundary(self.domain).getX()[self.domain.getDim()-1]
        b = inf(z)
        for ip in range(self.numP):
            om = 2.*np.pi/self.periods[ip] 
            if self.doTE: 
                self.adjTEpde.setValue(D=1j*om*self.mu*self.sigma)
                DEFxy = Zxys[ip] - self.wDxys[ip]
                gEx = grad(Exs[ip], Function(self.domain))                
                hatDxy = self.scxy[ip]*safeDiv(DEFxy.conjugate(), gEx[self.domain.getDim()-1])
                YTE = -1j*om*self.mu*self.AF*hatDxy
                XTE = -self.AF*hatDxy*Zxys[ip]*kronecker(2)[1]
                self.adjTEpde.setValue(Y=YTE, X=XTE) 
                ExStar = self.adjTEpde.getSolution()
                Y[2] += (-1j*om*self.mu*self.sigma*interpolate(Exs[ip],FS)*interpolate(ExStar,FS)).real()
            if self.doTM:
                self.adjTMpde.setValue(D=1j*om*self.mu)               
                DEFyx = Zyxs[ip] - self.wDyxs[ip] 
                hatDyx =self.scyx[ip]*safeDiv(DEFyx.conjugate(), interpolate(Hxs[ip], DEFyx.getFunctionSpace()))
                YTM = -self.AF*hatDyx*Zyxs[ip]
                XTM = self.rho*self.AF*hatDyx*kronecker(2)[1]
                self.adjTMpde.setValue(Y=YTM, X=XTM) 
                HxStar = self.adjTMpde.getSolution()
                g1=grad(Hxs[ip], FS)
                g2=grad(HxStar, FS)
                #gg =   inner(g1,g2)
                #bob = self.AF*hatDyx*g1[1]  
                Y[2] += (self.rhoG*(inner(g1,g2) - self.AF*hatDyx*g1[1])).real()
        return ArithmeticTuple(Y,X,y) 

    def _getInverseHessianApproximation(self, M, R, *args):
        H = Data(0.,(3,),Solution(self.domain))
        for ind1 in range(3):
            self.Hpde.setValue(Y=R[0][ind1], X=R[1][ind1])
            H[ind1] = self.Hpde.getSolution()
        return H



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
w0 = 1.e-10
a = 1000.
b = 10.
eps = 0.05
mtol=1e-4

#############################################
#############################################
# Forward

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

model = MT2Dmodel(domain, sigmaBG, sensors, periods, mu=mu, fixBottom=True, airLayer = 0.0)
model.setSigmaRhoForward(trueSigma, trueRho, Ground = groundF)

sigma1, rho1 = model.getSigmaRho()
saveSilo("silos/InputModel", sigma=sigma1, rho = rho1)

Exs, Hxs, Zxy, Zyx = model.getExHxZs()

if True:
    print("inversion 2D MT fixed base")
    M0 = Data(0., (3,), Solution(domain))
    costFn = MTInv2(domain, sigmaBG, sensors, periods, M0, 
                       Dxys = Zxy, Dyxs = Zyx, w0 = w0,  a=a, b=b, mu = mu,
                       doTE = True, doTM = True, eps = eps)
    solver = MinimizerLBFGS(J=costFn, imax=1000, m_tol = m_tol)
    solver.run(M0)
    M=solver.getResult()
    print("m",M[2])
    saveSilo("silos/Commeni4_2ndOrder_w0_%f_a_%f_b_%f"%(log(w0),a,b), m=M[2], sigma=costFn.sigma,
                rho=costFn.rho)    
print("finished")




