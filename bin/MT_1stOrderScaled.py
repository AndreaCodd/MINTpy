from esys.escript import *
from esys.escript import length
import numpy as np
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions
from esys.escript.pdetools import Locator, ArithmeticTuple, MaskFromTag
from esys.finley import ReadMesh, ReadGmsh
from esys.weipa import saveSilo
import argparse
import importlib, sys, os
import os.path
import logging
#from time import time
from esys.downunder import MinimizerLBFGS, MeteredCostFunction
sys.path.append(os.getcwd())


'''
First order BFGS inversion.

options
1. BC fixed on lower boundary or 1D
2. given errors or percentage for errors

Input
1. frequencies
2. sensors
3. data
4. errors if existing
5. mesh

'''


class MTInv(MeteredCostFunction):
    provides_inverse_Hessian_approximation=True
    def __init__(self, domain, sigmaBG, sensors, freqs, m0, Dxys=[], Dyxs=[], 
                       Wxys=[], Wyxs=[], Exys=[],Eyxs=[], fixBC = True,
                       w0=0., w1=1., wa=1., mu=4*np.pi*1e-7):
                       
        super(MTInv, self).__init__()

        self.domain = domain
        self.GrCF = MaskFromTag(self.domain,"core","buffer")     # CF
        self.GrCF.expand()
        self.GrF = Scalar(0.,Function(domain))                   # F            
        self.GrF.setTaggedValue("core",1.)
        self.GrF.setTaggedValue("buffer",1.)
        self.GrF.expand()
        self.loc = Locator(ReducedFunction(self.domain),[s for s in sensors])
                
        self.mu = mu 
        self.w0 = w0
        self.w1 = w1 
        self.wa = wa
        self.sigmaBG = sigmaBG
        self.setSigmaRho(m0)
        
        self.om = 2*np.pi*freqs
        self.numP = len(freqs)
         

        # set up PDEs
        self.TEpde = self.createPDE()
        self.adjTEpde =self.createPDE()        

        self.TEpde.setValue(A=kronecker(self.domain.getDim())) 
        self.adjTEpde.setValue(A=kronecker(self.domain.getDim())) 

        self.TMpde = self.createPDE()
        self.adjTMpde =self.createPDE()          
                           
        # BC
        self.fixBC = fixBC
        z = self.domain.getX()[self.domain.getDim()-1]
        t = sup(z)
        b = inf(z)
        air=1-self.GrCF       
        if self.fixBC:
            self.TEpde.setValue(q=whereZero(z-t)+whereZero(z-b), r=(z-b)/(t-b))         
            self.TMpde.setValue(q=air+whereZero(z-b), r=air)              
        else:
            self.TEpde.setValue(q=whereZero(z-t), r=(z-b)/(t-b))        
            self.TMpde.setValue(q=wherePositive(z), r=wherePositive(z))        
                               
        self.adjTEpde.setValue(q=self.TEpde.getCoefficient('q'), r=Data())
        self.adjTMpde.setValue(q=self.TMpde.getCoefficient('q'), r=Data()) 
        
        # measurements
        rr = 0.001
        adfn=[]                  
        xyz=ReducedFunction(self.domain).getX()
        for ind1 in range(len(sensors)):
            locPr = Locator(ReducedFunction(self.domain),sensors[ind1])
            RR = (length(xyz-locPr.getValue(xyz)))**2
            af = interpolate(whereNegative(RR-rr), Function(self.domain))
            adfn.append(af)
                   
        self.wDxys=[]    # Array of complex Functions: Data for each frequency
        self.wDyxs=[]    
        self.scxy=[]     # Array of real Functions: scale for errors  
        self.scyx=[]
        self.AFxys = []  # Array of real Functions: 1 at measurements for each frequency  
        self.AFyxs = []   
        self.Wxys=Wxys   # Array of Vectors: 1 at measurement for each frequency
        self.Wyxs=Wyxs      
        self.totDataArea = 0. 
        self.totData = 0.        

        for ind1 in range(self.numP):
            # Vectors of data, existance of data, scale for data
            Dxy = Dxys[ind1]  
            Dyx = Dyxs[ind1]
            Exy = Exys[ind1]
            Eyx = Eyxs[ind1]
            wxy = Wxys[ind1]  
            wyx = Wyxs[ind1]
           
            # Zero Functions 
            wDxy = ComplexScalar(0.+0.j, Function(self.domain))
            wDyx = ComplexScalar(0.+0.j, Function(self.domain))
            scxy = Scalar(0., Function(self.domain))
            scyx = Scalar(0., Function(self.domain))
            AFxy = Scalar(0., Function(self.domain))
            AFyx = Scalar(0., Function(self.domain))  
             
            for ind2 in range(len(sensors)):
                if wxy[ind2]==1:
                    self.totData += 1                                    # add 1 to total of measurements
                    wDxy += Dxy[ind2]*adfn[ind2]                         # data to FUNCTION space                   
                    scxy += 1./abs(Exy[ind2])**2*adfn[ind2]   # scale to FUNCTION space
                    AFxy += adfn[ind2]                                   # where measurements FUNCTION 
                    self.totDataArea += integrate(adfn[ind2])            # add area to total measurement area 
                if wyx[ind2]==1:
                    self.totData += 1                                    # add 1 to total of measurements   
                    wDyx += Dyx[ind2]*adfn[ind2]                         # data to FUNCTION space                                  
                    scyx +=  1./abs(Eyx[ind2])**2*adfn[ind2]                  # scale to FUNCTION space            
                    AFyx += adfn[ind2]                                   # where measurements FUNCTION
                    self.totDataArea += integrate(adfn[ind2])            # add area to total measurement area    
                        
            self.wDxys.append(wDxy)      # DATA   complex FUNCTION space
            self.wDyxs.append(wDyx)       
            self.scxy.append(scxy)       # SCALE real FUNCTION space
            self.scyx.append(scyx)         
            self.AFxys.append(AFxy)      # real FUNCTION space: 1 where measurement (FOR GETTING COMPUTED Zs)
            self.AFyxs.append(AFyx)         
        print("Total number of data points",self.totData)
        print("Total area for applied data", self.totDataArea)            
            
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
        print("done init")
        #1/0
        

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
        return rmY+rmX+rmy 

    def _getNorm(self, x):
        return Lsup(x)

    def setSigmaRho(self,m):
        mC = clip(m, minval=-10., maxval = 10.)
        mF = interpolate(mC, Function(self.domain))    
        air = 1-self.GrF
        self.sigma = self.GrF*exp(mF)*self.sigmaBG
        self.rhoG = self.GrF*exp(-mF)/self.sigmaBG  
        self.rho = self.rhoG + air*1.e14
        sigma = self.GrCF*self.sigmaBG*exp(m)
        self.sigma_boundary=interpolate(sigma, FunctionOnBoundary(self.domain))
        rho = 1./(self.sigmaBG)*self.GrCF*exp(-m) + (1-self.GrCF)*1.e14 
        self.rho_boundary=interpolate(rho, FunctionOnBoundary(self.domain))                   
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
        zCF = self.domain.getX()[self.domain.getDim()-1]
        bCF = inf(zCF)
        tCF = sup(zCF)
        
        for ip in range(self.numP):
            cu=self.mu*self.om[ip] 
            self.TEpde.setValue(D=1j*cu*self.sigma)
            if not self.fixBC:
                ke = (1+1j)*sqrt(cu*self.sigma_boundary/2.)
                self.TEpde.setValue(d=ke*whereZero(z-b), r=(zCF-bCF)/(tCF-bCF))                  

            Ex = self.TEpde.getSolution() 
            gEx = grad(Ex, Function(self.domain))
            Hy = -1./(1j*cu)*gEx[1]
            Zxy = safeDiv(interpolate(Ex, Hy.getFunctionSpace()), Hy)
            Exs.append(Ex)
            Zxys.append(Zxy)
            
            self.TMpde.setValue(A=self.rho*kronecker(2))
            self.TMpde.setValue(D=1j*cu)
            if not self.fixBC:
                kh = (1+1j)*sqrt(cu*self.rho_boundary/2.)
                self.TMpde.setValue(d=kh*whereZero(z-b))            
            Hx = self.TMpde.getSolution()
            gHx = grad(Hx, Function(self.domain))
            Ey = self.rho*gHx[1] 
            Zyx = safeDiv(Ey,interpolate(Hx, Ey.getFunctionSpace()))
            Hxs.append(Hx)
            Zyxs.append(Zyx)
        return Exs, Hxs, Zxys, Zyxs        
        
    def _getValue(self, m, *args):
        print("get Value m",m)
        if len(args)==0:
            args=self.getArguments(m)
        Exs, Hxs, Zxys, Zyxs = args        
        normTE = 0.
        normTM = 0.
        for ip in range(self.numP):
            DEFxy = (length(Zxys[ip]*self.AFxys[ip] - self.wDxys[ip]))**2*self.scxy[ip]
            normTE += integrate(DEFxy)            
            DEFyx = (length(Zyxs[ip]*self.AFyxs[ip] - self.wDyxs[ip]))**2*self.scyx[ip]
            normTM += integrate(DEFyx)
        gm = grad(m) 
        A1 = self.wa*self.w1*integrate(gm[0]**2) +  self.w1*integrate(gm[0]**2)   
        #A1=self.w1*integrate((length(grad(m)))**2)
        A0=self.w0*integrate(interpolate(m, Function(self.domain))**2) 
        print ("reg, mf_xy, mf_yx, totmf =",A0/2.+A1/2.,normTE/2.,normTM/2.,normTE/2.+ normTM/2.)
        return (A0+A1+normTE+normTM)/2            
 
    def _getGradient(self, m, *args):
        #print("getGradient m", m)
        if len(args)==0: 
            args=self.getArguments(m)
        Exs, Hxs, Zxys, Zyxs = args
        gm = grad(m) 
        FS = gm.getFunctionSpace()
        X= Data(0.,(2,),FS)
        X[0] = gm[0]*self.w1*self.wa
        X[1] = gm[1]*self.w1 

        Y = self.w0*interpolate(m,FS)
        y = Scalar(0., FunctionOnBoundary(self.domain))

        self.adjTMpde.setValue(A=self.rho*kronecker(self.domain.getDim()))          
        
        z = FunctionOnBoundary(self.domain).getX()[self.domain.getDim()-1]
        b = inf(z)
        
        for ip in range(self.numP):
            cu = self.om[ip]*self.mu 
            
            # TE
            self.adjTEpde.setValue(D=1j*cu*self.sigma)
            DEFxy = Zxys[ip]*self.AFxys[ip] - self.wDxys[ip]      #DEFxy = Zxys[ip] - self.wDxys[ip]
            gEx = grad(Exs[ip], Function(self.domain))                
            hatDxy = safeDiv(self.scxy[ip]*DEFxy.conjugate(), gEx[self.domain.getDim()-1])
            YTE = -1j*cu*self.AFxys[ip]*hatDxy
            XTE = -self.AFxys[ip]*hatDxy*Zxys[ip]*kronecker(2)[1]
            self.adjTEpde.setValue(Y=YTE, X=XTE)
            if not self.fixBC:
                ke = (1.+1j)*sqrt(cu*self.sigma_boundary/2.)            
                self.adjTEpde.setValue(d=ke*whereZero(z-b))            
            ExStar = self.adjTEpde.getSolution()
            EstarE = Exs[ip]*ExStar    
            Y += (-1j*cu*self.sigma*interpolate(EstarE,FS)).real()
            if not self.fixBC:
                EstarE_b = interpolate(EstarE, FunctionOnBoundary(self.domain))
                y += (ke*EstarE_b).real()
            
            # TM
            self.adjTMpde.setValue(D=1j*cu)
            DEFyx = Zyxs[ip]*self.AFyxs[ip] - self.wDyxs[ip]              
            #DEFyx = Zyxs[ip] - self.wDyxs[ip]          
            hatDyx =safeDiv(self.scyx[ip]*DEFyx.conjugate(), interpolate(Hxs[ip], DEFyx.getFunctionSpace()))
            YTM = -self.AFyxs[ip]*hatDyx*Zyxs[ip]
            XTM = self.rho*self.AFyxs[ip]*hatDyx*kronecker(2)[1]
            self.adjTMpde.setValue(Y=YTM, X=XTM)
            if not self.fixBC:            
                kh = (1.+1j)*sqrt(cu*self.rho_boundary/2.) 
                self.adjTMpde.setValue(d=kh*whereZero(z-b))
            HxStar = self.adjTMpde.getSolution()
            g1=grad(Hxs[ip], FS)
            g2=grad(HxStar, FS)
            gg = inner(g1,g2)  
            bob = self.AFyxs[ip]*hatDyx*g1[1]  
            Y += (self.rhoG*(gg - bob)).real()
            if not self.fixBC:
                HstarH_b = interpolate(Hxs[ip]*HxStar, FunctionOnBoundary(self.domain)) 
                y += (kh*HstarH_b).real()     
                       
        return ArithmeticTuple(Y,X,y) 

    def _getInverseHessianApproximation(self, m, r, *args):
        self.Hpde.setValue(X=r[1], Y=r[0])
        return self.Hpde.getSolution() 
 
 
 
## -------------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='2DMT BFGS', epilog="a.codd@uq.edu.au")
parser.add_argument(dest='config', metavar='configfile', type=str, help='python configuration')
args = parser.parse_args()

print()
print("** 2D isotropic MT BFGS inversion Testing**")

config = importlib.import_module(args.config)

print("Configuration "+args.config+" imported.")
if config.fixBC:
    print("fixed bottom boundary")
else:
    print("1D wave equation bottom boundary")
print()
# domain, sensors, frequencies, background density
domain = ReadMesh(config.flyfile, numDim=2)
GrCF = MaskFromTag(domain,"core","buffer")     # CF
GrCF.expand()
sigmaBG = config.sigmaBG

# data 
dataScale = config.dataScale  
Dxys = np.load(config.dataZxyName)*dataScale*config.signTE
Dyxs = np.load(config.dataZyxName)*dataScale*config.signTM
wxys = np.load(config.wxyName)
wyxs = np.load(config.wxyName) 
freqs = np.load(config.freqsName)
sensors = np.load(config.SensorsName)

if config.withEps:
    print("eps errors")
    eps=config.eps
    Exys = abs(Dxys)*eps
    Eyxs = abs(Dyxs)*eps
else:
    print("given errors")
    Exys = np.load(config.ExyName)*dataScale
    Eyxs = np.load(config.EyxName)*dataScale    

print("number of frequencies", len(freqs))
print("min frequency", min(freqs))
print("max frequency", max(freqs))
print("number of sensors",len(sensors))

print()

m0 = Scalar(0., ContinuousFunction(domain))

if config.Reduced: # reduce sensors and frequencies

    print("w1, wa",config.w1, config.wa)
    print("inversion 2D MT fixed base reduced frequencies Reduced sensors and frequencies")
    
    firstS = config.firstS
    lastS = config.lastS
    stepS = config.stepS
    firstF = config.firstF 
    lastF = config.lastF
    stepF = config.stepF
    
    redSensors = sensors[firstS:lastS:stepS]
    redFreqs = freqs[firstF:lastF:stepF]
    
    print("number of reduced sensors",len(redSensors))
    print("number of reduced frequencies",len(redFreqs))
    print("min and max frequency", min(redFreqs),max(redFreqs))
    
    xyD = Dxys[firstF:lastF:stepF] 
    xyW = wxys[firstF:lastF:stepF]
    xyE = Exys[firstF:lastF:stepF]
    
    yxD = Dyxs[firstF:lastF:stepF]
    yxW = wyxs[firstF:lastF:stepF]
    yxE = Eyxs[firstF:lastF:stepF]
    
    redDxys=[]
    redWxys=[]
    redExys=[]
    
    redDyxs=[]
    redWyxs=[]
    redEyxs=[]
        
    for dxy in xyD:
        redDxys.append(dxy[firstS:lastS:stepS])
    for wxy in xyW:
        redWxys.append(wxy[firstS:lastS:stepS])
    for exy in xyE:
        redExys.append(exy[firstS:lastS:stepS])
                        
    for dyx in yxD:
        redDyxs.append(dyx[firstS:lastS:stepS])
    for wyx in yxW:
        redWyxs.append(wyx[firstS:lastS:stepS])
    for eyx in yxE:
        redEyxs.append(eyx[firstS:lastS:stepS])
                          
    costFnRed = MTInv(domain, sigmaBG, redSensors, redFreqs, m0, 
                       Dxys = redDxys, Dyxs = redDyxs,  Wxys = redWxys, Wyxs = redWyxs,
                       Exys = redExys, Eyxs = redEyxs, fixBC = config.fixBC,
                       w0 = config.w0, w1 = config.w1, wa=config.wa, mu = config.mu)
                       
    args0 = costFnRed.getArguments(m0)
    redExs0, redHxs0, redZxys0, redZyxs0 = args0
                     
    solver = MinimizerLBFGS(J=costFnRed, imax=1000)
    solver.run(m0)
    m=solver.getResult()
    print("m",m)
    print("sigmaBG", sigmaBG)
    sigma = GrCF*exp(m)*sigmaBG
    rhoG = GrCF*exp(-m)/sigmaBG  

    saveSilo(config.silofile, m=m, sigma=sigma, rhoG=rhoG)
   
    redargs = costFnRed.getArguments(m)
    redExs, redHxs, redZxys, redZyxs = redargs
    
    redwDxys = costFnRed.wDxys
    redwDyxs = costFnRed.wDyxs
    redAFxys = costFnRed.AFxys        
    redAFyxs = costFnRed.AFyxs   
    redscxy = costFnRed.scxy
    redscyx = costFnRed.scyx
    redtotData=costFnRed.totData       
    
    redloc = Locator(ReducedFunction(domain),[s for s in redSensors])      
    redfinalZxy=[]
    redfinalZyx=[]
    redtotSxy=0.
    redtotSyx=0.
    redxySum=0.
    redyxSum=0.
    redinitZxy=[]
    redinitZyx=[]
    redxySum0=0.
    redyxSum0=0.    
        
    for ip in range(len(redFreqs)):
        redDxy = redDxys[ip]
        redExy = redExys[ip]
        redZxy = redloc(redZxys[ip]*redAFxys[ip])
        redZxy0 = redloc(redZxys0[ip]*redAFxys[ip])
        redfinalZxy.append(redZxy)
        redinitZxy.append(redZxy0)   
        redWxy = redWxys[ip]
        
        redDyx = redDyxs[ip]
        redEyx = redEyxs[ip]
        redZyx = redloc(redZyxs[ip]*redAFyxs[ip])        
        redZyx0 = redloc(redZyxs0[ip]*redAFyxs[ip]) 
        redfinalZyx.append(redZyx)
        redinitZyx.append(redZyx0)        
        redWyx = redWyxs[ip]   
                
        for k in range(len(redZxy)):
            reddxy = redDxy[k]
            redexy = redExy[k]
            redzxy = redZxy[k]
            redzxy0 = redZxy0[k]
            redwxy = redWxy[k]
            print("reddxy  redzxy redzxy0 redwxy")
            print(reddxy,  redzxy, redzxy0, redwxy)
            if redwxy == 1.:
                redtotSxy +=1.
                redxySum += abs(reddxy-redzxy)**2/(redexy)**2  
                redxySum0 += abs(reddxy-redzxy0)**2/(redexy)**2
                               
            reddyx = redDyx[k]
            redeyx = redEyx[k]
            redzyx = redZyx[k]
            redzyx0 = redZyx0[k]
            redwyx = redWyx[k]
            print("reddyx  redzyx redzyx0 redwyx")
            print(reddyx,  redzyx, redzyx0, redwyx)
            if redwyx == 1.:
                redtotSyx += 1.
                redyxSum += abs(reddyx-redzyx)**2/(redeyx)**2                  
                redyxSum0 += abs(reddyx-redzyx0)**2/(redeyx)**2 
                     
    print("Sum(|dxy-rxy|**2/|exy|**2) final init",redxySum, redxySum0)
    print("Sum(|dyx-ryx|**2/|eyx|**2) final init",redyxSum, redyxSum0)
    print("Total data xy", redtotSxy)
    print("Total data yx", redtotSyx)
    print("errorRMS xy",np.sqrt(redxySum)/redtotSxy, np.sqrt(redxySum0)/redtotSxy)   
    print("errorRMS yx",np.sqrt(redyxSum)/redtotSyx, np.sqrt(redyxSum0)/redtotSyx)
    print("errorRMS ",np.sqrt(redxySum+redyxSum)/(redtotSxy+redtotSyx), np.sqrt(redxySum0+redyxSum0)/(redtotSxy+redtotSyx))   
    
    np.save(config.FinalZxyName, redfinalZxy)
    np.save(config.FinalZyxName, redfinalZyx)          
    np.save(config.InitZxyName, redZxy0)
    np.save(config.InitZyxName, redZyx0)       
     
else: #test ALL
    print("w1, wa",config.w1, config.wa)
    print("inversion 2D MT fixed base ALL frequencies and ALL sensors")
    print("min and max frequency", min(freqs), max(freqs))
               
    costFn = MTInv(domain, sigmaBG, sensors, freqs, m0, 
                       Dxys = Dxys, Dyxs = Dyxs,  Wxys = wxys, Wyxs = wyxs, 
                       Exys = Exys, Eyxs = Eyxs, fixBC = config.fixBC,
                       w0 = config.w0, w1 = config.w1, wa=config.wa, mu = config.mu)

    args0 = costFn.getArguments(m0)
    Exs0, Hxs0, Zxys0, Zyxs0 = args0
  
    solver = MinimizerLBFGS(J=costFn, imax=1000)
    solver.run(m0)
    m=solver.getResult()
    sigma = GrCF*exp(m)*sigmaBG
    rhoG = GrCF*exp(-m)/sigmaBG  
    print("m",m)
    print("sigmaBG", sigmaBG)

    saveSilo(config.silofile, m=m, sigma=sigma, rhoG=rhoG)

  
    #for nf in range(len(redFreqs)):
    #    loc = Locator(ReducedFunction(domain),[s for s in redSensors])    
    args = costFn.getArguments(m)
    Exs, Hxs, Zxys, Zyxs = args
    wDxys = costFn.wDxys
    wDyxs = costFn.wDyxs
    AFxys = costFn.AFxys      # Array of real Functions: 1 at measurements for each frequency  
    AFyxs = costFn.AFyxs   
    scxy = costFn.scxy
    scyx = costFn.scyx
    totData=costFn.totData
 
    loc = Locator(ReducedFunction(domain),[s for s in sensors]) 
           
    finalZxy=[]
    finalZyx=[]
    totSxy=0.
    totSyx=0.
    xySum=0.
    yxSum=0.
    initZxy=[]
    initZyx=[]
    xySum0=0.
    yxSum0=0.       
    
    for ip in range(len(freqs)):
        Dxy = Dxys[ip]
        Exy = Exys[ip]
        Zxy = loc(Zxys[ip]*AFxys[ip])
        Zxy0 = loc(Zxys0[ip]*AFxys[ip])        
        finalZxy.append(Zxy)
        initZxy.append(Zxy0)   
        Wxy = wxys[ip]
        
        Dyx = Dyxs[ip]
        Eyx = Eyxs[ip]
        Zyx = loc(Zyxs[ip]*AFyxs[ip])        
        Zyx0 = loc(Zyxs0[ip]*AFyxs[ip]) 
        finalZyx.append(Zyx)
        initZyx.append(Zyx0)        
        Wyx = wyxs[ip]   
                
        for k in range(len(Zxy)):
            dxy = Dxy[k]
            exy = Exy[k]
            zxy = Zxy[k]
            zxy0 = Zxy0[k]
            wxy = Wxy[k]
            if wxy == 1.:
                totSxy +=1.
                xySum += abs(dxy-zxy)**2/(exy)**2  
                xySum0 += abs(dxy-zxy0)**2/(exy)**2
            dyx = Dyx[k]
            eyx = Eyx[k]
            zyx = Zyx[k]
            zyx0 = Zyx0[k]
            wyx = Wyx[k]
            if wyx == 1.:
                totSyx += 1.
                yxSum += abs(dyx-zyx)**2/(eyx)**2                  
                yxSum0 += abs(dyx-zyx0)**2/(eyx)**2 
                     
    print("Sum(|dxy-rxy|**2/|0.05*dxy|**2) final init",xySum, xySum0)
    print("Sum(|dyx-ryx|**2/|0.05*dyx|**2) final init",yxSum, yxSum0)
    print("Total data xy", totSxy)
    print("Total data yx", totSyx)
    print("errorRMS xy",np.sqrt(xySum)/totSxy, np.sqrt(xySum0)/totSxy)   
    print("errorRMS yx",np.sqrt(yxSum)/totSyx, np.sqrt(yxSum0)/totSyx)
    print("errorRMS ",np.sqrt(xySum+yxSum)/(totSxy+totSyx), np.sqrt(xySum0+yxSum0)/(totSxy+totSyx))   
    
    np.save(config.FinalZxyName, finalZxy)
    np.save(config.FinalZyxName, finalZyx)          
    np.save(config.InitZxyName, Zxy0)
    np.save(config.InitZyxName, Zyx0)           
        

        
        
        
        
        
        
        

        



  
        

