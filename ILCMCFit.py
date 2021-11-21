#Written by Brendon Madison of Univ. of Kansas Fall 2021

import numpy as np
import matplotlib.pyplot as plt

#Need beta to create the beamstrahlung distributions
from scipy.stats import beta
#Need gamma to create the JWW ISR distribution
from scipy.special import gamma

#Import the ConfigParser function
#It imports on 'configparser' and 'os' and 'unicodedata' and 'shutil' and 'numpy'
#so make sure you have those installed!
from MCMacros import ImportConfig
#Import the function for calculating Chi2
from MCMacros import GetChi2

#Use config parser to import dictionaries with settings, values
#Most of the numerical and simulation values come from these dictionaries
st,mc,ft,sm,ex = ImportConfig('MCFit.ini',False)

#We define the Monte-Carlo fitting as a function so we can use it in other scripts.
def MCFit(setd,mcv,ftv,smv,exv):

    #Define the (1-sqrts/250)^(1/eta) axis
    tsqrts = np.linspace(0.05,0.99,smv["res"])

    #Calculate the beamstrahlung distribution.
    #Using scipy's function is better at replicating the shape than using the "brute force" method
    #Though I'm not 100% why this is the case...
    circe = ftv["n1"]*beta(ftv["a1"], ftv["g1"]).pdf(tsqrts**ftv["etafit"])
    #circe = ftv["n1"]*(tsqrts**(ftv["etafit"]))**(ftv["a1"]-1) * (1-tsqrts**(ftv["etafit"]))**(ftv["g1"]-1)

    #The asymptotic background in the fractional sqrts tail
    asym = ftv["n2"]/tsqrts**3 * 0.5*(1 + np.tanh(ftv["hslp"]*(tsqrts - ftv["hpos"])))

    #The gaussian distribution used for the events that undergo no beamstrahlung
    unbeams = np.exp(-0.5*(tsqrts-0)**2 / (0.05)**2)
    unbeams = unbeams / np.amax(unbeams) * ftv["unev"]

    beams = circe + asym
    beams = beams / np.amax(beams) * ftv["bmev"]
    #adding the beamstrahlung and unbeamstrahlunged. May not be correct to do this considering pbeams is being used...
    #and they essentially do the same thing...
    #beams = beams + unbeams

    #The following section is for ISR. Either use Kuraev-Fadin or Jadach-Ward-Was
    #In the software there is an energy cutoff for ISR. We don't yet know exactly what it is but lets consider 1 GeV
    #This would corresponds to tsqrt ~=0.4 but the below calculates it
    #I've also noticed that, in the ILCSOFT data there is a sharp peak at 0.4, maybe from the low energy cutoff. 
    #There is also a discontinuity around 0.9. Possible high energy cutoff?

    #Create the cutoff and then the "x axis"
    CUTOFFT = (1.0-(exv["sq"]-mcv["cutE"])/exv["sq"])**(1/ftv["eta"])
    isrT = np.linspace(CUTOFFT,0.99,smv["res"])
    sqrts = exv["sq"]*(1-isrT**ftv["eta"])

    #The "beta" parameter that is used in the ISR distributions
    kfb = 2.0*mcv["aem"]/np.pi * (np.log((sqrts**2)/(mcv["ml"]**2)) - 1)

    #Check to see which ISR is used
    if mcv["useJWW"] == True:
        #Use the Jadach-Ward-Was ISR distribution
        ISR = np.exp(kfb/4 + mcv["aem"]/3.1415 * (3.1415*3.1415/3 - 0.5)) * np.exp(0.5772*kfb)/(gamma(1+kfb)) * kfb*(1-sqrts/exv["sq"])**(kfb-1) * (1+kfb/2 - 0.5*(1-sqrts/exv["sq"] * sqrts/exv["sq"]) - kfb*((1-sqrts/exv["sq"])/2 + (1+3*sqrts*sqrts/(exv["sq"]**2))*np.log(sqrts/exv["sq"])/8))
    elif mcv["useJWW"] == False:
        #Use the Kuraev-Fadin ISR distribution
        ISR = kfb/16.0 * ((8.0 + 3*kfb)*(isrT)**(ftv["etafit"]*(kfb/2 - 1)) - 4*(2-isrT**ftv["etafit"]))

    #Beam spread distribution as a function of fractional sqrts
    bsprdt = np.linspace(-10,10,smv["res"])
    bsprd = np.exp(-0.5*(bsprdt - 0)**2/(mcv["bes"])**2)

    #Normalize the distributions for the purpose of random choice
    ISRnorm = ISR/np.sum(ISR)
    beamsnorm = beams/np.sum(beams)
    bsprd = bsprd / np.sum(bsprd)

    #Import the external data. We start with 100 bin data
    #As a reminder, use fracdat[0] to access the first row and fracdat[:,0] to access the first column.
    #The data is in format of x,y,yerr
    if smv["use200"] == True:
        fracdat = np.genfromtxt('ISRTailFit/sqrtsfrac200.csv', delimiter=',')
        sqdat = np.genfromtxt('ISRTailFit/sqrts200.csv', delimiter=',')
        binres = sqdat[:,0][1] - sqdat[:,0][0] 
        lorange = np.amin(sqdat[:,0]) - 0.5*binres
        hirange = np.amax(sqdat[:,0]) - 0.5*binres
        nbins = int(200)
    elif smv["use200"] == False:
        fracdat = np.genfromtxt('ISRTailFit/sqrtsfrac100.csv', delimiter=',')
        sqdat = np.genfromtxt('ISRTailFit/sqrts100.csv', delimiter=',')
        binres = sqdat[:,0][1] - sqdat[:,0][0] 
        lorange = np.amin(sqdat[:,0]) - 0.5*binres
        hirange = np.amax(sqdat[:,0]) - 0.5*binres
        nbins = int(100)

    #Randomly generate some energies/sqrts!
    ransqrt = []
    ranfrac = []

    #You should use the same number of events in the data as there is for the Monte-Carlo model
    #Then they have the same cumulative statistics and you don't need to normalize... 
    TOTEVENTS = int(np.sum(sqdat[:,1]))
    sqnormed = sqdat[:,1]/TOTEVENTS * smv["simev"]
    fracnormed = fracdat[:,1]/TOTEVENTS * smv["simev"]
    TOTEVENTS = int(smv["simev"])

    num = 0
    #Check counter for checking if the simulated region has an acceptable number of events 
    checkcount = 0

    while num < TOTEVENTS:
        ranbeam = np.random.choice(tsqrts,1,p=beamsnorm)
        ranisr = np.random.choice(isrT,1,p=ISRnorm)
        
        ranbeam = np.random.choice([ranbeam[0],0.0],p=[mcv["pbeam"],1.0-mcv["pbeam"]])
        ranisr = np.random.choice([ranisr[0],0.0],p=[mcv["pisr"],1.0-mcv["pisr"]])
        
        ranbs = np.random.choice(bsprdt,1,p=bsprd)
        
        beamsqrts = 1.0 - np.sign(ranbeam)*(ranbeam)**ftv["eta"]
        isrsqrts = 1.0 - np.sign(ranisr)*(ranisr)**ftv["eta"]
        
        #our "monte-carlo" generated values of sqrts and fractional sqrts
        if smv["useadd"] == True:
            mcsqrts = exv["sq"]*(1 - (1-beamsqrts) - (1-isrsqrts)) + ranbs
        elif smv["useadd"] == False:
            mcsqrts = exv["sq"]*(beamsqrts*isrsqrts) + ranbs
        mcfrac = (1.0-mcsqrts/exv["sq"])**(1.0/ftv["eta"])
        
        ranfrac.append(mcfrac)
        ransqrt.append(mcsqrts)
        
        checkcount += 1
        if checkcount == smv["check"]:
            hists, binedges = np.histogram(ransqrt,bins=sqdat[:,0].size,range=[lorange,hirange])
            num = np.sum(hists)
            checkcount = 0

    #array-ify
    ransqrt = np.array(ransqrt)
    ranfrac = np.array(ranfrac)

    #Generate histograms for plotting and chi2 computation
    hists, binedges = np.histogram(ransqrt,bins=sqdat[:,0].size,range=[lorange,hirange])
    histf, bef = np.histogram(ranfrac,bins=fracdat[:,0].size,range=[np.amin(fracdat[:,0]),np.amax(fracdat[:,0])])

    #computing the error bars using Poisson error bars
    sqmodstd = np.sqrt(hists)
    sdaterr = np.sqrt(sqnormed)
    fracmodstd = np.sqrt(histf)
    fracdatstd = np.sqrt(fracnormed)

    #when there is no data in a bin it gives nan or zero. It should be zero. nan is bad.
    sqmodstd[np.isnan(sqmodstd)] = 0.0
    sdaterr[np.isnan(sdaterr)] = 0.0
    fracmodstd[np.isnan(fracmodstd)] = 0.0
    fracdatstd[np.isnan(fracdatstd)] = 0.0

    #Compute the Chi2 given the model and data as well as their estimated error values
    Chi2Val = GetChi2(sqnormed,hists,sdaterr,sqmodstd) 
    Chi2Frac = GetChi2(fracnormed,histf,fracdatstd,fracmodstd)

    #Put the output of the simulation into a csv file

    if setd["verbose"] == True:

        print("Total number of simulated events: ",TOTEVENTS)
        print("Degrees of Freedom: ", int(sqdat[:,0].size - 1 - mcv["dof"]))
        print("Model with sqrts has Reduced Chi2 of ",Chi2Val)
        print("Model with frac. sqrts has Reduced Chi2 of ",Chi2Frac)


    if setd["plots"] == True:
    
        #sqrts plot
        fig = plt.figure(figsize=(6,9),dpi=100)
        gs = fig.add_gridspec(2,1,hspace=0.5,wspace=0.75)
        axs = gs.subplots()
        axs[0].errorbar(sqdat[:,0],hists,yerr=sqmodstd,alpha=0.5,label='model',markersize=1, capsize=2)
        axs[0].errorbar(sqdat[:,0],sqnormed,yerr=sdaterr,alpha=0.5,label='data',markersize=1, capsize=2)
        axs[0].set(xlabel=r'$\sqrt{s}$ (GeV)',ylabel='Arbitrary Count',title=r'MC $\sqrt{s}$ Fit')
        axs[0].legend(loc=1, bbox_to_anchor=(0.35,0.95))
        axs[0].annotate('Chi2 = '+str(int(Chi2Val*100.0)/100.0), xy=(0.65, 0.95), xycoords='axes fraction')

        #Residual plot
        axs[1].errorbar(sqdat[:,0],sqnormed-hists,np.sqrt(sqmodstd**2 + sdaterr**2),alpha=0.5,label='Residual',markersize=1, capsize=2)
        axs[1].set(xlabel=r'$\sqrt{s}$ (GeV)',ylabel='Residual',title=r'MC $\sqrt{s}$ Residual')
        plt.savefig(str(setd["name"])+"/"+str(setd["runname"])+"_sqrts_"+str(setd["runnum"])+".png")
        plt.close()
        
        plt.errorbar(fracdat[:,0],histf,yerr=fracmodstd,alpha=0.5,label='model')
        plt.errorbar(fracdat[:,0],fracnormed,yerr=fracdatstd,alpha=0.5,label='data')
        
        #Fractional plot
        fig = plt.figure(figsize=(6,9),dpi=100)
        gs = fig.add_gridspec(2,1,hspace=0.5,wspace=0.75)
        axs = gs.subplots()
        axs[0].errorbar(fracdat[:,0],histf,yerr=fracmodstd,alpha=0.5,label='model',markersize=1, capsize=2)
        axs[0].errorbar(fracdat[:,0],fracnormed,yerr=fracdatstd,alpha=0.5,label='data',markersize=1, capsize=2)
        axs[0].set(xlabel=r'$\sqrt{s}$ (GeV)',ylabel='Arbitrary Count',title=r'MC $\sqrt{s}$ Fit')
        axs[0].legend(loc=1, bbox_to_anchor=(0.35,0.95))
        axs[0].annotate('Chi2 = '+str(int(Chi2Frac*100.0)/100.0), xy=(0.65, 0.95), xycoords='axes fraction')
        
        #Residual plot
        axs[1].errorbar(fracdat[:,0],fracnormed-histf,yerr=np.sqrt(fracmodstd**2 + fracdatstd**2),alpha=0.5,label='Residual',markersize=1, capsize=2)
        axs[1].set(xlabel=r'$(1-x_{\sqrt{s}})^{1/\eta}$ (unitless)',ylabel='Residual',title=r'MC Fractional $\sqrt{s}$ Residual')
        plt.savefig(str(setd["name"])+"/"+str(setd["runname"])+"_fracsq_"+str(setd["runnum"])+".png")
        plt.close()
    
    #return the Chi2 values
    return Chi2Val,Chi2Frac

XX1,XX2 = MCFit(st,mc,ft,sm,ex)