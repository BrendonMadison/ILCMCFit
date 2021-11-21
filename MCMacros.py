#Written by Brendon Madison of Univ. of Kansas Fall 2021

import configparser
import os
import unicodedata as unc

from numpy import sqrt
from numpy import isnan
from numpy.random import uniform

from shutil import copyfile

def ImportConfig(FileName,Verb):
    cfg = configparser.ConfigParser()
    cfg.read(FileName)

    #The settings dictionary reads from the .ini file (cfg)
    setd = {
      "name": cfg['settings']['name'],
      "runname": cfg['settings']['runname'],
      "runnum": int(cfg['settings']['runnum']),
      "plots": cfg.getboolean('settings','plots'),
      "verbose": cfg.getboolean('settings','verbose')
    }

    #Tries to create a new directory for this 'name'
    #Then tries to copy the .ini file you just used into 'name' directory
    #So that, later on after you've done a bunch of runs, you have a record of the settings of each run
    try:
        os.mkdir(setd["name"])
        if Verb == True:
            print("Created directory with name: ",setd["name"])
        try:
            copyfile(FileName,str(setd["name"]+"/"+setd["runname"]+".ini"))
            if Verb == True:
                print("Created copy of .ini file in ",setd["name"]," named ",str(setd['runname']+".ini"))
        except:
            if Verb == True:
                print(setd["name"]," already has copy of ",setd["runname"],".ini in it.")
    except:
        if Verb == True:
            print("Directory of following name already exists: ",setd["name"])
        try:
            copyfile(FileName,str(setd["name"]+"/"+setd["runname"]+".ini"))
            if Verb == True:
                print("Created copy of .ini file in ",setd["name"]," named ",str(setd['runname']+".ini"))
        except:
            if Verb == True:
                print(setd["name"]," already has copy of ",setd["runname"],".ini in it.")

    #The Monte-Carlo (MC) Values (mcv aka mc) dictionary
    #These are what you change as variables for the MC model generation
    #So these are what change if you are trying to fine the best Chi2
    mcv = {
      "cutE": float(cfg['MCValues']['cutE']),
      "aem": float(cfg['MCValues']['aem']),
      "ml": float(cfg['MCValues']['ml']),
      "useJWW": cfg.getboolean('MCValues','useJWW'),
      "bes": float(cfg['MCValues']['bes']),
      "pisr": float(cfg['MCValues']['pisr']),
      "pbeam": float(cfg['MCValues']['pbeam']),
      "dof": int(cfg['MCValues']['dof'])
    }

    #The Fit Values (ftv aka ft) dictionary
    #Imports the outcome and settings of the Guinea Pig (GP) fitting
    ftv = {
      "n1": float(cfg['FitValues']['n1']),
      "n2": float(cfg['FitValues']['n2']),
      "a1": float(cfg['FitValues']['a1']),
      "g1": float(cfg['FitValues']['g1']),
      "hslp": float(cfg['FitValues']['hslp']),
      "hpos": float(cfg['FitValues']['hpos']),
      "etafit": float(cfg['FitValues']['etafit']),
      "eta": float(cfg['FitValues']['eta']),
      "totev": float(cfg['FitValues']['totev']),
      "bmev": float(cfg['FitValues']['bmev']),
      "unev": float(float(cfg['FitValues']['totev']) - float(cfg['FitValues']['bmev'])),
      "fitlo": float(cfg['FitValues']['fitlo']),
      "fithi": float(cfg['FitValues']['fithi'])
    }

    #The Simulation (sim) Values (smv aka sm) dictionary
    #These are important for determining how slow and how long the Monte-Carlo is
    #Also how precise, and how high resolution, it will be in exchange.
    smv = {
      "res": int(cfg['SimValues']['res']),
      "simev": int(cfg['SimValues']['simev']),
      "check": int(cfg['SimValues']['check']),
      "use200": cfg.getboolean('SimValues','use200'),
      "scan": float(cfg['SimValues']['scan']),
      "useadd": cfg.getboolean('SimValues','useadd')
    }

    #The experimental values dictionary
    #This is for actual experimental values not simulation or fitting values that relate to experiment
    exv = {
      "bep": float(cfg['ExpValues']['bep']),
      "bem": float(cfg['ExpValues']['bem']),
      "sq": float(2.0*sqrt(float(cfg['ExpValues']['bep']) * float(cfg['ExpValues']['bem']))),
      "Ne": int(float(cfg['ExpValues']['Ne'])),
      "Np": int(float(cfg['ExpValues']['Np'])),
      "ml": float(cfg['ExpValues']['ml']),
      "Nc": int(float(cfg['ExpValues']['Nc'])),
      "sigz": float(cfg['ExpValues']['sigz']),
      "sigy": float(cfg['ExpValues']['sigy']),
      "sigx": float(cfg['ExpValues']['sigx']),
      "Ca": float(cfg['ExpValues']['Ca'])
    }

    #prints the settings dictionary to the user
    if Verb == True:
        print(f'\nConfiguration for {setd["name"]}:\n')
        print(f'\t ---Settings---\nName\t\t\t=\t{setd["name"]}\nRun Name\t\t=\t{setd["runname"]}')
        print(f'Make & Save Plots\t=\t{setd["plots"]}\nVerbosity\t\t=\t{setd["verbose"]} ')
        print(f'\n\t ---Monte-Carlo Values---\nISR Energy Cutoff (GeV)\t=\t{mcv["cutE"]}\nFine Structure \t\t=\t{mcv["aem"]}')
        print(f'ISR Lepton Mass (GeV)\t=\t{mcv["ml"]}\nUse JWW ISR? \t\t=\t{mcv["useJWW"]}')
        print(f'Energy Spread (GeV)\t=\t{mcv["bes"]}\nPercent ISR\t\t=\t{mcv["pisr"]}\nPercent Beams.\t\t=\t{mcv["pbeam"]}')
        print(f'\n\t ---Guinea Pig Fit Values---\nBeta Dist. Norm.\t=\t{ftv["n1"]}\nx^(-3) Background Norm.\t=\t{ftv["n2"]}')
        print(f'Beta Dist. {unc.lookup("greek small letter alpha")}_1\t\t=\t{ftv["a1"]}\nBeta Dist. {unc.lookup("greek small letter beta")}\t\t=\t{ftv["g1"]}')
        print(f'Heaviside Slope\t\t=\t{ftv["hslp"]}\nHeaviside Position\t=\t{ftv["hpos"]}\n(1-x) Power, {unc.lookup("greek small letter eta")}\t\t=\t{ftv["eta"]}')
        print(f'Fit Power, {unc.lookup("greek small letter eta")}\t\t=\t{ftv["etafit"]}\nFit Lower Cutoff\t=\t{ftv["fitlo"]}\nFit Higher Cutoff\t=\t{ftv["fithi"]}')
        print(f'Total GP Events\t\t=\t{int(ftv["totev"])}\nGP Beams. Events\t=\t{int(ftv["bmev"])}')
        print(f'\n\t ---Simulation Values---\nPoints in arrays\t=\t{smv["res"]}\nEvents per Sim\t\t=\t{smv["simev"]}')
        print(f'Events per check \t=\t{smv["check"]}\nUse 200 Bins?\t\t=\t{smv["use200"]}')
        print(f'Scanning Window \t=\t+-{100*smv["scan"]}%')
        if smv["useadd"] == True:
            print('Energy Summing \t\t=\tAdd')
        elif smv["useadd"] == False:
            print('Energy Summing \t\t=\tProduct')
        print(f'\n\t ---Experiment Values---\nPositron Beam Energy (GeV)\t=\t{exv["bep"]}\nElectron Beam Energy (GeV)\t=\t{exv["bem"]}')
        print(f'Center Mass Energy (GeV)\t=\t{exv["sq"]}\nPositrons per bunch\t\t=\t{exv["Np"]}\nElectrons per bunch\t\t=\t{exv["Ne"]}')
        print(f'Pairs per crossing \t\t=\t{exv["Nc"]}\nLepton mass (GeV) \t\t=\t{exv["ml"]}\nX Bunch Spread (m)\t\t=\t{exv["sigx"]}')
        print(f'Y Bunch Spread (m) \t\t=\t{exv["sigy"]}\nZ Bunch Spread (m) \t\t=\t{exv["sigz"]}\nCrossing Angle (radians)\t=\t{exv["Ca"]}')

    return setd,mcv,ftv,smv,exv
    
def GetChi2(data,model,daterr,moderr):
    #Computes the Chi2 statistic given a data and model histogram that both have error values
    chi2 = 0.0
    for i in range(len(data)):
        #Normally chi2 is (observed - calculated)^2 / (variance)
        #Here we have errors on the model/calculated and data/observed. 
        #So the total variance is the quadratic sum of those errors.
        tempchi = (data[i]-model[i])**2 / (daterr[i]**2 + moderr[i]**2)
        if isnan(tempchi) == True:
            tempchi = 0
        chi2 += tempchi
    
    return chi2
    
def MCRandom(mcv,smv):
    #Takes the mcv dictionary and changes its values randomly according to the smv dictionary variable "scan"
    #So a "scan" value of 0.5 means change the variables uniformly on the range of +- 50%
    
    #If you want to force set the seed then you can do so with this:
    #np.random.seed(5) where 5 is your seed number
    
    mcv["cutE"] *= uniform(1-smv["scan"],1+smv["scan"])
    mcv["aem"] *= uniform(1-smv["scan"],1+smv["scan"])
    mcv["ml"] *= uniform(1-smv["scan"],1+smv["scan"])
    mcv["bes"] *= uniform(1-smv["scan"],1+smv["scan"])
    mcv["pisr"] *= uniform(1-smv["scan"],1+smv["scan"])
    mcv["pbeam"] *= uniform(1-smv["scan"],1+smv["scan"])
    
    #Must keep some values positive and within a window that works
    #e.g. the probabilities have to be greater than 0 and less than 1
    if mcv["pisr"] < 0:
        mcv["pisr"] = 1e-3
    if mcv["pisr"] > 1:
        mcv["pisr"] = 1-1e-3
    if mcv["pbeam"] < 0:
        mcv["pbeam"] = 1e-3
    if mcv["pbeam"] > 1:
        mcv["pbeam"] = 1-1e-3
        
    if mcv["cutE"] < 0:
        mcv["cutE"] = 1e-3
    if mcv["aem"] < 0:
        mcv["aem"] = 1e-9
    if mcv["ml"] < 0:
        mcv["ml"] = 1e-9

    #bes doesn't need a check because it gets squared -- bes^2
        
    #returns the mcv dictionary that has the changed values
    return mcv
    
def MCAverage(mcv,bmc1,bmc2):
    #Takes the mcv dictionary and the "best" mcv dictionaries for fitting sqrts (bmc1) and fractional sqrts (bmc2)
    #Then rewrites mcv to be the average of bmc1 and bmc2
    
    mcv["cutE"] = 0.5*(bmc1["cutE"] + bmc2["cutE"])
    mcv["aem"] = 0.5*(bmc1["aem"] + bmc2["aem"])
    mcv["ml"] = 0.5*(bmc1["ml"] + bmc2["ml"])
    mcv["bes"] = 0.5*(bmc1["bes"] + bmc2["bes"])
    mcv["pisr"] = 0.5*(bmc1["pisr"] + bmc2["pisr"])
    mcv["pbeam"] = 0.5*(bmc1["pbeam"] + bmc2["pbeam"])
    
    #returns the mcv dictionary that has the averaged values
    return mcv
    
def MCWrite(X1,X2,setd,mcv,smv):
    #Takes the Chi2 values, MCvalues and appends them to a csv file.
    
    filelog = open(str(setd["name"])+"/"+str(setd["runname"])+".txt","a")
    #The format is {chi2 of sqrts, chi2 of frac sqrts, cutE, aem, ml, bes, pisr, pbeam, simev
    filelog.write(f'{setd["runnum"]},{X1},{X2},{mcv["cutE"]},{mcv["aem"]},{mcv["ml"]},{mcv["bes"]},{mcv["pisr"]},{mcv["pbeam"]},{smv["simev"]}\n')
    filelog.close()