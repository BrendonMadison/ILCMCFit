#Written by Brendon Madison of Univ. of Kansas Fall 2021

import numpy as np

#Import MCFit as we will do a bunch of fits with random values relative to the values given in the .ini file
from ILCMCFit import MCFit

#Import the settings
from MCMacros import ImportConfig
#Macro to change Monte-Carlo values during the scan
from MCMacros import MCRandom
#Macro to write simulation output to a file
from MCMacros import MCWrite
#Macro to average the values of the MC values dictionary
from MCMacros import MCAverage

#Set the number of times to randomly walk through the parameter space
nwalk = 500

st,mc,ft,sm,ex = ImportConfig('MCFit.ini',False)

#Keep the best chi2 and mc values.
#best1 is for sqrts and best2 is for fractional sqrts
best1 = 999999
best2 = 999999
run1 = 0
run2 = 0
bestmc1 = mc.copy()
bestmc2 = mc.copy()

for i in range(nwalk):
    #Change the MC values randomly in the range specified in the .ini file
    mc = MCRandom(mc,sm)
    #Get the Chi2 of the fit. May also print things to console and save plots
    X1,X2 = MCFit(st,mc,ft,sm,ex)
    #Write the simulation values to the csv file
    MCWrite(X1,X2,st,mc,sm)
    if X1 < best1:
        best1 = X1
        bestmc1 = mc.copy()
        run1 = st["runnum"]
    if X2 < best2:
        best2 = X2
        bestmc2 = mc.copy()
        run2 = st["runnum"]
    #increment the run number
    st["runnum"] += 1
    
    #Set the values of mc to the average of bestmc1 and bestmc2
    #So the hope is that we randomly integrate to better values over time
    #We do the average because each fit is better for fitting certain parts
    #Fitting for fractional sqrt(s) is good for fitting the lower tail and peak lower side/spread
    #Fitting for sqrt(s) is better for fitting the peak upper side/spread and shape
    mc = MCAverage(mc,bestmc1,bestmc2)
    
    #This is optional as I'm not sure its better but its slowly makes the step size smaller
    #With this the step size decreases by 0.5% so, after 100 steps, the step size is ~60% the original via (1-0.5%)^(100)
    sm["scan"] *= 0.995
    
#After random walking through the parameter space print out the "best fit" values
print("The lowest Chi2 value for sqrts: ",best1)
print("This occured on run number: ",run1)
print("The settings for lowest sqrts: ")
for j in bestmc1:
    print (j,':',bestmc1[j])
    
print("The lowest Chi2 value for fractional sqrts: ",best2)
print("This occured on run number: ",run2)
print("The settings for lowest fractional sqrts: ")
for j in bestmc2:
    print (j,':',bestmc2[j])