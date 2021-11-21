# ILCMCFit
Includes various files on analyzing Guinea Pig (GP) simulation data to characterize beamstrahlung, analyzing ILCSOFT dimuon data to characterize ISR and sqrt(s), and then a Monte-Carlo (MC) fitting program that uses the previous fits and analysis to fit for sqrt(s) in a way that allows us to back out the energy spread and contributions from beamstrahlung and ISR.

Results:<br />
0.) Fitted GP data for various fits. Found in GPFracS directory. Focused on fitting (1-sqrt(s)/250)^(1/eta) with eta=6. This was so that sqrt(s) could be pulled out of the GP data fitting and used in the Monte-Carlo. <br />
1.) The GP best fit was 1 beta distribution with a 1/x^3 "background" function. Similar to this was the 2 beta distribution fit, which used more parameters for similar Chi2.<br />
2.) Fitted ILCSOFT data for three kinds of ISR distributions. Found in ISRTailFit directory. Kuraev-Fadin (KF) and Jadach-Ward-Was (JWW) performed similarly. JWW seems to be slightly better but unsure. So allowed the use of both in future code.<br />
3.) Created a Monte-Carlo algorithm, in python, to fit for both sqrt(s) and (1-sqrt(s)/250)^(1/eta) to roughly 1.5 reduced Chi2 in both.<br />
4.) The Chi2 minimization algorithm used was an adaptive random walk. So the minimization took a while and, technically, could be minimized further.<br />
4.5.) Used both sqrt(s) and (1-sqrt(s)/250) in fitting and tried to minimize both Chi2. It is argued that this is good as sqrt(s) fits the the peak region better while (1-sqrt(s)/250) is good for tail region.<br />
5.) Result indicates ~90% of events have 1 significant beamstrahlung event and ~24% of events have an ISR photon above a minimum energy threshold of ~3.5 GeV<br />
6.) The energy spread was fitted to ~0.88 GeV using an underlying Gaussian distribution.<br />

To do?:<br />
1.) Get parameter correlation matrix and parameter variance<br />
2.) Improve fit minimization<br />
3.) Better characterize GP data for different GP settings<br />
4.) Try to determine what the artifact/resolution looking feature in the upper tail region of the (1-sqrt(s)/250)^(1/eta) distribution in the ILCSOFT data is.<br />
5.) Try to edit code to be feasibly to run on a computing cluster. This could be a solution to improving the minimization... just brute force on a cluster haha<br />
