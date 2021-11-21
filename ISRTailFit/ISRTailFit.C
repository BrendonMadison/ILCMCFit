  //Reminder to include the macros file by running ".L Macros.C" before you use this file else FitPlot will not work

void ISRTailFit(){
  //
  // Primary purpose: Uses the Kuraev-Fadin (KF), Nicrosini-Trentadue (NT) and Jadach-Ward-Was (JWW) Initial State Radiation (ISR) distributions to fit for DiMuon sqrt(s) distribution of ILC 250 GeV e-e+ collider.
  //
  // Secondary purpose: Since this likely poorly fits the sqrt(s) distribution, as ISR is only 1 component that contributes to the sqrt(s), this gives an idea of how the parameters need to be changed to better fit the sqrt(s) as simply using lepton mass and fine structure constant, while make sense in theory, don't work perfectly in implementation. 
  //
  // This version does not keep the fine structure or lepton mass constant
  //
  // To dump the text output to a file type ".> mylog.txt" in the prompt before the script
  gSystem->Load("libRooFit");
  using namespace RooFit; 

  //The main variable in this case is the sqrt(s) energy
  RooRealVar sqrts("sqrts","#sqrt{s} [GeV]",150,249.9,"GeV");

  //Fit parameter fine structure and lepton mass
  RooRealVar alpha("alpha","#alpha_{EM}",1.0/137,0.1/137,10/137);
  alpha.setConstant(kTRUE);
  alpha.removeError();
  RooRealVar me("Mu","#m_{l}",105.7e-3,5e-6,250e-3);
  //Mean beam energy. Fit for but should be ~250 GeV
  RooRealVar be("Epsilon","#langle #sqrt{s} #rangle",250,240,260);
  //Delta_l is a fundamental parameter but it is fit for here because the calculation is not trivial and depends on the above values as well as the ER zeta function.
  //It is also supposed to be ~=2*alpha*log(sqrts/m_l) ~= 1
  //RooGenericPdf dl("dl","#Delta_{l}","1 + (alpha/3.1415) * (1.5*log((sqrts*sqrts) / (me*me)) - 1/3 * 3.1415*3.1415 - 2)",RooArgSet(alpha,sqrts,me));
  RooRealVar dl("dl","#Delta_{l}",1,0,500);
  //define the beta parameter of the isr fit according to the Kuraev&Fadin and other ISR papers
  RooGenericPdf beta("beta","beta","2*alpha/3.1415926 * (log(sqrts*sqrts/(Mu*Mu)) - 1)",RooArgSet(alpha,sqrts,me));
  //The first PDF to be fit for is KF ISR
  RooGenericPdf KF("KF","KF","beta/16 * ((8 + 3*beta)*(1-sqrts/Epsilon)^(beta/2 - 1) - 4*(1+sqrts/Epsilon))",RooArgSet(beta,sqrts,be));
  //The second PDF to be fit for is NT ISR (to second order in beta)
  RooGenericPdf NT("NT","NT","(dl)*beta*(sqrts/Epsilon)^(beta-1) - 0.5*beta*(2-sqrts/Epsilon) + 0.125*beta*beta*((2-sqrts/Epsilon)*(3*log(1-sqrts/Epsilon) - 4*log(sqrts/Epsilon)) - (4*log(1-sqrts/Epsilon))/(sqrts/Epsilon) - 6 + sqrts/Epsilon)",RooArgSet(dl,beta,sqrts,be));
  //The third PDF to be fit for is JWW ISR (to what they call c or 3rd order correction)
  RooGenericPdf JWW("JWW","JWW","(exp(beta/4 + alpha/3.1415 * (3.1415*3.1415/3 - 0.5)) * exp(0.5772*beta)/(ROOT::Math::tgamma(1+beta)) * beta*(1-sqrts/Epsilon)^(beta-1))*(1+beta/2 - 0.5*(1-sqrts/Epsilon * sqrts/Epsilon) - beta*((1-sqrts/Epsilon)/2 + (1+3*sqrts*sqrts/(Epsilon*Epsilon))*log(sqrts/Epsilon)/8))",RooArgSet(beta,alpha,sqrts,be)); 

  // Access the histogram and tree stored in the root file
  TFile hfile("../../gdev/Test-DiMuon/DiMuons_LR.root");
  TTree* tree=(TTree*)hfile.Get("DiMuonTree");

  //Draws two histograms (twice for 100 and 200 bins) and saves them to csv files for the purpose of Monte-Carlo fitting later
  tree->Draw("(pfoSqrtsPlus)>>htemp(100,233,253)","","goff");
  TH1F *hs1 = (TH1F*) gDirectory->Get("htemp");
  Hist2CSV(hs1,std::string("sqrts100.csv"));
  
  tree->Draw("(pfoSqrtsPlus)>>htemp(200,233,253)","","goff");
  TH1F *hs2 = (TH1F*) gDirectory->Get("htemp");
  Hist2CSV(hs2,std::string("sqrts200.csv"));

  tree->Draw("(1-pfoSqrtsPlus/250)^(1/6)>>htemp(100,0.05,0.95)","","goff");
  TH1F *hf1 = (TH1F*) gDirectory->Get("htemp");
  Hist2CSV(hf1,std::string("sqrtsfrac100.csv"));
  
  tree->Draw("(1-pfoSqrtsPlus/250)^(1/6)>>htemp(200,0.05,0.95)","","goff");
  TH1F *hf2 = (TH1F*) gDirectory->Get("htemp");
  Hist2CSV(hf2,std::string("sqrtsfrac200.csv"));

  //Draw the histogram that is to be fitted for. You can change the options here. The range for x should be the same as used below in htemp.
  tree->Draw("(pfoSqrtsPlus)>>htemp(200,150.0,249.9)","","");

  //Grabs htemp off the directory
  TH1F *h1 = (TH1F*) gDirectory->Get("htemp");

  //Normalizes htemp and changes it name and draws with error bars.
  h1->SetTitle("#sqrt{s}");
  h1->Scale(1.0/(h1->Integral()));
  h1->Draw("E");
  
  Int_t UncutEvents = h1->GetEntries();
  Int_t TotEvents = UncutEvents;

  // Create a binned dataset that imports contents of h1 and associates its contents to observable 'sqrts' aka the beam energy
  RooDataHist db("db", "db", sqrts,Import(*h1));

  //Start with KF
  //create a RooAddPdf for passing the model to the plotting program.
  RooRealVar addfrac("% Signal","addfrac",1.0,0.0,1.0);
  addfrac.setConstant(kTRUE);
  addfrac.removeError();
  RooAddPdf sigAdd("sigAdd","Muon Free Kuraev-Fadin Tail Fit #sqrt{s}",RooArgList(KF),addfrac);

  //Temporary file to pass the PDF
  TFile f("temp.root","recreate");

  db.Write();
  sigAdd.Write();
  sqrts.Write();

  //Check that the PDF and data were saved.
  if(gDirectory->Get("db") == NULL) {
    cout << "Unable to return database." << endl;
  }

  //Calls FitPlot
  std::string TitleArg = "#sqrt{s} KF|KF(#sqrt{s},m_{l},#Epsilon,#alpha)";
  FitPlot(db.GetName(),sigAdd.GetName(),sqrts.GetName(),TitleArg.c_str(),TotEvents,UncutEvents);

  //Closes the temp file
  f.Close();

  //Removes the temp file
  gSystem->Exec("rm temp.root");

  //Next is NT
  //create a RooAddPdf for passing the model to the plotting program.
  RooRealVar addfrac2("% Signal","addfrac2",1.0,0.0,1.0);
  addfrac2.setConstant(kTRUE);
  addfrac2.removeError();
  RooAddPdf sigAdd2("sigAdd2","Muon Free Nicrosini-Trentadue Tail Fit #sqrt{s}",RooArgList(NT),addfrac2);

  //Temporary file to pass the PDF
  TFile f2("temp.root","recreate");

  db.Write();
  sigAdd2.Write();
  sqrts.Write();

  //Check that the PDF and data were saved.
  if(gDirectory->Get("db") == NULL) {
    cout << "Unable to return database." << endl;
  }

  TitleArg = "#sqrt{s} NT|NT(#sqrt{s},m_{l},#Epsilon,#alpha,#Delta_{l})";
  FitPlot(db.GetName(),sigAdd2.GetName(),sqrts.GetName(),TitleArg.c_str(),TotEvents,UncutEvents);

  //Closes the temp file
  f2.Close();

  //Removes the temp file
  gSystem->Exec("rm temp.root");

  //Now JWW
  //create a RooAddPdf for passing the model to the plotting program.
  RooRealVar addfrac3("% Signal","addfrac3",1.0,0.0,1.0);
  addfrac3.setConstant(kTRUE);
  addfrac3.removeError();
  RooAddPdf sigAdd3("sigAdd3","Muon Free Jadach-Ward-Was Tail Fit #sqrt{s}",RooArgList(JWW),addfrac3);

  //Temporary file to pass the PDF
  TFile f3("temp.root","recreate");

  db.Write();
  sigAdd3.Write();
  sqrts.Write();

  //Check that the PDF and data were saved.
  if(gDirectory->Get("db") == NULL) {
    cout << "Unable to return database." << endl;
  }

  TitleArg = "#sqrt{s} JWW|JWW(#sqrt{s},m_{l},#Epsilon,#alpha)";
  FitPlot(db.GetName(),sigAdd3.GetName(),sqrts.GetName(),TitleArg.c_str(),TotEvents,UncutEvents);

  //Closes the temp file
  f3.Close();

  //Removes the temp file
  gSystem->Exec("rm temp.root");

  //Close the histogram file
  hfile.Close();

}
