  //Reminder to include the macros file by running ".L Macros.C" before you use this file else FitPlot will not work

void GPFracSFit(){
  //
  // Fits Guinea Pig (GP) beamstrahlung data for the fractional sqrt(s) distributions given an ILC e+e- at sqrts~=250 GeV
  //
  // Two distributions are used. 
  // Inverse gaussian, which describes how a random ensemble randomly drifts towards a mean
  // Beta distribution, which describes random numbers that are constrained in the limit of 0->1
  // 
  // It is assumed the beta distributions will work well as the CIRCE parameterization, which is used to parameterize beamstrahlung, uses beta distributions.
  // See CIRCE paper here: https://arxiv.org/pdf/hep-ph/9607454.pdf 
  // and a more recent presentation here: https://indico.desy.de/event/10353/contributions/1414/attachments/1013/1143/circe2-handout.pdf
  //
  // Secondary purpose:
  // To get a paramterization that can allow for generating sets of sqrt(s) with only the beamstrahlung effects. The hope being that this can be included in a larger model that includes beam energy spread and initial state radiation. 
  // 
  //This specific file is used to correct the error of the double beta distribution fit that is caused by some issue in ROOT.
  //
  // To dump the text output to a file type ".> mylog.txt" in the prompt before the script
  gSystem->Load("libRooFit");
  using namespace RooFit; 

  //The variable below is for the x axis. Here we use the fractional sqrt(s) of the event.
  RooRealVar bsx("bsx","Fractional #sqrt{s} [1 - #sqrt{s}/(250 GeV)]^{#eta}",0.05,0.75,"GeV/250");

  //Keep in mind when declaring variable's names (the first string argument):
  //1.)A hash (#) will be put in front of it. So if you want your variable to be a greek beta symbol just call your variable "beta". However, if you want an uppercase B you will need to put "Beta" or use \thinspace (which is essentially nothing). So B would be "thinspace B" and b would be "thinspace b".
  //2.)For any superscript or subscript write as normal but do not use the squigly brackets, {}, as they will be added automatically by the plot macro. Additionally, only single digit numbers are supported for this operation. Furthermore, you should not have numbers in your variable name unless they are single digits AND inside the super or sub script(s).
  //i.e. Do not use "be5k" as it will be converted to "be{5}k" which will look odd on the plot. Furthermore, "be^45k" will be rewritten as "be^{4}{5}k" which will display more like "be^4 5k". 

  //Normalization/coefficient values. 1 and 2 are for addition. 3 is for the product sums if you want to do that.
  RooRealVar A1("Alpha_1","A_{1}",1.0,1.0e-9,10000);
  RooRealVar A2("Alpha_2","A_{2}",0.01,1.0e-9,10000);
  RooRealVar A3("Alpha_3","A_{3}",0.05,1.0e-9,10000);

  //Beta distribution fitting parameters
  RooRealVar alpha1("alpha_1","#alpha_{1}",1.24,0.5,2);
  RooRealVar gam1("gamma_1","#gamma_{1}",40,8,50);
  RooRealVar alpha2("alpha_2","#alpha_{2}",0.47,0.2,2);
  RooRealVar gam2("gamma_2","#gamma_{2}",40,0,50);

  //1/x^3 background fitting parameters
  RooRealVar Mu1("Mu_1","#Mu_{1}",1,1e-6,250);
  RooRealVar Kappa("Kappa","#Kappa",0.5,-10,10);

  //Inverse gaussian fitting parameters
  RooRealVar lam("lambda_1","#lambda_{1}",1,0,20);
  RooRealVar mu("mu_1","#mu_{1}",5,0,50);
  RooRealVar lam2("lambda_2","#lambda_{2}",2,0,20);
  RooRealVar mu2("mu_2","#mu_{2}",2,0,50);

  //Power law from the data plotting. Should be the same as what is used in the DrawOption string below. You can set it constant but its variable here.
  RooRealVar eta("eta","eta",6,3,7);
  eta.setConstant(kTRUE);
  eta.removeError();
  //(1-x)^(1/eta) looking distributions to fit with. Multiple functions so that sum can be used.
  //Inverse Gaussian
  RooGenericPdf invg1("invg1","invg1","sqrt(lambda_1/(2*3.14159*(bsx)^3)) * exp(-0.5*((lambda_1*(bsx-mu_1)^2)/(mu_1*mu_1*bsx)))",RooArgSet(lam,bsx,mu));
  RooGenericPdf invg2("invg2","invg2","sqrt(lambda_2/(2*3.14159*(bsx)^3)) * exp(-0.5*((lambda_2*(bsx-mu_2)^2)/(mu_2*mu_2*bsx)))",RooArgSet(lam2,bsx,mu2));
  //Beta distributions
  RooGenericPdf beta1("beta1","beta1","(((bsx^eta)^(alpha_1-1)) * ((1-bsx^eta)^(gamma_1-1)))",RooArgSet(bsx,eta,alpha1,gam1));
  RooGenericPdf beta2("beta2","beta2","(((bsx^eta)^(alpha_2-1)) * ((1-bsx^eta)^(gamma_2-1)))",RooArgSet(bsx,eta,alpha2,gam2));
  //Asymptotic background for the beta + linear fit that turns on after 0.5.
  //Its called lin because it started as a linear background and then lots of madness turned it into this
  RooGenericPdf lin("lin","lin","Alpha_2 / (bsx^3) * 0.5*(1+tanh(10*(bsx-0.6)))",RooArgSet(A2,bsx));

  //Some function combinations to fit
  RooGenericPdf fun1("fun1","fun1","Alpha_1 * beta1 + lin",RooArgSet(A1,beta1,lin));
  RooGenericPdf fun2("fun2","fun2","Alpha_1 * beta1",RooArgSet(A1,beta1));
  RooGenericPdf fun3("fun3","fun3","Alpha_1 * invg1 + Alpha_2 * invg2",RooArgSet(A1,invg1,A2,invg2));
  RooGenericPdf fun4("fun4","fun4","Alpha_1 * invg1",RooArgSet(A1,invg1));
  RooGenericPdf fun5("fun5","fun5","Alpha_1 * beta1 + Alpha_2 * beta2",RooArgSet(A1,beta1,A2,beta2));

  // Access the histogram and tree stored in the root file
  TFile hfile("GP/gp_ILC250_V17.root");
  TTree* tree=(TTree*)hfile.Get("GPtree");

  //Draw the histogram that is to be fitted for. You can change the DrawOptions and DrawCuts strings to satisfy this.
  TString DrawOption = "(1-sqrt(E1*E2)/125)^(1/5)>>htemp(100,0.05,0.75)";
  TCut DrawCuts = "2*sqrt(E1*E2) < 250";
  
  //Draw without cuts to get the total events
  tree->Draw(DrawOption,"","goff");
  TH1F *h1 = (TH1F*) gDirectory->Get("htemp");
  //We need this later for the stats display
  Int_t TotEvents = h1->GetEntries();
  h1->Scale(1.0/(h1->Integral()));
  //Repeat with our cuts
  tree->Draw(DrawOption,DrawCuts,"goff");
  TH1F *h2 = (TH1F*) gDirectory->Get("htemp");
  //We need this later for the stats display
  Int_t UncutEvents = h2->GetEntries();
  h2->Scale(1.0/(h2->Integral()));
  //h1->Draw("E");

  // Create a binned dataset that imports contents of h1 and associates its contents to observable 'sqrts' aka the beam energy
  RooDataHist db("db", "db", bsx,Import(*h2));

  //create a RooAddPdf for passing the model #1
  RooRealVar addfrac1("% Signal","addfrac1",1.0,0.0,1.0);
  addfrac1.setConstant(kTRUE);
  addfrac1.removeError();
  RooAddPdf sigAdd1("sigAdd1","GP Beams. Beta Inv",RooArgList(fun1),addfrac1);

  //Temporary file to pass the PDF
  TFile f("temp.root","recreate");

  db.Write();
  sigAdd1.Write();
  bsx.Write();

  //Check that the PDF and data were saved.
  if(gDirectory->Get("db") == NULL) {
    cout << "Unable to return database." << endl;
  }

  //Create the title string and then call the fitplot macro
  std::string TitleArg = "(1- #sqrt{s} )^{#eta}|#splitline{#Alpha_{1} Beta(x^{#eta},#alpha_{1} , #gamma_{2})}{ + #Alpha_{2} #Theta(x-0.5) / x^{3}}";
  FitPlot(db.GetName(),sigAdd1.GetName(),bsx.GetName(),TitleArg.c_str(),TotEvents,UncutEvents);

  //Cloce the file and remove the temp file
  f.Close();
  gSystem->Exec("rm temp.root");

  //create a RooAddPdf for passing the model #2
  RooRealVar addfrac2("% Signal","addfrac2",1.0,0.0,1.0);
  addfrac2.setConstant(kTRUE);
  addfrac2.removeError();
  RooAddPdf sigAdd2("sigAdd2","GP Beams. Circe",RooArgList(fun2),addfrac2);

  //Temporary file to pass the PDF
  TFile f2("temp.root","recreate");

  db.Write();
  sigAdd2.Write();
  bsx.Write();

  //Check that the PDF and data were saved.
  if(gDirectory->Get("db") == NULL) {
    cout << "Unable to return database." << endl;
  }

  TitleArg = "(1- #sqrt{s} )^{#eta}|#Alpha_{1} Beta(x^{#eta},#alpha_{1} , #gamma_{1})";
  FitPlot(db.GetName(),sigAdd2.GetName(),bsx.GetName(),TitleArg.c_str(),TotEvents,UncutEvents);

  //Cloce the file and remove the temp file
  f2.Close();
  gSystem->Exec("rm temp.root");

  //create a RooAddPdf for passing the model #3
  RooRealVar addfrac3("% Signal","addfrac3",1.0,0.0,1.0);
  addfrac3.setConstant(kTRUE);
  addfrac3.removeError();
  RooAddPdf sigAdd3("sigAdd3","GP Beams. 2 Inv. Gauss",RooArgList(fun3),addfrac3);

  //Temporary file to pass the PDF
  TFile f3("temp.root","recreate");

  db.Write();
  sigAdd3.Write();
  bsx.Write();

  //Check that the PDF and data were saved.
  if(gDirectory->Get("db") == NULL) {
    cout << "Unable to return database." << endl;
  }

  TitleArg = "(1- #sqrt{s} )^{#eta}|#splitline{#Alpha_{1} IG(x^{#eta},#mu_{1} , #lambda_{1})}{ + #Alpha_{2} IG(x^{#eta},#mu_{2},#lambda_{2})}";
  FitPlot(db.GetName(),sigAdd3.GetName(),bsx.GetName(),TitleArg.c_str(),TotEvents,UncutEvents);

  //Cloce the file and remove the temp file
  f3.Close();
  gSystem->Exec("rm temp.root");

  //create a RooAddPdf for passing the model #4
  RooRealVar addfrac4("% Signal","addfrac4",1.0,0.0,1.0);
  addfrac4.setConstant(kTRUE);
  addfrac4.removeError();
  RooAddPdf sigAdd4("sigAdd4","GP Beams. 1 Inv. Gauss",RooArgList(fun4),addfrac4);

  //Temporary file to pass the PDF
  TFile f4("temp.root","recreate");

  db.Write();
  sigAdd4.Write();
  bsx.Write();

  //Check that the PDF and data were saved.
  if(gDirectory->Get("db") == NULL) {
    cout << "Unable to return database." << endl;
  }

  TitleArg = "(1- #sqrt{s} )^{#eta}|#Alpha_{1} IG(x^{#eta},#mu_{1} , #lambda_{1})";
  FitPlot(db.GetName(),sigAdd4.GetName(),bsx.GetName(),TitleArg.c_str(),TotEvents,UncutEvents);

  //Cloce the file and remove the temp file
  f4.Close();
  gSystem->Exec("rm temp.root");

  //Wow, you wrote this five times now, maybe you should make a function?
  //Well, you know that copy and paste and replace are all functions, right? And they are pretty fast... hehe

  //create a RooAddPdf for passing the model #5
  RooRealVar addfrac5("% Signal","addfrac5",1.0,0.0,1.0);
  addfrac5.setConstant(kTRUE);
  addfrac5.removeError();
  RooAddPdf sigAdd5("sigAdd5","GP Beams. 2 Beta",RooArgList(fun5),addfrac5);

  //Temporary file to pass the PDF
  TFile f5("temp.root","recreate");

  db.Write();
  sigAdd5.Write();
  bsx.Write();

  //Check that the PDF and data were saved.
  if(gDirectory->Get("db") == NULL) {
    cout << "Unable to return database." << endl;
  }

  TitleArg = "(1- #sqrt{s} )^{#eta}|#splitline{#Alpha_{1} Beta(x^{#eta},#alpha_{1} , #gamma_{1})}{ + #Alpha_{2} Beta(x^{#eta},#alpha_{2},#gamma_{2})}";
  FitPlot(db.GetName(),sigAdd5.GetName(),bsx.GetName(),TitleArg.c_str(),TotEvents,UncutEvents);

  //Cloce the file and remove the temp file
  f5.Close();
  gSystem->Exec("rm temp.root");

  hfile.Close();

}
