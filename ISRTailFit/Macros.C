  //Contains the macros used. Namely FitPlot, which fits with RooFit and plots, and Hist2CSV, which converts a histogram's content to a csv text file.

#include "Macros.h"

void Hist2CSV(TH1* h, std::string FileName){
  //
  // Primary purpose: Takes a ILC Soft ROOT file and writes the contents of a given histogram(s) to csv.
  //                : Mainly for the purpose of quickly opening said data in python to run the Monte-Carlo fitting 



  //Takes a th1 (histogram) and then prints the bin centers, bin content and bin error
  //In general, there should be a file extension included on FileName. 
  //Technically it doesn't matter as it will be written in comma separated form regardless but be consistent and use .csv

   Int_t n = h->GetNbinsX();
   std::ofstream fo(FileName.c_str());
   Double_t hcon = 0.0;
   
   for (Int_t i=1; i<=n; i++) 
   {
      hcon = h->GetBinContent(i);
      fo << Form("%f,%f,%f",h->GetBinLowEdge(i)+h->GetBinWidth(i)/2, hcon,sqrt(hcon)) << endl;
   }
   fo.close();
   cout << "Dumped histogram contents to " << FileName << endl;
}

void FitPlot(const string argdb,const string argmodel, const string argx, const string argtit, const Int_t argev, const Int_t argcut){
  //Called in fitting functions. Takes the data and model and then fits with RooFit and plots it.
  //This streamlines and standardizes the plot process so that you don't have to write this each time you use RooFit...

  gSystem->Load("libRooFit");
  gSystem->Load("RooCrystalBall.cxx");
  using namespace RooFit;

  cout << "Using the database " << argdb << " and model " << argmodel <<  endl;

  cout << "Number events: " << argev << "\nNumber events after cuts: " << argcut << endl;

  gDirectory->ls();
  RooDataHist *ldb = (RooDataHist*)gDirectory->Get(argdb.c_str());
  
  //gDirectory->GetObject(argdb.c_str(),ldb);
  //  auto ldb = gDirectory->Get(argdb.c_str());
  RooAddPdf *lmodel = (RooAddPdf*)gDirectory->Get(argmodel.c_str());
  //gDirectory->GetObject(argmodel.c_str(),lmodel);
  RooRealVar *lx = (RooRealVar*)gDirectory->Get(argx.c_str());
  // gDirectory->GetObject(argx.c_str(),lx);
  //RooAbsPdf  lmodel = gDirectory->Get(argmodel.c_str());
  //auto lx = gDirectory->Get(argx.c_str());
  //Fit the data (db) to the model (crysball)


  //char *cmd;
  //cmd = Form("%s lmodel = gDirectory->Get(\"%s\");",argtype.c_str(),argmodel.c_str());
  //gROOT->ProcessLine(cmd);

  //char *cmd2;
  //cmd2 = Form("lmodel->fitTo(*ldb);");
  //gROOT->ProcessLine(cmd2);

  //check if the pointers are null
  if(ldb == NULL){
    cout << "Database returned null. Check that you are passing correctly, saving locally." << endl;
    return;
  }

  if(lmodel == NULL){
    cout << "Model returned null. Check that you are passing correctly, saving locally." << endl;
    return;
  }

  if(lx == NULL){
    cout << "X axis returned null. Check that you are passing correctly, saving locally." << endl;
    return;
  }

  lmodel->chi2FitTo(*ldb);

  // Plot data and results
  std::string argtitle = argtit.substr(0, argtit.find("|"));
  std::string argeqn = argtit.substr(argtit.find("|")+1);
  RooPlot *xframe=lx->frame(Title(Form("%s fit of \%s",lmodel->GetTitle(),argtitle.c_str())));
  ldb->plotOn(xframe);

  //char *cmd3;
  //cmd3 = Form("lmodel->plotOn(xframe,LineColor(kBlue));");
  //gROOT->ProcessLine(cmd3);

  lmodel->plotOn(xframe,LineColor(kBlue));
  
  //  lmodel->paramOn(xframe, Layout(0.25,0.25,0.75));
  //ldb->statOn(xframe, Layout(0.25,0.25,0.30));
  //  xframe->SetTitle("Title Argument");
  //  xframe->Draw();

  //Construct a histogram with pulls
  RooHist* hpull = xframe->pullHist();
  Double_t pmean = hpull->GetMean(2);
  Double_t prms = hpull->GetRMS(2);
  RooPlot* pframe = lx->frame(Title(Form("Pull of %s fit",argtitle.c_str())));
  pframe->addPlotable(hpull,"P");

  Int_t pt = 15;

  // Check goodness of fit ...
  //For some reason the below returns ChiSquare per bin...
  Double_t pChi = xframe->chiSquare();
  auto hdat = ldb->createHistogram(argx.c_str());
  RooArgSet asx(*lx);
  auto fmod = lmodel->asTF(asx,RooArgSet(),asx);
  fmod->SetNpx(hdat->GetNbinsX());

  //KS test statistic
  Double_t ksTest = hdat->KolmogorovTest(fmod->GetHistogram());
  
  //Double_t pChi = fmod->GetChisquare();

  Double_t nPars = lmodel->getParameters(*ldb)->getSize();
  Int_t ndof = hdat->GetNbinsX() - nPars;

  pChi = pChi * hdat->GetNbinsX();

  cout << "Reduced ChiSq: " << pChi << " / " << ndof << endl;
  cout << "Kolmogorov-Smirnov test result:" << ksTest << endl; 

  //Get the parameter set
  RooArgSet *params = lmodel->getParameters(*lx);

  //Write it to a text file
  params->writeToFile("Params.txt");

  //reads the parameter text file and puts it into a string so it can be used in the text box below
  std::ifstream pfile("Params.txt");
  std::string linestr;
  TString namestr;
  std::string meanstr;
  std::string errstr;
  std::string lineholder;
  Int_t strcounter = 0;
  std::vector<std::string> truncated;
  while (std::getline(pfile, linestr))
    {
      if(strcounter != 0)
	{
	  //We want to format the line string into scientific notation otherwise the text may go off the plot screen. Which looks bad.
          lineholder = linestr.substr(0, linestr.find("L("));
	  //Grabs the name of the parameter and then loops such that superscripts and subscripts will properly display
	  //Normally its impossible to pass a squigly bracket {} as it doesn't work with RooFit's parser
	  namestr = lineholder.substr(0,lineholder.find("=")+1);
	  //namestr.ReplaceAll("_","_{");
	  //namestr.ReplaceAll("^","^{");
	  for (int i = 0; i < 10; i++) {
	    namestr.ReplaceAll(Form("%i",i),Form("{%i}",i));
	  }
          meanstr = lineholder.substr(lineholder.find("=")+1,lineholder.find("+/-")-1);
          errstr = lineholder.substr(lineholder.find("-")+1);
	  //So we cut out the values we want and rewrite in scientific notation
          truncated.push_back(Form("#%s %.2e +- %.2e",std::string(namestr).c_str(),std::atof(meanstr.c_str()),std::atof(errstr.c_str())));
	}
      strcounter += 1;
    }

  gSystem->Exec("rm Params.txt");

  //lmodel->paramOn(xframe, Label(Form("#chi ^{2} / nDoF = %.2e / %i",pChi,ndof)),Layout(0.15,0.35,0.8));

  TCanvas* c = new TCanvas("c",argtit.c_str(),1200,1000);
  c->Divide(1,2);

  gStyle->SetTickLength(0.01,"y"); gStyle->SetTickLength(0.02,"x"); gStyle->SetLabelSize(0.05,"x y"); gStyle->SetLabelOffset(0.02,"x"); gStyle->SetTitleFontSize(0.1); gStyle->SetTitleOffset(0.9); gStyle->SetTitleXOffset(1); gStyle->SetTitleYOffset(1); gStyle->SetTitleXSize(0.075); gStyle->SetTitleYSize(0.05);

  c->cd(1); gPad->SetBottomMargin(0.2); gPad->SetTopMargin(0.2); gPad->SetRightMargin(0.2); xframe->GetXaxis()->SetTitle(lx->GetTitle()); xframe->GetXaxis()->SetRangeUser(xframe->GetXaxis()->GetXmin()*0.5,xframe->GetXaxis()->GetXmax()); xframe->GetYaxis()->SetTitle("Normalized Units"); xframe->Draw();

  //Add latex text box for the equation
  TLatex peq(0.82,0.75,Form("--Fit Equation--"));
  peq.SetNDC(kTRUE);
  peq.SetTextSize(0.05);
  peq.Draw();

  cout << argeqn << endl;
  TLatex ppdf;
  ppdf.SetNDC(kTRUE);
  ppdf.SetTextSize(0.035);
  ppdf.DrawLatex(0.82,0.65,TString(argeqn.c_str()));

  //Add latex text box for the parameters
  TLatex ptt(0.82,0.55,Form("--Fit Parameters--"));
  ptt.SetNDC(kTRUE);
  ptt.SetTextSize(0.05);
  ptt.Draw();

  //loop through the parameter vector
  TLatex pmt;
  pmt.SetTextSize(0.035);
  pmt.SetNDC(kTRUE);
  strcounter = 0;
  while (strcounter < truncated.size())
    {
      pmt.DrawLatex(0.82,0.50-strcounter*0.05,TString(truncated.at(strcounter)));
      cout << truncated.at(strcounter) << endl;
      strcounter += 1;
    }

  c->cd(2); gPad->SetBottomMargin(0.2); gPad->SetTopMargin(0.2); gPad->SetRightMargin(0.2); pframe->GetXaxis()->SetTitle(lx->GetTitle()); pframe->Draw();

  //Add latex text boxes
  TLatex tt(0.82,0.75,Form("--Fit Stats--"));
  //TLatex dmt(0.82,1.5,Form("Data Mean\t: %.2e",hdat->GetMean(2)));
  TLatex mt(0.82,0.7,Form("Pull Mean  : %.2e",pmean));
  TLatex rt(0.82,0.65,Form("Pull #sigma^{2}        : %.2e",prms));
  TLatex kt(0.82,0.6,Form("KS Test      : %.2e",ksTest));
  TLatex ct(0.82,0.55,Form("Red. #chi ^{2}      : %.2e",pChi/ndof)); 
  TLatex cnt(0.82,0.5,Form("#chi ^{2} / nDoF  : %.2e / %i ",pChi,ndof));
  tt.SetNDC(kTRUE);
  //dmt.SetNDC(kTRUE);
  mt.SetNDC(kTRUE);
  rt.SetNDC(kTRUE);
  kt.SetNDC(kTRUE);
  ct.SetNDC(kTRUE);
  cnt.SetNDC(kTRUE);
  tt.SetTextSize(0.05);
  mt.SetTextSize(0.035);
  rt.SetTextSize(0.035);
  kt.SetTextSize(0.035);
  ct.SetTextSize(0.035);
  cnt.SetTextSize(0.035);
  tt.Draw();
  //dmt.Draw();
  mt.Draw();
  rt.Draw();
  kt.Draw();
  ct.Draw();
  cnt.Draw();
  cnt.DrawLatex(0.82,0.45,Form("Events         : %i",argev));
  cnt.DrawLatex(0.82,0.4,Form("Post-cut Events       : %i",argcut));
  //TText *mtxt = new TText(hpull->GetXaxis()->GetXmax() * 0.965,pmean+3.0*prms,Form("Mean: %.2f",pmean));
  //TText *rtxt = new TText(hpull->GetXaxis()->GetXmax() * 0.965,pmean+1.0*prms,Form("RMS : %.2f",prms));
  //TText *ktxt = new TText(hpull->GetXaxis()->GetXmax() * 0.965,pmean-1.0*prms,Form("KS Test: %.2f",ksTest));
  //TText *ctxt = new TText(hpull->GetXaxis()->GetXmax() * 0.965,pmean-3.0*prms,Form("Red. Chi2: %.2f",pChi/ndof));
  //pframe->addObject(mtxt);
  //pframe->addObject(rtxt);
  //pframe->addObject(ktxt);
  //pframe->addObject(ctxt);

  TString filename = Form("%s_%s.pdf",argtitle.c_str(),lmodel->GetTitle());

  filename.ReplaceAll(" ","_");
  filename.ReplaceAll("#","");

  c->SaveAs(filename);

  delete c;

}
