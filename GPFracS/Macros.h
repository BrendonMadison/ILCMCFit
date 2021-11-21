#ifndef MACROS_H
#define MACROS_H

//These are needed to use the string and open+write the csv file and the temp files for FitPlot
#include <fstream>
#include <string>

void Hist2CSV(TH1* h, std::string FileName);

void FitPlot(const string argdb,const string argmodel, const string argx, const string argtit, const Int_t argev, const Int_t argcut);

#endif