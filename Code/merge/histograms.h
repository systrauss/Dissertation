#ifndef HISTOGRAMS_H
#define HISTOGRAMS_H

#ifndef _HISTOGRAMS_CXX
#define EXTERNAL extern
#else
#define EXTERNAL
#endif

//C++ libraries
#include <vector>

//ROOT libraries
#include <TH1F.h>
#include <TString.h>
#include <TFile.h>

void writeHistToFile(char* fileName, char* sOut, int nRunStart, int nRunStop);

#endif