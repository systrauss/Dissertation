/*
Please check the README for more information about this code and it's purpose.
*/

//C++ libraries
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

//ROOT libraries
#include <TROOT.h>
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TObject.h"
#include "TObjArray.h"

//Object files
#include "Coefficients.h"
#include "Constraints.h"
#include "histograms.h"
#include "Filelist.h"

using namespace std;

//File path info
extern std::string sFilepath;
extern std::string sTree;
extern std::string sRun;
extern std::string sFType;

//from histograms.cxx
extern std::vector<std::vector<TH1F*> > ge_en_ge_cut;
extern std::vector<std::vector<TH1F*> > sili_en_ge_cut;
extern std::vector<std::vector<TH1F*> > ge_en_sili_cut;
extern std::vector<std::vector<TH1F*> > sili_en_sili_cut;

//From coefficients.cxx
extern int nGeOrder; //Order of calibration i.e. 1 = linear.
extern int nGeDets; //Total number of signals from Germanium detectors
extern int nGeSegments; //number of segments in a single Germanium crystal, for adding purposes

extern int nSiLiOrder; //Order of calibration i.e. 1 = linear.
extern int nSiLiDets; //Total number of signals from SiLi
extern int nSiLiPlace; //Start of SiLis in generalized array detectors

std::vector<std::vector<double> > dGeBounds; //bounds for cuts
std::vector<std::vector<double> > dSiLiBounds; //bounds for cuts
int nGeConstraints;
int nSiLiConstraints;

std::vector<std::vector<double> > dGeGeBounds; //bounds for cuts
std::vector<std::vector<double> > dGeSiLiBounds; //bounds for cuts
std::vector<std::vector<double> > dSiLiSiLiBounds; //bounds for cuts
int nGeGeConstraints;
int nGeSiLiConstraints;
int nSiLiSiLiConstraints;

int main(int argc, char* argv[]) //Order of arguments: first run, last run, output file name, cut file name
{
	char buffer[50];
	int nRunStart = atoi(argv[1]);
	int nRunEnd = atoi(argv[2]);
	char* sOut = argv[3];
	char* sCut1 = argv[4];
	char* sCut2 = argv[5];
	readPaths();
	defineGeCoeff(); //From Coefficients.cxx
	defineSiLiCoeff(); //From Coefficients.cxx
	if(strcmp(sCut2,"0")==0)
	{
		sprintf(buffer,"GeCut_%s.dat",sCut1); //File name to input
		nGeConstraints = defineConstraints(buffer,dGeBounds); //From constraints.cxx
		sprintf(buffer,"SiLiCut_%s.dat",sCut1); //File name to input
		nSiLiConstraints = defineConstraints(buffer,dSiLiBounds);
		std::cout << "Regular Coincidence" << std::endl;
		writeHistToFile(Form("/scratch365/sstrauss/%s.root",sOut),sOut,nRunStart,nRunEnd, false); //fOut is in Filelist.h
	}
	else
	{
		sprintf(buffer,"Cut_Ge_%s_Ge_%s.dat",sCut1,sCut2); //File name to input
		nGeGeConstraints = defineTripleConstraints(buffer,dGeGeBounds); //From constraints.cxx
		sprintf(buffer,"Cut_Ge_%s_SiLi_%s.dat",sCut1,sCut2); //File name to input
		nGeSiLiConstraints = defineTripleConstraints(buffer,dGeSiLiBounds); //From constraints.cxx
		sprintf(buffer,"Cut_SiLi_%s_SiLi_%s.dat",sCut1,sCut2); //File name to input
		nSiLiSiLiConstraints = defineTripleConstraints(buffer,dSiLiSiLiBounds); //From constraints.cxx
		std::cout << "Triple Coincidence" << std::endl;
		writeHistToFile(Form("/scratch365/sstrauss/%s.root",sOut),sOut,nRunStart,nRunEnd, true); //fOut is in Filelist.h
	}
	return 0;
}