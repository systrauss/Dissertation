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
#include "TSystem.h"
#include "TROOT.h"
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
#include "TMath.h"

//Object files
#include "analysis.h"
#include "Constraints.h"
#include "histograms.h"
#include "Coefficients.h"
#include "Filelist.h"
#include "timing.h"

//Data tree files
// #include "eventData.h"
// #include "eventScaler.h"
// #include "ExpEventDictionary.h"
// #include "ExpEventLinkDef.h"

using namespace std;

//File path info
extern std::string sFilepath;
extern std::string sTree;
extern std::string sRun;
extern std::string sFType;

extern TChain* chain;
extern std::vector<std::vector<TH1F*> > ge_en_ge_cut;
extern std::vector<std::vector<TH1F*> > sili_en_ge_cut;
extern std::vector<std::vector<TH1F*> > ge_en_sili_cut;
extern std::vector<std::vector<TH1F*> > sili_en_sili_cut;


extern int nGeOrder; //Order of calibration i.e. 1 = linear.
extern int nGeDets; //Total number of signals from Germanium detectors
extern int nGeSegments; //number of segments in a single Germanium crystal, for adding purposes
extern std::vector<std::vector<double> > dGeCoefficients; //Ge Coefficients
extern std::vector<std::vector<double> > dGeCoeffRes; // Ge Residual Coefficients
extern std::vector<std::vector<double> > dGeRunCorr; //Run Correction Coefficients

extern int nSiLiOrder; //Order of calibration i.e. 1 = linear.
extern int nSiLiDets; //Total number of signals from SiLi
extern int nSiLiPlace; //Start of SiLis in generalized array detectors
extern std::vector<std::vector<double> > dSiLiCoefficients; //Coefficients
extern std::vector<std::vector<double> > dSiLiRunCorr; //Run Correction Coefficients

std::vector<std::vector<double> > dGeBounds; //bounds for cuts
std::vector<std::vector<double> > dSiLiBounds; //bounds for cuts
int nGeConstraints;
int nSiLiConstraints;

Timing tGates;

int main(int argc, char* argv[]) //Order of arguments: run #, output filename, cut filename, time low, time high
{
	char buffer[50];
	int nRunNum;
    nRunNum = atoi(argv[1]); //Run to do the cuts on
	char* sOut = argv[2]; //file title to write to
	char* sCut1 = argv[3]; //Cut file name indicator
	char* sTime = argv[4]; //Cut file name indicator
	gSystem->Load("$Directory$/libExpEvent.so","", kTRUE);
	readPaths(); //From Filelist.cxx
	makeChain(nRunNum); //From Filelist.cxx
	defineGeCoeff(); //From Coefficients.cxx
	defineGeCoeff(nRunNum); //From Coefficients.cxx, correction terms
	defineSiLiCoeff(); //From Coefficients.cxx
	defineSiLiCoeff(nRunNum);  //From Coefficients.cxx
	defineBGO(); //From constraints.cxx
	tGates = Timing(nGeDets/nGeSegments, nSiLiDets, sTime);
	//Debugging for the Timing class
	std::cout << "Germanium pairs: " <<  tGates.geGe.size() << std::endl;
	//for (int i = 0; i < tGates.geGe.size(); i++) std::cout << tGates.geGe[i][0] <<"\t"  << tGates.geGe[i][1] <<std::endl; 
	std::cout << "Germanium-SiLi pairs: " << tGates.geSiLi.size() << std::endl;
	//for (int i = 0; i < tGates.geSiLi.size(); i++) std::cout << tGates.geSiLi[i][0] <<"\t"  << tGates.geSiLi[i][1] <<std::endl; 
	std::cout << "SiLi pairs: " << tGates.siliSiLi.size() << std::endl;
	//for (int i = 0; i < tGates.siliSiLi.size(); i++) std::cout << tGates.siliSiLi[i][0] <<"\t"  << tGates.siliSiLi[i][1] <<std::endl;
	analysis ana(chain); //analysis class. Main part of code.
	sprintf(buffer,"GeCut_%s.dat",sCut1); //File name to input
	nGeConstraints = defineConstraints(buffer,dGeBounds); //From constraints.cxx
	sprintf(buffer,"SiLiCut_%s.dat",sCut1); //File name to input
	nSiLiConstraints = defineConstraints(buffer,dSiLiBounds);
	makeHistograms(nGeDets/nGeSegments,nGeConstraints,nSiLiDets,nSiLiConstraints); //from histograms.cxx
	ana.Loop(Form("$OutputDirectory$/%s_run_00%i.root",sOut,nRunNum),nRunNum, false); //fOut is in Filelist.h
}