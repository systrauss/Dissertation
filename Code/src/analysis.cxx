#define analysis_cxx
#include "analysis.h"

//C++ libraries
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

//ROOT libraries
#include <TSystem.h>
#include <TROOT.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TMath.h>

//object files
#include "histograms.h"
#include "Constraints.h"

//Data tree files
// #include "eventData.h"
// #include "eventScaler.h"
// #include "ExpEventDictionary.h"
// #include "ExpEventLinkDef.h"

//From Coefficients.cxx
extern int nGeOrder; //Order of calibration i.e. 1 = linear.
extern int nGeDets; //Total number of signals from Germanium detectors
extern int nGeSegments; //number of segments in a single Germanium crystal, for adding purposes
extern std::vector<std::vector<double> > dGeCoefficients; //Ge Coefficients
extern std::vector<std::vector<double> > dGeCoeffRes; // Ge Residual Coefficients
extern std::vector<std::vector<double> > dGeRunCorr; //Run Correction Coefficients

extern TChain* chain;

extern int nSiLiOrder; //Order of calibration i.e. 1 = linear.
extern int nSiLiDets; //Total number of signals from SiLi
extern int nSiLiPlace; //Start of SiLis in generalized array detectors
extern std::vector<std::vector<double> > dSiLiCoefficients; //Coefficients
extern std::vector<std::vector<double> > dSiLiRunCorr; //Run Correction Coefficients

//From Constraints.cxx
extern std::vector<std::vector<double> > dGeBounds; //bounds for cuts
extern std::vector<std::vector<double> > dSiLiBounds; //bounds for cuts
extern int nGeConstraints;
extern int nSiLiConstraints;

extern int nBGODets; //Total number of signals from BGO detectors
extern int nBGOPlace; //Start of BGOs in generalized array
extern std::vector<double> dBGOThreshold;

//From histograms.cxx
extern std::vector<std::vector<TH1F*> > ge_en_ge_cut;
extern std::vector<std::vector<TH1F*> > sili_en_ge_cut;
extern std::vector<std::vector<TH1F*> > ge_en_sili_cut;
extern std::vector<std::vector<TH1F*> > sili_en_sili_cut;

using namespace std;

void analysis::Loop(const char* fileName, int nRunNum, bool bTripleCoin)
{
//   In a ROOT session, you can do:
//      Root > .L analysis.C
//      Root > analysis t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   //Don't keep writing to fOut, so switch to the "user" space instead
   gROOT->cd();


   TRandom3 *randgen = new TRandom3(1); //For uniformity over bin, to remove artifacting
   std::vector<double> dGeEn;
   std::vector<double> dSiLiEn;
   std::vector<double> dBGO;

   //Times
   std::vector<double> dT_GeEnSeg;
   std::vector<double> dT_GeEn;
   std::vector<double> dT_SiLiEn;
   std::vector<double> dT_BGO;

   //Ge vector
   for(int i=0; i<nGeDets/nGeSegments;i++)
   {
       dGeEn.push_back(0);
       dT_GeEn.push_back(0);
   }
   for(int i=0; i<nGeDets;i++)
   {
       dT_GeEnSeg.push_back(0);
   }
   for(int i=0; i<nSiLiDets;i++)
   {
       dSiLiEn.push_back(0);
       dT_SiLiEn.push_back(0);
   }
   for(int i=0; i<nBGODets;i++)
   {
       dBGO.push_back(0);
       dT_BGO.push_back(0);
   }
   
   //main loop through the data
   Long64_t nentries = fChain->GetEntriesFast();
   
   //main loop through the data
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
    	if (jentry % 1000000 == 0) {
        	cout << "." << flush;
    	}
      // if (Cut(ientry) < 0) continue;
   //Germanium calibration
      for(int i=0 ; i < nGeDets ; i++)
      {
         if(i % nGeSegments == 0)
         {
            dGeEn[i/nGeSegments] = 0; //Reset the total clover energy
            dT_GeEnSeg[i] = 0; //Reset segment time
            dT_GeEn[i/nGeSegments] = 0; //Reset total time
         } 
         double eneC = 0, dADC = 0;
         dADC = ge_en_raw[i] + 0.5*(1 - 2*randgen->Rndm());
         if(ge_time[i] == 0 || buncher_time == 0)
         {
             dT_GeEnSeg[i] = -10000;
         }
         else
         {
             dT_GeEnSeg[i] = ge_time[i]-buncher_time;
         }
         for(int j=0 ; j<nGeOrder+1 ; j++)
         {
            eneC = eneC + dGeCoefficients[i][j]*pow(dADC,j);
         }
         eneC = dGeRunCorr[i][0]+dGeRunCorr[i][1]*eneC; //Run correction factor
         dGeEn[i] = eneC;
         dT_GeEn[i] = dT_GeEnSeg[i]; //This is the line that should fix the previous problem!
        //  if(i%nGeSegments!=0 && i!=14 && i!=15) 
        //  {
        //     dGeEn[i/nGeSegments] = dGeEn[i/nGeSegments]+eneC; //Add energies
        //     if(dT_GeEn[i/nGeSegments] == 0) dT_GeEn[i/nGeSegments] = dT_GeEnSeg[i]; //Timing
        //  }
        //  if(dADC > 0.5) std::cout << dADC << std::endl;
      }
      //SiLi Calibration
      for(int i=0 ; i < nSiLiDets ; i++)
      {
         dSiLiEn[i]=0;
         dT_SiLiEn[i] = 0;
         double dADC = 0;
         dADC = sili_en_raw[i] + 0.5*(1 - 2*randgen->Rndm());
         if(sili_time[i] == 0 || buncher_time == 0)
         {
             dT_SiLiEn[i] = -10000;
         }
         else
         {
             dT_SiLiEn[i] = sili_time[i]-buncher_time;
         }
         for(int j=0 ; j<nSiLiOrder+1 ; j++)
         {
            dSiLiEn[i] = dSiLiEn[i] + dSiLiCoefficients[i][j]*pow(dADC,j);
         }
         dSiLiEn[i] = dSiLiRunCorr[i][0]+dSiLiRunCorr[i][1]*dSiLiEn[i]; //Run-by-run correction
      }
      //BGOs, done dumb because of T. Battaglia not using arrays?!?
      dBGO[0] = bgo1_en;
      dBGO[1] = bgo4_en;
      dBGO[2] = bgo5_en;
      dBGO[3] = bgo6_en;
      dT_BGO[0] = bgo1_time-buncher_time;
      dT_BGO[1] = bgo4_time-buncher_time;
      dT_BGO[2] = bgo5_time-buncher_time;
      dT_BGO[3] = bgo6_time-buncher_time;

      //Run the constraints subroutine. At this time, it does use timing gates
      fillHistograms(nGeConstraints, dGeBounds, dGeEn, dGeEn, dSiLiEn, dBGO, dT_GeEn, dT_SiLiEn, dT_BGO, true);
      fillHistograms(nSiLiConstraints, dSiLiBounds, dSiLiEn, dGeEn, dSiLiEn, dBGO, dT_GeEn, dT_SiLiEn, dT_BGO, false);
   }

   //File to write out to
   TFile* fOut = new TFile(fileName,"RECREATE");

   //write to file
   fOut->cd();
   writeHistToFile(fOut);
   fOut->Close();
}
