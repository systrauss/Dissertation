
//+ Combined (simultaneous) fit of two histogram with separate functions
//  and some common parameters
//
// See http://root.cern.ch/phpBB3//viewtopic.php?f=3&t=11740#p50908
// for a modified version working with Fumili or GSLMultiFit
//
// N.B. this macro must be compiled with ACliC
//
//Author: L. Moneta - Dec 2010

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"


// definition of shared parameter
// background function, left
int iparBL[3] = { 0,      // Constant background value
                 1,    // Linear background value
                 2, //Shift of "0"
};

//background function, right
int iparBR[3] = { 0, // Constant background value
                  1, // Linear background value
                  2, //Shift of "0"
};

struct GlobalChi2 {
   GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
                ROOT::Math::IMultiGenFunction & f2) :
      fChi2_1(&f1), fChi2_2(&f2) {}

   // parameter vector is first background (in common 1 and 2)
   // and then is signal (only in 2)
   double operator() (const double *par) const {
      double p1[3];
      for (int i = 0; i < 3; ++i) p1[i] = par[iparBL[i] ];

      double p2[3];
      for (int i = 0; i < 3; ++i) p2[i] = par[iparBR[i] ];

      return (*fChi2_1)(p1) + (*fChi2_2)(p2);
   }

   const  ROOT::Math::IMultiGenFunction * fChi2_1;
   const  ROOT::Math::IMultiGenFunction * fChi2_2;
};

void combinedFit(TH1F* hist, double xMin, double xStart, double xEnd, double xMax) {

//   TH1D * hB = new TH1D("hB","histo B",100,0,100);
//   TH1D * hSB = new TH1D("hSB","histo S+B",100, 0,100);

//   TF1 * fB = new TF1("fB","expo",0,100);
//   fB->SetParameters(1,-0.05);
//   hB->FillRandom("fB");

//   TF1 * fS = new TF1("fS","gaus",0,100);
//   fS->SetParameters(1,30,5);

//   hSB->FillRandom("fB",2000);
//   hSB->FillRandom("fS",1000);

  // perform now global fit

  //TF1* fBL = new TF1("fBL","pol1",xMin,xStart); 
  //TF1* fBR = new TF1("fBR","pol1",xEnd,xMax); 

  TF1* fBL = new TF1("fBL","[0]+[1]*(x-[2])",xMin,xStart); 
  TF1* fBR = new TF1("fBR","[0]+[1]*(x-[2])",xEnd,xMax); 
  fBL->FixParameter(2,xMin);
  fBR->FixParameter(2,xMin);
 // TF1 * fSB = new TF1("fSB","expo + gaus(2)",0,100);

  ROOT::Math::WrappedMultiTF1 wfBL(*fBL,1);
  ROOT::Math::WrappedMultiTF1 wfBR(*fBR,1);

  ROOT::Fit::DataOptions opt;
  ROOT::Fit::DataRange rangeBL;
  // set the data range
  rangeBL.SetRange(xMin,xStart);
  ROOT::Fit::BinData dataBL(opt,rangeBL);
  ROOT::Fit::FillData(dataBL, hist);

  ROOT::Fit::DataRange rangeBR;
  rangeBR.SetRange(xEnd,xMax);
  ROOT::Fit::BinData dataBR(opt,rangeBR);
  ROOT::Fit::FillData(dataBR, hist);

  ROOT::Fit::Chi2Function chi2_BL(dataBL, wfBL);
  ROOT::Fit::Chi2Function chi2_BR(dataBR, wfBR);

  GlobalChi2 globalChi2(chi2_BL, chi2_BR);

  ROOT::Fit::Fitter fitter;

  const int Npar = 3;
  double par0[Npar] = {0,0,xMin};

  // create before the parameter settings in order to fix or set range on them
  fitter.Config().SetParamsSettings(3,par0);
//   // fix 5-th parameter
  fitter.Config().ParSettings(2).Fix();
//   // set limits on the third and 4-th parameter
//   fitter.Config().ParSettings(2).SetLimits(-10,-1.E-4);
//   fitter.Config().ParSettings(3).SetLimits(0,10000);
//   fitter.Config().ParSettings(3).SetStepSize(5);

  fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().SetMinimizer("Minuit2","Migrad");

  // fit FCN function directly
  // (specify optionally data size and flag to indicate that is a chi2 fit)
  fitter.FitFCN(3,globalChi2,0,dataBL.Size()+dataBR.Size(),true);
  ROOT::Fit::FitResult result = fitter.Result();
  result.Print(std::cout);
  hist->Draw();
  fBL->SetFitResult( result, iparBL);
  fBL->SetRange(rangeBL().first, rangeBL().second);
  fBL->Draw("same");
  fBR->SetFitResult( result, iparBR);
  fBR->SetRange(rangeBR().first, rangeBR().second);
  fBR->Draw("same");

  TF1* fB = new TF1("fB","[0]+[1]*(x-[2])");
  fB->SetFitResult( result, iparBL);
  double dBkdSigma = sqrt(TMath::Power((xEnd-xStart)*fB->GetParError(0),2)+TMath::Power(0.5*fB->GetParError(1)*(TMath::Power(xEnd,2)-TMath::Power(xStart,2)-2*xMin),2));
  //double dBkdArea = fB->GetParameter(0)*(xEnd-xStart)+fB->GetParameter(1)*(0.5*xEnd*xEnd - 0.5*xStart*xStart - xMin);
  double dBkdArea = fB->Integral(xStart,xEnd);

  double dArea = 0;
  for(int i=xStart; i <= xEnd; i++)
  {
      if(hist->GetBinContent(i)>0)
      { 
        dArea = dArea + hist->GetBinContent(i)-fB->Eval(i+0.5);
     // std::cout << i << "\t" << hist->GetBinContent(i) << "\t" << fB->Eval(i+0.5) << "\t" << hist->GetBinContent(i)-fB->Eval(i+0.5) << std::endl;
      }
  }
  std::cout << "Centroid: " << (xStart+xEnd)/2.0 << std::endl; 
  std::cout << "Area: " << dArea << " +/- " << sqrt(dArea) << " +/- " << sqrt(dBkdArea) << std::endl;
  std::cout << "Background Area: " << dBkdArea << std::endl;

//   TCanvas * c1 = new TCanvas("Simfit","Simultaneous fit of two histograms",
//                              10,10,700,700);
//   c1->Divide(1,2);
//   c1->cd(1);
//   gStyle->SetOptFit(1111);

//   fB->SetFitResult( result, iparB);
//   fB->SetRange(rangeB().first, rangeB().second);
//   fB->SetLineColor(kBlue);
//   hB->GetListOfFunctions()->Add(fB);
//   hB->Draw();

//   c1->cd(2);
//   fSB->SetFitResult( result, iparSB);
//   fSB->SetRange(rangeSB().first, rangeSB().second);
//   fSB->SetLineColor(kRed);
//   hSB->GetListOfFunctions()->Add(fSB);
//   hSB->Draw();


}

void Piecewise(TH1F* hist, double xMin, double xStart, double xEnd, double xMax)
{
    TF1* piecewise = new TF1("piecewise","([0]+[1]*x)*(TMath::Erfc((x-[2])*1000)+TMath::Erf((x-[3])*1000))", xMin, xMax);
    piecewise->FixParameter(2,xStart);
    piecewise->FixParameter(3,xEnd);
    piecewise->FixParameter(1,0);
    piecewise->FixParameter(0,36223);

    hist->Fit("piecewise","BLS","", xMin, xMax);
}
