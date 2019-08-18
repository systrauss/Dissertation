#include<TFile.h>
#include<TH1F.h>
#include<TF1.h>
#include<TMath.h>
#include<TRandom3.h>
#include<TFitResultPtr.h>

#include<sstream>
#include<string>
#include<iostream>

//linear+gaussian+skew+step guess, one iteration
void FitterSk1(TH1D* h1, double cent, double width, double rrr, double beta, double step, int xmin, int xmax)
{
	double pi = 3.14159265359;
//	TH1F* h1 = (TH1F*)inFile->Get(hist);
	TF1* lingaus = new TF1("lingaus","[0]+[1]*x+[2]*((1-[5]/100)*exp(-((x-[3])/(sqrt(2.0)*[4]))**2)+TMath::Erfc((x-[3])/(sqrt(2.0)*[4])+[4]/(sqrt(2.0)*[6]))*[5]/100*exp((x-[3])/[6])+TMath::Erfc((x-[3])/(sqrt(2.0)*[4]))*[7]/100)");
	lingaus->SetParameters(0,0,0,cent,width,rrr,beta,step);
	lingaus->FixParameter(3,cent);
	lingaus->FixParameter(4,width);
	lingaus->FixParameter(5,rrr);
	lingaus->FixParameter(6,beta);
	lingaus->FixParameter(7,step);
//	lingaus->SetParLimits(5,0,100.0); //R
//	lingaus->SetParLimits(7,0,100.0); //Step
	h1->Fit("lingaus","BLS","",xmin,xmax);
	lingaus->SetParLimits(3,cent-10,cent+10);
	lingaus->SetParLimits(4,0,width+10);

/*	for(int i=1; i<=500; i++)
	{
		h1->Fit("lingaus","BLQM","",xmin,xmax);
	}
	h1->Fit("lingaus","BLS","",xmin,xmax);
*/
}

//linear+gaussian+skew+step guess
void FitterSk(TH1D* h1, double cent, double width, double rrr, double beta, double step, int xMin, int xMax)
{
	double pi = 3.14159265359;
//	TH1F* h1 = (TH1F*)inFile->Get(hist);
	TF1* lin = new TF1("lin","[0]+[1]*x+[2]*TMath::Erfc((x-[3])/(sqrt(2.0)*[4]))*[5]/100");
	TF1* lingaus = new TF1("lingaus","[0]+[1]*x+[2]*((1-[5]/100)*exp(-((x-[3])/(sqrt(2.0)*[4]))**2)+TMath::Erfc((x-[3])/(sqrt(2.0)*[4])+[4]/(sqrt(2.0)*[6]))*[5]/100*exp((x-[3])/[6])+TMath::Erfc((x-[3])/(sqrt(2.0)*[4]))*[7]/100)");
	lingaus->SetParameters(0,0,0,cent,width,rrr,beta,step);
	lingaus->FixParameter(3,cent);
	lingaus->FixParameter(4,width);
	lingaus->FixParameter(5,rrr);
	lingaus->FixParameter(6,beta);
	lingaus->FixParameter(7,step);
	h1->Fit("lingaus","BLS","",xMin,xMax);	
	lingaus->SetParLimits(5,0,100.0); //R
	lingaus->SetParLimits(6,0,1000.0); //beta
	lingaus->SetParLimits(7,0,100.0); //Step
	lingaus->SetParLimits(3,cent-10,cent+10);
	lingaus->SetParLimits(4,0,width+10);

	for(int i=1; i<=100; i++)
	{
		h1->Fit("lingaus","BLQM","",xMin,xMax);
	}
	h1->Fit("lingaus","BLS","",xMin,xMax);

	TFitResultPtr r = h1->Fit("lingaus","LLS","",xMin,xMax);
	lin->SetParameter(0, lingaus->GetParameter(0));
	lin->SetParameter(1, lingaus->GetParameter(1));
	lin->SetParameter(2, lingaus->GetParameter(2));
	lin->SetParameter(3, lingaus->GetParameter(3));
	lin->SetParameter(4, lingaus->GetParameter(4));
	lin->SetParameter(5, lingaus->GetParameter(7));

	lin->SetParError(0, lingaus->GetParError(0));
	lin->SetParError(1, lingaus->GetParError(1));
	lin->SetParError(2, lingaus->GetParError(2));
	lin->SetParError(3, lingaus->GetParError(3));
	lin->SetParError(4, lingaus->GetParError(4));
	lin->SetParError(5, lingaus->GetParError(7));
	Double_t c0_err = lingaus->GetParError(0);
	Double_t c1_err = lingaus->GetParError(1);
	Double_t bkgd_integral_error = TMath::Sqrt( TMath::Power((xMax-xMin)*c0_err,2.0) + TMath::Power(0.5*(xMax-xMin)*(xMax-xMin)*c1_err,2.0) );
	Double_t bkgd_integral = lin->Integral(xMin,xMax);
	Double_t fn_integral = lingaus->Integral(xMin,xMax);
	Double_t fn_integral_error = lingaus->IntegralError( xMin, xMax, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray() );
	Double_t net_area = fn_integral - bkgd_integral;
	Double_t net_area_error = TMath::Sqrt( TMath::Power(fn_integral_error,2.0) + TMath::Power(bkgd_integral_error,2.0) );
	r->GetCovarianceMatrix().Print();
	printf("\n\nFit area = %10.6f+/-%10.6f\nBkgd area = %10.6f+/-%10.6f\nNet area = %10.6f+/-%10.6f\ \n",fn_integral,fn_integral_error,bkgd_integral,bkgd_integral_error,net_area,net_area_error);

}

//linear+gaussian guess
void FitterG(TFile* inFile, char* hist, double cent, double width, int xmin, int xmax)
{
	double pi = 3.14159265359;
	TH1F* h1 = (TH1F*)inFile->Get(hist);
	TF1* lingaus = new TF1("lingaus","[0]+[1]*x+[2]/(sqrt(2.0*pi)*[4])*exp(-((x-[3])/(sqrt(2.0)*[4]))**2)");
	lingaus->SetParNames("intercept","slope","area","centroid", "sigma");
	lingaus->SetParameters(0,0,0,cent,width);
	lingaus->FixParameter(3,cent);
	lingaus->FixParameter(4,width);
	h1->Fit("lingaus","BLS","",xmin,xmax);
	lingaus->SetParLimits(3,cent-5,cent+5);
	lingaus->SetParLimits(4,width-5,width+5);

	for(int i=1; i<=500; i++)
	{
		h1->Fit("lingaus","BLQM","",xmin,xmax);
	}
	h1->Fit("lingaus","BLS","",xmin,xmax);
}

//linear+gaussian*2 guess
void FitterG(TFile* inFile, char* hist, double cent, double width, double cent2, double width2, int xmin, int xmax)
{
	double pi = 3.14159265359;
	TH1F* h1 = (TH1F*)inFile->Get(hist);
	TF1* lingaus = new TF1("lingaus","[0]+[1]*x+[2]/(sqrt(2.0*pi)*[4])*exp(-((x-[3])/(sqrt(2.0)*[4]))**2)+[5]/(sqrt(2.0*pi)*[7])*exp(-((x-[6])/(sqrt(2.0)*[7]))**2)");
	lingaus->SetParNames("intercept","slope","area 1","centroid 1", "sigma 1","area 2","centroid 2", "sigma 2");
	lingaus->SetParameters(0,0,0,cent,width);
	lingaus->FixParameter(3,cent);
	lingaus->FixParameter(4,width);
	lingaus->FixParameter(6,cent2);
	lingaus->FixParameter(7,width2);
	h1->Fit("lingaus","BLS","",xmin,xmax);
	lingaus->SetParLimits(3,cent-5,cent+5);
	lingaus->SetParLimits(4,width-5,width+5);
	lingaus->SetParLimits(6,cent2-5,cent2+5);
	lingaus->SetParLimits(7,width2-5,width2+5);

	for(int i=1; i<=500; i++)
	{
		h1->Fit("lingaus","BLQM","",xmin,xmax);
	}
	h1->Fit("lingaus","BLS","",xmin,xmax);
}

//linear+gaussian*3 guess
void FitterG(TFile* inFile, char* hist, double cent, double width, double cent2, double width2, double cent3, double width3, int xmin, int xmax)
{
	double pi = 3.14159265359;
	TH1F* h1 = (TH1F*)inFile->Get(hist);
	TF1* lingaus = new TF1("lingaus","[0]+[1]*x+[2]/(sqrt(2.0*pi)*[4])*exp(-((x-[3])/(sqrt(2.0)*[4]))**2)+[5]/(sqrt(2.0*pi)*[7])*exp(-((x-[6])/(sqrt(2.0)*[7]))**2)+[8]/(sqrt(2.0*pi)*[10])*exp(-((x-[9])/(sqrt(2.0)*[10]))**2)");
	lingaus->SetParNames("intercept","slope","area 1","centroid 1", "sigma 1","area 2","centroid 2", "sigma 2","area 3","centroid 3", "sigma 3");
	lingaus->SetParameters(0,0,0,cent,width);
	lingaus->FixParameter(3,cent);
	lingaus->FixParameter(4,width);
	lingaus->FixParameter(6,cent2);
	lingaus->FixParameter(7,width2);
	lingaus->FixParameter(9,cent3);
	lingaus->FixParameter(10,width3);
	h1->Fit("lingaus","BLS","",xmin,xmax);
	lingaus->SetParLimits(3,cent-5,cent+5);
	lingaus->SetParLimits(4,width-5,width+5);
	lingaus->SetParLimits(6,cent2-5,cent2+5);
	lingaus->SetParLimits(7,width2-5,width2+5);
	lingaus->SetParLimits(9,cent3-5,cent3+5);
	lingaus->SetParLimits(10,width3-5,width3+5);

	for(int i=1; i<=500; i++)
	{
		h1->Fit("lingaus","BLQM","",xmin,xmax);
	}
	h1->Fit("lingaus","BLS","",xmin,xmax);
}

//guess for gamma singles, one old, one new peak
void FitterGS(TFile* inFile, char* hist, double cent, double width, double cent2, double width2, int xmin, int xmax)
{
	double pi = 3.14159265359;
	TH1F* h1 = (TH1F*)inFile->Get(hist);
	TF1* lingaus = new TF1("lingaus","[0]+[1]*x+[2]/(sqrt(2.0*pi)*[4])*exp(-((x-[3])/(sqrt(2.0)*[4]))**2)+[5]/(sqrt(2.0*pi)*[7])*exp(-((x-[6])/(sqrt(2.0)*[7]))**2)");
	lingaus->SetParameters(0,0,0,cent,width);
	lingaus->FixParameter(3,cent);
	lingaus->FixParameter(4,width);
	lingaus->FixParameter(6,cent2);
	lingaus->FixParameter(7,width2);
	h1->Fit("lingaus","BLS","",xmin,xmax);
	lingaus->SetParLimits(6,cent2-5,cent2+5);
	lingaus->SetParLimits(7,width2-5,width2+5);

	for(int i=1; i<=500; i++)
	{
		h1->Fit("lingaus","BLQM","",xmin,xmax);
	}
	h1->Fit("lingaus","BLS","",xmin,xmax);
}

//guess for gamma singles, one old, two new peaks
void FitterGS(TFile* inFile, char* hist, double cent, double width, double cent2, double width2, double cent3, double width3, int xmin, int xmax)
{
	double pi = 3.14159265359;
	TH1F* h1 = (TH1F*)inFile->Get(hist);
	TF1* lingaus = new TF1("lingaus","[0]+[1]*x+[2]/(sqrt(2.0*pi)*[4])*exp(-((x-[3])/(sqrt(2.0)*[4]))**2)+[5]/(sqrt(2.0*pi)*[7])*exp(-((x-[6])/(sqrt(2.0)*[7]))**2)+[8]/(sqrt(2.0*pi)*[10])*exp(-((x-[9])/(sqrt(2.0)*[10]))**2)");
	lingaus->SetParameters(0,0,0,cent,width);
	lingaus->FixParameter(3,cent);
	lingaus->FixParameter(4,width);
	lingaus->FixParameter(6,cent2);
	lingaus->FixParameter(7,width2);
	lingaus->FixParameter(9,cent3);
	lingaus->FixParameter(10,width3);
	h1->Fit("lingaus","BLS","",xmin,xmax);
	lingaus->SetParLimits(6,cent2-5,cent2+5);
	lingaus->SetParLimits(7,width2-5,width2+5);
	lingaus->SetParLimits(9,cent3-5,cent3+5);
	lingaus->SetParLimits(10,width3-5,width3+5);

	for(int i=1; i<=500; i++)
	{
		h1->Fit("lingaus","BLQM","",xmin,xmax);
	}
	h1->Fit("lingaus","BLS","",xmin,xmax);
}

//guess for gamma singles, two old, one new peaks
void FitterGS2(TFile* inFile, char* hist, double cent, double width, double cent2, double width2, double cent3, double width3, int xmin, int xmax)
{
	double pi = 3.14159265359;
	TH1F* h1 = (TH1F*)inFile->Get(hist);
	TF1* lingaus = new TF1("lingaus","[0]+[1]*x+[2]/(sqrt(2.0*pi)*[4])*exp(-((x-[3])/(sqrt(2.0)*[4]))**2)+[5]/(sqrt(2.0*pi)*[7])*exp(-((x-[6])/(sqrt(2.0)*[7]))**2)+[8]/(sqrt(2.0*pi)*[10])*exp(-((x-[9])/(sqrt(2.0)*[10]))**2)");
	lingaus->SetParameters(0,0,0,cent,width);
	lingaus->FixParameter(3,cent);
	lingaus->FixParameter(4,width);
	lingaus->FixParameter(6,cent2);
	lingaus->FixParameter(7,width2);
	lingaus->FixParameter(9,cent3);
	lingaus->FixParameter(10,width3);
	h1->Fit("lingaus","BLS","",xmin,xmax);
	lingaus->SetParLimits(9,cent3-5,cent3+5);
	lingaus->SetParLimits(10,width3-5,width3+5);

	for(int i=1; i<=500; i++)
	{
		h1->Fit("lingaus","BLQM","",xmin,xmax);
	}
	h1->Fit("lingaus","BLS","",xmin,xmax);
}

//linear+sine+gaussian, sine guess
void FitterGTS(TFile* inFile, char* hist, double cent, double width, double scale, double freq, double phase, int xmin, int xmax)
{
	double pi = 3.14159265359;
	TH1F* h1 = (TH1F*)inFile->Get(hist);
	TF1* lingaus = new TF1("lingaus","([0]+[1]*x+[5]/(sqrt(2.0*pi)*[7])*exp(-((x-[6])/(sqrt(2.0)*[7]))**2))*(1+[2]*sin([3]*x+[4]))");
	lingaus->SetParameters(0,0,0,cent,width);
	lingaus->FixParameter(6,cent);
	lingaus->FixParameter(7,width);
	lingaus->FixParameter(2,scale);
	lingaus->FixParameter(3,freq);
	lingaus->FixParameter(4,phase);
	h1->Fit("lingaus","BLS","",xmin,xmax);
	lingaus->SetParLimits(2,scale-scale,2);
	lingaus->SetParLimits(3,0.25,0.45);
	lingaus->SetParLimits(4,0,2*pi);

	for(int i=1; i<=500; i++)
	{
		h1->Fit("lingaus","BLQM","",xmin,xmax);
	}
	h1->Fit("lingaus","BLS","",xmin,xmax);
}

//linear+sine+gaussian*2, sine guess
void FitterGTS(TFile* inFile, char* hist, double cent, double width, double cent2, double width2, double scale, double freq, double phase, int xmin, int xmax)
{
	double pi = 3.14159265359;
	TH1F* h1 = (TH1F*)inFile->Get(hist);
	TF1* lingaus = new TF1("lingaus","([0]+[1]*x+[5]/(sqrt(2.0*pi)*[7])*exp(-((x-[6])/(sqrt(2.0)*[7]))**2)+[8]/(sqrt(2.0*pi)*[10])*exp(-((x-[9])/(sqrt(2.0)*[10]))**2))*(1+[2]*sin([3]*x+[4]))");
	lingaus->SetParameters(0,0,0,cent,width);
	lingaus->FixParameter(6,cent);
	lingaus->FixParameter(7,width);
	lingaus->FixParameter(9,cent2);
	lingaus->FixParameter(10,width2);
	lingaus->FixParameter(2,scale);
	lingaus->FixParameter(3,freq);
	lingaus->FixParameter(4,phase);
	h1->Fit("lingaus","BLS","",xmin,xmax);
	lingaus->SetParLimits(2,scale-scale,2);
	lingaus->SetParLimits(3,0.25,0.45);
	lingaus->SetParLimits(4,0,2*pi);

	for(int i=1; i<=500; i++)
	{
		h1->Fit("lingaus","BLQM","",xmin,xmax);
	}
	h1->Fit("lingaus","BLS","",xmin,xmax);
}

//linear+sine+gaussian*3, sine guess
void FitterGTS(TFile* inFile, char* hist, double cent, double width, double cent2, double width2, double cent3, double width3, double scale, double freq, double phase, int xmin, int xmax)
{
	double pi = 3.14159265359;
	TH1F* h1 = (TH1F*)inFile->Get(hist);
	TF1* lingaus = new TF1("lingaus","([0]+[1]*x+[5]/(sqrt(2.0*pi)*[7])*exp(-((x-[6])/(sqrt(2.0)*[7]))**2)+[8]/(sqrt(2.0*pi)*[10])*exp(-((x-[9])/(sqrt(2.0)*[10]))**2)+[11]/(sqrt(2.0*pi)*[13])*exp(-((x-[12])/(sqrt(2.0)*[13]))**2))*(1+[2]*sin([3]*x+[4]))");
	lingaus->SetParameters(0,0,0,cent,width);
	lingaus->FixParameter(6,cent);
	lingaus->FixParameter(7,width);
	lingaus->FixParameter(9,cent2);
	lingaus->FixParameter(10,width2);
	lingaus->FixParameter(12,cent3);
	lingaus->FixParameter(13,width3);
	lingaus->FixParameter(2,scale);
	lingaus->FixParameter(3,freq);
	lingaus->FixParameter(4,phase);
	h1->Fit("lingaus","BLS","",xmin,xmax);
	lingaus->SetParLimits(2,scale-scale,2);
	lingaus->SetParLimits(3,0.25,0.45);
	lingaus->SetParLimits(4,0,2*pi);

	for(int i=1; i<=500; i++)
	{
		h1->Fit("lingaus","BLQM","",xmin,xmax);
	}
	h1->Fit("lingaus","BLS","",xmin,xmax);
}

//linear+sine+gaussian, gaussian guess
void FitterGTG(TFile* inFile, char* hist, double cent, double width, double scale, double freq, double phase, int xmin, int xmax)
{
	double pi = 3.14159265359;
	TH1F* h1 = (TH1F*)inFile->Get(hist);
	TF1* lingaus = new TF1("lingaus","([0]+[1]*x+[5]/(sqrt(2.0*pi)*[7])*exp(-((x-[6])/(sqrt(2.0)*[7]))**2))*(1+[2]*sin([3]*x+[4]))");
	lingaus->SetParameters(0,0,0,cent,width);
	lingaus->FixParameter(6,cent);
	lingaus->FixParameter(7,width);
	lingaus->FixParameter(2,scale);
	lingaus->FixParameter(3,freq);
	lingaus->FixParameter(4,phase);
	h1->Fit("lingaus","BLS","",xmin,xmax);
	lingaus->SetParLimits(6,cent-5,cent+5);
	lingaus->SetParLimits(7,width-5,width+5);

	for(int i=1; i<=500; i++)
	{
		h1->Fit("lingaus","BLQM","",xmin,xmax);
	}
	h1->Fit("lingaus","BLS","",xmin,xmax);
}

//linear+sine+gaussian*2, gaussian*2 guess
void FitterGTG(TFile* inFile, char* hist, double cent, double width, double cent2, double width2, double scale, double freq, double phase, int xmin, int xmax)
{
	double pi = 3.14159265359;
	TH1F* h1 = (TH1F*)inFile->Get(hist);
	TF1* lingaus = new TF1("lingaus","([0]+[1]*x+[5]/(sqrt(2.0*pi)*[7])*exp(-((x-[6])/(sqrt(2.0)*[7]))**2)+[8]/(sqrt(2.0*pi)*[10])*exp(-((x-[9])/(sqrt(2.0)*[10]))**2))*(1+[2]*sin([3]*x+[4]))");
	lingaus->SetParameters(0,0,0,cent,width);
	lingaus->FixParameter(6,cent);
	lingaus->FixParameter(7,width);
	lingaus->FixParameter(9,cent2);
	lingaus->FixParameter(10,width2);
	lingaus->FixParameter(2,scale);
	lingaus->FixParameter(3,freq);
	lingaus->FixParameter(4,phase);
	h1->Fit("lingaus","BLS","",xmin,xmax);
	lingaus->SetParLimits(6,cent-5,cent+5);
	lingaus->SetParLimits(7,width-5,width+5);
	lingaus->SetParLimits(9,cent2-5,cent2+5);
	lingaus->SetParLimits(10,width2-5,width2+5);

	for(int i=1; i<=500; i++)
	{
		h1->Fit("lingaus","BLQM","",xmin,xmax);
	}
	h1->Fit("lingaus","BLS","",xmin,xmax);
}

//linear+sine+gaussian*3, gaussian*3 guess
void FitterGTG(TFile* inFile, char* hist, double cent, double width, double cent2, double width2,  double cent3, double width3, double scale, double freq, double phase, int xmin, int xmax)
{
	double pi = 3.14159265359;
	TH1F* h1 = (TH1F*)inFile->Get(hist);
	TF1* lingaus = new TF1("lingaus","([0]+[1]*x+[5]/(sqrt(2.0*pi)*[7])*exp(-((x-[6])/(sqrt(2.0)*[7]))**2)+[8]/(sqrt(2.0*pi)*[10])*exp(-((x-[9])/(sqrt(2.0)*[10]))**2)+[11]/(sqrt(2.0*pi)*[13])*exp(-((x-[12])/(sqrt(2.0)*[13]))**2))*(1+[2]*sin([3]*x+[4]))");
	lingaus->SetParameters(0,0,0,cent,width);
	lingaus->FixParameter(6,cent);
	lingaus->FixParameter(7,width);
	lingaus->FixParameter(9,cent2);
	lingaus->FixParameter(10,width2);
	lingaus->FixParameter(12,cent3);
	lingaus->FixParameter(13,width3);
	lingaus->FixParameter(2,scale);
	lingaus->FixParameter(3,freq);
	lingaus->FixParameter(4,phase);
	h1->Fit("lingaus","BLS","",xmin,xmax);
	lingaus->SetParLimits(6,cent-5,cent+5);
	lingaus->SetParLimits(7,width-5,width+5);
	lingaus->SetParLimits(9,cent2-5,cent2+5);
	lingaus->SetParLimits(10,width2-5,width2+5);
	lingaus->SetParLimits(12,cent3-5,cent3+5);
	lingaus->SetParLimits(13,width3-5,width3+5);

	for(int i=1; i<=500; i++)
	{
		h1->Fit("lingaus","BLQM","",xmin,xmax);
	}
	h1->Fit("lingaus","BLS","",xmin,xmax);
}

//linear+gaussian
void FitterI(TFile* inFile, char* hist, double cent, double width, int xMin, int xMax)
{
	double pi = 3.14159265359;
	TH1F* h1 = (TH1F*)inFile->Get(hist);
//	TF1* lingaus = new TF1("lingaus","pol1(0)+gaus(2)");
	TF1* lin = new TF1("lin","[0]+[1]*x");
	TF1* lingaus = new TF1("lingaus","[0]+[1]*x+[2]/(sqrt(2.0*pi)*[4])*exp(-((x-[3])/(sqrt(2.0)*[4]))**2)");
	lingaus->SetParameters(-20,0.022,800,0,0);
	lingaus->FixParameter(3,cent);
	lingaus->FixParameter(4,width);
//	h1->Fit("lingaus","BLQM","",xMin,xMax);
//	lingaus->SetParLimits(3,cent-5,cent+5);
//	lingaus->SetParLimits(4,width-5,width+5);
//	h1->Fit("lingaus","BLM","",xMin,xMax);

	for(int i=1; i<=500; i++)
	{
		h1->Fit("lingaus","LLQM","",xMin,xMax);
	}
	TFitResultPtr r = h1->Fit("lingaus","LLS","",xMin,xMax);
	lin->SetParameter(0, lingaus->GetParameter(0));
	lin->SetParameter(1, lingaus->GetParameter(1));
	Double_t bkgd_integral = lin->Integral(xMin,xMax);
	Double_t c0_err = lingaus->GetParError(0);
	Double_t c1_err = lingaus->GetParError(1);
	Double_t bkgd_integral_error = TMath::Sqrt( TMath::Power((xMax-xMin)*c0_err,2.0) + TMath::Power(0.5*(xMax-xMin)*(xMax-xMin)*c1_err,2.0) );
	
	Double_t fn_integral = lingaus->Integral(xMin,xMax);
	Double_t fn_integral_error = lingaus->IntegralError( xMin, xMax, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray() );
	Double_t net_area = fn_integral - bkgd_integral;
	Double_t net_area_error = TMath::Sqrt( TMath::Power(fn_integral_error,2.0) + TMath::Power(bkgd_integral_error,2.0) );
	r->GetCovarianceMatrix().Print();
	printf("\n\nFit area = %10.6f+/-%10.6f\nBkgd area = %10.6f+/-%10.6f\nNet area = %10.6f+/-%10.6f\n\n",fn_integral,fn_integral_error,bkgd_integral,bkgd_integral_error,net_area,net_area_error);
	
}

//linear+gaussian*2
void FitterI(TFile* inFile, char* hist, double cent, double width, double cent2, double width2, int xMin, int xMax)
{
	double pi = 3.14159265359;
	TH1F* h1 = (TH1F*)inFile->Get(hist);
//	TF1* lingaus = new TF1("lingaus","pol1(0)+gaus(2)");
	TF1* lin = new TF1("lin","[0]+[1]*x");
	TF1* lingaus = new TF1("lingaus","[0]+[1]*x+[2]/(sqrt(2.0*pi)*[4])*exp(-((x-[3])/(sqrt(2.0)*[4]))**2)+[5]/(sqrt(2.0*pi)*[7])*exp(-((x-[6])/(sqrt(2.0)*[7]))**2)");
	lingaus->SetParameters(-20,0.022,800,0,0);
	lingaus->FixParameter(3,cent);
	lingaus->FixParameter(4,width);
	lingaus->FixParameter(6,cent2);
	lingaus->FixParameter(7,width2);
//	h1->Fit("lingaus","BLQM","",xMin,xMax);
//	lingaus->SetParLimits(3,cent-5,cent+5);
//	lingaus->SetParLimits(4,width-5,width+5);
//	h1->Fit("lingaus","BLM","",xMin,xMax);

	for(int i=1; i<=500; i++)
	{
		h1->Fit("lingaus","LLQM","",xMin,xMax);
	}
	TFitResultPtr r = h1->Fit("lingaus","LLS","",xMin,xMax);
	lin->SetParameter(0, lingaus->GetParameter(0));
	lin->SetParameter(1, lingaus->GetParameter(1));
	Double_t bkgd_integral = lin->Integral(xMin,xMax);
	Double_t c0_err = lingaus->GetParError(0);
	Double_t c1_err = lingaus->GetParError(1);
	Double_t bkgd_integral_error = TMath::Sqrt( TMath::Power((xMax-xMin)*c0_err,2.0) + TMath::Power(0.5*(xMax-xMin)*(xMax-xMin)*c1_err,2.0) );
	
	Double_t fn_integral = lingaus->Integral(xMin,xMax);
	Double_t fn_integral_error = lingaus->IntegralError( xMin, xMax, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray() );
	Double_t net_area = fn_integral - bkgd_integral;
	Double_t net_area_error = TMath::Sqrt( TMath::Power(fn_integral_error,2.0) + TMath::Power(bkgd_integral_error,2.0) );
	r->GetCovarianceMatrix().Print();
	printf("\n\nFit area = %10.6f+/-%10.6f\nBkgd area = %10.6f+/-%10.6f\nNet area = %10.6f+/-%10.6f\n\n",fn_integral,fn_integral_error,bkgd_integral,bkgd_integral_error,net_area,net_area_error);
	
}

//linear+gaussian*3
void FitterI(TFile* inFile, char* hist, double cent, double width, double cent2, double width2, double cent3, double width3, int xMin, int xMax)
{
	double pi = 3.14159265359;
	TH1F* h1 = (TH1F*)inFile->Get(hist);
//	TF1* lingaus = new TF1("lingaus","pol1(0)+gaus(2)");
	TF1* lin = new TF1("lin","[0]+[1]*x");
	TF1* lingaus = new TF1("lingaus","[0]+[1]*x+[2]/(sqrt(2.0*pi)*[4])*exp(-((x-[3])/(sqrt(2.0)*[4]))**2)+[5]/(sqrt(2.0*pi)*[7])*exp(-((x-[6])/(sqrt(2.0)*[7]))**2)+[8]/(sqrt(2.0*pi)*[10])*exp(-((x-[9])/(sqrt(2.0)*[10]))**2)");
	lingaus->SetParameters(-20,0.022,800,0,0);
	lingaus->FixParameter(3,cent);
	lingaus->FixParameter(4,width);
	lingaus->FixParameter(6,cent2);
	lingaus->FixParameter(7,width2);
	lingaus->FixParameter(9,cent3);
	lingaus->FixParameter(10,width3);
//	h1->Fit("lingaus","BLQM","",xMin,xMax);
//	lingaus->SetParLimits(3,cent-5,cent+5);
//	lingaus->SetParLimits(4,width-5,width+5);
//	h1->Fit("lingaus","BLM","",xMin,xMax);

	for(int i=1; i<=500; i++)
	{
		h1->Fit("lingaus","LLQM","",xMin,xMax);
	}
	TFitResultPtr r = h1->Fit("lingaus","LLS","",xMin,xMax);
	lin->SetParameter(0, lingaus->GetParameter(0));
	lin->SetParameter(1, lingaus->GetParameter(1));
	Double_t bkgd_integral = lin->Integral(xMin,xMax);
	Double_t c0_err = lingaus->GetParError(0);
	Double_t c1_err = lingaus->GetParError(1);
	Double_t bkgd_integral_error = TMath::Sqrt( TMath::Power((xMax-xMin)*c0_err,2.0) + TMath::Power(0.5*(xMax-xMin)*(xMax-xMin)*c1_err,2.0) );
	
	Double_t fn_integral = lingaus->Integral(xMin,xMax);
	Double_t fn_integral_error = lingaus->IntegralError( xMin, xMax, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray() );
	Double_t net_area = fn_integral - bkgd_integral;
	Double_t net_area_error = TMath::Sqrt( TMath::Power(fn_integral_error,2.0) + TMath::Power(bkgd_integral_error,2.0) );
	r->GetCovarianceMatrix().Print();
	printf("\n\nFit area = %10.6f+/-%10.6f\nBkgd area = %10.6f+/-%10.6f\nNet area = %10.6f+/-%10.6f\n\n",fn_integral,fn_integral_error,bkgd_integral,bkgd_integral_error,net_area,net_area_error);
	
}

//linear+gaussian*4
void FitterI(TFile* inFile, char* hist, double cent, double width, double cent2, double width2, double cent3, double width3, double cent4, double width4, int xMin, int xMax)
{
	double pi = 3.14159265359;
	TH1F* h1 = (TH1F*)inFile->Get(hist);
//	TF1* lingaus = new TF1("lingaus","pol1(0)+gaus(2)");
	TF1* lin = new TF1("lin","[0]+[1]*x");
	TF1* lingaus = new TF1("lingaus","[0]+[1]*x+[2]/(sqrt(2.0*pi)*[4])*exp(-((x-[3])/(sqrt(2.0)*[4]))**2)+[5]/(sqrt(2.0*pi)*[7])*exp(-((x-[6])/(sqrt(2.0)*[7]))**2)+[8]/(sqrt(2.0*pi)*[10])*exp(-((x-[9])/(sqrt(2.0)*[10]))**2)+[11]/(sqrt(2.0*pi)*[13])*exp(-((x-[12])/(sqrt(2.0)*[13]))**2)");
	lingaus->SetParameters(-20,0.022,800,0,0);
	lingaus->FixParameter(3,cent);
	lingaus->FixParameter(4,width);
	lingaus->FixParameter(6,cent2);
	lingaus->FixParameter(7,width2);
	lingaus->FixParameter(9,cent3);
	lingaus->FixParameter(10,width3);
	lingaus->FixParameter(12,cent4);
	lingaus->FixParameter(13,width4);
//	h1->Fit("lingaus","BLQM","",xMin,xMax);
//	lingaus->SetParLimits(3,cent-5,cent+5);
//	lingaus->SetParLimits(4,width-5,width+5);
//	h1->Fit("lingaus","BLM","",xMin,xMax);

	for(int i=1; i<=500; i++)
	{
		h1->Fit("lingaus","LLQM","",xMin,xMax);
	}
	TFitResultPtr r = h1->Fit("lingaus","LLS","",xMin,xMax);
	lin->SetParameter(0, lingaus->GetParameter(0));
	lin->SetParameter(1, lingaus->GetParameter(1));
	Double_t bkgd_integral = lin->Integral(xMin,xMax);
	Double_t c0_err = lingaus->GetParError(0);
	Double_t c1_err = lingaus->GetParError(1);
	Double_t bkgd_integral_error = TMath::Sqrt( TMath::Power((xMax-xMin)*c0_err,2.0) + TMath::Power(0.5*(xMax-xMin)*(xMax-xMin)*c1_err,2.0) );
	
	Double_t fn_integral = lingaus->Integral(xMin,xMax);
	Double_t fn_integral_error = lingaus->IntegralError( xMin, xMax, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray() );
	Double_t net_area = fn_integral - bkgd_integral;
	Double_t net_area_error = TMath::Sqrt( TMath::Power(fn_integral_error,2.0) + TMath::Power(bkgd_integral_error,2.0) );
	r->GetCovarianceMatrix().Print();
	printf("\n\nFit area = %10.6f+/-%10.6f\nBkgd area = %10.6f+/-%10.6f\nNet area = %10.6f+/-%10.6f\n\n",fn_integral,fn_integral_error,bkgd_integral,bkgd_integral_error,net_area,net_area_error);
	
}

//linear+gaussian+sine
void FitterT(TFile* inFile, char* hist, double cent, double width, double scale, double freq, double phase, int xMin, int xMax)
{
	double pi = 3.14159265359;
	TH1F* h1 = (TH1F*)inFile->Get(hist);
	TF1* linsine = new TF1("linsine","([0]+[1]*x)*(1+[2]*sin([3]*x+[4]))");
	TF1* lingaus = new TF1("lingaus","([0]+[1]*x+[5]/(sqrt(2.0*pi)*[7])*exp(-((x-[6])/(sqrt(2.0)*[7]))**2))*(1+[2]*sin([3]*x+[4]))");
	lingaus->SetParameters(0,0,0,0,0,0,0,0);
	lingaus->FixParameter(6,cent);
	lingaus->FixParameter(7,width);
	lingaus->FixParameter(2,scale);
	lingaus->FixParameter(3,freq);
	lingaus->FixParameter(4,phase);

	for(int i=1; i<=500; i++)
	{
		h1->Fit("lingaus","BLQM","",xMin,xMax);
	}
	TFitResultPtr r = h1->Fit("lingaus","LLS","",xMin,xMax);
	linsine->SetParameter(0, lingaus->GetParameter(0));
	linsine->SetParameter(1, lingaus->GetParameter(1));
	linsine->SetParameter(2,scale);
	linsine->SetParameter(3,freq);
	linsine->SetParameter(4,phase);
	Double_t bkgd_integral = linsine->Integral(xMin,xMax);
//	Double_t c0_err = lingaus->GetParError(0);
//	Double_t c1_err = lingaus->GetParError(1);
//	Double_t bkgd_integral_error = TMath::Sqrt( TMath::Power((xMax-xMin)*c0_err,2.0) + TMath::Power(0.5*(xMax-xMin)*(xMax-xMin)*c1_err,2.0) );
	Double_t bkgd_integral_error = TMath::Sqrt(bkgd_integral);
	
	Double_t fn_integral = lingaus->Integral(xMin,xMax);
	Double_t fn_integral_error = lingaus->IntegralError( xMin, xMax, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray() );
	Double_t net_area = fn_integral - bkgd_integral;
	Double_t net_area_error = TMath::Sqrt( TMath::Power(fn_integral_error,2.0) + TMath::Power(bkgd_integral_error,2.0) );
	r->GetCovarianceMatrix().Print();
	printf("\n\nFit area = %10.6f+/-%10.6f\nBkgd area = %10.6f+/-%10.6f\nNet area = %10.6f+/-%10.6f\n\n",fn_integral,fn_integral_error,bkgd_integral,bkgd_integral_error,net_area,net_area_error);
}

//linear+gaussian*2+sine
void FitterT(TFile* inFile, char* hist, double cent, double width, double cent2, double width2, double scale, double freq, double phase, int xMin, int xMax)
{
	double pi = 3.14159265359;
	TH1F* h1 = (TH1F*)inFile->Get(hist);
	TF1* linsine = new TF1("linsine","([0]+[1]*x)*(1+[2]*sin([3]*x+[4]))");
	TF1* lingaus = new TF1("lingaus","([0]+[1]*x+[5]/(sqrt(2.0*pi)*[7])*exp(-((x-[6])/(sqrt(2.0)*[7]))**2)+[8]/(sqrt(2.0*pi)*[10])*exp(-((x-[9])/(sqrt(2.0)*[10]))**2))*(1+[2]*sin([3]*x+[4]))");
	lingaus->SetParameters(0,0,0,0,0,0,0,0);
	lingaus->FixParameter(6,cent);
	lingaus->FixParameter(7,width);
	lingaus->FixParameter(9,cent2);
	lingaus->FixParameter(10,width2);
	lingaus->FixParameter(2,scale);
	lingaus->FixParameter(3,freq);
	lingaus->FixParameter(4,phase);

	for(int i=1; i<=500; i++)
	{
		h1->Fit("lingaus","BLQM","",xMin,xMax);
	}
	TFitResultPtr r = h1->Fit("lingaus","LLS","",xMin,xMax);
	linsine->SetParameter(0, lingaus->GetParameter(0));
	linsine->SetParameter(1, lingaus->GetParameter(1));
	linsine->SetParameter(2,scale);
	linsine->SetParameter(3,freq);
	linsine->SetParameter(4,phase);
	Double_t bkgd_integral = linsine->Integral(xMin,xMax);
//	Double_t c0_err = lingaus->GetParError(0);
//	Double_t c1_err = lingaus->GetParError(1);
//	Double_t bkgd_integral_error = TMath::Sqrt( TMath::Power((xMax-xMin)*c0_err,2.0) + TMath::Power(0.5*(xMax-xMin)*(xMax-xMin)*c1_err,2.0) );
	Double_t bkgd_integral_error = TMath::Sqrt(bkgd_integral);
	
	Double_t fn_integral = lingaus->Integral(xMin,xMax);
	Double_t fn_integral_error = lingaus->IntegralError( xMin, xMax, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray() );
	Double_t net_area = fn_integral - bkgd_integral;
	Double_t net_area_error = TMath::Sqrt( TMath::Power(fn_integral_error,2.0) + TMath::Power(bkgd_integral_error,2.0) );
	r->GetCovarianceMatrix().Print();
	printf("\n\nFit area = %10.6f+/-%10.6f\nBkgd area = %10.6f+/-%10.6f\nNet area = %10.6f+/-%10.6f\n\n",fn_integral,fn_integral_error,bkgd_integral,bkgd_integral_error,net_area,net_area_error);
}

//linear+gaussian*3+sine
void FitterT(TFile* inFile, char* hist, double cent, double width, double cent2, double width2, double cent3, double width3, double scale, double freq, double phase, int xMin, int xMax)
{
	double pi = 3.14159265359;
	TH1F* h1 = (TH1F*)inFile->Get(hist);
	TF1* linsine = new TF1("linsine","([0]+[1]*x)*(1+[2]*sin([3]*x+[4]))");
	TF1* lingaus = new TF1("lingaus","([0]+[1]*x+[5]/(sqrt(2.0*pi)*[7])*exp(-((x-[6])/(sqrt(2.0)*[7]))**2)+[8]/(sqrt(2.0*pi)*[10])*exp(-((x-[9])/(sqrt(2.0)*[10]))**2)+[11]/(sqrt(2.0*pi)*[13])*exp(-((x-[12])/(sqrt(2.0)*[13]))**2))*(1+[2]*sin([3]*x+[4]))");
	lingaus->SetParameters(0,0,0,0,0,0,0,0);
	lingaus->FixParameter(6,cent);
	lingaus->FixParameter(7,width);
	lingaus->FixParameter(9,cent2);
	lingaus->FixParameter(10,width2);
	lingaus->FixParameter(12,cent3);
	lingaus->FixParameter(13,width3);
	lingaus->FixParameter(2,scale);
	lingaus->FixParameter(3,freq);
	lingaus->FixParameter(4,phase);

	for(int i=1; i<=500; i++)
	{
		h1->Fit("lingaus","BLQM","",xMin,xMax);
	}
	TFitResultPtr r = h1->Fit("lingaus","LLS","",xMin,xMax);
	linsine->SetParameter(0, lingaus->GetParameter(0));
	linsine->SetParameter(1, lingaus->GetParameter(1));
	linsine->SetParameter(2,scale);
	linsine->SetParameter(3,freq);
	linsine->SetParameter(4,phase);
	Double_t bkgd_integral = linsine->Integral(xMin,xMax);
//	Double_t c0_err = lingaus->GetParError(0);
//	Double_t c1_err = lingaus->GetParError(1);
//	Double_t bkgd_integral_error = TMath::Sqrt( TMath::Power((xMax-xMin)*c0_err,2.0) + TMath::Power(0.5*(xMax-xMin)*(xMax-xMin)*c1_err,2.0) );
	Double_t bkgd_integral_error = TMath::Sqrt(bkgd_integral);
	
	Double_t fn_integral = lingaus->Integral(xMin,xMax);
	Double_t fn_integral_error = lingaus->IntegralError( xMin, xMax, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray() );
	Double_t net_area = fn_integral - bkgd_integral;
	Double_t net_area_error = TMath::Sqrt( TMath::Power(fn_integral_error,2.0) + TMath::Power(bkgd_integral_error,2.0) );
	r->GetCovarianceMatrix().Print();
	printf("\n\nFit area = %10.6f+/-%10.6f\nBkgd area = %10.6f+/-%10.6f\nNet area = %10.6f+/-%10.6f\n\n",fn_integral,fn_integral_error,bkgd_integral,bkgd_integral_error,net_area,net_area_error);
}

//linear+gaussian*4+sine
void FitterT(TFile* inFile, char* hist, double cent, double width, double cent2, double width2, double cent3, double width3, double cent4, double width4, double scale, double freq, double phase, int xMin, int xMax)
{
	double pi = 3.14159265359;
	TH1F* h1 = (TH1F*)inFile->Get(hist);
	TF1* linsine = new TF1("linsine","([0]+[1]*x)*(1+[2]*sin([3]*x+[4]))");
	TF1* lingaus = new TF1("lingaus","([0]+[1]*x+[5]/(sqrt(2.0*pi)*[7])*exp(-((x-[6])/(sqrt(2.0)*[7]))**2)+[8]/(sqrt(2.0*pi)*[10])*exp(-((x-[9])/(sqrt(2.0)*[10]))**2)+[11]/(sqrt(2.0*pi)*[13])*exp(-((x-[12])/(sqrt(2.0)*[13]))**2)+[14]/(sqrt(2.0*pi)*[16])*exp(-((x-[15])/(sqrt(2.0)*[16]))**2))*(1+[2]*sin([3]*x+[4]))");
	lingaus->SetParameters(0,0,0,0,0,0,0,0);
	lingaus->FixParameter(6,cent);
	lingaus->FixParameter(7,width);
	lingaus->FixParameter(9,cent2);
	lingaus->FixParameter(10,width2);
	lingaus->FixParameter(12,cent3);
	lingaus->FixParameter(13,width3);
	lingaus->FixParameter(15,cent4);
	lingaus->FixParameter(16,width4);
	lingaus->FixParameter(2,scale);
	lingaus->FixParameter(3,freq);
	lingaus->FixParameter(4,phase);

	for(int i=1; i<=500; i++)
	{
		h1->Fit("lingaus","BLQM","",xMin,xMax);
	}
	TFitResultPtr r = h1->Fit("lingaus","LLS","",xMin,xMax);
	linsine->SetParameter(0, lingaus->GetParameter(0));
	linsine->SetParameter(1, lingaus->GetParameter(1));
	linsine->SetParameter(2,scale);
	linsine->SetParameter(3,freq);
	linsine->SetParameter(4,phase);
	Double_t bkgd_integral = linsine->Integral(xMin,xMax);
//	Double_t c0_err = lingaus->GetParError(0);
//	Double_t c1_err = lingaus->GetParError(1);
//	Double_t bkgd_integral_error = TMath::Sqrt( TMath::Power((xMax-xMin)*c0_err,2.0) + TMath::Power(0.5*(xMax-xMin)*(xMax-xMin)*c1_err,2.0) );
	Double_t bkgd_integral_error = TMath::Sqrt(bkgd_integral);
	
	Double_t fn_integral = lingaus->Integral(xMin,xMax);
	Double_t fn_integral_error = lingaus->IntegralError( xMin, xMax, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray() );
	Double_t net_area = fn_integral - bkgd_integral;
	Double_t net_area_error = TMath::Sqrt( TMath::Power(fn_integral_error,2.0) + TMath::Power(bkgd_integral_error,2.0) );
	r->GetCovarianceMatrix().Print();
	printf("\n\nFit area = %10.6f+/-%10.6f\nBkgd area = %10.6f+/-%10.6f\nNet area = %10.6f+/-%10.6f\n\n",fn_integral,fn_integral_error,bkgd_integral,bkgd_integral_error,net_area,net_area_error);
}
