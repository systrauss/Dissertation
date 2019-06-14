//This file is for the creation of the arrays of histograms to be used with the cuts code. It also contains the fill histograms command.

#include <histograms.h>

//C++ libraries
#include <vector>
#include <iostream>

//ROOT libraries
#include <TH1F.h>
#include <TString.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TObject.h>

//Libraries within Code
#include "timing.h"

//From Coefficients.cxx
extern int nGeOrder; //Order of calibration i.e. 1 = linear.
extern int nGeDets; //Total number of signals from Germanium detectors
extern int nGeSegments; //number of segments in a single Germanium crystal, for adding purposes

extern int nSiLiOrder; //Order of calibration i.e. 1 = linear.
extern int nSiLiDets; //Total number of signals from SiLi
extern int nSiLiPlace; //Start of SiLis in generalized array detectors

//From Constraints.cxx
extern std::vector<std::vector<double> > dGeBounds; //bounds for cuts
extern std::vector<std::vector<double> > dSiLiBounds; //bounds for cuts
extern int nGeConstraints;
extern int nSiLiConstraints;

extern int nBGODets; //Total number of signals from BGO detectors
extern int nBGOPlace; //Start of BGOs in generalized array
extern std::vector<double> dBGOThreshold;

//From main.cxx
extern const Timing tGates;

//From this file
std::vector<std::vector<TH1F*> > ge_en_ge_cut;
std::vector<std::vector<TH1F*> > sili_en_ge_cut;
std::vector<std::vector<TH1F*> > ge_en_sili_cut;
std::vector<std::vector<TH1F*> > sili_en_sili_cut;


//For regular coincidence
void makeHistograms(int nGeDets, int nGeCuts, int nSiLiDets, int nSiLiCuts)
{
	//Vector to make a new row
	std::vector<TH1F*> row;
	//Loop through and construct them all.
	for(int j=0 ; j < nGeDets ; j++)
	{
		ge_en_ge_cut.push_back(row);
		ge_en_sili_cut.push_back(row);
		for(int i=0 ; i < nGeCuts ; i++)
		{
			ge_en_ge_cut[j].push_back( new TH1F(Form("ge_en_ge_cut[%i][%i]",j,i),Form("Clover_%i",j),4001,0.5,4000.5));
		}
		for(int i=0 ; i < nSiLiCuts ; i++)
		{
			ge_en_sili_cut[j].push_back(new TH1F(Form("sili_en_ge_cut[%i][%i]",j,i),Form("Clover_%i",j+1),4001,0.5,4000.5));
		}
	}
	for(int j=0 ; j < nSiLiDets ; j++)
	{
		sili_en_ge_cut.push_back(row);
		sili_en_sili_cut.push_back(row);
		for(int i=0 ; i < nGeCuts ; i++)
		{
			sili_en_ge_cut[j].push_back(new TH1F(Form("ge_en_sili_cut[%i][%i]",j,i),Form("SiLi_%i",j),4001,0.5,4000.5));
		}
		for(int i=0 ; i < nSiLiCuts ; i++)
		{
			sili_en_sili_cut[j].push_back(new TH1F(Form("sili_en_sili_cut[%i][%i]",j,i),Form("SiLi_%i",j+1),4001,0.5,4000.5));
		}
	}
}

//Regular, singular coincidence
void fillHistograms(int nConstraints, std::vector<std::vector<double> >dConstraints, std::vector<double> dEnCut, std::vector<double> dGeDetectors, std::vector<double> dSiLiDetectors, std::vector<double> dBGO, std::vector<double> dGeT, std::vector<double> dSiLiT, std::vector<double> dBGOT, bool GeorSiLi/*true  = ge gate, false = sili gate*/)
{
	for(int i=0; i<nConstraints ; i++)
	{
		if((dEnCut[(int)dConstraints[i][0]] > dConstraints[i][1] && //gate detector lower limit
			dEnCut[(int)dConstraints[i][0]] < dConstraints[i][2] && //gate detector upper limit
			GeorSiLi                                             //&& //is germanium
			//dBGO[(int)dConstraints[i][0]] < dBGOThreshold[(int)dConstraints[i][0]]) //BGO Threshold met
		)|| //OR (above is Ge constraint, below is SiLi)
		   (dEnCut[(int)dConstraints[i][0]] > dConstraints[i][1] && //gate detector lower limit
		   	dEnCut[(int)dConstraints[i][0]] < dConstraints[i][2] && //gate detector upper limit
		   	!GeorSiLi)) //is SiLi
		{
			//Fill Germanium histograms, includes BGO thresholds and timing
			for(int j = 0 ; j < nGeDets/nGeSegments ; j++)
			{
				if(GeorSiLi && (int)dConstraints[i][0] != j) // if Ge gate AND gate detector is not this detector AND BGO is under threshold)
				{
					//Timing check - ge-ge
					for(int k = 0; k < tGates.geGe.size(); k++)
					{
						//if( (((dGeT[(int)dConstraints[i][0]]-dGeT[j])>=tGates.geGe[k][2] && (dGeT[(int)dConstraints[i][0]]-dGeT[j])<=tGates.geGe[k][3]) || ((dGeT[j]-dGeT[(int)dConstraints[i][0]])>=tGates.geGe[k][2] && (dGeT[j]-dGeT[(int)dConstraints[i][0]])<=tGates.geGe[k][3]))) std::cout << "Within timing gate for Ge-Ge" << std::endl;
						if((int)dConstraints[i][0]==tGates.geGe[k][0] && j==tGates.geGe[k][1])
						{
							if( (dGeT[(int)dConstraints[i][0]]-dGeT[j])>=tGates.geGe[k][2] && (dGeT[(int)dConstraints[i][0]]-dGeT[j])<=tGates.geGe[k][3] && dGeDetectors[j] > 5) ge_en_ge_cut[j][i]->Fill(dGeDetectors[j]);
						}
						else if(j==tGates.geGe[k][0] && (int)dConstraints[i][0]==tGates.geGe[k][1])
						{
							if( (dGeT[j]-dGeT[(int)dConstraints[i][0]])>=tGates.geGe[k][2] && (dGeT[j]-dGeT[(int)dConstraints[i][0]])<=tGates.geGe[k][3] && dGeDetectors[j] > 5) ge_en_ge_cut[j][i]->Fill(dGeDetectors[j]);
						}
					}
				}
				else if (!GeorSiLi)
				{
					//Timing check ge-sili
					for(int k = 0; k < tGates.geSiLi.size(); k++)
					{
						if((int)dConstraints[i][0]==tGates.geSiLi[k][1] && j==tGates.geSiLi[k][0])
						{
							if( (dGeT[j]-dSiLiT[(int)dConstraints[i][0]])>=tGates.geSiLi[k][2] && (dGeT[j]-dSiLiT[(int)dConstraints[i][0]])<=tGates.geSiLi[k][3] && dGeDetectors[j] > 5) ge_en_sili_cut[j][i]->Fill(dGeDetectors[j]);
						}
					}
				}
			}
			//Fill SiLi histograms, includes timing
			for(int j = 0 ; j < nSiLiDets ; j++)
			{
				if(GeorSiLi)
				{
					//Timing check ge-sili
					for(int k = 0; k < tGates.geSiLi.size(); k++)
					{
						if(j==tGates.geSiLi[k][1] && (int)dConstraints[i][0]==tGates.geSiLi[k][0])
						{
							if( (dGeT[(int)dConstraints[i][0]]-dSiLiT[j])>=tGates.geSiLi[k][2] && (dGeT[(int)dConstraints[i][0]]-dSiLiT[j])<=tGates.geSiLi[k][3] && dSiLiDetectors[j] > 20) sili_en_ge_cut[j][i]->Fill(dSiLiDetectors[j]);
						}
					}
				}
				else if (!GeorSiLi && (int)dConstraints[i][0] != j)
				{
					//timing check sili-sili
					for(int k = 0; k < tGates.siliSiLi.size(); k++)
					{
						if((int)dConstraints[i][0]==tGates.siliSiLi[k][0] && j==tGates.siliSiLi[k][1])
						{
							if( (dSiLiT[(int)dConstraints[i][0]]-dSiLiT[j])>=tGates.siliSiLi[k][2] && (dSiLiT[(int)dConstraints[i][0]]-dSiLiT[j])<=tGates.siliSiLi[k][3] && dSiLiDetectors[j] > 20) sili_en_sili_cut[j][i]->Fill(dSiLiDetectors[j]);
						}
						else if(j==tGates.siliSiLi[k][0] && (int)dConstraints[i][0]==tGates.siliSiLi[k][1])
						{
							if( (dSiLiT[j]-dSiLiT[(int)dConstraints[i][0]])>=tGates.siliSiLi[k][2] && (dSiLiT[j]-dSiLiT[(int)dConstraints[i][0]])<=tGates.siliSiLi[k][3] && dSiLiDetectors[j] > 20) sili_en_sili_cut[j][i]->Fill(dSiLiDetectors[j]);
						}
					}
				}
			}
		}
	}
}


void writeHistToFile(TFile* fOut)
{
	for(int i=0; i<nGeConstraints ; i++)
	{
		fOut->cd();
		gDirectory->mkdir(Form("Ge Cut Detector %i, Peak %f",(int)dGeBounds[i][0],(dGeBounds[i][1]+dGeBounds[i][2])/2));
		fOut->cd(Form("Ge Cut Detector %i, Peak %f",(int)dGeBounds[i][0],(dGeBounds[i][1]+dGeBounds[i][2])/2));
		for(int j=0; j<nGeDets/nGeSegments ; j++)
		{
			ge_en_ge_cut[j][i]->Write(Form("Clover_%i",j),TObject::kOverwrite,1E8);
		}
		for(int j=0; j<nSiLiDets ; j++)
		{
			sili_en_ge_cut[j][i]->Write(Form("SiLi_%i",j+1),TObject::kOverwrite,1E8);
		}
	}
	for(int i=0; i<nSiLiConstraints ; i++)
	{
		fOut->cd();
		gDirectory->mkdir(Form("SiLi Cut Detector %i, Peak %f",(int)dSiLiBounds[i][0]+1,(dSiLiBounds[i][1]+dSiLiBounds[i][2])/2));
		fOut->cd(Form("SiLi Cut Detector %i, Peak %f",(int)dSiLiBounds[i][0]+1,(dSiLiBounds[i][1]+dSiLiBounds[i][2])/2));
		for(int j=0; j<nGeDets/nGeSegments ; j++)
		{
			ge_en_sili_cut[j][i]->Write(Form("Clover_%i",j),TObject::kOverwrite,1E8);
		}
		for(int j=0; j<nSiLiDets ; j++)
		{
			sili_en_sili_cut[j][i]->Write(Form("SiLi_%i",j+1),TObject::kOverwrite,1E8);
		}
	}
	fOut->Write();

}