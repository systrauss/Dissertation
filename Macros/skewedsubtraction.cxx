/* This macro is meant to subtract off a skewed gaussian of the given 
parameters. 

The conversion coefficient, area of the gamma, efficiency of the 
gamma detector, and efficiency of the electron detector are used to calculate 
the area of the skewed gaussian, and give an input to the height.

A best guess based on other data must be used for the width/sigma. R and beta
are taken from calibration data. The centroid is based off two things: the gamma
centroid and the electron orbital given.*/

TH1F* Subtraction(TH1F* hist, double dICC, double dGeEff, double dSiLiEff, double dGeArea, double dCorrRatio, double dGeCentroid, int nOrbital, double dR, double dBeta, double dSigma, int xMin, int xMax)
{
    double dSiLiCentroid = 0;
    double dSkewCentroid = 0;
    double dSiLiArea = 0;
    double dHeight = 0;
    TH1F* hSubtracted = hist->Clone();
    //Getting sili centroid
    switch(nOrbital)
    {
        case 0:
            dSiLiCentroid = dGeCentroid - 50.239;
            break;
        case 1:
            dSiLiCentroid = dGeCentroid - 7.9303;
            break;
        case 2:
            dSiLiCentroid = dGeCentroid - 1.3;
            break
        default:
            cout << "Invalid electron orbital" << endl;
            return nullptr;
    }
    //calculating area under sili peak
    dSiLiArea = dICC*dGeArea/dGeEff*dSiLiEff*1/dCorrRatio;
    cout << "SiLi Area: " << dSiLiArea << endl;
    //calculating the height from the area.
    dHeight = dSiLiArea*100/(2*exp(-dSigma**2/(2*dBeta**2))*dR*dBeta-sqrt(2*3.1415926)*(dR-100)*dSigma);
    cout << "SiLi Height: " << dHeight << endl;
    //calculating the skewed shift in the centroid
    dSkewCentroid = (2*dBeta*dBeta*exp(-dSigma*dSigma/(2*dBeta*dBeta)))/((dR-100)*sqrt(2*3.14159)*dSigma+2*dR*dBeta*exp(-dSigma*dSigma/(2*dBeta*dBeta)));
    std::cout << "Centroid shift: " << dSkewCentroid << std::endl;
    dSkewCentroid = dSiLiCentroid-dSkewCentroid;
    //Skewed gaussian function being set
    TF1* skewed = new TF1("skewed","[0]*(1-[3]/100)*exp(-((x-[1])/(sqrt(2.0)*[2]))**2)+[0]*[3]/100*exp((x-[5])/[4])*TMath::Erfc((x-[5])/(sqrt(2.0)*[2])+[2]/(sqrt(2.0)*[4]))", xMin,xMax);
    skewed->SetParNames("height","peak","sigma","R","beta","centroid");
    skewed->SetParameters(dHeight,dSiLiCentroid,dSigma,dR,dBeta,dSiLiCentroid);

    //begin the subtraction
    for (int i=xMin; i < xMax; i++)
    {
        hSubtracted->SetBinContent(i,hSubtracted->GetBinContent(i)-skewed->Eval(hSubtracted->GetBinCenter(i)));
    }

    TCanvas* c1 = new TCanvas();
    c1->Divide(1,2);
    c1->cd(1);
    hist->Draw();
    skewed->SetLineColor(1);
    skewed->Draw("same");
    hSubtracted->SetLineColor(2);
    hSubtracted->Draw("same");
    c1->cd(2);
    hSubtracted->Draw();
    return hSubtracted;
}