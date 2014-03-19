#define myutils_cxx
#include "myutils.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "TH1F.h"
#include <string>

TH1F * myutils::GetFitSlicesY(TH2F* h2d, int param){
    
    TObjArray aSlices;
    h2d->FitSlicesY(0, 0, -1, 0, "QNR", &aSlices);
    TH1F* h1 = (TH1F*) aSlices[1];
    TH1F* h2 = (TH1F*) aSlices[2];
    // have to clone because aslices goes away after de-scoping
    TH1F* hcl1 = (TH1F*) h1->Clone();
    TH1F* hcl2 = (TH1F*) h2->Clone();
    
    hcl1->SetLineColor( h2d->GetLineColor() );
    hcl1->SetMarkerColor( h2d->GetMarkerColor() );
    hcl1->SetFillColor( h2d->GetFillColor() );
    hcl1->SetLineStyle( h2d->GetLineStyle() );
    hcl1->SetMarkerStyle( h2d->GetMarkerStyle() );
    hcl1->SetFillStyle( h2d->GetFillStyle() );
    hcl2->SetLineColor( h2d->GetLineColor() );
    hcl2->SetMarkerColor( h2d->GetMarkerColor() );
    hcl2->SetFillColor( h2d->GetFillColor() );
    hcl2->SetLineStyle( h2d->GetLineStyle() );
    hcl2->SetMarkerStyle( h2d->GetMarkerStyle() );
    hcl2->SetFillStyle( h2d->GetFillStyle() );
    
    if (param == 1){ return hcl1; }
    else if (param == 2){ return hcl2; }
    else{ std::cout << "WARNING INVALID PARAM!!!!" << std::endl; }
}

TGraph * myutils::ROC(TH1F *&S,TH1F *&B,int np){
	//Usage:
	//TH1F *S = new TH1F();
	//TH1F *B = new TH1F();
	//...Fill S and B with signal (S) and background (B)
	//Need S and B to have the same number of bins!
	//TH1F *curve = new TH1F("","",M,0,1); where M is however fine you want the ROC curve binning to be.
	//ROC(S,B,curve);
    
    TH1F* curve = new TH1F("curve","curve",np,0,1);
    
	const int n = S->GetNbinsX();
	
	vector<Likelihood> mydata;
	double maxratio=-1;
	for (int i=0;i<n;i++){
		Likelihood onepoint;
		onepoint.S=S->GetBinContent(i+1);
		onepoint.B=B->GetBinContent(i+1);
		if (B->GetBinContent(i+1)>0){
			onepoint.ratio=S->GetBinContent(i+1)/B->GetBinContent(i+1);
			if (onepoint.ratio>maxratio){
				maxratio=onepoint.ratio;
			}
		}
		else if (S->GetBinContent(i+1)>0){
			onepoint.ratio=-2;
		}
		else{
			onepoint.ratio=-1;
		}
		mydata.push_back(onepoint);
	}
	for (int i=0; i<n; i++){
		if (mydata[i].ratio==-2){
			mydata[i].ratio=maxratio+1;
		}
	}
    
	sort(mydata.begin(), mydata.end(), by_ratio());
	double totalB = B->Integral();
	double totalS = S->Integral();
    
    for (int i=0; i<n; i++){
		double myS = 0.;
		double myB = 0.;
		for (int j=i; j<n; j++){
			myS+=mydata[j].S/totalS;
			myB+=mydata[j].B/totalB;
		}
		//std::cout << " " << i << " " << myS << " " << myB << " " << mydata[i].ratio << " " << mydata[i].B  << std::endl;
        //curve->SetBinContent(curve->FindBin(myS),1-myB);
        //std::cout << "curve->FindBin(myS),myB = " << curve->FindBin(myS) << ", " << myB << std::endl;
        curve->SetBinContent(curve->FindBin(myS),myB);
		//curve->SetBinContent(curve->FindBin(myS),1/max(0.00001,myB));
	}

    std::vector<float> x;
    std::vector<float> y;
    for (int i = 0; i < curve->GetNbinsX(); i++){
        if (curve->GetBinContent(i+1) != 0){
            x.push_back(curve->GetBinCenter(i+1));
            y.push_back(curve->GetBinContent(i+1));
        }
    }
    TGraph* tg = new TGraph(x.size(),&(x[0]),&(y[0]));

    delete curve;
    
    return tg;
}

