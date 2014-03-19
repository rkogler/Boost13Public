#ifndef myutils_h
#define myutils_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TVectorF.h"

class myutils {
    public :
    
    myutils();
    ~myutils();
    TH1F* GetFitSlicesY(TH2F* h2d, int param);
    TGraph * ROC(TH1F *&S,TH1F *&B,int np);

    struct QM {
        double q;
        double m;
    };
    
    struct by_charge {
        bool operator()(QM const &a, QM const &b) {
            return a.q < b.q;
        }
    };
    
    struct Likelihood {
        double S;
        double B;
        double ratio;
    };
    
    struct by_ratio {
        bool operator()(Likelihood const &a, Likelihood const &b) {
            return a.ratio < b.ratio;
        }
    };
    
};

#endif

#ifdef myutils_cxx
myutils::myutils()
{
    std::cout << "building utilities" << std::endl;
}

myutils::~myutils()
{
    
}


#endif // #ifdef vJetSubstructureAnalysis_cxx