// -*- C++ -*-
//System includes
#include <map>
#include <algorithm>

// BOOST 2012 Substructure methods
#include "BOOSTFastJets.h"

//Generator Interfaces
#include "HepMC/GenParticle.h"
#include "HepMC/GenEvent.h"

//Rivet framework
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

//Projections
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/FastJets.hh"
//#include "Rivet/Projections/MissingMomentum.hh"

// Histogram booking
#include "LWH/Histogram1D.h"
//#include "LWH/Histogram2D.h"

typedef std::map<std::string,AIDA::IHistogram1D*> BookedHistos;
namespace Rivet
{

/// Generic analysis looking at various distributions of final state particles
class MC_GENSTUDY_JETCHARGE : public Analysis
{
public:
    /// Constructor
    MC_GENSTUDY_JETCHARGE()
        : Analysis("MC_GENSTUDY_JETCHARGE")
    {
        for(unsigned int i=0; i < 4; i++) _nPassing[i]=0;
    }

    /// @name Analysis methods
    //@{
    /// Book histograms and initialise projections before the run
    void init()
    {
        // Projections
        const FinalState fs;
        addProjection(fs, "FS");
        std::vector<std::pair<double, double > > muonEtaRanges;
        muonEtaRanges.push_back(make_pair(-2.4,2.4));
        WFinder muWFinder(fs, muonEtaRanges, 25*GeV, MUON,
                          40*GeV,1000*GeV,25*GeV,0.6,true,false,80.4,true);
        // Tag a W in the event, focus on jets that don't come from the W decay.
        addProjection(muWFinder,"muWFinder");
        FastJets JetProjection(muWFinder.remainingFinalState(),FastJets::ANTIKT, 0.6); //FastJets::KT,0.7
        addProjection(JetProjection,"Jets");
        ///////////////
        // Histograms
        ///////////////
        _histograms["JetMult"]		= bookHistogram1D("JetMult"		, 6, -0.5, 5.5);
        //Jet Kinematics
        _histograms["JetPt"]		= bookHistogram1D("JetPt"		, 50, 33, 300);
        _histograms["JetE"]		= bookHistogram1D("JetE"		, 25, 20, 300);
        _histograms["JetEta"]		= bookHistogram1D("JetEta"		, 25, -2, 2);
        _histograms["JetRapidity"]	= bookHistogram1D("JetRapidity"		, 25, -2, 2);
        //_histograms["JetPhi"]		= bookHistogram1D("JetPhi"		, 25, 0, TWOPI);
        _histograms["JetMass"]		= bookHistogram1D("JetMass"		, 100, 0, 40);
        //_histograms["Jet2Mass"]		= bookHistogram1D("Jet2Mass"		, 100, 0, 100);
        //_histograms["Jet3Mass"]		= bookHistogram1D("Jet3Mass"		, 100, 0, 100);
        _histograms["SubJetMult"]		= bookHistogram1D("SubJetMult"		, 15, -0.5, 29.5);
        _histograms["SubJet2Mass"]	= bookHistogram1D("SubJet2Mass"		, 100, 0, 35);
        _histograms["SubJet3Mass"]	= bookHistogram1D("SubJet3Mass"		, 100, 0, 45);
        _histograms["SubJetDeltaR"]	= bookHistogram1D("SubJetDeltaR"	, 50, 0, 1.0);
        _histograms["SubJetMass"]		= bookHistogram1D("SubJetMass"		, 100, 0, 12);
        _histograms["SubJetSumEt"]	= bookHistogram1D("SubJetSumEt"		, 30, 0, 175);

        //Jet Charge Histos
        _histograms["WCharge"]		= bookHistogram1D("WCharge"		, 3, -1.5, 1.5);

        _histograms["WJetChargeK5"]	= bookHistogram1D("WJetChargeK5"	, 50, -3, 3);
        _histograms["QuarkJetChargeK5"]	= bookHistogram1D("QuarkJetChargeK5"	, 50, -3, 3);
        _histograms["GluonJetChargeK5"]	= bookHistogram1D("GluonJetChargeK5"	, 50, -3, 3);

        _histograms["WJetChargeK3"]	= bookHistogram1D("WJetChargeK3"	, 50, -3, 3);
        _histograms["QuarkJetChargeK3"]	= bookHistogram1D("QuarkJetChargeK3"	, 50, -3, 3);
        _histograms["GluonJetChargeK3"]	= bookHistogram1D("GluonJetChargeK3"	, 50, -3, 3);

        _histograms["QuarkNegTwoThirdsK5"]= bookHistogram1D("QuarkNegTwoThirdsK5"	, 50, -3, 3);
        _histograms["QuarkNegOneThirdK5"]	= bookHistogram1D("QuarkNegOneThirdK5"	, 50, -3, 3);
        _histograms["QuarkOneThirdK5"]	= bookHistogram1D("QuarkOneThirdK5"	, 50, -3, 3);
        _histograms["QuarkTwoThirdsK5"]	= bookHistogram1D("QuarkTwoThirdsK5"	, 50, -3, 3);

        _histograms["QuarkNegTwoThirdsK3"]= bookHistogram1D("QuarkNegTwoThirdsK3"	, 50, -3, 3);
        _histograms["QuarkNegOneThirdK3"]	= bookHistogram1D("QuarkNegOneThirdK3"	, 50, -3, 3);
        _histograms["QuarkOneThirdK3"]	= bookHistogram1D("QuarkOneThirdK3"	, 50, -3, 3);
        _histograms["QuarkTwoThirdsK3"]	= bookHistogram1D("QuarkTwoThirdsK3"	, 50, -3, 3);

        _histograms["ChargeSignPurity"]   = bookHistogram1D("ChargeSignPurity"    ,50,33,300);
        _histograms["QuarkJetEta"]	= bookHistogram1D("QuarkJetEta"		, 25, -2, 2);
        _histograms["GluonJetEta"]	= bookHistogram1D("GluonJetEta"		, 25, -2, 2);
        _histograms["QuarkJetPt"]         = bookHistogram1D("QuarkJetPt"          ,50,33,300);
        _histograms["GluonJetPt"]         = bookHistogram1D("GluonJetPt"          ,50,33,300);

        _histograms["JetPullTheta"]       = bookHistogram1D("JetPullTheta"	,50,-PI,PI);
        _histograms["JetPullMag"]         = bookHistogram1D("JetPullMag"          ,50,0,0.04);
        //_hist2DJetChargeWPt		= bookHistogram2D("JetChargeVsWPt"	,50,-0.3,0.3,50,24,300);
        _histograms["TruthDeltaR"]        = bookHistogram1D("TruthDeltaR"         ,50,0,0.7);
        _histograms["TruthPdgID"]         = bookHistogram1D("TruthPdgID"          ,7,-0.5,6.5);
        //Dipolarity
        _histograms["Dipolarity"]         =  bookHistogram1D("Dipolarity"          ,50,0.0,1.5);

        //N-subjettiness histos
        _histograms["JetMassFilt"]	= bookHistogram1D("JetMassFilt"		, 60, 0, 50);
        _histograms["JetMassTrim"]	= bookHistogram1D("JetMassTrim"		, 60, 0, 50);
        _histograms["JetMassPrune"]	= bookHistogram1D("JetMassPrune"	, 60, 0, 20);
        _histograms["NSubJettiness"]	= bookHistogram1D("NSubJettiness"	, 40, -0.005, 1.005);
        _histograms["NSubJettiness1Iter"]	= bookHistogram1D("NSubJettiness1Iter"	, 40, -0.005, 1.005);
        _histograms["NSubJettiness2Iter"]	= bookHistogram1D("NSubJettiness2Iter"	, 40, -0.005, 1.005);
    }
    /// quickly calculate standard deviation of pt distribution in jets
    virtual void pt_stddev(const PseudoJets& jets, double& mean,double& stddev,const double N)
    {
        foreach(const fastjet::PseudoJet& jet, jets)
        mean+=jet.pt();
        mean=mean/N;
        foreach(const fastjet::PseudoJet& jet, jets)
        stddev+=((jet.pt()-mean)*(jet.pt()-mean));
        stddev=stddev/N;
    }
    virtual void fillChargeHistograms(const fastjet::PseudoJet& jet, const FastJets& JetProjection,
                                      const double k, const int wCharge,
                                      const double weight, const int pdgId)
    {
        stringstream kStr;
        kStr<<"K"<<static_cast<int>(k*10);
        const double jetCharge = wCharge*JetCharge(JetProjection,jet,k,1*GeV);
        _histograms["WJetCharge"+kStr.str()]->fill(jetCharge,weight);
        if(abs(pdgId) < 7)
        {
            _histograms["QuarkJetCharge"+kStr.str()]->fill(jetCharge,weight);
            switch(wCharge*PID::threeCharge(pdgId))
            {
            case -2:
                _histograms["QuarkNegTwoThirds"+kStr.str()]->fill(jetCharge,weight);
                break;
            case -1:
                _histograms["QuarkNegOneThird"+kStr.str()]->fill(jetCharge,weight);
                break;
            case 1:
                _histograms["QuarkOneThird"+kStr.str()]->fill(jetCharge,weight);
                break;
            case 2:
                _histograms["QuarkTwoThirds"+kStr.str()]->fill(jetCharge,weight);
                break;
            }
        }
        else if(abs(pdgId)  == 21)
        {
            _histograms["GluonJetCharge"+kStr.str()]->fill(jetCharge,weight);
        }

    }

    virtual void analyzeSubJets(const fastjet::PseudoJet& jet,const double weight)
    {
        const double ptmin=0.5*GeV;
        double sumEt=0.0;
        fastjet::ClusterSequence clusterSeq(jet.validated_cs()->constituents(jet),fastjet::JetDefinition(fastjet::kt_algorithm,0.6));

        PseudoJets subJets=clusterSeq.exclusive_jets_up_to(3);

        fastjet::ClusterSequence antiKTClusterSeq(jet.validated_cs()->constituents(jet),fastjet::JetDefinition(fastjet::antikt_algorithm,0.1));
        PseudoJets smallSubJets=antiKTClusterSeq.inclusive_jets(ptmin);
        int smallJetMult = smallSubJets.size();
        _histograms["SubJetMult"]->fill(smallJetMult,weight);
        unsigned int nSubJets=subJets.size();

        if(nSubJets==3)
            _histograms["SubJet3Mass"]->fill((subJets.at(0)+subJets.at(1)+subJets.at(2)).m(),weight);

        for(unsigned int j=0; j!=nSubJets; ++j)
        {
            sumEt+=subJets.at(j).Et();
            _histograms["SubJetMass"]->fill(subJets.at(j).m());
            for(unsigned int k=(j+1); k!=nSubJets; ++k)
            {
                _histograms["SubJetDeltaR"]->fill(subJets.at(j).delta_R(subJets.at(k)),weight);
                _histograms["SubJet2Mass"]->fill((subJets.at(j)+subJets.at(k)).m(),weight);
            }
        }
        _histograms["SubJetSumEt"]->fill(sumEt,weight);
    }
    /// Perform the per-event analysis
    void analyze(const Event& event)
    {
        _nPassing[0]++;
        const WFinder& muWFinder = applyProjection<WFinder>(event,"muWFinder");
        if (muWFinder.bosons().size() != 1)
            vetoEvent;
        _nPassing[1]++;
        //Dipolarity(0.,0.);
        const double weight = event.weight();
        const FastJets& JetProjection=applyProjection<FastJets>(event, "Jets");
        const PseudoJets& jets = JetProjection.pseudoJetsByPt(35.0*GeV);
        if (jets.size() > 0)
        {
            _nPassing[2]++;
            const unsigned int jetMult=jets.size();
            _histograms["JetMult"]->fill(jetMult);
            /// Rather than loop over all jets, just take the first hard
            /// one, Make sure entire jet is within fiducial volume
            if(jets.front().eta() > -(2.5-0.6) && jets.front().eta() < (2.5-0.6))
            {
                if(jets.front().has_valid_cs())
                    analyzeSubJets(jets.front(),weight);

                foreach (const fastjet::PseudoJet& jet, jets)
                {
                    _histograms["JetMassFilt"]->fill(Filter(JetProjection.clusterSeq(),jet, FastJets::CAM, 3, 0.3).m(), weight);
                    _histograms["JetMassTrim"]->fill(Trimmer(JetProjection.clusterSeq(),jet, FastJets::CAM, 0.03, 0.3).m(), weight);
                    _histograms["JetMassPrune"]->fill(Pruner(JetProjection.clusterSeq(),jet, FastJets::CAM, 0.4, 0.1).m(), weight);
                    PseudoJets constituents = jet.constituents();
                    if (constituents.size() > 10)
                    {
                        PseudoJets axes(GetAxes(JetProjection.clusterSeq(), 2, constituents, FastJets::CAM, 0.5));
                        _histograms["NSubJettiness"]->fill(TauValue(2, 1, constituents, axes), weight);
                        UpdateAxes(2, constituents, axes);
                        _histograms["NSubJettiness1Iter"]->fill(TauValue(2, 1, constituents, axes), weight);
                        UpdateAxes(2, constituents, axes);
                        _histograms["NSubJettiness2Iter"]->fill(TauValue(2, 1, constituents, axes), weight);
                    }
                }
                _nPassing[3]++;
                const double wCharge=PID::charge(muWFinder.bosons().front().pdgId());
                //const double jetCharge=wCharge*JetProjection.JetCharge(jets.front(),0.5,1*GeV);
                const std::pair<double,double> tvec=JetPull(JetProjection,jets.front());
                _histograms["Dipolarity"]->fill(Dipolarity(jets.front()),weight);
                _histograms["JetMass"]->fill(jets.front().m(),weight);
                _histograms["JetPt"]->fill(jets.front().pt(),weight);
                _histograms["JetE"]->fill(jets.front().E(),weight);
                _histograms["JetEta"]->fill(jets.front().eta(),weight);
                _histograms["JetRapidity"]->fill(jets.front().rapidity(),weight);
                //histograms["JetPhi"]->fill(jets.front().phi(),weight);
                //_hist2DJetChargeWPt->fill(jetCharge,muWFinder.bosons().front().momentum().pT(),weight);
                //_histograms["WJetCharge"]->fill(jetCharge,weight);
                _histograms["WCharge"]->fill(wCharge,weight);
                _histograms["JetPullMag"]->fill(tvec.first,weight);
                if(tvec.first > 0)
                {
                    _histograms["JetPullTheta"]->fill(tvec.second,weight);
                }
                HepMC::GenParticle* truthParton=NULL;
                double truthDelR(0);
                foreach (HepMC::GenParticle* const p, particles(event.genEvent()))
                {
                    if((p->pdg_id() != 21) and (abs(p->pdg_id()) > 6)) continue;
                    //This may be slow, but its the path of minimal obfuscation
                    const double delR = jets.front().delta_R(fastjet::PseudoJet(p->momentum().px(),
                                        p->momentum().py(),
                                        p->momentum().pz(),
                                        p->momentum().e()));
                    if(truthParton==NULL)
                    {
                        truthDelR = delR;
                        truthParton = p;
                    }
                    else if(delR < 0.6 && truthParton->momentum().perp() < p->momentum().perp())
                    {
                        //if(null) assign truthParton, updateDelR
                        truthDelR = delR;
                        truthParton = p;
                    }
                }
                _histograms["TruthDeltaR"]->fill(truthDelR,weight);
                truthParton = NULL;
                foreach (HepMC::GenParticle* const p, particles(event.genEvent()))
                {
                    //if(null) assign truthParton, updateDelR
                    if((p->pdg_id() != 21) and (abs(p->pdg_id()) > 6)) continue;
                    //This may be slow, but its the path of minimal obfuscation
                    const double delR = jets.front().delta_R(fastjet::PseudoJet(p->momentum().px(),
                                        p->momentum().py(),
                                        p->momentum().pz(),
                                        p->momentum().e()));
                    if(truthParton==NULL)
                    {
                        truthDelR = delR;
                        truthParton = p;
                    }
                    else if(delR < 0.4 && truthParton->momentum().perp() < p->momentum().perp())
                    {
                        //if(null) assign truthParton, updateDelR
                        truthDelR = delR;
                        truthParton = p;
                    }
                }
                const int pdgId = truthParton->pdg_id();
                _histograms["TruthPdgID"]->fill((abs(pdgId)==21) ? 0 :abs(pdgId), weight);
                fillChargeHistograms(jets.front(), JetProjection, 0.3, static_cast<int>(wCharge), weight, pdgId);
                fillChargeHistograms(jets.front(), JetProjection, 0.5, static_cast<int>(wCharge), weight, pdgId);
                if(abs(pdgId) < 7)
                {
                    _histograms["QuarkJetPt"]->fill(jets.front().pt(),weight);
                    _histograms["QuarkJetEta"]->fill(jets.front().eta(),weight);
                    if(wCharge*PID::charge(pdgId) < 0.0)
                    {
                        _histograms["ChargeSignPurity"]->fill(jets.front().pt(),weight);
                    }
                }
                else if(pdgId == 21)
                {
                    _histograms["GluonJetPt"]->fill(jets.front().pt(),weight);
                    _histograms["GluonJetEta"]->fill(jets.front().eta(),weight);
                }
            }
        }
        else
        {
            vetoEvent;
        }
    }
    /// Finalize
    void finalize()
    {
        cout<<"Cut summary: "<<endl;
        cout<<"| Inclusive | "<<_nPassing[0]<< " | "<<endl;
        cout<<"| Found W   | "<<_nPassing[1]<< " | "<<endl;
        cout<<"| >1 Jet    | "<<_nPassing[2]<< " | "<<endl;
        cout<<"| Fiducial  | "<<_nPassing[3]<< " | "<<endl;
        cout<<"Mean Jet Charge (k=0.3): "<<_histograms["WJetChargeK3"]->mean()<<" +/- "<<_histograms["WJetChargeK3"]->rms()<<endl;
        cout<<"Mean Jet Charge (k=0.5): "<<_histograms["WJetChargeK5"]->mean()<<" +/- "<<_histograms["WJetChargeK5"]->rms()<<endl;

        // foreach(BookedHistos::value_type H,_histograms){
        // 	normalize(H.second);
        // }
        // normalize(_hist2DJetChargeWPt);

    }
    //@}
private:
    ///@param _histograms Indexed by histogram name for easy management
    ///until Rivet Autobooking takes over, allows any number of
    ///histograms to be added "on the fly" in the init() method.
    //@{
    BookedHistos _histograms;
    //AIDA::IHistogram2D *_hist2DJetChargeWPt;
    //@}
    /// @param _nPassing Event count for efficiency studies
    //@{
    int _nPassing[4];
    //@}
};
// The hook for the plugin system
DECLARE_RIVET_PLUGIN(MC_GENSTUDY_JETCHARGE);
}
