// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


    /// @brief Add a short analysis description here
    class WpWm_lvlv : public Analysis {
    public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(WpWm_lvlv);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

        _docut = 0;
        if (getOption("DOCUT") =="YES") _docut = 1;
        if (getOption("DOCUT") =="NO") _docut = 0;
        std::cout << "received docut " << _docut <<"\n";

        // The basic final-state projection:
        // all final-state particles within
        // the given eta acceptance
        const FinalState fs(Cuts::abseta < 4.5);

        // The final-state particles declared above are clustered using FastJet with
        // the anti-kT algorithm and a jet-radius parameter 0.4
        // muons and neutrinos are excluded from the clustering
        FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
        declare(jetfs, "jets");
        
        // FinalState of direct photons and bare muons and electrons in the event - ignore taus but if want to include use TauFinder
        DirectFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);
        DirectFinalState photons(Cuts::abspid == PID::PHOTON);
        // Dress the bare direct leptons with direct photons within dR < 0.1,
        // and apply some fiducial cuts on the dressed leptons depending on param passed
        Cut lepton_cuts;
        if (_docut==1){lepton_cuts= Cuts::abseta < 2.5 && Cuts::pT > 27.0*GeV;} 
        else{lepton_cuts= Cuts::abseta < 10.0 && Cuts::pT > 0.001*GeV;}
        
        DressedLeptons dressed_leps(photons, bare_leps, 0.1, lepton_cuts);
        declare(dressed_leps, "leptons_stable");

        declare(MissingMomentum(), "METFinder");

        // Book histograms
        int n_nbins = 10;
        int n_pt = 200;
        double max_pt = 3000.0;
        int n_rap = 60;
        double max_rap = 6.0;
        int n_pid = 15;
        //jet plots
        book(_h["n_jet"], "n_jet", n_nbins, 0.0, n_nbins);
        book(_h["n_bjet"], "n_bjet", n_nbins, 0.0, n_nbins);
        book(_h["tagjet1_index"], "tagjet1_index", n_nbins, 0.0, n_nbins);
        book(_h["tagjet2_index"], "tagjet2_index", n_nbins, 0.0, n_nbins);
        book(_h2["tagjets_index"], "tagjets_index", int(n_nbins/2), 0.0, int(n_nbins/2),int(n_nbins/2), 0.0, int(n_nbins/2));
        book(_h["jet3_index"], "jet3_index", n_nbins, 0.0, n_nbins);
        book(_h["pt_tagjet1"], "pt_tagjet1", n_pt, 0.0, max_pt); 
        book(_h["pt_tagjet2"], "pt_tagjet2", n_pt, 0.0, max_pt);
        book(_h["m_tagjets"], "m_tagjets", n_pt, 0.0, max_pt);
        book(_h["dy_tagjets"], "dy_tagjets", n_rap, 0, max_rap);
        book(_h2["m_dy_tagjets"], "m_dy_tagjets", n_pt, 0.0, max_pt, n_rap, 0, max_rap);
        book(_h["eta_tagjet1"], "eta_tagjet1", n_rap, -1*max_rap, max_rap); 
        book(_h["eta_tagjet2"], "eta_tagjet2", n_rap, -1*max_rap, max_rap);
        book(_h["eta_tagjets"], "eta_tagjets", n_rap, -1*max_rap, max_rap);
        book(_h["deta_tagjets"], "deta_tagjets", n_rap, 0, max_rap);
        book(_h["phi_tagjet1"], "phi_tagjet1", n_rap, 0, max_rap);
        book(_h["phi_tagjet2"], "phi_tagjet2", n_rap, 0, max_rap);
        book(_h["dphi_tagjets"], "dphi_tagjets", n_rap, 0, max_rap);
        // //lepton plots
        book(_h2["leptons_pids"], "leptons_pids", 2*n_pid, -1*n_pid, n_pid, 2*n_pid, -1*n_pid, n_pid);
        book(_h["n_lepton_stable"], "n_lepton_stable", n_nbins, 0.0, n_nbins);
        book(_h["lepton_pt"], "lepton_pt", int(n_pt/3), 0.0, max_pt/3);
        book(_h["lepton_eta"], "lepton_eta", n_rap, -1*max_rap, max_rap);
        book(_h["m_ll"], "m_ll", int(n_pt), 0.0, max_pt);
        //other
        book(_h["MET"], "MET", int(n_pt), 0.0, max_pt);
        book(_h["m_T"], "m_T", int(n_pt*2), 0.0, max_pt*2);
        book(_h["centrality"], "centrality", int(n_nbins*2), 0.0, int(n_nbins/2)); 
        book(_h["jet3_centrality"], "jet3_centrality", int(n_nbins*2), 0.0, int(n_nbins/2)); 
        book(_c["found_VBS_pair"],"found_VBS_pair");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

        // Retrieve dressed leptons, sorted by pT
        Particles leptons_stable = apply<FinalState>(event, "leptons_stable").particles();
        int nlep_stable = leptons_stable.size();
        if (nlep_stable!=2)  vetoEvent; // meaning both are e,mu and not tau

        const Particle& lep1 = leptons_stable[0];
        const Particle& lep2 = leptons_stable[1]; 
        if (lep1.charge() == lep2.charge()) vetoEvent; // want opposite charge leptons

        const FourMomentum fourvec_ll = lep1.mom() + lep2.mom(); 
        const double m_ll = fourvec_ll.mass()/GeV;
        if (_docut==1 && m_ll<80.0) vetoEvent; 

        const double lep1_pid = lep1.pid();
        const double lep2_pid = lep2.pid();
        if (fabs(lep1_pid)==fabs(lep2_pid)) vetoEvent;

        // Retrieve clustered jets, sorted by pT, with a minimum pT cut
        Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV);
        Jets btagging_jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::absrap < 2.5 && Cuts::pT > 20*GeV);
        // Remove all jets within dR < 0.2 of a dressed lepton
        idiscardIfAnyDeltaRLess(jets, leptons_stable, 0.2);
        idiscardIfAnyDeltaRLess(btagging_jets, leptons_stable, 0.2);

        int njets = jets.size();
        if (njets < 2 || njets > 3)  vetoEvent;  

        int nbtags = count(btagging_jets, hasBTag());
        if (_docut==1 && nbtags>0) vetoEvent;

        bool foundVBSJetPair = false; // look in opposite hemispheres and pair should have highest mjj
        double max_mjj = 0;
        int tag1_jet_index = -1 ,tag2_jet_index = -1;
        for (int i = 0; i < njets; i++) {
        const Jet& i_jet = jets[i];
            for (int j = 0; j < njets; j++) {
                if (i!=j){
                const Jet& j_jet = jets[j];
                const double mjj = (i_jet.mom() + j_jet.mom()).mass()/GeV;
                const double eta_prod = i_jet.eta()*j_jet.eta();
                if  (eta_prod < 0.0 && mjj>max_mjj){
                    max_mjj = mjj;
                    foundVBSJetPair = true;
                    tag1_jet_index = i;
                    tag2_jet_index = j;
                    }
                }
            }
        }
        if (tag2_jet_index < tag1_jet_index) swap(tag1_jet_index, tag2_jet_index); // organize tag jets by pt  
        if (!foundVBSJetPair)  vetoEvent;

        const FourMomentum tag1_jet = jets[tag1_jet_index].mom();
        const FourMomentum tag2_jet = jets[tag2_jet_index].mom();
        if (_docut==1 && (tag1_jet.pT()<25.0 || tag2_jet.pT()<25.0)) vetoEvent; 

        const double m_tagjets = (tag1_jet + tag2_jet).mass()/GeV;
        if (_docut==1 && m_tagjets<500.0) vetoEvent;

        const double dy_tagjets = fabs(tag1_jet.rap() - tag2_jet.rap());
        
        const double centrality_piece_1 = std::min(lep1.eta(),lep2.eta()) - std::min(tag1_jet.eta(),tag2_jet.eta()); 
        const double centrality_piece_2 = std::max(tag1_jet.eta(),tag2_jet.eta()) - std::max(lep1.eta(),lep2.eta());
        const double centrality = std::min(centrality_piece_1, centrality_piece_2);
        if (_docut==1 && centrality<0.5) vetoEvent;

        double jet3_centrality = -1.0;
        int jet3_index = -1;
        if (njets>=3) {
            // take as 3rd jet highest pt one which is not a tag jet
            for (int k = 0; k < njets; k++) {
                if (k!=tag1_jet_index && k!=tag2_jet_index){
                    jet3_index = k;
                    break;
                }
            }
            if (jet3_index!=-1){
                const FourMomentum jet3 = jets[jet3_index].mom();
                jet3_centrality = fabs(jet3.rap() - 1/2*(tag1_jet.rap()+tag2_jet.rap())/(tag1_jet.rap()-tag2_jet.rap()));
            }
        } 

        const MissingMomentum& METfinder = apply<MissingMomentum>(event, "METFinder");
        const double scalar_MET = METfinder.missingPt()/GeV;
        if (_docut==1 && scalar_MET<15.0) vetoEvent;

        const FourMomentum fourvec_MET = METfinder.missingMomentum();
        const double m_T = (fourvec_MET + fourvec_ll).mass()/GeV;

        //jet plots
        _h["n_jet"]->fill(njets);
        _h["n_bjet"]->fill(nbtags);
        _h["tagjet1_index"]->fill(tag1_jet_index);
        _h["tagjet2_index"]->fill(tag2_jet_index);
        _h2["tagjets_index"]->fill(tag1_jet_index,tag2_jet_index);
        _h["jet3_index"]->fill(jet3_index);
        _h["pt_tagjet1"]->fill(tag1_jet.pt());
        _h["pt_tagjet2"]->fill(tag2_jet.pt());
        // above this worked before
        _h["m_tagjets"]->fill(m_tagjets);
        _h["dy_tagjets"]->fill(dy_tagjets);
        _h2["m_dy_tagjets"]->fill(m_tagjets,dy_tagjets);
        _h["eta_tagjet1"]->fill(tag1_jet.eta());
        _h["eta_tagjet2"]->fill(tag2_jet.eta());
        _h["eta_tagjets"]->fill(tag1_jet.eta()); _h["eta_tagjets"]->fill(tag2_jet.eta()); // fill both to the same hists
        _h["deta_tagjets"]->fill(deltaEta(tag1_jet,tag2_jet));
        _h["phi_tagjet1"]->fill(tag1_jet.phi());
        _h["phi_tagjet2"]->fill(tag2_jet.phi());
        _h["dphi_tagjets"]->fill(deltaPhi(tag1_jet,tag2_jet));
        //lepton plots
        _h2["leptons_pids"]->fill(lep1_pid,lep2_pid); // to check that only have e and mu and of same charge
        _h["n_lepton_stable"]->fill(nlep_stable);
        _h["lepton_pt"]->fill(lep1.pT()); _h["lepton_pt"]->fill(lep2.pT()); // fill both to the same hists
        _h["lepton_eta"]->fill(lep1.eta()); _h["lepton_eta"]->fill(lep2.eta());
        _h["m_ll"]->fill(m_ll);
        // other
        _h["MET"]->fill(scalar_MET);
        _h["m_T"]->fill(m_T);
        _h["centrality"]->fill(centrality);
        _h["jet3_centrality"]->fill(jet3_centrality);
        _c["found_VBS_pair"]->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {

        double veto_survive_sumW = dbl(*_c["found_VBS_pair"]);
        double veto_survive_frac = veto_survive_sumW / sumW();
        std::cout << "survived veto, will norm to this: " << veto_survive_frac << "\n";
        double norm_to = veto_survive_frac*crossSection()/picobarn; // norm to generated cross-section in pb (after cuts)
        
        std::vector<std::string> hist_names_1d = {"n_jet","pt_tagjet1","pt_tagjet2","m_tagjets",
        "dy_tagjets","eta_tagjet1","eta_tagjet2","eta_tagjets", "deta_tagjets", 
        "phi_tagjet1","phi_tagjet2","dphi_tagjets",
        "n_lepton_stable","lepton_pt","lepton_eta",
        "m_ll","MET","m_T","centrality","jet3_centrality","tagjet1_index","tagjet2_index","jet3_index"};       
        for (auto&& i_name : hist_names_1d){ normalize(_h[i_name], norm_to);}
        // also norm few 2d
        normalize(_h2["m_dy_tagjets"], norm_to);
        normalize(_h2["leptons_pids"],norm_to);
        normalize(_h2["tagjets_index"],norm_to);

    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, Histo2DPtr> _h2;
    // map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    int _docut;
    /// @}


    };


    RIVET_DECLARE_PLUGIN(WpWm_lvlv);

}
