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
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Tools/Cutflow.hh"
#include <fstream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace Rivet {


    /// @brief Add a short analysis description here
    class ZZ_llvv : public Analysis {
    public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ZZ_llvv);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
        _cut_mode = getOption("cut");
        
        std::ifstream json_file("ZZ_llvv_cuts.json");
        _jcuts = json::parse(json_file);
        std::cout << "++++++ to check json 1 var got photon MET_by_HT" << _jcuts["MET_by_HT"] << "\n";

        // The basic final-state projection:
        // all final-state particles within
        const FinalState fs;

        // FinalState of direct photons and bare muons and electrons in the event - ignore taus but if want to include use TauFinder
        DirectFinalState bare_e(Cuts::abspid == PID::ELECTRON);
        DirectFinalState bare_mu(Cuts::abspid == PID::MUON);
        DirectFinalState photons_for_dressing(Cuts::abspid == PID::PHOTON);
        // Dress the bare direct leptons with direct photons within dR < 0.1,
        // and apply some fiducial cuts on the dressed leptons depending on param passed
        DressedLeptons dressed_e(photons_for_dressing, bare_e, 0.1);
        DressedLeptons dressed_mu(photons_for_dressing, bare_mu, 0.1);
        declare(dressed_e, "e_stable");
        declare(dressed_mu, "mu_stable");

        // The final-state particles declared above are clustered using FastJet with
        // the anti-kT algorithm and a jet-radius parameter 0.4
        // muons and neutrinos are excluded from the clustering, also veto electrons(+muons but this is redundant) there

        VetoedFinalState hadrons(FinalState(Cuts::absetaIn(0.0, 10.0)));
        hadrons.addVetoOnThisFinalState(dressed_e);
        hadrons.addVetoOnThisFinalState(dressed_mu);
        declare(hadrons, "hadrons");
        FastJets jetsfs(hadrons, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
        declare(jetsfs, "jets");

        declare(MissingMomentum(), "METFinder");

        // plots common with others
        std::vector<std::string> hists_files {"jet_hists.json", "lepton_hists.json", "ZZ_llvv_hists.json"};
        for (auto & i_hists_file_name : hists_files) {
            std::ifstream i_hist_file(i_hists_file_name);
            json i_hist = json::parse(i_hist_file);
            for (json::iterator it = i_hist.begin(); it != i_hist.end(); ++it) {
                book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
                _hist_names.push_back(it.key());
            }
        }
        //counter for efficiency
        book(_c["pos_w_initial"],"pos_w_initial");
        book(_c["neg_w_initial"],"neg_w_initial");
        for (std::string& i_clip : _clips){ 
            std::string i_pos_c_name = "pos_w_final_clip_" + i_clip;
            std::string i_neg_c_name = "neg_w_final_clip_" + i_clip;
            book(_c[i_pos_c_name], i_pos_c_name);
            book(_c[i_neg_c_name], i_neg_c_name);
        }

        // Cut-flows
        _cutflows.addCutflow("ZZ_llvv_selections", {"two_lep", "lep_pt",
                                                    "m_ll", "dR_ll", "MET", "dphi_MET_pt_ll", "MET_by_HT", "n_jets"});

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
        // save weights before cuts
        double ev_nominal_weight =  event.weights()[0];
        if (ev_nominal_weight>=0){_c["pos_w_initial"]->fill();} // dont need anything in bracket as this will be weight on weight
        else {_c["neg_w_initial"]->fill();}

        _cutflows.fillinit();

        // Retrieve dressed leptons, sorted by pT
        Particles e_stable;
        Particles mu_stable;

        e_stable = apply<FinalState>(event, "e_stable").particlesByPt(Cuts::absetaIn(0.0, _jcuts["eta_lepton"]));
        mu_stable = apply<FinalState>(event, "mu_stable").particlesByPt(Cuts::absetaIn(0.0, _jcuts["eta_lepton"]));
    

        Particles leptons = e_stable + mu_stable; 
        // sort by pt as not clear if e+m is in order invidually not but necessary together
        std::sort(leptons.begin(), leptons.end(), [](Particle const &a, Particle const &b) {
            return a.pT() > b.pT(); // biggest pT will be first in array
            });

        idiscardIfAnyDeltaRLess(e_stable, mu_stable, 0.4);

        int nlep = leptons.size();
        if (nlep!=_jcuts["n_lepton_stable"])  vetoEvent; 
        _cutflows.fillnext();

        if ((leptons[0].pT()<_jcuts["pt_lepton1"] || leptons[1].pT()<_jcuts["pt_lepton2"])) vetoEvent;
        _cutflows.fillnext();
        
        const FourMomentum fourvec_ll = leptons[0].mom() + leptons[1].mom();
        double m_ll = fourvec_ll.mass()/GeV;
        if ((m_ll<_jcuts["m_ll"][0] || m_ll>_jcuts["m_ll"][1])) vetoEvent;
        _cutflows.fillnext();

        double dR_ll = deltaR(leptons[0].rap(),leptons[0].phi(), leptons[1].rap(),leptons[1].phi());
        if (dR_ll>_jcuts["dR_ll"]) vetoEvent;
        _cutflows.fillnext();

        const MissingMomentum& METfinder = apply<MissingMomentum>(event, "METFinder");
        const double MET = METfinder.missingPt()/GeV;
        if (MET<_jcuts["pt_MET"]) vetoEvent;
        _cutflows.fillnext();

        const FourMomentum fourvec_MET = METfinder.missingMomentum();
        double dphi_MET_pt_ll = deltaPhi(fourvec_MET,fourvec_ll);
        if (dphi_MET_pt_ll<_jcuts["dphi_MET_pt_ll"]) vetoEvent;
        _cutflows.fillnext();

        double HT = METfinder.scalarEt();
        double MET_by_HT = MET / HT;
        if (MET_by_HT<_jcuts["MET_by_HT"]) vetoEvent;
        _cutflows.fillnext();

        // // Retrieve clustered jets, sorted by pT, with a minimum pT cut
        Jets jets;
        if (_cut_mode=="SR"){
            jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > dbl(_jcuts["pt_jet"])*GeV && Cuts::absetaIn(0.0, _jcuts["eta_jet"]));
        }
        else {
            jets = apply<FastJets>(event, "jets").jetsByPt();
        }

        int n_jets = jets.size();
        if (n_jets < _jcuts["n_jets"])  vetoEvent;  
        _cutflows.fillnext();

        const FourMomentum tag1_jet = jets[0].mom();
        const FourMomentum tag2_jet = jets[1].mom();
        const double m_tagjets = (tag1_jet + tag2_jet).mass()/GeV;
        const double dy_tagjets =  deltaRap(tag1_jet, tag2_jet); // they dont do the cut

        // do clipping - in case with 2 + z to avoid much work take two Z with highest pt
        std::vector<FourMomentum> hs_bosons_Z = {};
        const Particles all_particles = event.allParticles();
        for(const Particle& p_rivet : all_particles){
            ConstGenParticlePtr p_hepmc = p_rivet.genParticle();
            int status = p_hepmc->status();
            if (abs(status)==23 or abs(status)==22){
                int i_pid = p_hepmc->pid();
                FourMomentum i_mom = p_hepmc->momentum();
                if (abs(i_pid) == 23){hs_bosons_Z.push_back(i_mom);}
            }
        }
        std::sort(hs_bosons_Z.begin(), hs_bosons_Z.end(), [](FourMomentum const &a, FourMomentum const &b) {return a.pT() > b.pT(); }); // biggest pT will be first in array
        bool have_two_hs_bosons = false;
        double hs_diboson_mass = 0.0;
        int n_z_hs = hs_bosons_Z.size(); 
        if (n_z_hs > 1){
            hs_diboson_mass = (hs_bosons_Z[0] + hs_bosons_Z[1]).mass()/GeV;
            have_two_hs_bosons = true;
        }
        if (!have_two_hs_bosons) vetoEvent; // just in case reject events where dont have zz somehow
    
        //jet plots
        _h["n_jets"]->fill(n_jets);
        _h["pt_tagjet1"]->fill(tag1_jet.pt());
        _h["pt_tagjet2"]->fill(tag2_jet.pt());
        _h["eta_tagjets"]->fill(tag1_jet.eta()); _h["eta_tagjets"]->fill(tag2_jet.eta());
        _h["phi_tagjets"]->fill(tag1_jet.phi()); _h["phi_tagjets"]->fill(tag2_jet.phi());
        _h["m_tagjets"]->fill(m_tagjets);
        _h["dy_tagjets"]->fill(dy_tagjets);
        _h["dphi_tagjets"]->fill(deltaPhi(tag1_jet,tag2_jet));
        //lepton plots
        _h["n_lepton_stable"]->fill(nlep);
        _h["pt_lepton"]->fill(leptons[0].pT()); _h["pt_lepton"]->fill(leptons[0].pT()); 
        _h["eta_lepton"]->fill(leptons[1].eta()); _h["eta_lepton"]->fill(leptons[1].eta());    
        //ana-specific
        _h["n_Z_hs"]->fill(n_z_hs);
        _h["m_ZZ"]->fill(hs_diboson_mass);
        _h["m_ll"]->fill(m_ll);
        _h["dR_ll"]->fill(dR_ll);
        _h["dy_ll"]->fill(deltaRap(leptons[0].mom(), leptons[1].mom()));
        _h["dphi_ll"]->fill(deltaPhi(leptons[0].mom(),leptons[1].mom()));
        _h["MET_by_HT"]->fill(MET_by_HT);
        _h["pt_MET"]->fill(MET);
        _h["dphi_MET_pt_ll"]->fill(dphi_MET_pt_ll);
        // clipping-related
        _h["pt_Z_clip_inf"]->fill(fourvec_ll.pT());
        if (ev_nominal_weight>=0){_c["pos_w_final_clip_inf"]->fill();}
        else {_c["neg_w_final_clip_inf"]->fill();}
        for (std::string& i_clip : _clips){
            if (i_clip=="inf") continue; // as done above without cut
            int i_clip_num = std::stoi(i_clip);
            std::string i_pos_c_name = "pos_w_final_clip_" + i_clip;
            std::string i_neg_c_name = "neg_w_final_clip_" + i_clip;
            std::string i_hist_name = "pt_Z_clip_" + i_clip;
            if (hs_diboson_mass < i_clip_num) {
                _h[i_hist_name]->fill(fourvec_ll.pT());
                if (ev_nominal_weight>=0){_c[i_pos_c_name]->fill();}
                else {_c[i_neg_c_name]->fill();}
                }
            }

    }

    /// Normalise histograms etc., after the run
    void finalize() {
        std::string cut_str = _cutflows.str();
        std::string cutflow_file = getOption("OUTDIR") + "/cutflow.txt";
        std::ofstream ofs (cutflow_file, std::ofstream::out); 
        ofs << cut_str;
        ofs.close();

        double pos_w_sum_initial = dbl(*_c["pos_w_initial"]); // from which also number of entries can be obtained
        double neg_w_sum_initial = dbl(*_c["neg_w_initial"]);
        double pos_w_sum_final = dbl(*_c["pos_w_final_clip_inf"]);
        double neg_w_sum_final = dbl(*_c["neg_w_final_clip_inf"]);
        MSG_INFO("\n pos weights initial final ratio " << pos_w_sum_initial <<" " << pos_w_sum_final <<" "<< pos_w_sum_final/pos_w_sum_initial << "\n" );
        MSG_INFO("\n neg weights initial final ratio " << neg_w_sum_initial <<" " << neg_w_sum_final <<" "<< neg_w_sum_final/neg_w_sum_initial << "\n" );

        // normalize all to 1 since in case of mostly negative weights not clear what it will do
        for (auto & i_name : _hist_names){ 
            std::cout << "normalizeing hist " << i_name <<" to 1; " ;
            normalize(_h[i_name], 1.0);
            }
    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, Histo2DPtr> _h2;
    map<string, CounterPtr> _c;
    int _docut;
    json _jcuts;
    Cutflows _cutflows;
    std::vector<std::string> _hist_names;
    std::vector<std::string> _clips {"inf", "3000", "2000", "1500", "1000", "700"};
    Cut _eta_jet_cut;
    std::string _cut_mode;
    /// @}

    };


    RIVET_DECLARE_PLUGIN(ZZ_llvv);

}
