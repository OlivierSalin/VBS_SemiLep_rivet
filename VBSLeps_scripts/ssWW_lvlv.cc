// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Tools/Cutflow.hh"
#include "Rivet/Math/MathUtils.hh"
#include <fstream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace Rivet {


  /// @brief Add a short analysis description here
  class ssWW_lvlv : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ssWW_lvlv);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

        std::ifstream json_file("Cuts/ssWW_lvlv_cuts.json");
        _jcuts = json::parse(json_file);
        std::cout << "++++++ to check json 1 var got photon pt min" << _jcuts["abs_diff_m_z"] << "\n";

        _electron_eta_cut = (Cuts::absetaIn(_jcuts["eta_electron"][0][0], _jcuts["eta_electron"][0][1])) ||
                            (Cuts::absetaIn(_jcuts["eta_electron"][1][0], _jcuts["eta_electron"][1][1]));
        _muon_eta_cut = Cuts::absetaIn(0.0, _jcuts["eta_muon"]);

        // The basic final-state projection:
        // all final-state particles within
        // the given eta acceptance
        const FinalState fs;

        // The final-state particles declared above are clustered using FastJet with
        // the anti-kT algorithm and a jet-radius parameter 0.4
        // muons and neutrinos are excluded from the clustering
        FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
        declare(jetfs, "jets");

        // FinalState of direct photons and bare muons and electrons in the event - ignore taus but if want to include use TauFinder
        DirectFinalState bare_e(Cuts::abspid == PID::ELECTRON);
        DirectFinalState bare_mu(Cuts::abspid == PID::MUON);
        DirectFinalState photons_for_dressing(Cuts::abspid == PID::PHOTON);
        DressedLeptons dressed_e(photons_for_dressing, bare_e, 0.1);
        DressedLeptons dressed_mu(photons_for_dressing, bare_mu, 0.1);
        declare(dressed_e, "e_stable");
        declare(dressed_mu, "mu_stable");


        declare(MissingMomentum(), "METFinder");

        // Book histograms
        // plots common with others
        std::ifstream jet_hist_file("Hists_VBSLep/jet_hists.json");
        json jet_hist = json::parse(jet_hist_file);
        for (json::iterator it = jet_hist.begin(); it != jet_hist.end(); ++it) {
            book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
            _hist_names.push_back(it.key());
        }
        //lepton plots
        std::ifstream lep_hist_file("Hists_VBSLep/lepton_hists.json");
        json y_hist = json::parse(lep_hist_file);
        for (json::iterator it = y_hist.begin(); it != y_hist.end(); ++it) {
            book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
            _hist_names.push_back(it.key());
        }
        // plots that are not in other ana
        std::ifstream ana_hist_file("Hists_VBSLep/ssWW_lvlv_hists.json");
        json ana_hist = json::parse(ana_hist_file);
        for (json::iterator it = ana_hist.begin(); it != ana_hist.end(); ++it) {
            book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
            _hist_names.push_back(it.key());
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
        _cutflows.addCutflow("ssWW_lvlv_selections", {"pt_MET","n_lep","lep_charge_dR", "m_ll","n_jets_bjets",
                                                      "m_tagjets","dy_tagjets","two_HS_bosons"});

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
        // save weights before cuts
        double ev_nominal_weight =  event.weights()[0];
        if (ev_nominal_weight>=0){_c["pos_w_initial"]->fill();}
        else {_c["neg_w_initial"]->fill();}

        _cutflows.fillinit();

        const MissingMomentum& METfinder = apply<MissingMomentum>(event, "METFinder");
        const double scalar_MET = METfinder.missingPt()/GeV;
        if (scalar_MET<_jcuts["pt_MET"]) vetoEvent;
        _cutflows.fillnext();

        const FourMomentum fourvec_MET = METfinder.missingMomentum();

        // Retrieve clustered jets, sorted by pT, with a minimum pT cut
        Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 20*GeV && Cuts::absrap < 4.5);
        Jets btagging_jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::absrap < 2.5 && Cuts::pT > 20*GeV);

        // Retrieve dressed leptons, sorted by pT
        Particles e_stable;
        Particles mu_stable;
        e_stable = apply<FinalState>(event, "e_stable").particlesByPt(_electron_eta_cut);
        mu_stable = apply<FinalState>(event, "mu_stable").particlesByPt(_muon_eta_cut);

        // Remove jet if lepton too close
        idiscardIfAny(jets, e_stable, [](const Jet& jet, const Particle& lep) {
            return deltaR(lep, jet) < 0.4 && (lep.pT() / jet.pT() > 0.5);
        });
        // otherwise remove lepton
        idiscardIfAny(e_stable, jets, [](const Particle& lep, const Jet& jet) {
            return deltaR(lep, jet) < 0.4 && (lep.pT() / jet.pT() < 0.5);
        });

        Particles leptons = e_stable + mu_stable;
        idiscardIfAnyDeltaRLess(btagging_jets, leptons, 0.2);
        // sort by pt as not clear if e+m is in order invidually not but necessary together
        std::sort(leptons.begin(), leptons.end(), [](Particle const &a, Particle const &b) {
            return a.pT() > b.pT(); // biggest pT will be first in array
        });

        int nlep_stable = leptons.size();
        int cut_nlep = _jcuts["n_lepton_stable"];
        if (nlep_stable!=cut_nlep)  vetoEvent;
        _cutflows.fillnext();

        const Particle& lep1 = leptons[0];
        const Particle& lep2 = leptons[1];
        if (lep1.pT()<dbl(_jcuts["pt_lepton12"])*GeV) vetoEvent;
        if (lep2.pT()<dbl(_jcuts["pt_lepton12"])*GeV) vetoEvent;
        double eta_ee = _jcuts["eta_electron_ee"];
        if (fabs(lep1.pid())==11 && fabs(lep2.pid())==11 && (lep1.eta()>eta_ee || lep2.eta()>eta_ee)) vetoEvent;


        if (lep1.charge() * lep2.charge() < 0) vetoEvent; // want same charge leptons
        if (deltaR(lep1, lep2) < 0.2)  vetoEvent;
        _cutflows.fillnext();

        const FourMomentum fourvec_ll = lep1.mom() + lep2.mom();
        const double m_ll = fourvec_ll.mass()/GeV;
        if (m_ll<_jcuts["m_ll"]) vetoEvent;
        _cutflows.fillnext();



        const double lep1_abs_pid = lep1.pid();
        const double lep2_abs_pid = lep2.pid();
        double m_z = 91.18;
        if ((lep1_abs_pid == 11 && lep2_abs_pid == 11 && fabs(m_ll - m_z) < _jcuts["abs_diff_m_z"])) vetoEvent;

        int njets = jets.size();
        if (njets < _jcuts["n_jets"])  vetoEvent;
        int nbtags = count(btagging_jets, hasBTag());
        if (nbtags>_jcuts["n_b_jets"]) vetoEvent;
        _cutflows.fillnext();

        const FourMomentum tag1_jet = jets[0].mom();
        const FourMomentum tag2_jet = jets[1].mom();
        if ((tag1_jet.pT()<_jcuts["pt_tagjet1"] || tag2_jet.pT()<_jcuts["pt_tagjet2"])) vetoEvent;

        const double m_tagjets = (tag1_jet + tag2_jet).mass()/GeV;
        if (m_tagjets<_jcuts["m_tagjetsSR"]) vetoEvent;
        _cutflows.fillnext();

        const double dy_tagjets = fabs(tag1_jet.rap() - tag2_jet.rap());
        if (dy_tagjets<_jcuts["dy_tagjets"]) vetoEvent;
        _cutflows.fillnext();

        const double m_T = (fourvec_MET + fourvec_ll).mass()/GeV;

        // do clipping - in case with 2+ w to avoid much work take two w with highest pt
        std::vector<FourMomentum> hs_bosons = {}; // also can be WZ for WZ CR
        const Particles all_particles = event.allParticles();
        for(const Particle& p_rivet : all_particles){
            ConstGenParticlePtr p_hepmc = p_rivet.genParticle();
            int status = p_hepmc->status();
            if (abs(status)==23 or abs(status)==22){
                int i_pid = p_hepmc->pid();
                FourMomentum i_mom = p_hepmc->momentum();
                if (abs(i_pid) == 24 or abs(i_pid) == 23){
                    hs_bosons.push_back(i_mom);
                }
            }
        }
        std::sort(hs_bosons.begin(), hs_bosons.end(), [](FourMomentum const &a, FourMomentum const &b) {return a.pT() > b.pT(); }); // biggest pT will be first in array
        bool have_two_hs_bosons = false;
        double hs_diboson_mass = 0.0;
        int n_z_hs = hs_bosons.size();
        if (n_z_hs > 1){
              hs_diboson_mass = (hs_bosons[0] + hs_bosons[1]).mass() / GeV;
              have_two_hs_bosons = true;
        }
        if (!have_two_hs_bosons) vetoEvent; // just in case reject events where dont have ww somehow
        _cutflows.fillnext();

        //jet plots
        _h["n_jets"]->fill(njets);
        _h["pt_tagjet1"]->fill(tag1_jet.pt());
        _h["pt_tagjet2"]->fill(tag2_jet.pt());
        _h["m_tagjets"]->fill(m_tagjets);
        _h["dy_tagjets"]->fill(dy_tagjets);
        _h["eta_tagjets"]->fill(tag1_jet.eta()); _h["eta_tagjets"]->fill(tag2_jet.eta()); // fill both to the same hists
        _h["dphi_tagjets"]->fill(deltaPhi(tag1_jet,tag2_jet));
        //lepton plots
        _h["n_lepton_stable"]->fill(nlep_stable);
        _h["pt_lepton"]->fill(lep1.pT()); _h["pt_lepton"]->fill(lep2.pT()); // fill both to the same hists
        _h["eta_lepton"]->fill(lep1.eta()); _h["eta_lepton"]->fill(lep2.eta());
        // other
        if (lep1_abs_pid == 11 && lep2_abs_pid == 11){_h["m_ee"]->fill(m_ll);}
        _h["n_b_jet"]->fill(nbtags);
        _h["MET"]->fill(scalar_MET);
        _h["m_T"]->fill(m_T);

        // clipping-related
        _h["m_diboson"]->fill(hs_diboson_mass);
        _h["m_ll_clip_inf"]->fill(m_ll);
        _h["m_tagjets_clip_inf"]->fill(m_tagjets);
        if (ev_nominal_weight>=0){_c["pos_w_final_clip_inf"]->fill();}
        else {_c["neg_w_final_clip_inf"]->fill();}
        for (std::string& i_clip : _clips){
            if (i_clip=="inf") continue; // as done above without cut
            int i_clip_num = std::stoi(i_clip);
            std::string i_pos_c_name = "pos_w_final_clip_" + i_clip;
            std::string i_neg_c_name = "neg_w_final_clip_" + i_clip;
            std::string i_hist_name_m_ll = "m_ll_clip_" + i_clip;
            std::string i_hist_name_m_tagjets = "m_tagjets_clip_" + i_clip;
            if (hs_diboson_mass < i_clip_num) {
                _h[i_hist_name_m_ll]->fill(m_ll);
                _h[i_hist_name_m_tagjets]->fill(m_tagjets);
                if (ev_nominal_weight>=0){_c[i_pos_c_name]->fill();}
                else {_c[i_neg_c_name]->fill();}
            }
        }

    } // end of anylize


    /// Normalise histograms etc., after the run
    void finalize() {
        std::cout << _cutflows.str();

        double pos_w_sum_initial = dbl(*_c["pos_w_initial"]); // from which also number of entries can be obtained
        double neg_w_sum_initial = dbl(*_c["neg_w_initial"]);
        double pos_w_sum_final = dbl(*_c["pos_w_final_clip_inf"]);
        double neg_w_sum_final = dbl(*_c["neg_w_final_clip_inf"]);
        MSG_INFO("\n pos weights initial final ratio " << pos_w_sum_initial <<" " << pos_w_sum_final <<" "<< pos_w_sum_final/pos_w_sum_initial << "\n" );
        MSG_INFO("\n neg weights initial final ratio " << neg_w_sum_initial <<" " << neg_w_sum_final <<" "<< neg_w_sum_final/neg_w_sum_initial << "\n" );

        // normalize all to 1 since in case of mostly negative weights not clear what it will do
        for (auto & i_name : _hist_names){
            normalize(_h[i_name], 1.0);
        }

    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, Histo2DPtr> _h2;
    // map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    std::string _cut_mode;
    Cut _electron_eta_cut;
    Cut _muon_eta_cut;
    json _jcuts;
    Cutflows _cutflows;
    std::vector<std::string> _hist_names;
    std::vector<std::string> _clips {"inf", "3000", "2000", "1500", "1000", "700"};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ssWW_lvlv);

}
